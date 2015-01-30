////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License (version 2) as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program;if not, write to the
//
// Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

//
#include "birdseed-dev/PriorsReader.h"
//
#include "broadutil/BroadUtil.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cerrno>
#include <string>

#ifndef _MSC_VER
#include <sys/mman.h>
#endif

using namespace std;
using namespace birdseed::dev;

// Version 1 had male priors with "m" appended to SNP name.
// Version 2 has "-2" appended to SNP name for diploid priors, and "-1" appended for haploid.
static const uint32_t kVersion = 2;
static const uint32_t kMagicNumber = 71459;

struct BinaryPriorsHeader {
  enum {kChipTypeLen = 31};
  uint32_t magic;
  uint32_t version;
  char chipType[kChipTypeLen + 1];
  uint32_t numPriors;
  uint32_t priorsFileOffset; // Offset from beginning of file where the priors start.
};

const char *BinaryPriorsReader::kDefaultChipType = "BI_SNP";

bool BasicPriorsReader::nextPriorWithoutPloidySuffix(std::string* cur_name, Priors* cur_priors)
{
  if (!nextPrior(cur_name, cur_priors)) {
    return false;
  }
  assert(endswith(cur_name->c_str(), "-1") || endswith(cur_name->c_str(), "-2"));
  *cur_name = cur_name->substr(0, cur_name->length() - 2);
  return true;
}

void
BasicPriorsReader::writePriorsFile(std::string path, std::string chipType)
{
  if (Util::stringEndsWith(path, ".tsv")) {
    writeTsvPriorsFile(path, chipType);
  } else if (Util::stringEndsWith(path, ".txt")) {
    writeTextPriorsFile(path, chipType);
  } else if (Util::stringEndsWith(path, ".priors")) {
    writeBinaryPriorsFile(path, chipType);
  } else {
    throw BroadException("Cant determine priors file format from the filename.", __FILE__, __LINE__, path.c_str());
  }
}

void
BasicPriorsReader::writeTextPriorsFile(std::string path, std::string chipType)
{
  // We ignore chipType -- there isnt space for it in the header.
  ofstream strm(path.c_str());
  std::string cur_name;
  Priors      cur_priors(MAX_NUM_CLUSTERS);

  startPriorIteration();
  while (nextPriorWithoutPloidySuffix(&cur_name, &cur_priors)) {
    strm << cur_name;
    // Reverse order of priors for writing to file.
    for (int i = cur_priors.getNumPriors() - 1; i >= 0 ; --i) {
      strm << ";" << cur_priors.getPrior(i).tostring();
    }
    strm << endl;
  }
  strm.close();
}

void BasicPriorsReader::writeTsvPriorsFile(std::string path, std::string chipType)
{
  affx::TsvFile tsv;
  // this comes first...
  tsv.defineFile("probeset_id\tcopy_number\tBB\tAB\tAA");
  // add our meta info...
  tsv.addHeader("created-by", "birdseed/PriorsReader.cpp");
  tsv.addHeader("chip-type", chipType);
  // open for writing...
  tsv.writeTsv_v1(path);

  std::string cur_name;
  Priors      cur_priors(MAX_NUM_CLUSTERS);

  startPriorIteration();
  while (nextPriorWithoutPloidySuffix(&cur_name, &cur_priors)) {
    int num_priors = cur_priors.getNumPriors();
    assert((num_priors == 2) || (num_priors == 3));
    int copy_number = (num_priors == 2) ? 1 : 2; // change priors to copy number
    //
    tsv.set(0, "probeset_id", cur_name);
    tsv.set(0, "copy_number", copy_number);
    // Tsv uses whitespace to seperate columns.  Use ","s to join them.
    if (copy_number == 2) {
      tsv.set(0, "AA", cur_priors.getPrior(Priors::AA_INDEX).tostring(","));
      tsv.set(0, "AB", cur_priors.getPrior(Priors::AB_INDEX).tostring(","));
      tsv.set(0, "BB", cur_priors.getPrior(Priors::BB_INDEX).tostring(","));
    } else {
      tsv.set(0, "AA", cur_priors.getPrior(Priors::A_INDEX).tostring(","));
      tsv.set(0, "AB", "null");
      tsv.set(0, "BB", cur_priors.getPrior(Priors::B_INDEX).tostring(","));
    }
    tsv.writeLevel(0);
  }

  tsv.close();
}

void BasicPriorsReader::writeBinaryPriorsFile(std::string path, std::string chipType)
{
  // whoops! too big.
  if (strlen(chipType.c_str()) > BinaryPriorsHeader::kChipTypeLen) {
    throw BroadException("ChipType string too long.", __FILE__, __LINE__);
  }

  // Create the header
  BinaryPriorsHeader header;
  // zero out the header
  memset(&header, 0, sizeof(header));
  header.magic = kMagicNumber;
  header.version =  kVersion;
  strncpy(header.chipType, chipType.c_str(), sizeof(header.chipType));
  header.numPriors = numPriors();
  header.priorsFileOffset = sizeof(header);

  FILE *fp = fopen_check(path.c_str(), "wb");
  fwrite_check(&header, sizeof(header), 1, fp);

  //
  std::string cur_name;
  Priors      cur_priors(MAX_NUM_CLUSTERS);

  startPriorIteration();
  while (nextPrior(&cur_name, &cur_priors) == true) {
    //
    if (cur_name.length() > kMaxPriorNameLen) {
      throw BroadException("Prior name too long", __FILE__, __LINE__, cur_name.c_str());
    }
    NamedPriors namedPriors(MAX_NUM_CLUSTERS);
    // zero out the struct...
    memset(namedPriors.priorName, 0, sizeof(namedPriors));
    // ... and fill it
    strcpy(namedPriors.priorName, cur_name.c_str());
    namedPriors.priors = cur_priors;
    //
    fwrite_check(&namedPriors, sizeof(namedPriors), 1, fp);
  }
  fclose_check(fp);
}

///
/// BasicMapPriorsReader methods.
/// These are used by both Text and Tsv readers.
///
void BasicMapPriorsReader::startPriorIteration()
{
  priorIt = priorsMap.begin();
}

bool BasicMapPriorsReader::nextPrior(std::string* cur_name, Priors* cur_priors)
{
  assert((cur_name != NULL) && (cur_priors != NULL));

  if (priorIt == priorsMap.end()) {
    cur_name->erase();  // clear out to be clean.
    return false;
  }
  // set the output vars and say ok.
  (*cur_name) = (*priorIt).first;
  (*cur_priors) = (*priorIt).second;
  priorIt++;
  return true;
}

const Priors* BasicMapPriorsReader::getPrior(const std::string &priorName)
{
  PriorsMap::const_iterator it = priorsMap.find(priorName);
  if (it == priorsMap.end()) {
    Verbose::warn(2, string("Could not find prior in text priors file: ") + priorName);
    return NULL;
  }
  return &(*it).second;
}

const char* BasicMapPriorsReader::nextPriorName()
{
  if (priorIt == priorsMap.end()) {
    return NULL;
  }
  return (*priorIt++).first.c_str();
}

size_t BasicMapPriorsReader::numPriors()
{
  return priorsMap.size();
}

///
/// Text reader
///

TextPriorsReader::TextPriorsReader(std::string path)
{
  FILE *fp = fopen_check(path.c_str(), "r");
  char buf[10000];

  while (fgets(buf, sizeof(buf), fp) != NULL) {
    if (strlen(buf) >= sizeof(buf) - 1) {
      throw BroadException("Line too long in text priors file.", __FILE__, __LINE__, path.c_str());
    }
    const char *priorName = strtok(buf, ";");
    if (priorName == NULL || strlen(priorName) == 0) {
      throw BroadException("Error parsing text priors file.",  __FILE__, __LINE__, path.c_str());
    }
    Priors priors(MAX_NUM_CLUSTERS);
    // Finny produces priors in the reverse order of what Josh expects.
    for (int i = MAX_NUM_CLUSTERS - 1; i >= 0; --i) {
      const char *priorsString = strtok(NULL, ";");
      if (priorsString == NULL) {
        if (i == 0) {
          priors.setNumPriors(2);
          // JHG: "down" dont you mean "back"? And what is with the numeric constants?
          // Slide the priors down because we were populating them backwards.
          priors.setPrior(priors.getPrior(1), 0);
          priors.setPrior(priors.getPrior(2), 1);
          break;
        }
        throw BroadException("Error parsing text priors file.",  __FILE__, __LINE__, path.c_str());
      }

      priors.setPrior(Prior(priorsString), i);
    }
    assert(priors.getNumPriors() == 2 || priors.getNumPriors() == 3);
    // Mark the prior as diploid or haploid
    string strPriorName = priorName;
    if (strPriorName[strPriorName.length() - 1] == BIRDSEED_SUFFIX_CHAR) {
      strPriorName = strPriorName.substr(0, strPriorName.length() - 1);
    }
    if (priors.getNumPriors() == 2) {
      strPriorName += "-1";
    } else {
      strPriorName += "-2";
    }

    if (!priorsMap.insert(PriorsMap::value_type(strPriorName, priors)).second) {
      throw BroadException("Prior seen more than once text priors file.",  __FILE__, __LINE__, priorName);
    }
  }
  if (ferror(fp)) {
    throw BroadException("Error reading text priors file.", __FILE__, __LINE__, path.c_str(), errno);
  }
  fclose(fp);
}

///
/// Tsv reader functions.
///

///
/// @brief     A helper function to decode a prior string into a prior object
/// @param     prior     the prior to modify
/// @param     str
static void
prior_set_from_string(Prior* prior, std::string& str)
{
  assert(prior != NULL);
  // The string from the tsv file looks like:
  // "0.41626,0.20745,0.00303,0.00006,0.00066,210"
  // use an int to match "%u".
  int numObs;
  int rv = sscanf(str.c_str(), "%lf,%lf,%lf,%lf,%lf,%u",
                  &prior->m_mean[0], &prior->m_mean[1],
                  &prior->m_covarMatrix[0][0], &prior->m_covarMatrix[1][0], &prior->m_covarMatrix[1][1],
                  &numObs);
  // copy the int to the uint32_t
  prior->m_numObservations = numObs;
  // ok?
  if (rv != 6) {
    throw BroadException("Error parsing TSV prior", __FILE__, __LINE__, str.c_str());
  }

  // the diagonal is the same and not stored.
  prior->m_covarMatrix[0][1] = prior->m_covarMatrix[1][0];
}

const Priors* TsvPriorsReader::getPrior(const std::string &priorName)
{
  PriorsMap::const_iterator it = priorsMap.find(priorName);
  if (it == priorsMap.end()) {
    // If this prior name has 0 copies (i.e. suffix of "-0") then
    // return null as prior file doesn't have entries for 0 copy snps
    // (i.e. chrY snps in females)
    int length = priorName.length();
    if (length > 2 && priorName[length-2] == '-' && priorName[length-1] == '0') {
      return NULL;
    }
    // Didn't find a prior complain
    Verbose::warn(2, "No prior for name: '" + priorName + "'. Perhaps you have the wrong model file?");
    return NULL;
  }
  return &(*it).second;
}

TsvPriorsReader::TsvPriorsReader(std::string path, std::set<const char *, Util::ltstr> *probeSetsToLoad)
{
  affx::TsvFile tsv;
  std::string f_probeset_id;
  int f_copy_number;
  string chip;
  if (tsv.open(path) != affx::TSV_OK) {
    string msg = "Couldn't open file: '" + path + "' to read.";
    throw BroadException(msg.c_str(), __FILE__, __LINE__);
  }

  // check the fields we what are there.
  // we allow "id" as a alias for "probeset_id"
  int probeset_id_cidx = tsv.cname2cidx(0, "probeset_id", "id");
  if (probeset_id_cidx < 0) {
    throw BroadException("Missing required column: ", __FILE__, __LINE__, "probeset_id");
  }
  int copy_number_cidx = tsv.cname2cidx(0, "copy_number");
  if (copy_number_cidx < 0) {
    throw BroadException("Missing required column: ", __FILE__, __LINE__, "copy_number");
  }
  // allow for case in column names.
  int aa_cidx = tsv.cname2cidx(0, "AA", "aa");
  if (aa_cidx < 0) {
    throw BroadException("Missing required column: ", __FILE__, __LINE__, "AA");
  }
  int ab_cidx = tsv.cname2cidx(0, "AB", "ab");
  if (ab_cidx < 0) {
    throw BroadException("Missing required column: ", __FILE__, __LINE__, "AB");
  }
  int bb_cidx = tsv.cname2cidx(0, "BB", "bb");
  if (bb_cidx < 0) {
    throw BroadException("Missing required column: ", __FILE__, __LINE__, "BB");
  }

  // fill the map from the TsvFile...
  while (tsv.nextLevel(0) == affx::TSV_OK) {
    Priors priors(MAX_NUM_CLUSTERS);

    tsv.get(0, probeset_id_cidx, f_probeset_id);
    tsv.get(0, copy_number_cidx, f_copy_number);
    // Skip probesets we aren't requested to load. If probeSetsToLoad is NULL read everything.
    if (probeSetsToLoad != NULL) {
      string name = f_probeset_id;
      if (name[name.length() - 1] == BIRDSEED_SUFFIX_CHAR) {
        name = name.substr(0, name.length()  - 1);
      }
      if (probeSetsToLoad != NULL && probeSetsToLoad->find(name.c_str()) == probeSetsToLoad->end())
        continue;
    }

    // only two values (1 or 2) allowed.
    if (f_copy_number == 1) {
      priors.setNumPriors(2);
    } else if (f_copy_number == 2) {
      priors.setNumPriors(3);
    } else { // opps!
      std::string err_str;
      tsv.get(0, copy_number_cidx, err_str);
      err_str = "'" + err_str + "' at " + path + ":" + ToStr(tsv.lineNumber());
      throw BroadException("Bad value for copynumber : ", __FILE__, __LINE__, err_str.c_str());
    }

    //
    std::string f_aa;
    Prior p_aa;
    tsv.get(0, aa_cidx, f_aa);
    prior_set_from_string(&p_aa, f_aa);
    priors.setPrior(p_aa, Priors::AA_INDEX); // AA_INDEX = 0, A_INDEX = 0
    //
    if (f_copy_number == 2) {
      std::string f_ab;
      Prior p_ab;
      tsv.get(0, ab_cidx, f_ab);
      prior_set_from_string(&p_ab, f_ab);
      priors.setPrior(p_ab, Priors::AB_INDEX); // AB_INDEX = 1
    }
    //
    std::string f_bb;
    Prior p_bb;
    tsv.get(0, bb_cidx, f_bb);
    prior_set_from_string(&p_bb, f_bb);
    priors.setPrior(p_bb, ((f_copy_number == 2) ? Priors::BB_INDEX : Priors::DiploidIndex(Priors::B_INDEX))); // BB_INDEX = 2, B_INDEX = 1;

    // to be compatible with legacy files chop off the 'm' if it exists.
    string name = f_probeset_id;
    if (name[name.length() - 1] == BIRDSEED_SUFFIX_CHAR) {
      name = name.substr(0, name.length()  - 1);
    }
    // Add to the map
    if (!priorsMap.insert(PriorsMap::value_type(name + "-" + ToStr(f_copy_number), priors)).second) {
      throw BroadException("Prior seen more than once text priors file.",  __FILE__, __LINE__, f_probeset_id.c_str());
    }
  }
  //
  tsv.close();
}

///
/// Binary prior reader funcs...
///

BinaryPriorsReader::BinaryPriorsReader(std::string path, std::string expectedChipType):
    pMMap(NULL)
{
  // We cant handle Binary format files on the sparc.
#ifdef __sparc__
  printf("Binary format files cant be handled on a sparc at this time.\n"
         "Please ask a programmer to read 'birdseed/Prior.h' to find out why.\n");
  assert(0);
#endif

  FILE *fp = fopen_check(path.c_str(), "rb");
  BinaryPriorsHeader header;
  fread_check(&header, sizeof(header), 1, fp);
  if (header.magic != kMagicNumber) {
    throw BroadException("Unexpected magic number in binary priors file.", __FILE__, __LINE__, path.c_str());
  }
  if (header.version != kVersion) {
    throw BroadException("Unexpected version number in binary priors file.", __FILE__, __LINE__, path.c_str());
  }

  // dont test if expectedChipType is blank.  (Due to force)
  if (strcmp(expectedChipType.c_str(), "") != 0) {
    if (strcmp(header.chipType, expectedChipType.c_str()) != 0) {
      std::stringstream strm;
      strm << "Actual chip type " << header.chipType << " different from expected chip type " << expectedChipType;
      throw BroadException(strm.str().c_str(), __FILE__, __LINE__);
    }
  }

  countPriors = header.numPriors;
  // Just mmap the whole file.  If the header gets too big we might change this
  mappedSize = header.priorsFileOffset + countPriors * sizeof(NamedPriors);

#ifndef _MSC_VER
  // NYI -- Linux-specific
  pMMap = mmap(NULL, mappedSize, PROT_READ, MAP_SHARED, fileno(fp), 0);
  if (pMMap == (void *) - 1) {
    throw BroadException("mmap error.", __FILE__, __LINE__, path.c_str(), errno);
  }
  allPriors = (const NamedPriors *)((char *)pMMap + header.priorsFileOffset);
#endif

  // File isn't closed in event of exception, but that is a fatal error anyway.
  // Does the file need to stay open for the duration of the mmap?
  fclose(fp);
}

BinaryPriorsReader::~BinaryPriorsReader()
{
#ifndef _MSC_VER
  // NYI -- Linux-specific
#if !NDEBUG
  int ret =
#endif
    munmap(pMMap, mappedSize);
  assert(ret == 0);
#endif
}

static int namedPriorsCompare(const void *v1, const void *v2)
{
  BinaryPriorsReader::NamedPriors *p1 = (BinaryPriorsReader::NamedPriors *)v1;
  BinaryPriorsReader::NamedPriors *p2 = (BinaryPriorsReader::NamedPriors *)v2;

  return strcmp(p1->priorName, p2->priorName);
}


const Priors *BinaryPriorsReader::getPrior(const std::string &priorName)
{
  NamedPriors key(MAX_NUM_CLUSTERS);
  strncpy(key.priorName, priorName.c_str(), sizeof(key.priorName));
  // Shouldn't be necessary, but just to make sure
  key.priorName[sizeof(key.priorName) - 1] = '\0';

  NamedPriors *prior = (NamedPriors *)bsearch(&key, allPriors, countPriors, sizeof(NamedPriors), namedPriorsCompare);
  if (prior == NULL) {
    Verbose::warn(2, string("Could not find prior in binary priors file: ") + priorName);
    return NULL;
  }
  return &prior->priors;
}

// These methods override the Basic ones.
// The binary uses a simple index to step through the contents.
void BinaryPriorsReader::startPriorIteration()
{
  priorIndex = 0;
}

bool BinaryPriorsReader::nextPrior(std::string* cur_name, Priors* cur_priors)
{
  assert((cur_name != NULL) && (cur_priors != NULL));

  if (priorIndex >= countPriors) {
    cur_name->erase();
    return false;
  }
  // found one.
  (*cur_name) = allPriors[priorIndex].priorName;
  (*cur_priors) = allPriors[priorIndex].priors;
  priorIndex++;
  return true;
}

const char *BinaryPriorsReader::nextPriorName()
{
  if (priorIndex >= countPriors) {
    return NULL;
  }
  return allPriors[priorIndex++].priorName;
}

size_t BinaryPriorsReader::numPriors()
{
  return countPriors;
}

bool SnpPriorKey::ploidyDiffersByGender(const char *snpName) const
{
  SpecialSnpMap::const_iterator iter = specialSnps.find(snpName);
  // If not in specialSnpMap, assumed to be diploid for both males & females
  if (iter == specialSnps.end()) {
    return false;
  }
  const std::pair<int, int> &p = iter->second;
  return (p.first != p.second);
}
// Return a string based on the number of chromosomes expected given the gender.
// for example a chrX snp in males has 1 copy so expected is: 'snpname-1' and a
// chrY snp in females has 0 copies so expected is: 'snpname-0'. regular autosomal
// snps would be 'snpname-2'
std::string SnpPriorKey::priorKeyForSnp(const std::string &snpName, Gender gender) const
{
  std::string key;
  SpecialSnpMap::const_iterator iter = specialSnps.find(snpName);
  if (iter == specialSnps.end()) {
    key = snpName + "-2"; // "-2" for diploid by default.
    return key;
  } else {
    const std::pair<int, int> &p = iter->second;
    // same for men and women so just append it (i.e. mitochondiral snps)
    if (p.first == p.second) {
      key = snpName + "-" + ToStr(p.first);
      return key;
    }
    // give female key if unknown or female
    else if (gender == GENDER_FEMALE || gender == GENDER_UNKNOWN) {
      key = snpName + "-" + ToStr(p.second);
    }
    // give male key if male
    else if (gender == GENDER_MALE) {
      key = snpName + "-" + ToStr(p.first);
    } else {
      std::string msg1 = "Don't recognize gender: '" + ToStr(gender) + "'";
      std::string msg2 = " for snp: '" + snpName + "'";
      throw BroadException(msg1.c_str(), __FILE__, __LINE__, msg2.c_str());
    }
  }
  if (key.empty()) {
    std::string msg1 = " for snp: '" + snpName + "'";
    throw BroadException("Error. Can't determine key.", __FILE__, __LINE__, msg1.c_str());
  }
  return key;
}
/******************************************************************/
/**************************[END OF PriorsReader.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
