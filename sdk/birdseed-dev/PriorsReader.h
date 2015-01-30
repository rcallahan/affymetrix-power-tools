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

#ifndef _PRIORSREADER_H_
#define _PRIORSREADER_H_

//
#ifndef WIN32
#include <dirent.h>
#endif
#include "birdseed-dev/Clusters.h"
#include "birdseed-dev/Prior.h"
//
#include "broadutil/BroadException.h"
#include "chipstream/SpecialSnps.h"
#include "util/Util.h"
//
#include <map>
#include <memory>
#include <set>
#include <string>
//

// Classes:
// + BasicPriorsReader (virtual)
//    + BasicMapPriorsReader (keeps data in a std::map)
//      + TextPriorsReader   (the broad text   format)
//      + TsvPriorsReader    (the affy  text   format)
//    + BinaryPriorsReader   (the broad binary format)
//
static const char BIRDSEED_SUFFIX_CHAR = 'm';

namespace birdseed
{
namespace dev {
class BasicPriorsReader
{
public:
  // Includes optional gender character.  Is this big enough?
  enum {kMaxPriorNameLen = 31};

  // Let
  struct NamedPriors {
    NamedPriors(size_t numPriors):
        priorName(),
        priors(numPriors) {}

    char priorName[kMaxPriorNameLen+1];
    Priors priors;
  };

  virtual ~BasicPriorsReader() {};

  // The caller should not delete the returned object.
  // The lifetime of the returned object is only guaranteed
  // until the next call to getPrior() or until
  // the PriorsReader is destroyed.  If a caller wants the information
  // for longer than that, a copy must be made.
  virtual const Priors* getPrior(const std::string &priorName) = 0;

  // Set prior iteration back to the beginning
  virtual void startPriorIteration() = 0;
  // Returns a different prior name with each call, in arbitrary order.
  // Returns NULL when iteration is done.  prior name may have "m" appended
  // to it if it is the male version of the prior.
  /// @remarks This interacts with nextPrior().  (It shares the iterator.)
  virtual const char *nextPriorName() = 0;

  /// @brief     Get the next name and prior info
  /// @param     cur_name   pointer to the std::string to modify.
  /// @param     cur_priors pointer to the Priors object to modify
  /// @return    false when at the end of the data.
  /// @remarks   must call startPriorIteration() to get to the beginning.
  virtual bool nextPrior(std::string* cur_name, Priors* cur_priors) = 0;

  /// @brief     Get the next name (without suffix indicating ploidy) and prior info
  /// @param     cur_name   pointer to the std::string to modify.
  /// @param     cur_priors pointer to the Priors object to modify
  /// @return    false when at the end of the data.
  /// @remarks   must call startPriorIteration() to get to the beginning.
  bool nextPriorWithoutPloidySuffix(std::string* cur_name, Priors* cur_priors);

  // This is only guaranteed to return a valid answer if called after startPriorIteration()
  virtual size_t numPriors() = 0;

  /// @brief     Write the collection of priors to a file.
  /// @param     outputFile   path of the output file.
  /// @param     chipType     chipType to write into the file (not used by all formats.)
  /// @remarks   The format of the file will be taken from the ending of the outputFile.
  ///            The endings allowed by this function are:
  ///            ".priors" => Broad Binary.
  ///            ".tsv" => Affy Tsv.
  ///            ".txt" => Broad Text.
  void writePriorsFile(std::string outputFile, std::string chipType);
  /// @brief     Write the collection of priors to a file.
  /// @param     outputFile   path of the output file.
  /// @param     chipType     ignored by this file format.
  void writeTextPriorsFile(std::string outputFile, std::string chipType);
  /// @brief     Write the collection of priors to a file.
  /// @param     outputFile   path of the output file.
  /// @param     chipType     chipType to write into the file.
  void writeTsvPriorsFile(std::string outputFile, std::string chipType);
  /// @brief     Write the collection of priors to a file.
  /// @param     outputFile   path of the output file.
  /// @param     chipType     chipType to write into the file.
  void writeBinaryPriorsFile(std::string outputFile, std::string chipType);
};

/** A class which contains the map used by the Tsv and Text format which subclass it. */
class BasicMapPriorsReader: public BasicPriorsReader
{
protected:
  /// A handy typedef
  typedef std::map<std::string, Priors> PriorsMap;
  /// Both Text and Tsv formats keep their data in a std::map.
  PriorsMap priorsMap;
  /// An iterator over the map data.
  PriorsMap::const_iterator priorIt;
public:
  virtual void startPriorIteration();
  virtual const char *nextPriorName();
  virtual bool nextPrior(std::string* cur_name, Priors* cur_priors);
  virtual const Priors* getPrior(const std::string &priorName);
  virtual size_t numPriors();
};

/** Used to read the Broad Text Priors format. */
class TextPriorsReader: public BasicMapPriorsReader
{
public:
  TextPriorsReader(std::string path);
};

/** Used to read the Affy Tsv Priors format. */
class TsvPriorsReader: public BasicMapPriorsReader
{
public:
  TsvPriorsReader(std::string path, std::set<const char *, Util::ltstr> *probeSetsToLoad = NULL);
  const Priors* getPrior(const std::string &priorName);
};

/** Used to read the Broad Binary Priors format.  Rather
    than using a std::map like Tsv and Text, a counter is
    used to walk the Priors.  This means that all the
    iterator methods have different implementations than
    the other two.
 */
class BinaryPriorsReader: public BasicPriorsReader
{
public:
  static const char *kDefaultChipType;

private:
  // Start of memory-mapped area
  void *pMMap;
  size_t mappedSize;
  size_t countPriors;
  // These are sorted by prior name
  const NamedPriors *allPriors;
  size_t priorIndex;

public:
  BinaryPriorsReader(std::string path, std::string expectedChipType = kDefaultChipType);
  ~BinaryPriorsReader();
  const Priors *getPrior(const std::string &priorName);

  void startPriorIteration();
  bool nextPrior(std::string* cur_name, Priors* cur_priors);
  const char *nextPriorName();
  virtual size_t numPriors();
};

// Object that encapsulates our naming conventions for special (chrX, chrY, chrM, etc.) snps
class SnpPriorKey
{

public:
  // sMap is assumed to have lifetime longer than SnpPriorKey, since we hang onto a reference
  // rather than copying it over and over again.
  SnpPriorKey(const SpecialSnpMap &sMap):
      specialSnps(sMap) {
  }

  bool ploidyDiffersByGender(const char *snpName) const;

  // Return a string based on the number of chromosomes expected given the gender.
  // for example a chrX snp in males has 1 copy so expected is: 'snpname-1' and a
  // chrY snp in females has 0 copies so expected is: 'snpname-0'. regular autosomal
  // snps would be 'snpname-2'
  std::string priorKeyForSnp(const std::string &snpName, Gender gender) const;

private:
  /// Map of snp names to the copies in males and females.
  const SpecialSnpMap &specialSnps;

};

template<class BASICPRIORSREADER, class PRIORS>
class PriorsReaderTemplate
{
public:
  // Note that it is the responsibility of this class to destroy the basicReader
  // when it is done with it.
  // specialSnps is assumed to have lifetime longer than this class, since we hang onto a reference
  // rather than copying it over and over again.
  PriorsReaderTemplate(const SpecialSnpMap &specialSnps, BASICPRIORSREADER *basicReader):
      snpKey(specialSnps),
      basicPriorsReader(basicReader) {
  }

  bool ploidyDiffersByGender(const char *snpName) {
    return snpKey.ploidyDiffersByGender(snpName);
  }

  // The caller should not delete the returned object.
  // The lifetime of the returned object is only guaranteed
  // until the next call to getPriorsForSNP() or until
  // the PriorsReader is destroyed.  If a caller wants the information
  // for longer than that, a copy must be made.
  const PRIORS *getPriorsForSNP(const char *snpName, Gender gender = GENDER_UNKNOWN) {
    return basicPriorsReader->getPrior(getPriorName(snpName, gender).c_str());
  }

  // Look up the key using our "special" (i.e. chrY, chrX, chrM) snps aware object
  std::string getPriorName(const char *snpName, Gender gender = GENDER_UNKNOWN) {
    return snpKey.priorKeyForSnp(snpName, gender);
  }

  size_t numPriors() {
    return basicPriorsReader->numPriors();
  }
private:
  // Object that encapsulates oure naming conventions for special snps.
  SnpPriorKey snpKey;
  std::auto_ptr<BASICPRIORSREADER> basicPriorsReader;
};

typedef PriorsReaderTemplate<BasicPriorsReader, Priors> PriorsReader;

};
};
#endif /* _PRIORSREADER_H_ */

/******************************************************************/
/**************************[END OF PriorsReader.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
