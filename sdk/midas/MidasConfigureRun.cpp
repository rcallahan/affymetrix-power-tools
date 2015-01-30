////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

/// @file   MidasConfigureRun.cpp
/// @brief  Class for configuring and running the MidasEngine 
///         which computes the MidasSpliceDetector statistical 
///         analysis of exon data.

//
#include "midas/MidasConfigureRun.h"
//
#include "midas/MidasEngine.h"
//
#include "util/Convert.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/Util.h"
//
#include <algorithm>
#include <cerrno>
#include <iostream>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>

//
#ifndef WIN32
#include  <unistd.h>
#endif

using namespace std;
using namespace affx;

typedef map<string, string>::iterator mapSSIter_t;
typedef map<string, int>::iterator mapSIIter_t;
typedef const vector<string>::iterator cvecSIter_t;
typedef const vector<string>::const_iterator cvecSCIter_t;

std::string 
MidasFileRoot(const std::string& str)
{
  //const std::string unixSeparator = "/";
  //const std::string windowsSeparator = "\\";
  std::string::size_type pos = 0;
  pos = str.rfind (Fs::pathSep());
  if (pos != std::string::npos)
    // return root based on unix separators if any were found
    return str.substr (pos+1);
  
  pos = str.rfind (Fs::osPathSep());
  if (pos != std::string::npos)
    // else return root based on windows separators if any were found
    return str.substr (pos+1);
  
  // return full name if no separators were found
  return str;
}


/**
 * @brief Constructor.
 *
 * Copy input file names, output directory name, flags,
 * and log stabilization factor to local storage.
 * option checks.
 *
 * @param celsFileName Cels file name.
 * @param geneDataFileName Gene data file name.
 * @param exonDataFileName Exon data file name.
 * @param metaFileName meta probeset file name.
 * @param wantPvalues Caller requests p values.
 * @param wantFstats Caller requests F statistics.
 * @param wantNormalized Caller requests normalized exon signal.
 * @param logStabilize Amount to add to data before taking logarithm.
 * @param commandLine Copy of the command line, for output tag.
 * @param execVersion Version, cvs id string.
 * @param noLogTransform Do not log transform data.
 * @param keepPath Keep cel file path.
 *
 *  Errors: exception if no output chosen
 */
midasConfigureRun::midasConfigureRun (std::string celsFileName, std::string geneDataFileName,
                                      std::string exonDataFileName, std::string metaFileName, 
                                      std::string outputDirectory, const bool wantPvalues, 
                                      const bool wantFstats, const bool wantNormalized,
                                      const float& logStabilize, const std::string& commandLine,
                                      const std::string& execVersion,const bool noLogTransform,
                                      const bool keepPath)
    : celsFileName (celsFileName), geneDataFileName (geneDataFileName),
      exonDataFileName (exonDataFileName),
      metaFileName (metaFileName), outputDirectory (outputDirectory),
      wantPvalues (wantPvalues), wantFstats (wantFstats), wantNormalized (wantNormalized),
      logStabilize (logStabilize), commandLine (commandLine), noLogTransform (noLogTransform),
      comment ("############################################################"),
      execVersion (execVersion), keepPath (keepPath)
{
  // require at least one output option
  if (! wantPvalues && ! wantFstats && ! wantNormalized)
    Err::errAbort("At least one output option must be selected");
  // no warning yet
  pWarningMessage = 0;

  // initialize vector of tags to propagate from gene, exon data files
  const string dataFileTagsIn[] = { "guid", "chip_type", "lib_set_name", "lib_set_version" };
  const int dataSize = sizeof (dataFileTagsIn) / sizeof (dataFileTagsIn[0]);
  const vector<string> dataFileTagsString (dataFileTagsIn, dataFileTagsIn + dataSize);
  dataFileTags = dataFileTagsString;

  //
  groupIdMax=0;
  groupIdsUsed=0;
}

// validate meta probeset file - must contain probeset_id, probeset_list columns
void midasConfigureRun::checkMetaFile()
{
  if (metaTsv.open(metaFileName)!=TSV_OK) {
      Err::errAbort("Open failed: "+string(metaFileName));
  }
  // gotta have probeset_id and probeset_list
  if (metaTsv.cname2cidx(0,PROBESET_ID)<0) {
      Err::errAbort("No '" PROBESET_ID "' in '"+string(metaFileName)+"'");
  }
  if (metaTsv.cname2cidx(0,PROBESET_LIST)<0) {
      Err::errAbort("No '" PROBESET_LIST "' in '"+string (metaFileName)+"'");
  }
  // finished with the meta file
  metaTsv.close();
  
  // propagate all tags, prepend 
  propagateMetaTags(metatagBucketTsv,"meta_probeset_",metaTsv);
}

// read cels.txt file, generate map between cel_files name and group_id string
void midasConfigureRun::readCelsFile()
{
#ifdef WIN32
  celsTsv.m_optEscapeOk = false;
#endif

  std::string cel_file_str;
  std::string group_id_str;

  celsTsv.bind(0, "cel_files", &cel_file_str, TSV_BIND_REQUIRED);
  celsTsv.bind(0, "group_id", &group_id_str, TSV_BIND_REQUIRED);

  //printf("\nopen: celsFileName: %s\n",celsFileName.c_str());
  if (celsTsv.open(celsFileName)!=TSV_OK) {
      Err::errAbort(celsFileName);
  }

  //
  //celsTsv.dump();

  int row=0;

  //
  while (celsTsv.nextLevel(0)==TSV_OK) {
    /// @todo  This adds a "C:"
    if (!keepPath) {
      cel_file_str=MidasFileRoot(cel_file_str);
    }
    //
    if (!celFileGroupIdStringMap.insert(make_pair(cel_file_str,group_id_str)).second) {
      string msg = "Duplicate cel_files entry in '"+string(celsFileName)+"' at row "+ToStr(row);
      msg += ": " + cel_file_str + ", " + group_id_str;
      Err::errAbort(msg);
    }
    row++;
  }
  celsTsv.close();

  //
  propagateMetaTags(metatagBucketTsv,"cel_group_",celsTsv);
}


// Allocate a groupId for the string.
// Remember to set groupIdMax=0.
int
midasConfigureRun::allocGroupId(std::string groupIdStr)
{
  int groupId;

  map<string,int>::iterator i;
  i=groupIdStringToIntMap.find(groupIdStr);
  if (i==groupIdStringToIntMap.end()) {
    groupId=groupIdMax;
    // printf("allocGroupId: g_str='%s' g_id=%d\n",groupIdStr.c_str(),groupIdMax); //dbg
    groupIdStringToIntMap.insert(make_pair(groupIdStr,groupIdMax));
    groupIdMax++;
  } else {
    groupId=(*i).second;
  }
  return groupId;
}

// read genedata file - for each cel_file name found in both the genedata.txt and
// cels.txt files, assign an (int) group_id index, write the cel_file name to
// the geneDataCelFiles vector and the (int) group_id index to the groupId Ints
// vector.  celFilesUsed is the size of both of these vectors.
void midasConfigureRun::readGeneDataFile()
{
  int groupIdInt;
  string groupIdStr;
  string celFile;

  if (geneDataTsv.open(geneDataFileName)!=TSV_OK) {
      Err::errAbort(geneDataFileName);
  }
  // have to have "probeset_id" as 0th column
  if (geneDataTsv.cname2cidx(0,PROBESET_ID)<0) {
      Err::errAbort("No probeset_id");
  }
  // backquotes "\" are allowed as windows seps -- they arent escape chars.
  geneDataTsv.m_optEscapeOk=false;

  //
  int colcnt=geneDataTsv.getColumnCount(0);
  
  for (int col=1;col<colcnt;col++) {
    geneDataTsv.cidx2cname(0,col,celFile);

    if (!keepPath) {
      celFile=MidasFileRoot(celFile);
    }

    // printf("col=%4d %s\n",col,celFile.c_str()); //dbg
    std::map<std::string,std::string>::iterator cfi;
    cfi=celFileGroupIdStringMap.find(celFile);
    if (cfi==celFileGroupIdStringMap.end()) {
      // printf("continue\n");
      continue;
    }
    groupIdInt=allocGroupId(cfi->second);
    celFileGroupIdIntMap.insert(make_pair(celFile,groupIdInt));

    // append this cel_file name and group_id int to the appropriate vectors
    geneDataCelFiles.push_back(celFile);
    groupIdInts.push_back(groupIdInt);
    ++celFilesUsedInGeneData;
    ++groupIdsUsed;
  }
  // done with the file
  geneDataTsv.close();

  //
  propagateMetaTags(metatagBucketTsv,"gene_summary_",geneDataTsv);

  // Check for vaild inputs
  if (groupIdsUsed<2)  {
      Err::errAbort("Need more than 2 groups; "+ToStr(groupIdMax)+" given.");
  }
  if (celFilesUsedInGeneData==groupIdMax) {
      Err::errAbort("need at least one group with more than one data column.");
  }

  //printf("groupIdMax=%d   celFilesUsed=%d\n",groupIdMax,celFilesUsedInGeneData); //dbg
}

// read exon data file header - here we're only validating that the list
// of cel_file names matches those found in the gene data file
void midasConfigureRun::readExonDataFile()
{
  // there?
  if (exonDataTsv.open(exonDataFileName)!=TSV_OK) {
      Err::errAbort("Cant open '"+string(exonDataFileName)+"'");    
  }
  // correct header?
  if (exonDataTsv.cname2cidx(0,PROBESET_ID)<0) {
      Err::errAbort("'" PROBESET_ID "' not found in '"+string(exonDataFileName)+"'");
  }
  exonDataTsv.m_optEscapeOk=false;

  // count of cel_files used (data columns)
  int celFilesUsedInExonData = 0;
  
  string colname;
  string msg;
  int cidx;
  string celFile;
  
  // First column must be PROBESET_ID
  exonDataTsv.cidx2cname(0,0,colname);
  if (colname!=PROBESET_ID) {
    msg="zeroth column is " + colname + " not " PROBESET_ID;
    Err::errAbort(msg);
  }

  int exonColumns=exonDataTsv.getColumnCount(0);

  std::map<int,int> groupIdsUsedMap;
  mapSSIter_t celFileGroupIdStringIter;
  mapSIIter_t celFileGroupIdIntIter;
  //mapSSIter_t celFileIter;

  while (exonDataTsv.nextLevel(0)==affx::TSV_OK) {
    for (cidx=1;cidx<exonColumns;cidx++) {
      exonDataTsv.cidx2cname(0,cidx,celFile);
      // check if the exon data column name is a cel_file name in cels.txt
      celFileGroupIdStringIter.operator = (celFileGroupIdStringMap.find(celFile));
      // if a cel_file name is found in both the cels.txt and exon data files
      // but not the genedata file, warn the user and continue without it
      celFileGroupIdIntIter = celFileGroupIdIntMap.find(celFile);
      if (celFileGroupIdIntIter == celFileGroupIdIntMap.end()) {
        string msg = "Warning: cel_file name " + celFile +
          " is found in the " + string (exonDataFileName) +
          " file but not in the " + string (geneDataFileName) + " file";
        addWarning (msg);
        continue;
      }

      groupIdsUsedMap[celFileGroupIdIntIter->second]=1;
      ++celFilesUsedInExonData;
    }
  }
  // check if any cel_file name was found in the genedata file, not in the exon data file
  if (celFilesUsedInGeneData != celFilesUsedInExonData) {
    // if a cel_file name has been found in genedata file but not the exon data file,
    // warn the user, delete the file name from the geneDataCelFiles vector and
    // the group index from the groupIdInts vector
    const size_t geneDataCelFilesCount = geneDataCelFiles.size();
    msg="";
    for (unsigned int i = 0; i < geneDataCelFilesCount; ++i) {
      // check each cel_file name already found in the genedata file
      if (exonDataTsv.cname2cidx(0,geneDataCelFiles[i])<0) {
        msg+=geneDataCelFiles[i]+" ";
      }
    }
    if (msg!="") {
      msg=string("Warning: mismatch between cel_files in genedata, exon data files.\n"
                 "The following cel_files were found in the genedata file, not in the exon data file:\n")
        +msg+string("\n");
      addWarning(msg);
    }
  }

  // 
  groupIdsUsed=groupIdsUsedMap.size();

  // need to recheck if too few groups or data columns remain
  if (groupIdsUsed < 2) {
      msg = ToStr (groupIdsUsed) + " experimental group(s) identified; need at least two to process.\n";
      Err::errAbort(msg);
  }
  // @todo: add this back in
  if (celFilesUsedInGeneData == groupIdsUsed) {
    msg = "Unable to process this data - need at least one group with more than one data column.\n";
    Err::errAbort(msg);
  }
  
  // all done
  exonDataTsv.close();
  metatagBucketTsv.addHeader("exon_summary_filename",exonDataFileName);
  propagateMetaTags(metatagBucketTsv,"exon_summary_",exonDataTsv);
}

void
midasConfigureRun::addMetaTags(affx::TsvFile& tsvfile)
{
  tsvfile.addHeader("guid",affxutil::Guid::GenerateNewGuid());
  tsvfile.addHeader("exec_guid",execGuid);
  tsvfile.addHeader("exec_version",execVersion);
  tsvfile.addHeader("create_date",timeString);
  tsvfile.addHeader("cmd=",commandLine);
  tsvfile.addHeaderComment("###################");
}

//////////

// If no warning message is already present, generate a new one
// else concatenate the input message to the existing message.
void midasConfigureRun::addWarning (const std::string& message)
{
  // cout << "warning: " << message << "\n"; // dbg

  if (! pWarningMessage)
    pWarningMessage = new string (message);
  else
    *pWarningMessage += "\n" + message;
}

void
midasConfigureRun::throwIfFileExists(std::string fileName)
{
  ifstream inFstream;
  string msg;

  inFstream.open (fileName.c_str(), ios_base::in);
  if (!inFstream.fail()) {
    msg = "Output file " + fileName + " already exists.  Please rename or delete it.";
    Err::errAbort(msg);
  }
}

// Call Err::errAbort on any error condition.
// Return pointer to non-fatal warning message, if any
std::string* 
midasConfigureRun::configure ()
{
  // check if requested output files already exist
  string outDirString (outputDirectory);
  if (outDirString == "")
    outDirString = ".";

  // open for input - abort if output file already exists
  if (wantPvalues) {
    pvaluesFileName = Fs::join(outDirString,PVALUES_OUTPUT);
    throwIfFileExists(pvaluesFileName);
  }
  if (wantFstats) {
    fstatsFileName = Fs::join(outDirString,FSTATS_OUTPUT);
    throwIfFileExists(pvaluesFileName);
  }
  if (wantNormalized) {
    normalizedFileName = Fs::join(outDirString,NORMALIZED_OUTPUT);
    throwIfFileExists(pvaluesFileName);
  }

  // validate meta probeset file
  checkMetaFile();

  // build map from cel_file name to group_id string
  readCelsFile();

  celFilesUsedInGeneData = 0;
  groupIdsUsed = 0;

  // read, validate genedata file
   readGeneDataFile();

  // read, validate exon data file
  readExonDataFile();

  //const unsigned int numGroups = geneDataCelFiles.size();

  //
  execGuid = affxutil::Guid::GenerateNewGuid();

  // set up headers for requested outputs
  timeString=Util::getTimeStamp();

  return pWarningMessage;
}

void
midasConfigureRun::propagateMetaTags(TsvFile& destTsv,std::string prefix,TsvFile& srcTsv)
{
  destTsv.addHeadersFrom(srcTsv,prefix,TSV_ADD_KEYS);
  destTsv.addHeaderComment("###################");
}

void 
midasConfigureRun::deleteOutputs_1(bool want,std::string filename)
{
  if (wantPvalues&&(filename!="")) {
    unlink(filename.c_str());
  }
}

// Delete output files - close if open
void midasConfigureRun::deleteOutputs()
{
  deleteOutputs_1(wantPvalues   ,pvaluesFileName   );
  deleteOutputs_1(wantFstats    ,fstatsFileName    );
  deleteOutputs_1(wantNormalized,normalizedFileName);
 }

// Construct, run the midas engine
// Engine calls Err::errAbort on any error condition.
void midasConfigureRun::run()
{
  // the engine object constructor requires user input file names, cel_file and group_id (int) lists,
  // output file names, flags, and log stabilization factor
  midasEngine engine (geneDataCelFiles, groupIdInts, metaFileName, geneDataFileName, exonDataFileName,
	pvaluesFileName.c_str(), fstatsFileName.c_str(), normalizedFileName.c_str(), wantPvalues,
	wantFstats, wantNormalized, logStabilize, noLogTransform);
  // Add metatags to the bucket...
  addMetaTags(metatagBucketTsv);
  /// @todo  why do we have to copy it to another bucket?
  // If the engine was made at the start we could be adding it directly.
  engine.metatagBucketTsv.addHeadersFrom(metatagBucketTsv,TSV_ADD_ALL);
  // now run 
  engine.analyze();
}
