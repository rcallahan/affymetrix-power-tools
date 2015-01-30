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

/// @file   MidasEngine.cpp
/// @brief  Iterator which reads input data files, calls
///         the MidasSpliceDetector for each probeset_id (gene) and
///         probeset_list_id (exon) data pair, writes output to file.

//
#include "midas/MidasEngine.h"
//
#include "file/TsvFile/TsvFile.h"
#include "portability/affy-base-types.h"
#include "util/Convert.h"
#include "util/PgOptions.h"
//
#include <cerrno>
#include <iostream>
#include <sstream>
#include <stdexcept>
//

using namespace std;
using namespace affx;

void define_midas_opts(PgOptions* opts)
{
  opts->setUsage("midas - Microarray Detection of Alternative Splicing.\n"
                 "Usage:\n"
                 "   midas --cel-files cels.txt -g gene.summary.txt -e exon.summary.txt -m metaprobeset.txt \\\n"
                 "            -o out -pvalues -fstats -normalized -stabilize 8.0 -no-log-transform");
  
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "",
                     "false");
  opts->defineOption("", "cel-files", PgOpt::STRING_OPT,
                     "File defining mapping between cel_files and experimental group_id.",
                     "");
  opts->defineOption("g", "genedata", PgOpt::STRING_OPT,
                     "File containing gene-level summary data.",
                     "");
  opts->defineOption("e", "exondata", PgOpt::STRING_OPT,
                     "File containing exon-level summary data.",
                     "");
  opts->defineOption("m", "metaprobeset", PgOpt::STRING_OPT,
                     "File containing mapping (meta-probeset file) between exon probesets and gene level probesets.",
                     "");
  opts->defineOption("o", "out-dir", PgOpt::STRING_OPT,
                     "Directory which will contain the midas output.",
                     "");
  opts->defineOption("pv", "pvalues", PgOpt::BOOL_OPT,
                     "Output p values.",
                     "true");
  opts->defineOption("f", "fstats", PgOpt::BOOL_OPT,
                     "Output F statistics.",
                     "false");
  opts->defineOption("n", "normalized", PgOpt::BOOL_OPT,
                     "Output normalized exon intensities.",
                     "false");
  opts->defineOption("s", "stabilize", PgOpt::DOUBLE_OPT,
                     "Stabilization factor for log data transforms.",
                     "8.0");
  opts->defineOption("nol", "no-logtrans", PgOpt::BOOL_OPT,
                     "Do not log transform gene and exon level summary data.",
                     "false");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
  opts->defineOption("", "keep-path", PgOpt::BOOL_OPT,
                     "Keep cel file path.",
                     "false");
}

/**
 *  \brief Constructor.
 *
 *  Input data files have already been opened by the caller -
 *  here we copy parameters, perform minor size and
 *  option checks.  Output file headers have also already
 *  been written - here we open to append.
 *
 *  @param celFiles List of cel_file names.
 *  @param groups List of group indices.
 *  @param metaFileName Meta probeset file name.
 *  @param geneDataFileName Gene data file name.
 *  @param exonDataFileName Exon data file name.
 *  @param pvaluesFileName Output file name for p values.
 *  @param fstatsFileName Output file name for F statistics.
 *  @param normalizedFileName Output file name for normalized exon signal.
 *  @param wantPvalues Caller requests p values.
 *  @param wantFstats Caller requests F statistics.
 *  @param wantNormalized Caller requests normalized exon signal
 *  @param logStabilize Amount to add to data before taking logarithm.
 *  @param noLogTransform Do not log transform data.
 *
 *  Errors: exception if parameter checks fail
 */
midasEngine::midasEngine(const std::vector<std::string>& celFiles,
                         const std::vector<int>& groups,
                         std::string metaFileName,
                         std::string geneDataFileName,
                         std::string exonDataFileName,
                         std::string pvaluesFileName,
                         std::string fstatsFileName,
                         std::string normalizedFileName,
                         int wantPvalues, 
                         int wantFstats,
                         int wantNormalized,
                         const float& logStabilize,
                         int noLogTransform)
    : celFiles (celFiles), groups (groups), metaFileName (metaFileName),
      geneDataFileName (geneDataFileName), exonDataFileName (exonDataFileName),
      pvaluesFileName (pvaluesFileName), fstatsFileName (fstatsFileName),
      normalizedFileName (normalizedFileName), wantPvalues (wantPvalues), wantFstats (wantFstats),
      wantNormalized (wantNormalized), logStabilize (logStabilize), noLogTransform (noLogTransform),
      spliceDetector (groups, wantPvalues, wantFstats, wantNormalized, logStabilize, noLogTransform),
      numGroups (groups.size())
{
  // require at least one output option
  if (! wantPvalues && ! wantFstats && ! wantNormalized)
    Err::errAbort("At least one output option must be selected");

  // celFile list, group list must be of the same size
  if (celFiles.size() != numGroups)
    Err::errAbort("Mismatch between size of cel_file names, group ids");
}

/**
 *  \brief Destructor.
 *
 *  Delete each vector<int> in the map between probesetListId
 *  and its vector of probesetIds.
*/
midasEngine::~midasEngine()
{
  typedef map<int, vector<int>* >::iterator mapVecIter_t;
  mapVecIter_t mapPos;
  for (mapPos = probesetListIdToProbesetIdVectorMap.begin();
      mapPos != probesetListIdToProbesetIdVectorMap.end(); ++mapPos)
    delete mapPos->second;
}

////////////////////

/**
 *  \brief Read meta
 *
 *  Read meta probeset file.  Generate map between probeset_list_ids
 *  (from probeset_list column) and vector of probeset_id
 *  (gene) ids (from probeset_id column).  The correspondence
 *  between probeset_list_id and probeset_id id is potentially
 *  one to many.
 *
 *  Errors: call errAbort if any type mismatch
 *  (require that all id strings can be converted to an int)
 *  or in the case of any error return from NextEntryConverted().
*/
void midasEngine::readMeta()
{
  // printf("readmeta: %s\n",metaFileName.c_str());

  // We open and close the file here
  if (metaTsv.open(metaFileName)!=TSV_OK) {
    string msg = "Problem opening file: " + ToStr(metaFileName);
    Err::errAbort(msg);
  }
  metaTsv.m_optEscapeOk=false;

  //
  typedef const map<int, vector<int>* >::iterator cmapVecIter_t;
  //string stringOfProbesetListIds;
  string probeset_list;
  vector<string> probeset_list_vec;
  int32_t probesetId;
  int row = 0;

  while (metaTsv.nextLevel(0)==TSV_OK) {
    ++row;
    //
    metaTsv.get(0,"probeset_id",probesetId);

    // This is a space seperated list
    metaTsv.get(0,"probeset_list",probeset_list);
    // printf("readmeta: probeset_list='%s'\n",probeset_list.c_str());
    affx::splitstr(probeset_list,' ',probeset_list_vec);
    //
    for (unsigned int i=0;i<probeset_list_vec.size();i++) {
      int probesetListId = -1;
      const char* s_start=probeset_list_vec[i].c_str();
      char* s_end=NULL;
      probesetListId=strtol(s_start,&s_end,0);
      // error converting
      if (s_start==s_end) {
        string msg = "Unable to convert '%s'"+probeset_list_vec[i];
        Err::errAbort(msg);
      }
 
      // has this probesetListId been previously seen?
      cmapVecIter_t probesetListIdToProbesetIdVecPos
        = probesetListIdToProbesetIdVectorMap.find (probesetListId);
      if (probesetListIdToProbesetIdVecPos == probesetListIdToProbesetIdVectorMap.end()) {
        // if not, store probesetId in a new vector
        vector<int>* pVectorInt = new vector<int>;
        (*pVectorInt).push_back (probesetId);
        probesetListIdToProbesetIdVectorMap.insert (make_pair (probesetListId, pVectorInt) );
      }
      else {
        // else a vector<int> already exists for this probesetListId - add to it
        // for now, do not bother checking if this probesetId has already been associated
        // with this probesetListId
        vector<int>* pVectorInt = probesetListIdToProbesetIdVecPos->second;
        (*pVectorInt).push_back (probesetId);
      }
    }
  }
  //
  metaTsv.close();
}

/**
 *  \brief Read, process gene data
 *
 *  Store cel_file data as a vector of vector of floats.
 *  Fill in map with key gene (probeset) id, value gene data row.
 *
 *  Errors: Call errAbort if the gene (probeset) label is a duplicate,
 *  in the case of any error in reading the input file or converting data.
*/
void midasEngine::readGeneData ()
{
  //int rv;
  // We open and close the file here
  if (geneDataTsv.open(geneDataFileName)!=TSV_OK) {
    string msg = "Problem opening file: " + ToStr(geneDataFileName);
    Err::errAbort(msg);
  }
  // printf("engine: readgenedata: %s\n",geneDataFileName.c_str());
  // The genedata file must have probeset_id as the first column. 
  if (geneDataTsv.cname2cidx(0,"probeset_id")!=0) {
    string msg="probeset_id is not the 0th column.";
    Err::errAbort(msg);
  }

  int row=0;
  int geneDataProbesetId;

  while (geneDataTsv.nextLevel(0)==TSV_OK) {
    if (geneDataTsv.get(0,"probeset_id",geneDataProbesetId)!=TSV_OK) {
        Err::errAbort("Bad Conversion: probeset_id");
    }
    //
    if (!geneIdToRowMap.insert(make_pair(geneDataProbesetId,row)).second) {
      // abort if a duplicate
      string msg = "Duplicate gene label " + ToStr (geneDataProbesetId)
       	+ "in file " + ToStr (geneDataFileName) + " at row " + ToStr (row + 1);
      Err::errAbort(msg);
    }
    // printf("gdpid=%d ",geneDataProbesetId); // dbg
    // the number of data fields in gene data is the number of groups
    geneData.push_back(vector<float>(numGroups));
    //
    for (int col=1;col<geneDataTsv.getColumnCount(0);col++) {
      if (geneDataTsv.get(0,col,geneData[row][col-1])!=TSV_OK) {
          Err::errAbort("Bad Conversion: geneData");
      }
      // printf(" %f",geneData[row][col-1]); //dbg
    }
    row++;
    // printf("\n"); //dbg
  }
  geneDataTsv.close();
}

////////////////////

void
midasEngine::openOutput_Pvalues(std::string filename)
{
  pvaluesTsv.defineFile(PROBESET_LIST_ID "\t" PROBESET_ID "\t" "pvalue");
  pvaluesTsv.addHeadersFrom(metatagBucketTsv,TSV_ADD_ALL);
  pvaluesTsv.writeTsv_v1(filename);
}
void
midasEngine::openOutput_Fstats(std::string filename)
{
  fstatsTsv.defineFile(PROBESET_LIST_ID "\t" PROBESET_ID "\t" "fstatistic");
  fstatsTsv.addHeadersFrom(metatagBucketTsv,TSV_ADD_ALL);
  fstatsTsv.writeTsv_v1(filename);
}
void
midasEngine::openOutput_Normalized(std::string filename)
{
  int cidx=0;

  // first two are fixed
  normalizedTsv.defineColumn(0,cidx++,PROBESET_LIST_ID);
  normalizedTsv.defineColumn(0,cidx++,PROBESET_ID);
  // the gene data columns
  for (unsigned int i = 0; i < numGroups; ++i) {
    normalizedTsv.defineColumn(0,cidx++,celFiles[i]);
  }
  normalizedTsv.addHeadersFrom(metatagBucketTsv,TSV_ADD_ALL);
  normalizedTsv.writeTsv_v1(filename);
}

////////////////////

/**
 *  \brief Read, process exon data
 *
 *  Read exon data; for each row, access associated gene data using
 *  the previously constructed maps and write the requested output to file.
 *
 *  Errors: call errAbort for any i/o or format problem in the input file,
 *  if there is no matching gene id from the map.
 *  We also check for an error in the gene id to gene data row map - "should" not happen.
 *  Ofstreams are configured to throw for any file write error - no explicit checking
 *  is done here.
*/
void midasEngine::processExonData ()
{
  typedef const map<int, vector<int>* >::iterator cmapVecIter_t;
  typedef const map<int, int>::iterator mapIter_t;

  int32_t probesetListId;
  float fstat;
  vector<float> normalizedSignal (numGroups);
  vector<float> probesetListData (numGroups);
  float pvalue;
  //int rv;

  if (exonDataTsv.open(exonDataFileName)!=TSV_OK) {
      Err::errAbort("Cant open exon data: "+exonDataFileName);
  }

  int col_cnt=exonDataTsv.getColumnCount(0);

  // open outputfiles
  if (wantPvalues)
    openOutput_Pvalues(pvaluesFileName);
  if (wantFstats)
    openOutput_Fstats(fstatsFileName);
  if (wantNormalized)
    openOutput_Normalized(normalizedFileName);

  while (exonDataTsv.nextLevel(0)==TSV_OK) {
    // suck in data
    exonDataTsv.get(0,"probeset_id",probesetListId);
    for (int col=1;col<col_cnt;col++) {
      /// @todo map names to columns
      if (exonDataTsv.get(0,col,probesetListData[col-1])!=TSV_OK) {
        std::string bad_str;
        exonDataTsv.get(0,col,bad_str);
        Err::errAbort("Bad float conversion: '"+bad_str+"'");
      }
    }
    
    // are any probeset (gene) ids associated with this probeset_list_id?
    cmapVecIter_t geneIdVecPos = probesetListIdToProbesetIdVectorMap.find (probesetListId);
    // skip it if not
    if (geneIdVecPos==probesetListIdToProbesetIdVectorMap.end()) {
      continue;
    }
    vector<int>* pGeneIds = geneIdVecPos->second;
    const size_t geneIdCount = (*pGeneIds).size();
    for (unsigned int i = 0; i < geneIdCount; ++i) {
      const int geneId = (*pGeneIds)[i];
      // map gene id to gene data row
      mapIter_t geneRowIter = geneIdToRowMap.find (geneId);
      // skip if not
      if (geneRowIter==geneIdToRowMap.end()) {
        continue;
      }
      const int geneRow = geneRowIter->second;

      // run the detector
      spliceDetector.runMidasSingle(probesetListData,geneData[geneRow],&fstat,&normalizedSignal,&pvalue);
      
      // dump the data
      if (wantPvalues) {
        pvaluesTsv.set(0,"probeset_list_id",probesetListId);
        pvaluesTsv.set(0,"probeset_id"     ,geneId);
        pvaluesTsv.set(0,"pvalue"          ,pvalue);
        pvaluesTsv.writeLevel(0);
      }
      if (wantFstats) {
        fstatsTsv.set(0,"probeset_list_id",probesetListId);
        fstatsTsv.set(0,"probeset_id"     ,geneId);
        fstatsTsv.set(0,"fstatistic"      ,fstat);
        fstatsTsv.writeLevel(0);
      }
      if (wantNormalized) {
        normalizedTsv.set(0,"probeset_list_id",probesetListId);
        normalizedTsv.set(0,"probeset_id"     ,geneId);
        for (unsigned int g=0;g<numGroups;g++) {
          normalizedTsv.set(0,g+2,normalizedSignal[g]);
        }
        normalizedTsv.writeLevel(0);
      }
    } // for
  } // while
  
  exonDataTsv.close();

  // shut
  if (wantPvalues)
    pvaluesTsv.close();
  if (wantFstats)
    fstatsTsv.close();
  if (wantNormalized)
    normalizedTsv.close();
}

/**
 *  \brief Process gene, exon data
 *
 *  Construct splice detector object, read in, process
 *  gene and exon data.
 *
 *  Errors: none detected here.
*/
void midasEngine::analyze()
{
  readMeta();
  readGeneData();			// read gene data
  processExonData ();			// process exon data
}
