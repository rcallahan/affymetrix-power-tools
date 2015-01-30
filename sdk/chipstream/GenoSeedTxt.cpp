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
/**
 * @file   GenoSeedTxt.cpp
 * @author Chuck Sugnet
 * @date   Fri Jun 30 17:31:22 2006
 * 
 * @brief Implements interface for a genotyping oracle used for
 * genotyping seeds in brlmm. In this case the genotype seed calls are
 * simply read from a text file.
 * 
 */

//
#include "chipstream/GenoSeedTxt.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Fs.h"

using namespace std;
using namespace affx;

/** 
 * Open up a text file and read the genotypes from it. The chrX
 * snps will be used to determine the het call rate on chrX. If
 * supplied.
 * 
 * @param fileName - Name of text file to read genotypes from.
 * @param layout - Physical description of chip and probesets.
 * @param chrXSnps - Which snps are on chrX.
 */
GenoSeedTxt::GenoSeedTxt(const std::string& fileName, ChipLayout &layout, 
                         std::vector<std::string> &colName,
                         std::map<std::string,bool> &chrXSnps) {
  m_ChipCount = 0;
  readInPrecompGenoTypes(fileName, colName, chrXSnps);

  m_GenderName = "supplied-genotypes-chrX-het-rate";
  m_GenderDescription = "chrX het rates from supplied genotypes";
  
}

/** 
 * Read in the precomputed genotypes from a tab delimited file. Must
 * be a column for each Celfile with same header name as celfile and an row
 * for every snp to be analyzed. There can be extra columns and order of columns
 * is not important as the resulting order will be the same as the colNames
 * vector.
 * 
 * @param fileName - path to text file to read genotypes from.
 * @param colNames - Name of the columns we're expecting to find in this file.
 * @param gtypes - Map to store the genotypes in.
 */
void GenoSeedTxt::readInPrecompGenoTypes(const std::string& fileName, 
                                         std::vector<std::string> &colNames,
                                         std::map<std::string,bool> &chrXSnps) {
  int i;
  affx::TsvFile tsv;
  m_ChipCount = colNames.size();
  vector<GType> calls(colNames.size());
  vector<int> intCalls(colNames.size(), -1);
//  vector<uint64_t> hetChrxCalls(colNames.size(), 0);
//  vector<uint64_t> hetChrxPossibleCalls(colNames.size(), 0);
//  vector<uint64_t> madeCalls(colNames.size(), 0);
//  vector<uint64_t> possibleCalls(colNames.size(), 0);
  //
  vector<uint32_t> hetChrxCalls(colNames.size(), 0);
  vector<uint32_t> hetChrxPossibleCalls(colNames.size(), 0);
  vector<uint32_t> madeCalls(colNames.size(), 0);
  vector<uint32_t> possibleCalls(colNames.size(), 0);
  
  string probeset;
  m_KnownCalls.clear();
  /* Expecting one column for the id and one for each cel file. */
  tsv.bind(0,"probeset_id", &probeset, TSV_BIND_REQUIRED);
  Verbose::out(1, "Reading in seed genotypes from " + ToStr(fileName));
  if(tsv.open(fileName) != TSV_OK) {
    Err::errAbort("Couldn't open file: " + ToStr(fileName) + " to read.");
  }
  // Read CEL/sample column names from hints file
  int columnCount = tsv.getColumnCount(0);
  std::set<std::string> hintSampleNames;
  for (i = 1; i < columnCount; i++) { // i = 0 : probeset_id
    hintSampleNames.insert(tsv.getColumnName(0, i));
  }
  // compare hints file column names to CELfile names
  std::string basename;
  for(i = 0; i < colNames.size(); i++) {
    basename = Fs::basename(colNames[i]);
    if (hintSampleNames.find(basename) != hintSampleNames.end()) {
      tsv.bind(0, basename.c_str(), &intCalls[i], TSV_BIND_REQUIRED);        
    }
    else {
      Verbose::out(3, "Couldn't find hints for " + basename + " in file " + fileName);
    }
  }

  m_KnownCalls.clear();
  while(tsv.nextLevel(0) == TSV_OK) {
    map<string,bool>::iterator iter = chrXSnps.find(probeset);
    bool isChrX = (iter != chrXSnps.end() && iter->second == true);
    if(m_KnownCalls.find(probeset) != m_KnownCalls.end()) {
      Err::errAbort("Found probeset: " + probeset + " specified multiple times in file: " + ToStr(fileName));
    }
    
    /* convert to char as tsv doesn't bind char. */
    for(uint32_t i = 0; i < intCalls.size(); i++) {
      // this converts -1 to NN, 0 to AA, 1 to AB, 2 to BB
      calls[i] = GType_from_int(intCalls[i]);
      /* keep track of call rate. */
      possibleCalls[i]++;
      if(calls[i] != NN) {
        madeCalls[i]++;
      /* keep track of het rate on chrx for gender determination. */
        if(isChrX) {
          hetChrxPossibleCalls[i]++;
          if(calls[i] == AB) {
            hetChrxCalls[i]++;
          }
        }
      }
    }
    m_KnownCalls[probeset] = calls;
  }
  /* Calculate and save our per chip call rates and het rate on chrX */
  m_CallRates.clear();
  m_HetChrXCallRates.clear();
  for(size_t chipIx = 0; chipIx < m_ChipCount; chipIx++) {
    float callRate = (float) madeCalls[chipIx]/possibleCalls[chipIx];
    m_CallRates.push_back(callRate);
    if(hetChrxPossibleCalls[chipIx] > 0) {
      float hetCalls = (float) hetChrxCalls[chipIx] / hetChrxPossibleCalls[chipIx];
      m_HetChrXCallRates.push_back(hetCalls);
    }
    else {
      m_HetChrXCallRates.push_back(-1);
    }
  }
  tsv.close();
}
