////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
 * @file   DMET2Check.cpp
 * @author Mybrid Spalding
 * @date   Wed Jun  4 11:23:13 PDT 2008
 * @brief  Allele translation table regression class for DMET2 files.
 */

#include "translation/regression/DMET2Check.h"
//
#include "translation/regression/DMET2Model.h"
//
#include <cassert>
#include <iostream>
#include <sstream>
//

///////////////////////////////////////////////////////////////////////////////
// DMET2Check
///////////////////////////////////////////////////////////////////////////////


/*****************************************************************************/
/**
 * DMET2Check::DMET2Check
 * Synopsis: Default constructor.
 *
 * @param rte - the runtime environment
 * @param dmet2MarkerFile - the "gold" file.
 * @param dmet3MarkerFile - the "generated" file.
 *
 * @return - constructor
 */
/*****************************************************************************/
DMET2Check::DMET2Check(RunTimeEnvironment & rte,
                       const string & dmet2File,
                       const string & dmet3File,
                       const string & genoFile)
    :  m_dmet2File(dmet2File), m_dmet3File(dmet3File),
    m_genoFile(genoFile), m_rte(&rte)
{


}
// end DMET2Check::DMET2Check
/*****************************************************************************/
///////////////////////////////////////////////////////////////////////////////
// DMET2MarkerCheck
/*****************************************************************************/
/**
 * DMET2MarkerCheck::check
 * Synopsis: virtual function from RegressionCheck. This routine compares
 * the DMET2 and DMET3 regresssion files.
 *
 * @param msg - to be returned, the
 *
 * @return - true, if both files are identical. There is no epsilon.
 */
/*****************************************************************************/
bool DMET2MarkerCheck::check(string & msg)
{

  bool pass = true;

  DMET2MarkerModel d2(m_rte, m_dmet2File, true);
  DMET2MarkerModel d3(m_rte, m_dmet3File, false);


  set< string >                  missingExperiments;
  set< string >                  extraExperiments;
  map< string, vector< string> > missingExperimentProbeSets;
  map< string, vector< string> > extraExperimentProbeSets;
  map< int, stringstream *>      diffs;
  map< string, int >             diffStats;
  int                            dmet2MissingMarkersCount  = 0;
  int                            dmet3ExtraMarkersCount    = 0;

  for (int col = 0; col < d2.m_rows[0].size(); col++) {
    diffStats[ d2.getColumnName(col)] = 0;
  }

  map< string, int > dmet2Stats;

  map< string, int > dmet3Stats;
  dmet3Stats["INVALID_LINES"] = 0;


  map< string, string > dmet3DuplicateMarkers;

  // Account for the missing experiments in DMET3 and ignore the
  // individual records as being part of the missing experiment.

  set< string >::iterator itSS;

  for (itSS = d2.m_experiments.begin(); itSS != d2.m_experiments.end();
       itSS++) {

    if (d3.m_experiments.find(*itSS) == d3.m_experiments.end()) {
      pass = false;
      missingExperiments.insert(*itSS);
    }
  }

  map< string, int>::iterator itSI;
  for (itSI = d2.m_egaIndex.begin(); itSI != d2.m_egaIndex.end(); itSI++) {

    string egaKey = itSI->first;
    int d2Row      = itSI->second;

    string experiment = d2.m_rows[d2Row][DMET2MarkerModel::EXPERIMENT];
    string probeSet    = d2.m_rows[d2Row][DMET2MarkerModel::PROBE_SET];

    // Skip over missing experiment records as the entire
    // experiment will be reported as missing.
    if (missingExperiments.find(experiment) != missingExperiments.end()) {
      dmet2MissingMarkersCount++;
      continue;
    }

    map< string, int>::iterator jtSI = d3.m_egaIndex.find(egaKey);
    // Missing line = missing probe set.
    if (jtSI == d3.m_egaIndex.end()) {
      dmet2MissingMarkersCount++;
      missingExperimentProbeSets[experiment].push_back(probeSet);
      continue;
    }

    int d3Row = jtSI->second;

    // The experiment / probe set exists, so compare the columns.
    for (int col = 0; col < d2.m_rows[d2Row].size(); col ++) {

      if (d2.m_rows[d2Row][col] != d3.m_rows[d3Row][col]) {
        int dmet3LineNum = d3.m_fileRowIndex[d3Row];

        pass = false;
        diffStats[ d2.getColumnName(col)] += 1;

        if (diffs.count(dmet3LineNum)  == 0) {
          diffs[dmet3LineNum] = new stringstream();
          (*diffs[dmet3LineNum]) << "### DMET3 Line " << dmet3LineNum << ":";
        }
        (*diffs[dmet3LineNum]) << " column(" << d2.getColumnName(col) << ")  DMET2=(" << d2.m_rows[d2Row][col] << ") DMET3=(" << d3.m_rows[d3Row][col] << ") |";

      }

    }


  }

  for (itSS = d3.m_experiments.begin() ; itSS != d3.m_experiments.end();
       itSS++) {
    if (d2.m_experiments.find(*itSS) == d2.m_experiments.end()) {
      if (extraExperiments.find(*itSS) == extraExperiments.end()) {
        extraExperiments.insert(*itSS);
      }
    }
  }

  for (itSI = d3.m_egaIndex.begin(); itSI != d3.m_egaIndex.end(); itSI++) {

    if (d2.m_egaIndex.find(itSI->first) == d2.m_egaIndex.end()) {

      dmet3ExtraMarkersCount++;
      int d3Row = itSI->second;
      string experiment = d3.m_rows[d3Row][DMET2MarkerModel::EXPERIMENT];
      string probeSet    = d3.m_rows[d3Row][DMET2MarkerModel::PROBE_SET];
      extraExperimentProbeSets[experiment].push_back(probeSet);
      pass = false;
    }
  }

  stringstream summarySStr;

  summarySStr << "--------------------------------------------------\n";
  summarySStr << "*** DMET2 MARKER REGRESSION *** \n";
  summarySStr << m_dmet2File << " <=> " << m_dmet3File << endl;
  summarySStr << "--------------------------------------------------\n";
  summarySStr << "--------------------------------------------------\n";
  summarySStr << "Diffs \n";
  summarySStr << "--------------------------------------------------\n";

  map< int, stringstream* >::iterator itISS;

  for (itISS = diffs.begin(); itISS != diffs.end(); itISS++) {
    summarySStr << (*itISS->second).str() << endl;
  }

  summarySStr << "--------------------------------------------------\n";
  summarySStr << "SUMMARY\n";
  summarySStr << "--------------------------------------------------\n";
  summarySStr << "DMET3 missing experiment total: " << missingExperiments.size()  << endl;
  for (itSS = missingExperiments.begin() ; itSS !=  missingExperiments.end();
       itSS++) {
    summarySStr << *itSS << endl;
  }
  summarySStr << "DMET2 marker             total: " << d2.m_rows.size() << endl;
  summarySStr << "DMET2 missing marker     total: " <<  dmet2MissingMarkersCount << endl;
  summarySStr << "DMET3 extra experiment   total: " << extraExperiments.size() << endl;
  for (itSS = extraExperiments.begin(); itSS !=  extraExperiments.end() ;
       itSS++) {
    summarySStr << *itSS << endl;
  }
  summarySStr << "DMET3 marker             total: " << d3.m_rows.size() << endl;
  summarySStr << "DMET3 extra marker       total: " << dmet3ExtraMarkersCount << endl;
  summarySStr << "DMET2 compared markers   total: " << (d2.m_rows.size()  - dmet2MissingMarkersCount) << endl;
  summarySStr << "Diff                     total: " << diffs.size() << endl;

  for (itSI = diffStats.begin(); itSI != diffStats.end(); itSI++) {
    summarySStr << itSI->first << " diff total: " << itSI->second << endl;
  }
  summarySStr << "--------------------------------------------------\n";

  if (! pass) {

    msg = summarySStr.str();
  }

  return pass;

}
// end DMET2MarkerCheck::check
/*****************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// DMET2HaplotypeCheck
///////////////////////////////////////////////////////////////////////////////

/*****************************************************************************/
/**
 * DMET2HaplotypeCheck::check
 * Synopsis: virtual function from RegressionCheck. This routine compares
 * the DMET2 and DMET3 regresssion files.
 *
 * @param msg - to be returned, the
 *
 * @return - true, if both files are identical. There is no epsilon.
 */
/*****************************************************************************/
bool DMET2HaplotypeCheck::check(string & msg)
{

  bool pass = true;

  DMET2HaplotypeModel d2(m_rte, m_dmet2File);
  DMET2HaplotypeModel d3(m_rte, m_dmet3File);

  set< string >                  missingExperiments;
  set< string >                  extraExperiments;
  map< int, stringstream *>      d2Diffs;
  map< int, stringstream *>      d3Diffs;
  map< string, int >             diffStats;
  set< string >::iterator itSS;
  int dmet2DifferentCallCount = 0;
  int dmet3DifferentCallCount = 0;

  stringstream summarySStr;

  for (int col = 0; col < d2.m_rows[0].size(); col++) {
    diffStats[ d2.getColumnName(col)] = 0;
  }

  // Account for the missing or extra experiments in DMET3 and ignore the
  // individual missing records as being part of the missing experiment.

  for (itSS = d2.m_experiments.begin(); itSS != d2.m_experiments.end();
       itSS++) {

    if (d3.m_experiments.find(*itSS) == d3.m_experiments.end()) {
      pass = false;
      missingExperiments.insert(*itSS);
    }
  }

  for (itSS = d3.m_experiments.begin() ; itSS != d3.m_experiments.end();
       itSS++) {
    map< string, int>::iterator jtSI;
    if (d2.m_experiments.find(*itSS) == d2.m_experiments.end()) {
      if (extraExperiments.find(*itSS) == extraExperiments.end()) {
        extraExperiments.insert(*itSS);
      }
    }
  }

  map< string, int>::iterator itSI;
  for (itSI = d2.m_egcIndex.begin(); itSI != d2.m_egcIndex.end(); itSI++) {

    string egcKey = itSI->first;
    int d2Row     = itSI->second;
    string experiment = d2.m_rows[d2Row][DMET2HaplotypeModel::EXPERIMENT];
    string gene = d2.m_rows[d2Row][DMET2HaplotypeModel::GENE];
    string call = d2.m_rows[d2Row][DMET2HaplotypeModel::CALL];

    // Skip over missing experiment records as the entire
    // experiment will be reported as missing.
    if (missingExperiments.find(experiment) != missingExperiments.end()) {
      continue;
    }

    map< string, int>::iterator jtSI;
    if ((jtSI = d3.m_egcIndex.find(egcKey)) == d3.m_egcIndex.end()) {
      int dmet2LineNum = d2.m_fileRowIndex[d2Row];
      diffStats[ d2.getColumnName(DMET2HaplotypeModel::CALL)] += 1;
      d2Diffs[dmet2LineNum] = new stringstream();
      (*d2Diffs[dmet2LineNum]) << "<<< DMET2 Line " << dmet2LineNum << ": ";
      (*d2Diffs[dmet2LineNum]) << "different Call: " << experiment << " ";
      (*d2Diffs[dmet2LineNum]) << gene << "  \"" << call << "\"";
      dmet2DifferentCallCount++;
      pass = false;
      continue;
    }

    int d3Row = jtSI->second;

    // The experiment / probe set exists, so compare the columns.
    for (int col = 0; col < d2.m_rows[d2Row].size(); col ++) {

      if (d2.m_rows[d2Row][col] != d3.m_rows[d3Row][col]) {
        int dmet3LineNum = d3.m_fileRowIndex[d3Row];

        pass = false;
        diffStats[ d2.getColumnName(col)] += 1;

        if (d3Diffs.count(dmet3LineNum)  == 0) {
          d3Diffs[dmet3LineNum] = new stringstream();
          (*d3Diffs[dmet3LineNum]) << "### DMET3 Line " << dmet3LineNum << ":";
        }
        (*d3Diffs[dmet3LineNum]) << " different " << d2.getColumnName(col) << ": " << experiment << " " << gene << " " <<   "DMET2=(" << d2.m_rows[d2Row][col] << ") DMET3=(" << d3.m_rows[d3Row][col] << ") |";

      }

    }
  }

  for (itSI = d3.m_egcIndex.begin(); itSI != d3.m_egcIndex.end(); itSI++) {

    string egcKey = itSI->first;
    int d3Row     = itSI->second;
    string experiment = d3.m_rows[d3Row][DMET2HaplotypeModel::EXPERIMENT];
    string gene = d3.m_rows[d3Row][DMET2HaplotypeModel::GENE];
    string call = d3.m_rows[d3Row][DMET2HaplotypeModel::CALL];

    // Skip over extra experiment records as the entire
    // experiment will be reported as extra.
    if (extraExperiments.find(experiment) != extraExperiments.end()) {
      continue;
    }

    map< string, int>::iterator jtSI;
    if ((jtSI = d2.m_egcIndex.find(egcKey)) == d2.m_egcIndex.end()) {
      int dmet3LineNum = d3.m_fileRowIndex[d3Row];
      diffStats[ d3.getColumnName(DMET2HaplotypeModel::CALL)] += 1;
      d3Diffs[dmet3LineNum] = new stringstream();
      (*d3Diffs[dmet3LineNum]) << ">>> DMET3 Line " << dmet3LineNum << ": ";
      (*d3Diffs[dmet3LineNum]) << "different Call: " << experiment << " ";
      (*d3Diffs[dmet3LineNum]) << gene << "  \"" << call << "\"";
      dmet3DifferentCallCount++;
      pass = false;
    }
  }

  summarySStr << "--------------------------------------------------\n";
  summarySStr << "*** DMET2 HAPLOTYPE REGRESSION *** \n";
  summarySStr << m_dmet2File << " <=> " << m_dmet3File << endl;
  summarySStr << "--------------------------------------------------\n";
  summarySStr << "--------------------------------------------------\n";
  summarySStr << "Diffs \n";
  summarySStr << "--------------------------------------------------\n";
  map< int, stringstream* >::iterator itISS;

  for (itISS = d2Diffs.begin(); itISS != d2Diffs.end(); itISS++) {
    summarySStr << (*itISS->second).str() << endl;
  }
  for (itISS = d3Diffs.begin(); itISS != d3Diffs.end(); itISS++) {
    summarySStr << (*itISS->second).str() << endl;
  }
  summarySStr << "--------------------------------------------------\n";
  summarySStr << "SUMMARY\n";
  summarySStr << "--------------------------------------------------\n";
  summarySStr << "DMET3 missing experiment total: " << missingExperiments.size() << endl;


  for (itSS = missingExperiments.begin() ; itSS !=  missingExperiments.end();
       itSS++) {
    summarySStr << *itSS << endl;
  }
  summarySStr << "DMET3 extra experiment   total: " << extraExperiments.size() << endl;
  for (itSS = extraExperiments.begin(); itSS !=  extraExperiments.end() ;
       itSS++) {
    summarySStr << *itSS << endl;
  }

  summarySStr << "DMET2 differing call     total: " << dmet2DifferentCallCount << endl;
  summarySStr << "DMET3 differing call     total: " << dmet3DifferentCallCount << endl;
  summarySStr << "Diff                     total: " << d2Diffs.size() + d3Diffs.size() << endl;

  for (itSI = diffStats.begin(); itSI != diffStats.end(); itSI++) {
    summarySStr << itSI->first << " diff total: " << itSI->second << endl;
  }
  summarySStr << "--------------------------------------------------\n";


  if (! pass) {
    msg = summarySStr.str();
  }

  return pass;

}
// end DMET2HaplotypeCheck::check
/*****************************************************************************/

