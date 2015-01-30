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
 * @file   DMET3TestCaseCheck.cpp
 * @author Mybrid Spalding
 * @date   Wed Jun  4 11:23:13 PDT 2008
 * @brief  DMET3 test case regression report check 
 */

#include "translation/regression/DMET3TestCaseCheck.h"
//
#include "translation/regression/ComprehensiveReportTableModel.h"
#include "translation/regression/DMET3TestCaseTableModel.h"
//
#include "util/Fs.h"
#include "util/Util.h"
//
#include <cassert>
#include <iostream>
#include <sstream>
//

///////////////////////////////////////////////////////////////////////////////
// DMET3TestCaseCheck.cpp
///////////////////////////////////////////////////////////////////////////////

using namespace std;

/*****************************************************************************/
/**
 * DMET3TestCaseCheck::DMET3TestCaseCheck
 * Synopsis: Default constructor.
 *
 * @param rte - the runtime environment
 * @param testCase - the DMET3TestCaseTableModel record being run
  * @param dmet3MarkerFile - the "generated" file.
 *
 */
/*****************************************************************************/
DMET3TestCaseCheck::DMET3TestCaseCheck(RunTimeEnvironment & rte,
                                       std::string testCase,
                                       std::vector< std::string > & testCaseRow,
                                       DMET3ReportTestCheckFiles & tcf)
  :  m_testCase(testCase) , m_testCaseRow(testCaseRow),  m_rte(rte), m_tcf(tcf)
{


}
// end DMET3TestCaseCheck::DMET3TestCaseCheck
/*****************************************************************************/
/*****************************************************************************/
/**
 * DMET3TestCaseCheck::check
 * Synopsis: virtual function from RegressionCheck. This routine compares
 * the results and expected results with test codes defined in the test case. 
 *
 * @param msg - to be returned the
 *
 * @return - true, if both files are identical. There is no epsilon.
 */
/*****************************************************************************/
bool DMET3TestCaseCheck::check(string & msg)
{

  bool pass = true;

  msg = msg + string("| TEST CASE | ") + m_testCase + string(" | DESCRIPTION | ") + m_testCaseRow[DMET3TestCaseTableModel::LOG_MESSAGE];

  if ( m_testCaseRow[DMET3TestCaseTableModel::TEST_CODES].empty() ) {
    pass = checkAllColumns( msg ) && pass;
  }
  

  msg = msg + string("| TEST CASE | ") + m_testCase + string(" | END TEST CASE |" );

  Verbose::out(ADT_VERBOSE_EXCEPTION, msg );
    
  return pass;

}
// end DMET3TestCaseCheck::check
/*****************************************************************************/
/*****************************************************************************/
/**
 * DMET3TestCaseCheck::checkAllColumns
 *
 * Helper function to check. The default test function where
 * all columns are checked for identical data. 
 * 
 *
 * @param msg - to be returned the
 *
 * @return - true, if both files are identical. There is no epsilon.
 */
/*****************************************************************************/
bool DMET3TestCaseCheck::checkAllColumns( string & msg ) {

  bool pass = true;
  bool okToContinue = true;
  
  map< int, stringstream *>      diffs;
  map< string, int >             diffStats;

  ComprehensiveReportTableModel expResults( m_rte, m_tcf.m_resultFiles["comprehensive"]);

  string testResultsFile = Fs::join(m_rte.m_adtOpts.m_outputDir,
                                    m_testCaseRow[DMET3TestCaseTableModel::BASE_NAME] + "_comprehensive.rpt");

  if (! Fs::isReadable(testResultsFile)) {
    APT_ERR_ABORT( testResultsFile + " file not found.");
  }

  ComprehensiveReportTableModel testResults( m_rte, testResultsFile );

  for ( int row = 0; (row < expResults.size() && okToContinue); row++ ) {

    // PROBE SET
    // If the probe sets don't match, we abandon ship.
    string expProbeSet = testResults.m_rows[row][ testResults.getColumnIndex("Probe Set ID")];
    string testProbeSet = testResults.m_rows[row][ testResults.getColumnIndex("Probe Set ID")];

      // NUMBER OF ROWS
      if ( row >= testResults.size() ) {
        diffs[row] = new stringstream();
        (*diffs[row]) << "### TEST STOPPED because test results end of table encountered at expected results row: " << row;
        okToContinue = false;
        pass = false;
        break;
      }

      if ( expProbeSet != testProbeSet ) {
      diffs[row] = new stringstream();
      (*diffs[row]) << "### TEST STOPPED because expected and test results have different probe sets: " << expProbeSet << " != " << testProbeSet;
      pass = false;
      okToContinue = false;
    }

    for ( int col = 0; (col < expResults.m_rows[row].size() && okToContinue);
                        col++ ) {

      if  (testResults.m_rows[row].size() <= col) {
        testResults.m_rows[row].push_back("");
      }
      if ( expResults.m_rows[row][col] != testResults.m_rows[row][col] ) {
        pass = false;
        if ( diffs.count(row) == 0 ) {

          diffs[row] = new stringstream();

          (*diffs[row]) << "### " << expResults.m_rows[row][ expResults.getColumnIndex("Gene")] << " | " << expProbeSet << " | ERRORS";

        }
        (*diffs[row]) << " | " << expResults.getColumnName(col) << ": expected \"" << expResults.m_rows[row][col] << "\" got \"" << testResults.m_rows[row][col] ;
        if ( testResults.m_rows[row][col].empty() ) {
          (*diffs[row]) << "EMPTY";
        }
        (*diffs[row]) << "\"";
      }
      
    }

  }

  if ( pass ) {
    msg = msg + " | PASSED";
    return pass;
  }

  msg = msg + string(" | FAILED Most likely reason : ") + m_testCaseRow[DMET3TestCaseTableModel::FAILED_MESSAGE] + string(" | ");

  stringstream summarySStr;
  summarySStr << endl;
  summarySStr << "--------------------------------------------------" << endl;
  summarySStr << "*** DMET3 REGRESSION ERROR REPORT *** \n";
  summarySStr << "Test Case: " << m_testCase << endl;
  summarySStr << "--------------------------------------------------" << endl;
  
  map< int, stringstream* >::iterator itISS;

  for (itISS = diffs.begin(); itISS != diffs.end(); itISS++) {
    summarySStr << (*itISS->second).str() << endl;
  }

  summarySStr << "--------------------------------------------------" << endl;

  msg = msg + summarySStr.str();
    
  return pass;

}
// DMET3TestCaseCheck::checkAllColumns
/*****************************************************************************/

