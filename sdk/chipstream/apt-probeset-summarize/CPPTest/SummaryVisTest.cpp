////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

/// @file   SummaryVisTest.cpp
/// @brief  Cppunit class for testing SummaryVis.

#include "chipstream/apt-probeset-summarize/CPPTest/SummaryVisTest.h"
#include "chipstream/apt-probeset-summarize/SummaryVis.h"
#include "util/Fs.h"
#include "util/TextFileCheck.h"
#include "util/Util.h"
//

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SummaryVisTest );

using namespace std;

void SummaryVisTest::testSummaryVis()
{
  const char* programName = "apt-summary-vis";
  const char* outputFile = "./output/summaryvis.txt";
  const char* genomePosFile = "./data/HuEx-1_0-st-v2.annot.hg16.subset.csv";
  const char* probesetIdFile1 = "./data/probeset-list10.txt";
  const char* probesetIdFile2 = "./data/probeset-list75.txt";
  const char* transcriptClusterIdFile1 = "./data/transcript-cluster-id-list10.txt";
  const char* transcriptClusterIdFile2 = "./data/transcript-cluster-id-list75.txt";
  const char* summaryFile = "./data/exon.summary.txt";
  const char* summary0File = "./data/exon.summary0.txt";
  const char* summary1File = "./data/exon.summary1.txt";
  const char* summary2File = "./data/exon.summary2.txt";
  const char* summary3File = "./data/exon.summary3.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName, 
    "-o", outputFile, 
    "-g", genomePosFile,
    "-probeset-ids", probesetIdFile1,
    "-probeset-ids", probesetIdFile2,
    "-transcript-cluster-ids", transcriptClusterIdFile1,
    "-transcript-cluster-ids", transcriptClusterIdFile2,
    "-prepend-filename",
    summaryFile, summary0File, summary1File,
    summary2File, summary3File,
    NULL};
  int argcTest = ( sizeof (argvTest) / sizeof (argvTest[0]) ) - 1;

  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  bool summaryOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    summaryVis summary (argcTest, argvTest, version);
    summary.run();
    summaryOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (summaryOk == true);

  // check output
  TextFileCheck check ("output/summaryvis.txt", "data/summaryvis.txt", 61);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void SummaryVisTest::testWiggleColName()
{
  const char* programName = "apt-summary-vis";
  const char* outputFile = "./output/summaryvis-col-name.wig";
  const char* genomePosFile = "./data/HuEx-1_0-st-v2.annot.hg16.subset.csv";
  const char* thyroidColName = "huex_wta_thyroid_A1.CEL";
  const char* summaryFile = "./data/exon.summary.txt";
  const char* summary0File = "./data/exon.summary0.txt";
  const char* summary1File = "./data/exon.summary1.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName,
    "-o", outputFile, 
    "-g", genomePosFile,
    "-wiggle-col-name", thyroidColName, 
    summaryFile, summary0File, summary1File, 
    NULL};
  int argcTest = ( sizeof (argvTest) / sizeof (argvTest[0]) ) - 1;
  
  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  bool summaryOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    summaryVis summary (argcTest, argvTest, version);
    summary.run();
    summaryOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (summaryOk == true);

  // check output
  TextFileCheck check ("output/summaryvis-col-name.wig", "data/summaryvis-col-name.wig", 41);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void SummaryVisTest::testWiggleColIndex()
{
  const char* programName = "apt-summary-vis";
  const char* outputFile = "./output/summaryvis-col-index.wig";
  const char* genomePosFile = "./data/HuEx-1_0-st-v2.annot.hg16.subset.csv";
  const char* summaryFile = "./data/exon.summary.txt";
  const char* summary0File = "./data/exon.summary0.txt";
  const char* summary1File = "./data/exon.summary1.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName,
    "-o", outputFile,
    "-g", genomePosFile,
    "-wiggle-col-index", "36",
    summaryFile, summary0File, summary1File, 
    NULL};
  int argcTest = ( sizeof (argvTest) / sizeof (argvTest[0]) ) - 1;

  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  bool summaryOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    summaryVis summary (argcTest, argvTest, version);
    summary.run();
    summaryOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (summaryOk == true);

  // check output
  TextFileCheck check ("output/summaryvis-col-index.wig", "data/summaryvis-col-index.wig", 41);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}
