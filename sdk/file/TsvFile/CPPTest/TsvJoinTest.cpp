////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////


/// @file   TsvJoinTest.cpp
/// @brief  Cppunit class for testing TsvJoin.

#include "file/TsvFile/CPPTest/TsvJoinTest.h"
#include "file/TsvFile/TsvJoin.h"
#include "util/Fs.h"
#include "util/TextFileCheck.h"
#include "util/Util.h"
//

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( TsvJoinTest );

using namespace std;

void TsvJoinTest::testTsvJoin()
{
  const char* programName = "apt-tsv-join";
  const char* outputFile = "./output/tsvjoin.txt";
  const char* probesetFile = "./data/probeset-annot.csv";
  const char* transcriptFile = "./data/transcript-annot.csv";
  const char* transcript1File = "./data/transcript-annot1.csv";
  const char* transcript2File = "./data/transcript-annot2.csv";
  const char* transcript3File = "./data/transcript-annot3.csv";
  const char* transcript4File = "./data/transcript-annot4.csv";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName, 
    "-k",
    "transcript_cluster_id", 
    "-o", 
    outputFile,
    probesetFile,
    transcriptFile,
    transcript1File,
    transcript2File,
    transcript3File,
    transcript4File, 
    NULL};

  int argcTest = ( sizeof (argvTest) / sizeof (argvTest[0]) ) - 1;

  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  bool joinOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    tsvJoin join (argcTest, argvTest, version);
    join.run();
    joinOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (joinOk == true);

  // check output
  TextFileCheck check ("output/tsvjoin.txt", "data/tsvjoin.txt", 4);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}
