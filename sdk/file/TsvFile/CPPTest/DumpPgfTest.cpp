////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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


/// @file   DumpPgfTest.cpp
/// @brief  Cppunit class for testing DumpPgf.

#include "file/TsvFile/CPPTest/DumpPgfTest.h"
#include "file/TsvFile/DumpPgf.h"
#include "util/Fs.h"
#include "util/TextFileCheck.h"
#include "util/Util.h"
//

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DumpPgfTest );

using namespace std;

void DumpPgfTest::testDumpEntirePgf()
{
  const char* programName = "apt-dump-pgf";
  const char* outputFile = "./output/dump-entire-pgf.txt";
  const char* pgfFile = "./data/HuEx-1_0-st-v1.ed.pgf";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName, 
    "-out-file", outputFile, 
    "-p", pgfFile,
    NULL};

  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  bool dumpOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    dumpPgf pgfDump (argvTest, version);
    pgfDump.run();
    dumpOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (dumpOk == true);

  // check output
  TextFileCheck check ("output/dump-entire-pgf.txt", "data/dump-entire-pgf.txt", 4);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void DumpPgfTest::testDumpProbesetIds()
{
  const char* programName = "apt-dump-pgf";
  const char* outputFile = "./output/dump-probeset-ids.txt";
  const char* pgfFile = "./data/HuEx-1_0-st-v1.ed.pgf";
  const char* probesetIdFileName = "./data/probeset-ids.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName,
    "-out-file", outputFile, 
    "-p", pgfFile,
    "-probeset-ids", probesetIdFileName, 
    NULL};

  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  bool dumpOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    dumpPgf pgfDump (argvTest, version);
    pgfDump.run();
    dumpOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (dumpOk == true);

  // check output
  TextFileCheck check ("output/dump-probeset-ids.txt", "data/dump-probeset-ids.txt", 4);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void DumpPgfTest::testDumpProbeIds()
{
  const char* programName = "apt-dump-pgf";
  const char* outputFile = "./output/dump-probe-ids.txt";
  const char* pgfFile = "./data/HuEx-1_0-st-v1.ed.pgf";
  const char* clfFile = "./data/HuEx-1_0-st-v1.ed.clf";
  const char* probeIdFileName = "./data/probe-ids.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName, 
    "-out-file", outputFile, 
    "-p", pgfFile,
    "-c", clfFile,
    "-probe-ids", probeIdFileName, 
    NULL};

  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  bool dumpOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    dumpPgf pgfDump (argvTest, version);
    pgfDump.run();
    dumpOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (dumpOk == true);

  // check output
  TextFileCheck check ("output/dump-probe-ids.txt", "data/dump-probe-ids.txt", 4);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void DumpPgfTest::testDumpSequentialClf()
{
  const char* programName = "apt-dump-pgf";
  const char* outputFile = "./output/dump-sequential-clf.txt";
  const char* pgfFile = "./data/HuEx-1_0-st-v1.ed.pgf";
  const char* clfFile = "./data/HuEx-1_0-st-v1.sequential.clf";
  const char* probeIdFileName = "./data/probe-ids.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName, 
    "-out-file", outputFile,
    "-p", pgfFile,
    "-c", clfFile,
    "-probe-ids", probeIdFileName, 
    NULL};

  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  bool dumpOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    dumpPgf pgfDump (argvTest, version);
    pgfDump.run();
    dumpOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (dumpOk == true);

  // check output
  TextFileCheck check ("output/dump-sequential-clf.txt", "data/dump-sequential-clf.txt", 4);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void DumpPgfTest::testDumpType()
{
  const char* programName = "apt-dump-pgf";
  const char* outputFile = "./output/dump-type.txt";
  const char* pgfFile = "./data/HuEx-1_0-st-v1.ed.pgf";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName,
    "-out-file", outputFile, 
    "-p", pgfFile,
    "-probeset-type", "redo-coverage",
    "-probeset-type", "HG-U133-2.0",
    NULL};

  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  bool dumpOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    dumpPgf pgfDump (argvTest, version);
    pgfDump.run();
    dumpOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (dumpOk == true);

  // check output
  TextFileCheck check ("output/dump-type.txt", "data/dump-type.txt", 4);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void DumpPgfTest::testDumpTypeOr()
{
  const char* programName = "apt-dump-pgf";
  const char* outputFile = "./output/dump-type-or.txt";
  const char* pgfFile = "./data/HuEx-1_0-st-v1.ed.pgf";

  // set up dummy argc, argv, run as if called from main()
  const char* argvTest[] = {
    programName, 
    "-out-file", outputFile, 
    "-p", pgfFile,
    "-or",
    "-probeset-type", "redo-coverage", 
    "-probeset-type", "HG-U133-2.0",
    NULL};

  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output");
  }

  bool dumpOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    dumpPgf pgfDump (argvTest, version);
    pgfDump.run();
    dumpOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (dumpOk == true);

  // check output
  TextFileCheck check ("output/dump-type-or.txt", "data/dump-type-or.txt", 4);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}
