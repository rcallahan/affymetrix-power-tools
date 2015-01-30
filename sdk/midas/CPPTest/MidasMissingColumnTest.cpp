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

/// @file   MidasMissingColumnTest.cpp
/// @brief  Test warning message generation from gene data missing a cel column

//
#include "midas/CPPTest/MidasMissingColumnTest.h"
//
#include "midas/MidasConfigureRun.h"
#include "midas/MidasCreateDirectory.h"
#include "midas/MidasEngine.h"
#include "util/Fs.h"
//
#include <climits>
#include <iostream>

#define PROGRAM_NAME "MidasMissingColumnTest"
#define CELS_FILE "./data/Cels.txt"
#define META_FILE "./data/Meta.txt"
#define MISSING_COLUMN_FILE "./data/GeneDataMissingColumn.txt"
#define EXON_DATA_FILE "./data/ExonData.txt"
#define OUTPUT_DIR "output"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MidasMissingColumnTest );

using namespace std;

void MidasMissingColumnTest::testMissingColumn()
{
  PgOptions opts;
  define_midas_opts(&opts);

  // delete output from previous run, if any
  string outDirectory (OUTPUT_DIR);
  string pvaluesFileName = Fs::join(outDirectory,PVALUES_OUTPUT);
  Fs::rmIfExists(pvaluesFileName);

  // set up dummy argc, argv, run as if called from main()
  const char* argvMissingColumn[] = {
    PROGRAM_NAME,
    "--cel-files", CELS_FILE, 
    "-g", MISSING_COLUMN_FILE,
    "-e", EXON_DATA_FILE, 
    "-m", META_FILE,
    "-o", OUTPUT_DIR, 
    NULL};
  opts.parseArgv(argvMissingColumn);
  // create output directory if not already present
  std::string msg = midasCreateDirectory (opts.get("out-dir"));
  CPPUNIT_ASSERT (msg == "");

  bool configureRunOk = false;
  bool caughtException = false;
  try
  {
    // create object to configure, run midas
    const string version ("NON-OFFICIAL-RELEASE");
    const string cvsId ("$Id: MidasMissingColumnTest.cpp,v 1.15 2009-09-18 03:37:29 mspald Exp $");
    const string execVersion = version + " " + cvsId;
    midasConfigureRun configureRun (opts.get("cel-files"), opts.get("genedata"),
                                    opts.get("exondata"), opts.get("metaprobeset"),  opts.get("out-dir"),
                                    opts.getBool("pvalues"), opts.getBool("fstats"), opts.getBool("normalized"),
                                    opts.getDouble("stabilize"), opts.commandLine(), execVersion);
    // the point of this test is to generate a warning message in configure()
    string* msg = configureRun.configure();
    CPPUNIT_ASSERT (msg != 0);
    // back out the pvalues output file
    configureRun.deleteOutputs();
    // check that the pvalues output file is no longer present
    ifstream inFstream;
    inFstream.open (pvaluesFileName.c_str(), ios_base::in);
    bool pvaluesFileNameNotFound = inFstream.fail();
    CPPUNIT_ASSERT (pvaluesFileNameNotFound == true);
    // make sure no exception has been thrown
    configureRunOk = true;
  }
  catch (exception& e)
  {
    // should not execute the following line
    caughtException = true;
  }
  CPPUNIT_ASSERT (configureRunOk == true);
  CPPUNIT_ASSERT (caughtException == false);
}
