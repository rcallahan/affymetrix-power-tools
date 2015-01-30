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

/// @file   MidasBadGeneDataTest.cpp
/// @brief  Test meta file with probeset_id mismatch

//
#include "midas/CPPTest/MidasBadGeneDataTest.h"
//
#include "midas/MidasConfigureRun.h"
#include "midas/MidasCreateDirectory.h"
#include "midas/MidasEngine.h"
#include "util/Fs.h"
#include "util/CPPTest/Setup.h"
//
#include <climits>
#include <iostream>

#define PROGRAM_NAME "MidasBadGeneDataTest"
#define CELS_FILE "./data/Cels.txt"
#define META_FILE "./data/Meta.txt"
#define BAD_GENE_DATA_FILE "./data/BadGeneData.txt"
#define EXON_DATA_FILE "./data/ExonData.txt"
#define OUTPUT_DIR "output"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MidasBadGeneDataTest );

using namespace std;

void MidasBadGeneDataTest::testBadGeneData()
{
  PgOptions opts;
  define_midas_opts(&opts);

  // delete output from previous run, if any
  string outDirectory (OUTPUT_DIR);
  string pvaluesFileName = Fs::join(outDirectory,PVALUES_OUTPUT);
  Fs::rmIfExists(pvaluesFileName);
  
  // set up dummy argc, argv, run as if called from main()
  const char* argvBadGeneData[] = {
    PROGRAM_NAME,
    "--cel-files", CELS_FILE,
    "-g", BAD_GENE_DATA_FILE,
    "-e", EXON_DATA_FILE,
    "-m", META_FILE,
    "-o", OUTPUT_DIR, 
    NULL};

  opts.parseArgv(argvBadGeneData);

  // create output directory if not already present
  std::string msg1=midasCreateDirectory (opts.get("out-dir"));
  CPPUNIT_ASSERT (msg1 == "");

  // create object to configure, run midas
  const string version ("NON-OFFICIAL-RELEASE");
  const string cvsId ("$Id: MidasBadGeneDataTest.cpp,v 1.16 2009-09-18 03:37:29 mspald Exp $");
  const string execVersion = version + " " + cvsId;
  midasConfigureRun configureRun (opts.get("cel-files"), opts.get("genedata"),
                                  opts.get("exondata"), opts.get("metaprobeset"),  opts.get("out-dir"),
                                  opts.getBool("pvalues"), opts.getBool("fstats"), opts.getBool("normalized"),
                                  opts.getDouble("stabilize"), opts.commandLine(), execVersion);
  string* msg2 = configureRun.configure();
  CPPUNIT_ASSERT (msg2 == 0);
  // run the midas engine
  NEGATIVE_TEST(configureRun.run(), Except);
}
