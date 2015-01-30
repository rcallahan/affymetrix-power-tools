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

/// @file   MidasCelFilesMissingTest.cpp
/// @brief  Test meta probeset file with probeset_id mismatch

//
#include "midas/CPPTest/MidasCelFilesMissingTest.h"
//
#include "midas/MidasConfigureRun.h"
#include "midas/MidasCreateDirectory.h"
#include "midas/MidasEngine.h"
#include "util/Fs.h"
#include "util/CPPTest/Setup.h"
//
#include <climits>
#include <iostream>


#define PROGRAM_NAME "MidasCelFilesMissingTest"
#define EXTRA_CELS_FILE "./data/ExtraCels.txt"
#define META_FILE "./data/Meta.txt"
#define GENE_DATA_FILE "./data/GeneData.txt"
#define ONE_GROUP_EXON_DATA_FILE "./data/OneGroupExonData.txt"
#define OUTPUT_DIR "output"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MidasCelFilesMissingTest );

using namespace std;

void MidasCelFilesMissingTest::testCelFilesMissing()
{
  PgOptions opts;
  define_midas_opts(&opts);

  // delete output from previous run, if any
  string outDirectory (OUTPUT_DIR);
  string pvaluesFileName = Fs::join(outDirectory,PVALUES_OUTPUT);
  Fs::rmIfExists(pvaluesFileName);

  // set up dummy argc, argv, run as if called from main()
  const char* argvCelFilesMissing[] = {
    PROGRAM_NAME,
    "--cel-files", EXTRA_CELS_FILE, 
    "-g", GENE_DATA_FILE,
    "-e", ONE_GROUP_EXON_DATA_FILE,
    "-m", META_FILE,
    "-o", OUTPUT_DIR,
    NULL};

  opts.parseArgv(argvCelFilesMissing);
  // create output directory if not already present
  std::string msg1 = midasCreateDirectory(opts.get("out-dir"));
  CPPUNIT_ASSERT (msg1 == "");

  // create object to configure, run midas
  const string version ("NON-OFFICIAL-RELEASE");
  const string cvsId ("$Id: MidasCelFilesMissingTest.cpp,v 1.16 2009-09-18 03:37:29 mspald Exp $");
  const string execVersion = version + " " + cvsId;
  midasConfigureRun configureRun (opts.get("cel-files"), opts.get("genedata"),
                                  opts.get("exondata"), opts.get("metaprobeset"),  opts.get("out-dir"),
                                  opts.getBool("pvalues"), opts.getBool("fstats"), opts.getBool("normalized"),
                                  opts.getDouble("stabilize"), opts.commandLine(), execVersion);
  NEGATIVE_TEST(configureRun.configure(), Except);
}
