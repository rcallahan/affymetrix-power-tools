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

/// @file   CelExtractTest.cpp
/// @brief  Cppunit class for testing CelExtract.

#include "chipstream/apt-cel-extract/CPPTest/CelExtractTest.h"
#include "chipstream/apt-cel-extract/CelExtract.h"
#include "util/Fs.h"
#include "util/MixedFileCheck.h"
#include "util/Util.h"
//

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( CelExtractTest );

using namespace std;

void CelExtractTest::testExtractEntirePgf()
{
  // create the output dir
  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  // throw exceptions rather than exit
  Err::setThrowStatus(true);
  bool extractOk = false;

  // try extraction
  try
  {
      CelExtractOptions o;
      o.pgfFile = "./data/HuEx-test.pgf";
      o.clfFile = "./data/HuEx-test.clf";
      o.celFiles.push_back("./data/huex_cerebellum.CEL");
      o.celFiles.push_back("./data/huex_heart.CEL");
      o.celFiles.push_back("./data/huex_muscle.CEL");
      o.outFile = "./output/extract-entire-pgf.txt";
      CelExtract ce(o);
      ce.extract();
      extractOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (extractOk == true);

  // back to regular exit err handler
  Err::setThrowStatus(false);

  // check output
  MixedFileCheck check ("output/extract-entire-pgf.txt", "data/extract-entire-pgf.txt", 0.1,0,0);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void CelExtractTest::testExtractProbesetIds()
{
  // create the output dir
  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  // throw exceptions rather than exit
  Err::setThrowStatus(true);
  bool extractOk = false;

  // try extraction
  try
  {
      CelExtractOptions o;
      o.pgfFile = "./data/HuEx-test.pgf";
      o.clfFile = "./data/HuEx-test.clf";
      o.celFiles.push_back("./data/huex_cerebellum.CEL");
      o.celFiles.push_back("./data/huex_heart.CEL");
      o.celFiles.push_back("./data/huex_muscle.CEL");
      o.outFile = "./output/extract-probeset-ids.txt";
      o.probesetIdsFiles.push_back("./data/probeset-ids.txt");
      CelExtract ce(o);
      ce.extract();
      extractOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (extractOk == true);

  // back to regular exit err handler
  Err::setThrowStatus(false);

  // check output
  MixedFileCheck check ("output/extract-probeset-ids.txt", "data/extract-probeset-ids.txt", 0.1,0,0);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void CelExtractTest::testExtractProbeIds()
{
  // create the output dir
  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  // throw exceptions rather than exit
  Err::setThrowStatus(true);
  bool extractOk = false;

  // try extraction
  try
  {
      CelExtractOptions o;
      o.pgfFile = "./data/HuEx-test.pgf";
      o.clfFile = "./data/HuEx-test.clf";
      o.celFiles.push_back("./data/huex_cerebellum.CEL");
      o.celFiles.push_back("./data/huex_heart.CEL");
      o.celFiles.push_back("./data/huex_muscle.CEL");
      o.outFile = "./output/extract-probe-ids.txt";
      o.probeIdsFiles.push_back("./data/probe-ids.txt");
      CelExtract ce(o);
      ce.extract();
      extractOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (extractOk == true);

  // back to regular exit err handler
  Err::setThrowStatus(false);

  // check output
  MixedFileCheck check ("output/extract-probe-ids.txt", "data/extract-probe-ids.txt", 0.1,0,0);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void CelExtractTest::testExtractPmGcbg()
{
  // create the output dir
  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  // throw exceptions rather than exit
  Err::setThrowStatus(true);
  bool extractOk = false;

  // try extraction
  try
  {
      CelExtractOptions o;
      o.pgfFile = "./data/HuEx-test.pgf";
      o.clfFile = "./data/HuEx-test.clf";
      o.celFiles.push_back("./data/huex_cerebellum.CEL");
      o.celFiles.push_back("./data/huex_heart.CEL");
      o.celFiles.push_back("./data/huex_muscle.CEL");
      o.outFile = "./output/extract-pm-gcbg.txt";
      o.bgpFile = "./data/antigenomic.bgp";
      o.probesetIdsFiles.push_back("./data/probeset-ids.txt");
      o.analysisString = "pm-gcbg";
      o.subtractBackground = true;
      o.pmOnly = true;
      CelExtract ce(o);
      ce.extract();
      extractOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (extractOk == true);

  // back to regular exit err handler
  Err::setThrowStatus(false);

  // check output
  MixedFileCheck check ("output/extract-pm-gcbg.txt", "data/extract-pm-gcbg.txt", 1.0,0,0);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void CelExtractTest::testExtractPmMm()
{
  // create the output dir
  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  // throw exceptions rather than exit
  Err::setThrowStatus(true);
  bool extractOk = false;

  // try extraction
  try
  {
      CelExtractOptions o;
      o.pgfFile = "./data/HuEx-test.pgf";
      o.clfFile = "./data/HuEx-test.clf";
      o.celFiles.push_back("./data/huex_cerebellum.CEL");
      o.celFiles.push_back("./data/huex_heart.CEL");
      o.celFiles.push_back("./data/huex_muscle.CEL");
      o.outFile = "./output/extract-pm-mm.txt";
      o.probeIdsFiles.push_back( "./data/probe-ids.txt");
      o.analysisString = "pm-mm";
      o.subtractBackground = true;
      o.ignoreProbesWithoutMM = true;

      CelExtract ce(o);
      ce.extract();
      extractOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (extractOk == true);

  // back to regular exit err handler
  Err::setThrowStatus(false);

  // check output
  MixedFileCheck check ("output/extract-pm-mm.txt", "data/extract-pm-mm.txt", 0.1,0,0);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}

void CelExtractTest::testExtractCdf()
{
  // create the output dir
  if ( !Fs::dirExists("output") ) {
    Fs::mkdir("output", false);
  }

  // throw exceptions rather than exit
  Err::setThrowStatus(true);
  bool extractOk = false;

  // try extraction
  try
  {
      CelExtractOptions o;
      o.cdfFile = "../../../rawq/test/data/Test3.CDF";
      o.celFiles.push_back("../../../rawq/test/data/Test3.CEL");
      o.outFile = "./output/extract-cdf.txt";
      o.probesetIdsFiles.push_back("./data/probeset-names.txt");
      o.force = true;
      o.analysisString = "pm-mm";
      o.subtractBackground = true;
      o.ignoreProbesWithoutMM = true;
      CelExtract ce(o);
      ce.extract();
      extractOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    cout << e.what() << endl;
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (extractOk == true);

  // back to regular exit err handler
  Err::setThrowStatus(false);

  // check output
  MixedFileCheck check ("output/extract-cdf.txt", "data/extract-cdf.txt", 1.0,0,0);
  string errorMsg;
  bool success = check.check (errorMsg);
  if (! success)
    cout << errorMsg << endl;
  CPPUNIT_ASSERT( success == true);
}
