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

/// @file   MidasNormalTest.cpp
/// @brief  Cppunit class for testing midas methods.

//
#include "midas/CPPTest/MidasNormalTest.h"
//
#include "midas/MidasConfigureRun.h"
#include "midas/MidasCreateDirectory.h"
#include "midas/MidasEngine.h"
#include "midas/MidasSpliceDetector.h"
//
#include "portability/affy-base-types.h"
#include "util/Fs.h"
//
#include <climits>
#include <cmath>
#include <iomanip>
#include <iostream>

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MidasNormalTest );

using namespace std;

bool MidasNormalTest::floatsCloseEnough (const float& f1, const float& f2, int digits)
{
  const float tolerance = 1 / pow (10.0, digits);
  if (fabs (f1 - f2) > tolerance)
    return false;
  return true;
}

bool MidasNormalTest::vectorsCloseEnough (const vector<float>& vf1, const vector<float>& vf2, int digits)
{
  const unsigned int size = vf1.size();
  if (vf2.size() != size)
    return false;
  const float tolerance = 1 / pow (10.0, digits);
  for (unsigned int i = 0; i < size; ++i)
    if (fabs (vf1[i] - vf2[i]) > tolerance)
      return false;
  return true;
}

bool MidasNormalTest::vectorsOfVectorsCloseEnough (const vector<vector<float> >& vvf1,
      const vector<vector<float> >& vvf2, int digits)
{
  const unsigned int size = vvf1.size();
  if (vvf2.size() != size)
    return false;
  for (unsigned int i = 0; i < size; ++i)
    if (!vectorsCloseEnough (vvf1[i], vvf2[i], digits))
	return false;
  return true;
}

void MidasNormalTest::resetOptionValues (PgOpt *options[])
{
  for (int optIx = 0; (options != NULL) && (options[optIx] != NULL); optIx++)
    options[optIx]->resetToDefault();
}

void MidasNormalTest::testSpliceDetector()
{
  // experimental groups
  const int groupsIn[] = { 1, 1, 1, 2, 2, 3, 3, 3 };
  const int groupsSize = sizeof (groupsIn) / sizeof (groupsIn[0]);
  vector<int> groups (groupsIn, groupsIn + groupsSize);

  // test 1 gene data
  const float gData1In[] = { 22.3, 24.5, 27.7, 29.1, 34.6, 26.6, 100.9, 95.8 };
  const int dataSize = sizeof (gData1In) / sizeof (gData1In[0]);
  const vector<float> gData1 (gData1In, gData1In + dataSize);

  // test 1 exon data
  vector<vector<float> > pData1;
  const float pData11In[] = { 1, 3, 4.5, 10.1, 15.2, 32.2, 47.5, 43.3 };
  const float pData12In[] = { 11.7, 12.5, 14.6, 14.6, 20.3, 40.1, 50.2, 48.8 };
  const float pData13In[] = { 6.6, 8.3, 7, 7.9, 16.9, 28.8, 25.5, 30.3 };
  const float pData14In[] = { 44.1, 56.1, 43.1, 66.3, 77.8, 155.6, 140.6, 183.8 };
  pData1.push_back (vector<float> (pData11In, pData11In + dataSize));
  pData1.push_back (vector<float> (pData12In, pData12In + dataSize));
  pData1.push_back (vector<float> (pData13In, pData13In + dataSize));
  pData1.push_back (vector<float> (pData14In, pData14In + dataSize));

  // test 1 expected fstat results

  const float fstatExpect1In[] = { 3.879666328, 0.177621007, 0.01964771189, 0.3792025149 };
  const int test1Count = sizeof (fstatExpect1In) / sizeof (fstatExpect1In[0]);
  const vector<float> fstatExpect1 (fstatExpect1In, fstatExpect1In + test1Count);

  // test 1 expected spliceIndex results
  vector<vector<float> > spliceIndexExpect1;
  const float spliceIndexExpect11In[] = { -0.7059602737, -0.7059602737, -0.7059602737, -0.253100425, -0.253100425, 0, 0, 0 };
  const float spliceIndexExpect12In[] = { -0.1494996548, -0.1494996548, -0.1494996548, -0.15231435, -0.15231435, 0, 0, 0 };
  const float spliceIndexExpect13In[] = { -0.07029172033, -0.07029172033, -0.07029172033, 0, 0, -0.01261007786,
    -0.01261007786, -0.01261007786 };
  const float spliceIndexExpect14In[] = { -0.29951033, -0.29951033, -0.29951033, -0.1287990808, -0.1287990808, 0, 0, 0 };
  spliceIndexExpect1.push_back (vector<float> (spliceIndexExpect11In, spliceIndexExpect11In + dataSize));
  spliceIndexExpect1.push_back (vector<float> (spliceIndexExpect12In, spliceIndexExpect12In + dataSize));
  spliceIndexExpect1.push_back (vector<float> (spliceIndexExpect13In, spliceIndexExpect13In + dataSize));
  spliceIndexExpect1.push_back (vector<float> (spliceIndexExpect14In, spliceIndexExpect14In + dataSize));

  // test 1 expected pvalue results
  const float pvalueExpect1In[] = { 0.1472325623, 0.8454658985, 0.9806691408, 0.7131428123 };
  const vector<float> pvalueExpect1 (pvalueExpect1In, pvalueExpect1In + test1Count);

  // test 2 gene data
  const float gData2In[] = { 26.9, 32.6, 38.6, 100.5, 160.4, 210.6, 240.6, 193.4 };
  const vector<float> gData2 (gData2In, gData2In + dataSize);

  // test 2 exon data
  vector<vector<float> > pData2;
  const float pData21In[] = { 11.1, 13.1, 14.1, 55.6, 50.6, 100.3, 110.4, 93.6 };
  const float pData22In[] = { 28.7, 35.8, 36.2, 110.5, 145.8, 339.9, 400.8, 267.8 };
  pData2.push_back (vector<float> (pData21In, pData21In + dataSize));
  pData2.push_back (vector<float> (pData22In, pData22In + dataSize));

  // test 2 expected fstat results
  const float fstatExpect2In[] = { 0.3296247125, 17.93391418 };
  const int test2Count = sizeof (fstatExpect2In) / sizeof (fstatExpect2In[0]);
  const vector<float> fstatExpect2 (fstatExpect2In, fstatExpect2In + test2Count);

  // test 2 expected spliceIndex results
  vector<vector<float> > spliceIndexExpect2;
  const float spliceIndexExpect21In[] = { 0, 0, 0, -0.1270999908, -0.1270999908, -0.04168279842, -0.04168279842, -0.04168279842 };
  const float spliceIndexExpect22In[] = { -0.4010516107, -0.4010516107, -0.4010516107, -0.426741451, -0.426741451, 0, 0, 0 };
  spliceIndexExpect2.push_back (vector<float> (spliceIndexExpect21In, spliceIndexExpect21In + dataSize));
  spliceIndexExpect2.push_back (vector<float> (spliceIndexExpect22In, spliceIndexExpect22In + dataSize));

  // test 2 expected pvalue results
  const float pvalueExpect2In[] = { 0.7423245907, 0.02144354023 };
  const vector<float> pvalueExpect2 (pvalueExpect2In, pvalueExpect2In + test2Count);

  // test 4 gene, exon data identical: gene, exon data
  const float gData4In[] = { 22.3, 24.5, 27.7, 29.1, 34.6, 26.6, 100.9, 95.8 };
  const vector<float> gData4 (gData4In, gData4In + dataSize);

  // test 4 gene, exon data identical: expected fstat, pvalue results
  const float fstatExpect4 = 0.0;
  const float pvalueExpect4 = 1.0;

  // test 4 gene, exon data identical: expected spliceIndex results
  const float spliceIndexExpect4In[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
  const vector<float> spliceIndexExpect4 (spliceIndexExpect4In, spliceIndexExpect4In + dataSize);

  // request all outputs
  bool wantPvalues = true;
  bool wantFstats = true;
  bool wantNormalized = true;

  const float logStabilize = 8.0;
  bool noLogTransform = false;
  midasSpliceDetector m (groups, wantPvalues, wantFstats, wantNormalized, logStabilize, noLogTransform);

  // test single exon call
  vector<float> fstatFound1 (test1Count);
  vector<vector<float> > spliceIndexFound1 (test1Count, vector<float> (dataSize));
  vector<float> pvalueFound1 (test1Count);
  for (int i = 0; i < test1Count; ++i)
    m.runMidasSingle (pData1[i], gData1, &fstatFound1[i], &spliceIndexFound1[i], &pvalueFound1[i]);

  vector<float> fstatFound2 (test2Count);
  vector<vector<float> > spliceIndexFound2 (test2Count, vector<float> (dataSize));
  vector<float> pvalueFound2 (test2Count);
  for (int i = 0; i < test2Count; ++i)
    m.runMidasSingle (pData2[i], gData2, &fstatFound2[i], &spliceIndexFound2[i], &pvalueFound2[i]);

  // test multiple exons call
  vector<float> fstatFound1a (test1Count);
  vector<vector<float> > spliceIndexFound1a (test1Count, vector<float> (dataSize));
  vector<float> pvalueFound1a (test1Count);
  m.runMidasMultiple (pData1, gData1, &fstatFound1a, &spliceIndexFound1a, &pvalueFound1a);

  // test multiple exons call, no fstat requested
  wantPvalues = true;
  wantFstats = false;
  wantNormalized = true;
  midasSpliceDetector m1 (groups, wantPvalues, wantFstats, wantNormalized, logStabilize, noLogTransform);
  vector<vector<float> > spliceIndexFound1b (test1Count, vector<float> (dataSize));
  vector<float> pvalueFound1b (test1Count);
  m1.runMidasMultiple (pData1, gData1, 0, &spliceIndexFound1b, &pvalueFound1b);

  // test multiple exons call, no spliceIndex requested
  wantPvalues = true;
  wantFstats = true;
  wantNormalized = false;
  midasSpliceDetector m2 (groups, wantPvalues, wantFstats, wantNormalized, logStabilize, noLogTransform);
  vector<float> fstatFound1c (test1Count);
  vector<float> pvalueFound1c (test1Count);
  m2.runMidasMultiple (pData1, gData1, &fstatFound1c, 0, &pvalueFound1c);

  // test multiple exons call, no pvalue requested
  wantPvalues = false;
  wantFstats = true;
  wantNormalized = true;
  midasSpliceDetector m3 (groups, wantPvalues, wantFstats, wantNormalized, logStabilize, noLogTransform);
  vector<float> fstatFound1d (test1Count);
  vector<vector<float> > spliceIndexFound1d (test1Count, vector<float> (dataSize));
  m3.runMidasMultiple (pData1, gData1, &fstatFound1d, &spliceIndexFound1d, 0);

  // test gene, exon data identical
  wantPvalues = true;
  wantFstats = true;
  wantNormalized = true;
  midasSpliceDetector m4 (groups, wantPvalues, wantFstats, wantNormalized, logStabilize, noLogTransform);
  float fstatFound4;
  float pvalueFound4;
  vector<float> spliceIndexFound4 (dataSize);
  m4.runMidasSingle (gData4, gData4, &fstatFound4, &spliceIndexFound4, &pvalueFound4);

  // the F statistic can be quite large; compare log values rather than
  // using an involved normalization scheme
  vector<float> logFstatExpect1 (test1Count);
  for (int i = 0; i < test1Count; ++i)
    logFstatExpect1 [i] = log (fstatExpect1 [i]);
  vector<float> logFstatFound1 (test1Count);
  for (int i = 0; i < test1Count; ++i)
    logFstatFound1 [i] = log (fstatFound1 [i]);
  CPPUNIT_ASSERT (vectorsCloseEnough (logFstatFound1, logFstatExpect1, 5));
  CPPUNIT_ASSERT (vectorsOfVectorsCloseEnough (spliceIndexFound1, spliceIndexExpect1));
  CPPUNIT_ASSERT (vectorsCloseEnough (pvalueFound1, pvalueExpect1));

  vector<float> logFstatExpect2 (test2Count);
  for (int i = 0; i < test2Count; ++i)
    logFstatExpect2 [i] = log (fstatExpect2 [i]);
  vector<float> logFstatFound2 (test2Count);
  for (int i = 0; i < test2Count; ++i)
    logFstatFound2 [i] = log (fstatFound2 [i]);
  CPPUNIT_ASSERT (vectorsCloseEnough (logFstatFound2, logFstatExpect2, 5));
  CPPUNIT_ASSERT (vectorsOfVectorsCloseEnough (spliceIndexFound2, spliceIndexExpect2));
  CPPUNIT_ASSERT (vectorsCloseEnough (pvalueFound2, pvalueExpect2));

  vector<float> logFstatFound1a (test1Count);
  for (int i = 0; i < test1Count; ++i)
    logFstatFound1a [i] = log (fstatFound1a [i]);
  CPPUNIT_ASSERT (vectorsCloseEnough (logFstatFound1a, logFstatExpect1, 5));
  CPPUNIT_ASSERT (vectorsOfVectorsCloseEnough (spliceIndexFound1a, spliceIndexExpect1));
  CPPUNIT_ASSERT (vectorsCloseEnough (pvalueFound1a, pvalueExpect1));

  CPPUNIT_ASSERT (vectorsOfVectorsCloseEnough (spliceIndexFound1b, spliceIndexExpect1));
  CPPUNIT_ASSERT (vectorsCloseEnough (pvalueFound1b, pvalueExpect1));

  vector<float> logFstatFound1c (test1Count);
  for (int i = 0; i < test1Count; ++i)
    logFstatFound1c [i] = log (fstatFound1c [i]);
  CPPUNIT_ASSERT (vectorsCloseEnough (logFstatFound1c, logFstatExpect1, 5));
  CPPUNIT_ASSERT (vectorsCloseEnough (pvalueFound1c, pvalueExpect1));

  vector<float> logFstatFound1d (test1Count);
  for (int i = 0; i < test1Count; ++i)
    logFstatFound1d [i] = log (fstatFound1d [i]);
  CPPUNIT_ASSERT (vectorsCloseEnough (logFstatFound1d, logFstatExpect1, 5));
  CPPUNIT_ASSERT (vectorsOfVectorsCloseEnough (spliceIndexFound1d, spliceIndexExpect1));

  // checks for test of gene, exon data identical
  CPPUNIT_ASSERT (floatsCloseEnough (fstatFound4, fstatExpect4, 5));
  CPPUNIT_ASSERT (floatsCloseEnough (pvalueFound4, pvalueExpect4, 5));
  CPPUNIT_ASSERT (vectorsCloseEnough (spliceIndexFound4, spliceIndexExpect4, 5));
}

void MidasNormalTest::testConfigureRun()
{
  PgOptions opts;
  define_midas_opts(&opts);

  // delete output from previous run, if any
  const char* outputDir = "output";
  string outDirectory (outputDir);
  string pvaluesFileName = Fs::join(outDirectory, PVALUES_OUTPUT);
  string fstatsFileName =    Fs::join(outDirectory, FSTATS_OUTPUT);
  string normalizedFileName = Fs::join(outDirectory, NORMALIZED_OUTPUT);

  unlink (pvaluesFileName.c_str());
  unlink (fstatsFileName.c_str());
  unlink (normalizedFileName.c_str());

  const char* programName = "MidasNormalTest";
  const char* celFiles = "./data/Cels.txt";
  const char* metaFile = "./data/Meta.txt";
  const char* geneDataFile = "./data/GeneData.txt";
  const char* exonDataFile = "./data/ExonData.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvNoErrorsTest[] = {programName, 
                                    "--cel-files", celFiles, 
                                    "-g", geneDataFile,
                                    "-e", exonDataFile, 
                                    "-m", metaFile,
                                    "-o", outputDir,
                                    "-f",
                                    "-n", 
                                    NULL};
  opts.parseArgv(argvNoErrorsTest);

  std::string msg = midasCreateDirectory(opts.get("out-dir"));
  CPPUNIT_ASSERT (msg == "");

  bool configureRunOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    const string cvsId ("$Id: MidasNormalTest.cpp,v 1.25 2009-09-18 03:37:29 mspald Exp $");
    const string execVersion = version + " " + cvsId;
    // create object to configure, run midas
    midasConfigureRun configureRun (opts.get("cel-files"), opts.get("genedata"),
                                    opts.get("exondata"), opts.get("metaprobeset"),  opts.get("out-dir"),
                                    opts.getBool("pvalues"), opts.getBool("fstats"), opts.getBool("normalized"),
                                    opts.getDouble("stabilize"), opts.commandLine(), execVersion);
    // configure step may return a non-fatal warning message
    std::string* msg = configureRun.configure();
    CPPUNIT_ASSERT (msg == 0);
    // run the midas engine
    configureRun.run();
    configureRunOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (configureRunOk == true);
  // check outputs
  const string probeset_list_id_string (PROBESET_LIST_ID);
  vector<string> null_vector;
  int32_t probeset_list_id;
  int rv;

  //
  const string pvalue_string ("pvalue");
  float pvalue;
  affx::TsvFile Pvalues_tsv;
  rv=Pvalues_tsv.open(pvaluesFileName);
  CPPUNIT_ASSERT(rv==affx::TSV_OK);
  //
  rv=Pvalues_tsv.nextLevel(0);
  CPPUNIT_ASSERT(rv==affx::TSV_OK);
  rv=Pvalues_tsv.get(0,probeset_list_id_string,probeset_list_id);
  CPPUNIT_ASSERT( probeset_list_id == 2590411);
  rv=Pvalues_tsv.get(0,pvalue_string,pvalue);
  CPPUNIT_ASSERT( floatsCloseEnough (pvalue, 0.147233, 5) == true);
  //
  Pvalues_tsv.close();

  //
  const string fstatistic_string ("fstatistic");
  float fstatistic;
  affx::TsvFile Fstats_tsv;
  rv=Fstats_tsv.open(fstatsFileName);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  //
  rv=Fstats_tsv.nextLevel(0);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  rv=Fstats_tsv.get(0,probeset_list_id_string,probeset_list_id);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  CPPUNIT_ASSERT( probeset_list_id == 2590411);
  rv=Fstats_tsv.get(0,fstatistic_string,fstatistic);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  CPPUNIT_ASSERT( floatsCloseEnough (fstatistic, 3.87967, 5) == true);
  //
  Pvalues_tsv.close();


  const char* columnNamesIn[] =
    { "Sample1-1.median.ccel", "Sample2-1.median.ccel",
      "Sample3-1.median.ccel", "Sample4-2.median.ccel",
      "Sample5-2.median.ccel", "Sample6-3.median.ccel",
      "Sample7-3.median.ccel", "Sample8-3.median.ccel"
    };
  const float normalizedIn[] = {
    -0.70596, -0.70596,
    -0.70596, -0.25310,
    -0.25310,  0.00000,
    -1.66533e-16, -1.66533e-16
  };
  const int normalizedDataSize = sizeof (normalizedIn) / sizeof (normalizedIn[0]);

  //
  affx::TsvFile Norm_tsv;
  rv=Norm_tsv.open(normalizedFileName);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  rv=Norm_tsv.nextLevel(0);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  rv=Norm_tsv.get(0,probeset_list_id_string,probeset_list_id);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  CPPUNIT_ASSERT( probeset_list_id == 2590411);
  // This just checks the first line.
  for (int i=0;i<normalizedDataSize;i++) {
    // @todo: check order of columns?
    float expected;
    Norm_tsv.get(0,columnNamesIn[i],expected);
    CPPUNIT_ASSERT( floatsCloseEnough (expected,normalizedIn[i],4) == true);
  }
  Norm_tsv.close();
}

void MidasNormalTest::testKeepPath()
{
  PgOptions opts;
  define_midas_opts(&opts);

  // delete output from previous run, if any
  const char* outputDir = "output";
  string outDirectory (outputDir);
  string pvaluesFileName = Fs::join(outDirectory,PVALUES_OUTPUT);

  unlink (pvaluesFileName.c_str());

  // test cels file containing prefixes
  const char* programName = "MidasNormalTest";
  const char* celsFile = "./data/PathCels.txt";
  const char* metaFile = "./data/Meta.txt";
  const char* geneDataFile = "./data/GeneData.txt";
  const char* exonDataFile = "./data/ExonData.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvKeepPathTest[] = {programName,
                                    "--cel-files", celsFile, 
                                    "-g", geneDataFile,
                                    "-e", exonDataFile, 
                                    "-m", metaFile, 
                                    "-o", outputDir, 
                                    NULL};
  opts.parseArgv(argvKeepPathTest);
  
  std::string msg = midasCreateDirectory(opts.get("out-dir"));
  CPPUNIT_ASSERT (msg == "");

  bool configureRunOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    const string cvsId ("$Id: MidasNormalTest.cpp,v 1.25 2009-09-18 03:37:29 mspald Exp $");
    const string execVersion = version + " " + cvsId;
    // create object to configure, run midas
    midasConfigureRun configureRun (opts.get("cel-files"), opts.get("genedata"),
                                    opts.get("exondata"), opts.get("metaprobeset"),  opts.get("out-dir"),
                                    opts.getBool("pvalues"), opts.getBool("fstats"), opts.getBool("normalized"),
                                    opts.getDouble("stabilize"), opts.commandLine(), execVersion);
    // configure step may return a non-fatal warning message
    std::string* msg = configureRun.configure();
    CPPUNIT_ASSERT (msg == 0);
    // run the midas engine
    configureRun.run();
    configureRunOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (configureRunOk == true);

  
//  // check outputs
//  const string probeset_list_id_string (PROBESET_LIST_ID);
//  const string pvalue_string ("pvalue");
//  vector<string> null_vector;
//  CTSVFileData tsv;
//  int32_t probeset_list_id;
//
//  CTSVFileConvertOutputs pvaluesOutputs ("", null_vector, probeset_list_id_string, null_vector,
//                                        pvalue_string);
//  CPPUNIT_ASSERT( tsv.Open(pvaluesFileName.c_str(), CTSVFileData::TSV_STREAM_MODE,
//	                   "tab", NULL, &pvaluesOutputs) == true );
//  float pvalue;
//  CTSVFileConvertData pvaluesData (0, 0, &probeset_list_id, 0, &pvalue);
//  CPPUNIT_ASSERT( tsv.NextEntryConverted (pvaluesData) == true);
//  CPPUNIT_ASSERT( probeset_list_id == 2590411);
//  CPPUNIT_ASSERT( floatsCloseEnough (pvalue, 0.147233, 5) == true);
//  tsv.Close();
}

void MidasNormalTest::testKeepWinPath()
{
  PgOptions opts;
  define_midas_opts(&opts);

  // delete output from previous run, if any
  const char* outputDir = "output";
  string outDirectory (outputDir);
  string pvaluesFileName = Fs::join(outDirectory,PVALUES_OUTPUT);

  unlink (pvaluesFileName.c_str());

  // test cels file containing windows prefixes
  const char* programName = "MidasNormalTest";
  const char* celsFile = "./data/WinPathCels.txt";
  const char* metaFile = "./data/Meta.txt";
  const char* geneDataFile = "./data/GeneData.txt";
  const char* exonDataFile = "./data/ExonData.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvKeepWinPathTest[] = {programName, 
                                       "--cel-files", celsFile, 
                                       "-g", geneDataFile,
                                       "-e", exonDataFile, 
                                       "-m", metaFile,
                                       "-o", outputDir,
                                       NULL};

  opts.parseArgv(argvKeepWinPathTest);
  std::string msg = midasCreateDirectory(opts.get("out-dir"));
  CPPUNIT_ASSERT (msg == "");

  bool configureRunOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    const string cvsId ("$Id: MidasNormalTest.cpp,v 1.25 2009-09-18 03:37:29 mspald Exp $");
    const string execVersion = version + " " + cvsId;
    // create object to configure, run midas
    midasConfigureRun configureRun (opts.get("cel-files"), opts.get("genedata"),
                                    opts.get("exondata"), opts.get("metaprobeset"),  opts.get("out-dir"),
                                    opts.getBool("pvalues"), opts.getBool("fstats"), opts.getBool("normalized"),
                                    opts.getDouble("stabilize"), opts.commandLine(), execVersion);
    // configure step may return a non-fatal warning message
    std::string* msg = configureRun.configure();
    CPPUNIT_ASSERT (msg == 0);
    // run the midas engine
    configureRun.run();
    configureRunOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    bool caughtException = true;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (configureRunOk == true);

  // check outputs
  const string probeset_list_id_string (PROBESET_LIST_ID);
  const string pvalue_string ("pvalue");
  vector<string> null_vector;

  //
  affx::TsvFile Pvalues_tsv;
  Pvalues_tsv.open(pvaluesFileName);
  int32_t probeset_list_id;
  float pvalue;
  int rv;
  //
  rv=Pvalues_tsv.nextLevel(0);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  rv=Pvalues_tsv.get(0,probeset_list_id_string,probeset_list_id);
  CPPUNIT_ASSERT( probeset_list_id == 2590411);
  Pvalues_tsv.get(0,pvalue_string,pvalue);
  CPPUNIT_ASSERT( floatsCloseEnough (pvalue, 0.147233, 5) == true);
  //
  Pvalues_tsv.close();
}

void MidasNormalTest::testKeepPathKeep()
{
  PgOptions opts;
  define_midas_opts(&opts);

  // delete output from previous run, if any
  const char* outputDir = "output";
  string outDirectory (outputDir);
  string pvaluesFileName = Fs::join(outDirectory,PVALUES_OUTPUT);

  unlink (pvaluesFileName.c_str());

  // test cels file, gene, exon data containing prefixes, --keep-path
  const char* programName = "MidasNormalTest";
  const char* celsFile = "./data/PathCels.txt";
  const char* metaFile = "./data/Meta.txt";
  const char* geneDataFile = "./data/PathGeneData.txt";
  const char* exonDataFile = "./data/PathExonData.txt";

  // set up dummy argc, argv, run as if called from main()
  const char* argvKeepPathKeepTest[] = {programName,
                                        "--cel-files", celsFile, 
                                        "-g", geneDataFile,
                                        "-e", exonDataFile, 
                                        "-m", metaFile, 
                                        "-o", outputDir, 
                                        "--keep-path", 
                                        NULL};

  opts.parseArgv(argvKeepPathKeepTest);
  std::string msg = midasCreateDirectory(opts.get("out-dir"));
  CPPUNIT_ASSERT (msg == "");

  bool configureRunOk = false;
  try
  {
    const string version ("NON-OFFICIAL-RELEASE");
    const string cvsId ("$Id: MidasNormalTest.cpp,v 1.25 2009-09-18 03:37:29 mspald Exp $");
    const string execVersion = version + " " + cvsId;
    // create object to configure, run midas
    midasConfigureRun configureRun (opts.get("cel-files"), opts.get("genedata"),
                                    opts.get("exondata"), opts.get("metaprobeset"),  opts.get("out-dir"),
                                    opts.getBool("pvalues"), opts.getBool("fstats"), opts.getBool("normalized"),
                                    opts.getDouble("stabilize"), opts.commandLine(), execVersion, opts.getBool("no-logtrans"),
                                    opts.getBool("keep-path"));
    // configure step may return a non-fatal warning message
    
    std::string* msg = configureRun.configure();
    CPPUNIT_ASSERT (msg == 0);
    // run the midas engine
    configureRun.run();
    configureRunOk = true;
  }
  catch (exception& e)
  {
    // should not execute this code
    bool caughtException = true;
    cout << e.what() << endl;
    CPPUNIT_ASSERT (caughtException == false);
  }
  CPPUNIT_ASSERT (configureRunOk == true);

  // check outputs
  const string probeset_list_id_string (PROBESET_LIST_ID);
  const string pvalue_string ("pvalue");
  vector<string> null_vector;

  //
  affx::TsvFile Pvalues_tsv;
  Pvalues_tsv.open(pvaluesFileName);
  int32_t probeset_list_id;
  float pvalue;
  int rv;
  //
  rv=Pvalues_tsv.nextLevel(0);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  rv=Pvalues_tsv.get(0,probeset_list_id_string,probeset_list_id);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  CPPUNIT_ASSERT( probeset_list_id == 2590411);
  Pvalues_tsv.get(0,pvalue_string,pvalue);
  CPPUNIT_ASSERT(rv == affx::TSV_OK);
  CPPUNIT_ASSERT( floatsCloseEnough (pvalue, 0.147233, 5) == true);
  //
  Pvalues_tsv.close();
}
