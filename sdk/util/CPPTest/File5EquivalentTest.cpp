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

// You should have received a copy of the GNU General Public License
// along with this program;if not, write to the
// Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

////////////////////////////////////////////////////////////////

/**

 * @file   File5EquivalentTest.cpp
 * @author vliber
 * last change by rsatin on 10/13/09
 *
 @brief  Testing the File5_File::equivalent() function.
 *
 */


#include "util/AffxString.h"
#include "util/Convert.h"
#include "util/Util.h"
#include "util/Verbose.h"
#include "file5/File5_File.cpp"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
//
#include "util/CPPTest/Setup.h"

using namespace std;

class File5EquivalentTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( File5EquivalentTest );
  CPPUNIT_TEST(testDifferencesPositive);
  CPPUNIT_TEST(testEpsilonNegative);
//CPPUNIT_TEST(testIntegerNegative);       //rsatin TODO fix ClearQuest AFFY00024379 
  CPPUNIT_TEST(testStringNegative);
//CPPUNIT_TEST(testNonFiniteNegative1);    //rsatin TODO fix ClearQuest AFFY00024379 
//CPPUNIT_TEST(testNonFiniteNegative2);    //rsatin TODO fix ClearQuest AFFY00024379 
  CPPUNIT_TEST(testParametersNegative);
  CPPUNIT_TEST(testCorrelationPositive);
  CPPUNIT_TEST(testCorrelationNegative);
//CPPUNIT_TEST(testIgnoreListPositive);    //rsatin TODO fix ClearQuest AFFY00024379 
  CPPUNIT_TEST_SUITE_END();

  void testDifferencesPositive();
  void testEpsilonNegative();
  void testIntegerNegative();
  void testStringNegative();
  void testNonFiniteNegative1();
  void testNonFiniteNegative2();
  void testParametersNegative();
  void testCorrelationPositive();
  void testCorrelationNegative();
  void testIgnoreListPositive();
};
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(File5EquivalentTest );

void File5EquivalentTest::testDifferencesPositive()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  double dEpsilon; 

  Verbose::out(1, "***File5EquivalentTest testcases***");
  Verbose::out(1, "**File5EquivalentTest::testDifferencesPositive**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_false_diff2.ref.a5";
  CPPUNIT_ASSERT( affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "ProbeEffects",   setIgnore, dEpsilon=0.005 ) );
  CPPUNIT_ASSERT( affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "SNPPosteriors",  setIgnore, dEpsilon=0.005 ) );
  CPPUNIT_ASSERT( affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "SketchSNP",      setIgnore, dEpsilon=0.005 ) );
  CPPUNIT_ASSERT( affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "WaveCorrection", setIgnore, dEpsilon=0.005 ) );
}

// datasets SketchCN, AntigenomicProbes, MedianSignals and header parameters have non-allowed differences
void File5EquivalentTest::testEpsilonNegative()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  double dEpsilon; 

  Verbose::out(1, "**File5EquivalentTest::testEpsilonNegative**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_false_diff2.ref.a5";
  CPPUNIT_ASSERT(  affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "ProbeEffects",   setIgnore, dEpsilon=0.00000001 ) );
  CPPUNIT_ASSERT(  affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "SNPPosteriors",  setIgnore, dEpsilon=0.00000001 ) );
  CPPUNIT_ASSERT( !affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "SketchSNP",      setIgnore, dEpsilon=0.00000001 ) );
  CPPUNIT_ASSERT( !affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "WaveCorrection", setIgnore, dEpsilon=0.00000001 ) );
}

// dataset AntigenomicProbes.ProbeID int32 values have non-allowed single digit differences: 1037413 vs 1037412, etc.
void File5EquivalentTest::testIntegerNegative()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  setIgnore.insert("Cyto2.AntigenomicProbes.XXXXXXXX");
  double dEpsilon; 
  double dCorrelationCutoff;

  Verbose::out(1, "**File5EquivalentTest::testIntegerNegative**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_false_diff2.ref.a5";
  CPPUNIT_ASSERT( !affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "AntigenomicProbes", setIgnore, dEpsilon=0.005, dCorrelationCutoff=0.99 ) );
}

// dataset MedianSignals.probeset_id has string value has non-allowed difference: 'C-10XYZ' vs 'C-10T0Q'
void File5EquivalentTest::testStringNegative()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  setIgnore.insert("Cyto2.MedianSignals.MedianSignal");
  setIgnore.insert("Cyto2.MedianSignals.Position");
  double dEpsilon; 
  double dCorrelationCutoff;

  Verbose::out(1, "**File5EquivalentTest::testStringNegative**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_false_diff2.ref.a5";
  CPPUNIT_ASSERT( !affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "MedianSignals", setIgnore, dEpsilon=0.005, dCorrelationCutoff=0.99 ) );
}

// dataset SketchCN.Sketch has non-allowed difference: NaN vs 1782.38155470
void File5EquivalentTest::testNonFiniteNegative1()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  double dEpsilon; 

  Verbose::out(1, "**File5EquivalentTest::testNonFiniteNegative1**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_false_diff2.ref.a5";
  CPPUNIT_ASSERT( !affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "SketchCN", setIgnore, dEpsilon=0.005 ) );
}

// dataset SketchCN.Sketch has non-allowed difference: Inf vs 1782.38155470
void File5EquivalentTest::testNonFiniteNegative2()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  double dEpsilon;

  Verbose::out(1, "**File5EquivalentTest::testNonFiniteNegative2**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_true_diff2.ref.a5";
  CPPUNIT_ASSERT( !affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "SketchCN", setIgnore, dEpsilon=0.005 ) );
}

// some parameters are different with all values treated as text strings
void File5EquivalentTest::testParametersNegative()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  setIgnore.insert("MultiData.DmetBiAllelic.Signal");
  double dEpsilon;

  Verbose::out(1, "**File5EquivalentTest::testParametersNegative**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_false_diff2.ref.a5";
  CPPUNIT_ASSERT( !affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "Parameters", setIgnore, dEpsilon=0.1 ) );
}

// datasets WaveCorrection has corr=0.99999999820667 above threshold
void File5EquivalentTest::testCorrelationPositive()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  setIgnore.insert("Cyto2.WaveCorrection.X12");
  double dEpsilon;
  double dCorrelationCutoff;

  Verbose::out(1, "**File5EquivalentTest::testCorrelationPositive**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_false_diff2.ref.a5";
  CPPUNIT_ASSERT( affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "WaveCorrection", setIgnore, dEpsilon=0.00000001, dCorrelationCutoff=0.999 ) );
}

// datasets WaveCorrection has corr=0.1946 below threshold
void File5EquivalentTest::testCorrelationNegative()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  setIgnore.insert("Cyto2.WaveCorrection.X12");
  double dEpsilon;
  double dCorrelationCutoff;

  Verbose::out(1, "**File5EquivalentTest::testCorrelationNegative**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_true_diff2.ref.a5";
  CPPUNIT_ASSERT( !affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "WaveCorrection", setIgnore, dEpsilon=0.00000001, dCorrelationCutoff=0.90 ) );
}

// ignore column with integer differences Cyto2.AntigenomicProbes.ProbeID  
void File5EquivalentTest::testIgnoreListPositive()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::set<std::string> setIgnore;
  setIgnore.insert("Cyto2.Parameters.Parameter");
  setIgnore.insert("Cyto2.AntigenomicProbes.ProbeID");
  double dEpsilon;

  Verbose::out(1, "**File5EquivalentTest::testIgnoreListPositive**");
  std::string strFileName1 = "./input/a5/ref/demean_false_diff1.ref.a5";
  std::string strFileName2 = "./input/a5/ref/demean_false_diff2.ref.a5";
  CPPUNIT_ASSERT( affx::File5_File::equivalent( strFileName1, strFileName2, "Cyto2", "AntigenomicProbes", setIgnore, dEpsilon=0.001 ) );
}
