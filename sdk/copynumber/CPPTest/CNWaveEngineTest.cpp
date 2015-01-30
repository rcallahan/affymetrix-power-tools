////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
#include "copynumber/CNWaveEngine.h"
#include "copynumber/CPPTest/Setup.h"
//
#include "util/Util.h"
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
using namespace std;
/**
 * @class CNWaveEngineTest
 * @brief cppunit class for testing CNWaveEngine functions.
 * last change by vliber on 11/02/09
 */

class CNWaveEngineTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNWaveEngineTest); 
  CPPUNIT_TEST(defaultDefineOptionsTest);
  CPPUNIT_TEST(addCychpTest);
  CPPUNIT_TEST(checkOptionsTest);
  CPPUNIT_TEST(getAnnotationParameterTest);
  CPPUNIT_TEST(isCyto2ReferenceTest);
  CPPUNIT_TEST(isCopyNumberReferenceTest);
  

 
  
  CPPUNIT_TEST_SUITE_END();

public:  
  void defaultDefineOptionsTest();
  void addCychpTest();
  void checkOptionsTest();
  void getAnnotationParameterTest();
  void isCyto2ReferenceTest();
  void isCopyNumberReferenceTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNWaveEngineTest );



void CNWaveEngineTest::defaultDefineOptionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNWaveEngineTest::defaultDefineOptionsTest****");	
	CNWaveEngine cnWave;
	CPPUNIT_ASSERT(cnWave.getOpt("cn-reference-input")==""); 
	CPPUNIT_ASSERT(cnWave.getOpt("config-file")==""); 
	CPPUNIT_ASSERT(cnWave.getOpt("cychp-files")=="");
	CPPUNIT_ASSERT(cnWave.getOpt("cn-reference-output")==""); 
	CPPUNIT_ASSERT(cnWave.getOpt("analysis")=="additional-waves-reference-method"); 
	CPPUNIT_ASSERT(cnWave.getOpt("explain")==""); 
	CPPUNIT_ASSERT(cnWave.getOptInt("xChromosome")==24);
	CPPUNIT_ASSERT(cnWave.getOptInt("yChromosome")==25);
    CPPUNIT_ASSERT(cnWave.getOpt("cychps")=="");
}

void CNWaveEngineTest::addCychpTest()
{
	Verbose::out(1, "****CNWaveEngineTest::addCychpTest****");	
	CNWaveEngine cnWave;
	POSITIVE_TEST(cnWave.addCychp(INPUT+ "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel"));
	POSITIVE_TEST(cnWave.addCychp(INPUT+ "/Cyto/HapMap-As_NA18540_A01_01_CytoF_NN_20090121.cel"));
	POSITIVE_TEST(cnWave.addCychp(INPUT+ "/Cyto/HapMap-As_NA18547_A08_01_CytoF_NN_20090121.cel"));
	std::vector< std::string > values;
	values=cnWave.getOptVector("cychps");
	CPPUNIT_ASSERT(values.size()==3);
}

void CNWaveEngineTest::checkOptionsTest()
{
	Verbose::out(1, "****CNWaveEngineTest::checkOptionsTest****");

	//Message: Option temp-reference-file cannot be found in the options for this engine.
	CNWaveEngine cnWave1;
	NEGATIVE_TEST(cnWave1.getOpt("temp-reference-file"),Except);

	//Must specify an cn-reference-output.
	CNWaveEngine cnWave2;
	NEGATIVE_TEST(cnWave2.checkOptions(),Except);

	//The cn-reference-input specified either does not exist, or is not the correct type for the CNWaveEngine.
	CNWaveEngine cnWave3;
	cnWave3.setOpt("cn-reference-output",OUTPUT + "/Cyto/ref.a5");
	cnWave3.setOpt("cn-reference-input",INPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	NEGATIVE_TEST(cnWave3.checkOptions(),Except);

    //The cn-reference-output must not be the same as the cn-reference-input
	CNWaveEngine cnWave4;
	cnWave4.setOpt("cn-reference-output",INPUT + "/Cyto/ref.a5");
	cnWave4.setOpt("cn-reference-input",INPUT + "/Cyto/ref.a5");
	NEGATIVE_TEST(cnWave4.checkOptions(),Except);
	
	
	//No CYCHP files specified
	CNWaveEngine cnWave5;
	cnWave5.setOpt("cn-reference-output",OUTPUT + "/Cyto/ref.a5");
    cnWave5.setOpt("cn-reference-input",INPUT + "/Cyto/ref.a5");
	NEGATIVE_TEST(cnWave5.checkOptions(),Except);

	
	//positive. todo should be negative
	//Message: A file specified as a CYCHP input does not exist:
	//../../regression-data/data/copynumber-cyto/cppunit/input/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121_test.cel
	CNWaveEngine cnWave6;
	POSITIVE_TEST(cnWave6.addCychp(INPUT+ "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121_test.cel"));
	cnWave6.setOpt("cn-reference-output",OUTPUT + "/Cyto/ref.a5");
    cnWave6.setOpt("cn-reference-input",INPUT+ "/Cyto/ref.a5");
	NEGATIVE_TEST(cnWave6.checkOptions(),Except);

    //positive
	CNWaveEngine cnWave7;
	POSITIVE_TEST(cnWave7.addCychp(INPUT+ "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel"));
	cnWave7.setOpt("cn-reference-output",OUTPUT + "/Cyto/ref.a5");
    cnWave7.setOpt("cn-reference-input",INPUT+ "/Cyto/ref.a5");
	POSITIVE_TEST(cnWave7.checkOptions());
}


void CNWaveEngineTest::getAnnotationParameterTest()
{
    Verbose::out(1, "****CNWaveEngineTest::getAnnotationParameterTest****");
   	CPPUNIT_ASSERT(CNWaveEngine::getAnnotationParameter(INPUT+ "/Cyto/ref.a5","affymetrix-array-type")=="CytogeneticsFocused_Array");
	CPPUNIT_ASSERT(CNWaveEngine::getAnnotationParameter(INPUT+ "/Cyto/ref.a5","affymetrix-algorithm-param-option-verbose")=="4");
    //todo vliber no indication that file is bad
	POSITIVE_TEST(CNWaveEngine::getAnnotationParameter(INPUT+ "/Cyto/ref1.a5","affymetrix-array-type"));

}
void CNWaveEngineTest::isCyto2ReferenceTest()
{
    Verbose::out(1, "****CNWaveEngineTest::isCyto2ReferenceTest****");
   	CPPUNIT_ASSERT(CNWaveEngine::isCyto2Reference(INPUT+ "/Cyto/ref.a5")==true);
	CPPUNIT_ASSERT(CNWaveEngine::isCyto2Reference(INPUT+ "/Cyto/ref1.a5")==false);
}

void CNWaveEngineTest::isCopyNumberReferenceTest()
{
    Verbose::out(1, "****CNWaveEngineTest::isCopyNumberReferenceTest****");
   	CPPUNIT_ASSERT(CNWaveEngine::isCopyNumberReference(INPUT+ "/Cyto/ref.a5")==false);	
}
	

