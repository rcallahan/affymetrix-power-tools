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
#include "copynumber/CNTrisomyEngine.h"
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
 * @brief cppunit class for testing CNTrisomyEngine functions.
 * last change by vliber on 11/02/09
 */

class CNTrisomyEngineTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNTrisomyEngineTest); 
  CPPUNIT_TEST(defaultDefineOptionsTest);
  CPPUNIT_TEST(checkOptionsTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void defaultDefineOptionsTest();
  void checkOptionsTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNTrisomyEngineTest );



void CNTrisomyEngineTest::defaultDefineOptionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNTrisomyEngineTest::defaultDefineOptionsTest****");	
	CNTrisomyEngine cnTrisomy;
	CPPUNIT_ASSERT(cnTrisomy.getOpt("in-file")==""); 
	CPPUNIT_ASSERT(cnTrisomy.getOpt("out-file")==""); 
	CPPUNIT_ASSERT(cnTrisomy.getOpt("chr-numbers")=="");
	CPPUNIT_ASSERT(cnTrisomy.getOpt("chr-display")==""); 
	CPPUNIT_ASSERT(cnTrisomy.getOpt("syndrome-tests")==""); 
	CPPUNIT_ASSERT(cnTrisomy.getOpt("syndrome-names")==""); 
	CPPUNIT_ASSERT(cnTrisomy.getOptDouble("confidence-threshold")==0.0);
	CPPUNIT_ASSERT(cnTrisomy.getOpt("gender-tag")=="Y-gender-call");
    CPPUNIT_ASSERT(cnTrisomy.getOpt("qc-names")=="MAPD,snp-qc");
	CPPUNIT_ASSERT(cnTrisomy.getOpt("qc-ops")=="le,ge");
	CPPUNIT_ASSERT(cnTrisomy.getOpt("qc-thrs")=="0.27,1.1");
	CPPUNIT_ASSERT(cnTrisomy.getOptBool("test")==false);    
}

void CNTrisomyEngineTest::checkOptionsTest()
{
	Verbose::out(1, "****CNTrisomyEngineTest::checkOptionsTest****");	
	CNTrisomyEngine cnTrisomy1;
	//Must specify an input file.
	NEGATIVE_TEST(cnTrisomy1.checkOptions(),Except);

	CNTrisomyEngine cnTrisomy2;
    cnTrisomy2.setOpt("in-file",INPUT+"/Cyto/ref.a5");
	//Must specify an output file.
	NEGATIVE_TEST(cnTrisomy2.checkOptions(),Except);

	CNTrisomyEngine cnTrisomy3;
    cnTrisomy3.setOpt("in-file",INPUT+"/Cyto/ref.a5");
	cnTrisomy3.setOpt("out-file",INPUT+"/Cyto/ref.a5");
	//Must specify chromosomes to analyze.
	NEGATIVE_TEST(cnTrisomy3.checkOptions(),Except);

	CNTrisomyEngine cnTrisomy4;
    cnTrisomy4.setOpt("in-file",INPUT+"/Cyto/ref.a5");
	cnTrisomy4.setOpt("out-file",INPUT+"/Cyto/ref.a5");
	cnTrisomy4.setOpt("chr-display",INPUT+"/Cyto/ref.a5");
	//Must specify the display values for each chromosome.
	NEGATIVE_TEST(cnTrisomy4.checkOptions(),Except);

	CNTrisomyEngine cnTrisomy5;
    cnTrisomy5.setOpt("in-file",INPUT+"/Cyto/ref.a5");
	cnTrisomy5.setOpt("out-file",INPUT+"/Cyto/ref.a5");
	cnTrisomy5.setOpt("chr-display",INPUT+"/Cyto/ref.a5");
	cnTrisomy5.setOpt("chr-numbers",INPUT+"/Cyto/ref.a5");
	//Must specify the syndrome names.
	NEGATIVE_TEST(cnTrisomy5.checkOptions(),Except);

	CNTrisomyEngine cnTrisomy6;
    cnTrisomy6.setOpt("in-file",INPUT+"/Cyto/ref.a5");
	cnTrisomy6.setOpt("out-file",INPUT+"/Cyto/ref.a5");
	cnTrisomy6.setOpt("chr-display",INPUT+"/Cyto/ref.a5");
	cnTrisomy6.setOpt("chr-numbers",INPUT+"/Cyto/ref.a5");
	cnTrisomy6.setOpt("syndrome-names",INPUT+"/Cyto/ref.a5");
	//Message: Must specify the tests for each syndrome.
	NEGATIVE_TEST(cnTrisomy6.checkOptions(),Except);

	CNTrisomyEngine cnTrisomy7;
    cnTrisomy7.setOpt("in-file",INPUT+"/Cyto/ref.a5");
	cnTrisomy7.setOpt("out-file",INPUT+"/Cyto/ref.a5");
	cnTrisomy7.setOpt("chr-display",INPUT+"/Cyto/ref.a5");
	cnTrisomy7.setOpt("chr-numbers",INPUT+"/Cyto/ref.a5");
	cnTrisomy7.setOpt("syndrome-names",INPUT+"/Cyto/ref.a5");
	cnTrisomy7.setOpt("syndrome-tests",INPUT+"/Cyto/ref.a5");
	POSITIVE_TEST(cnTrisomy7.checkOptions());
}