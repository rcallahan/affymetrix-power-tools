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
#include "copynumber/CNCorrelationEngine.h"
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
 * @class CNCorrelationEngineTest
 * @brief cppunit class for testing CNCorrelationEngine functions.
 * last change by vliber on 11/20/09
 */

class CNCorrelationEngineTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNCorrelationEngineTest); 
  CPPUNIT_TEST(defaultDefineOptionsTest);
  CPPUNIT_TEST(checkOptionsTest);
  CPPUNIT_TEST(runTest);
  
  CPPUNIT_TEST_SUITE_END();

public:  
  void defaultDefineOptionsTest();
  void checkOptionsTest();
  void runTest();  
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNCorrelationEngineTest );



void CNCorrelationEngineTest::defaultDefineOptionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNCorrelationEngineTest::defaultDefineOptionsTest****");	
	CNCorrelationEngine cnCyto;
	CPPUNIT_ASSERT(cnCyto.getOpt("cychp-files")==""); 
	CPPUNIT_ASSERT(cnCyto.getOptInt("start-file-number")==-1);
	CPPUNIT_ASSERT(cnCyto.getOptInt("end-file-number")==-1); 
	CPPUNIT_ASSERT(cnCyto.getOpt("a5-output")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("text-output")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("explain")=="");
	CPPUNIT_ASSERT(cnCyto.getOptInt("xChromosome")==24);
	CPPUNIT_ASSERT(cnCyto.getOptInt("yChromosome")==25);
    CPPUNIT_ASSERT(cnCyto.getOpt("cychps")=="");
}

void CNCorrelationEngineTest::checkOptionsTest()
{
	Verbose::out(1, "****CNCorrelationEngineTest::checkOptionsTest****");
	//Fatal Error: Please specify the cychp-files parameter
	CNCorrelationEngine cnCyto1;
	POSITIVE_TEST(cnCyto1.setOpt("cychp-files",""));
	NEGATIVE_TEST(cnCyto1.checkOptions(), Except);

	//Fatal Error: No CYCHP files specified.
	CNCorrelationEngine cnCyto2;
   	POSITIVE_TEST(cnCyto2.setOpt("cychp-files",INPUT + "/Cyto/cychpFiles_list_bad1.txt"));
	NEGATIVE_TEST(cnCyto2.checkOptions(), Except);

	//Fatal Error: A file specified as a CYCHP input does not exist: ....HapMap-As_NA18603_B03_01_NN_20081218.ca_test.cychp
	CNCorrelationEngine cnCyto3;
   	POSITIVE_TEST(cnCyto3.setOpt("cychp-files",INPUT + "/Cyto/cychpFiles_list_bad2.txt"));
	NEGATIVE_TEST(cnCyto3.checkOptions(), Except);

	CNCorrelationEngine cnCyto4;
   	POSITIVE_TEST(cnCyto4.setOpt("cychp-files",INPUT + "/Cyto/cychpFiles_list.txt"));
	CPPUNIT_ASSERT(cnCyto4.getOptVector("cychp-files").size()==1);
	CPPUNIT_ASSERT(cnCyto4.getOptVector("cychps").size()==0);
	POSITIVE_TEST(cnCyto4.checkOptions());
	CPPUNIT_ASSERT(cnCyto4.getOptVector("cychps").size()==3);

}

void CNCorrelationEngineTest::runTest()
{
	Verbose::out(1, "****CNCorrelationEngineTest::runTest****");
	//Fatal Error: a5-output must be specified, text-output is optional.
	CNCorrelationEngine cnCyto1;
	POSITIVE_TEST(cnCyto1.setOpt("cychp-files",INPUT + "/Cyto/cychpFiles_list.txt"));
	POSITIVE_TEST(cnCyto1.setOpt("a5-output",""));
	NEGATIVE_TEST(cnCyto1.run(), Except);

	//Fatal Error: start-file-number must be greater than or equal to 1.
	CNCorrelationEngine cnCyto2;
	POSITIVE_TEST(cnCyto2.setOpt("cychp-files",INPUT + "/Cyto/cychpFiles_list.txt"));
	POSITIVE_TEST(cnCyto2.setOpt("a5-output",OUTPUT+"/test/a5"));
	POSITIVE_TEST(cnCyto2.setOpt("start-file-number","0"));
	NEGATIVE_TEST(cnCyto2.run(), Except);

	//Fatal Error: end-file-number must be greater than or equal to 1.
	CNCorrelationEngine cnCyto3;
	POSITIVE_TEST(cnCyto3.setOpt("cychp-files",INPUT + "/Cyto/cychpFiles_list.txt"));
	POSITIVE_TEST(cnCyto3.setOpt("a5-output",OUTPUT+"/test/a5"));
	POSITIVE_TEST(cnCyto3.setOpt("start-file-number","1"));
	POSITIVE_TEST(cnCyto3.setOpt("end-file-number","0"));
	NEGATIVE_TEST(cnCyto3.run(), Except);

	//Fatal Error: start-file-number must be less than or equal to 3.
	CNCorrelationEngine cnCyto4;
	POSITIVE_TEST(cnCyto4.setOpt("a5-output",OUTPUT+"/test/a5"));
	POSITIVE_TEST(cnCyto4.setOpt("start-file-number","4"));
	POSITIVE_TEST(cnCyto4.setOpt("end-file-number","1"));
	POSITIVE_TEST(cnCyto4.setOpt("cychp-files",INPUT + "/Cyto/cychpFiles_list.txt"));
	NEGATIVE_TEST(cnCyto4.run(), Except);

	//Fatal Error: end-file-number must be less than or equal to 3. 
	CNCorrelationEngine cnCyto5;
	POSITIVE_TEST(cnCyto5.setOpt("a5-output",OUTPUT+"/test/a5"));
	POSITIVE_TEST(cnCyto5.setOpt("start-file-number","1"));
	POSITIVE_TEST(cnCyto5.setOpt("end-file-number","4"));
	POSITIVE_TEST(cnCyto5.setOpt("cychp-files",INPUT + "/Cyto/cychpFiles_list.txt"));
	NEGATIVE_TEST(cnCyto5.run(), Except);

	//Fatal Error: end-file-number must be greater than or equal to start-file-number.
	CNCorrelationEngine cnCyto6;
	POSITIVE_TEST(cnCyto6.setOpt("a5-output",OUTPUT+"/test/a5"));
	POSITIVE_TEST(cnCyto6.setOpt("cychp-files",INPUT + "/Cyto/cychpFiles_list.txt"));
	POSITIVE_TEST(cnCyto6.setOpt("start-file-number","2"));
	POSITIVE_TEST(cnCyto6.setOpt("end-file-number","1"));
	NEGATIVE_TEST(cnCyto6.run(), Except);	
}
	

	

