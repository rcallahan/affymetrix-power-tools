////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNAnalysisMethodChipstream.h"
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
#include <map>
#include <string>
//
using namespace std;
/**
 * @class CNAnalysisMethodChipstreamTest
 * @brief cppunit class for testing CNAnalysisMethodChipstream functions.
 * last change by vliber on 10/28/09
 */

class CNAnalysisMethodChipstreamTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodChipstreamTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodChipstreamTest );

void CNAnalysisMethodChipstreamTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodChipstreamTest::functionsTest****");
	
	CNAnalysisMethodChipstream cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="chipstream");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber Chipstream");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="chipstream");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodChipstream::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==1);
	CPPUNIT_ASSERT(sv[0].name=="NDBandwidth");
	CPPUNIT_ASSERT(sv[0].type==2);
	CPPUNIT_ASSERT(sv[0].value=="10.0");
	CPPUNIT_ASSERT(sv[0].defaultVal=="10.0");
	CPPUNIT_ASSERT(sv[0].minVal=="10.0");
	CPPUNIT_ASSERT(sv[0].maxVal=="20.0");
	CPPUNIT_ASSERT(sv[0].descript=="Normal-Diploid Bandwidth");

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodChipstream::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="chipstream.NDBandwidth=10.0");
    CPPUNIT_ASSERT(sd.getDocName()=="chipstream");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber Chipstream");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==1);
	CPPUNIT_ASSERT(v[0].asString()=="10.0");
	CPPUNIT_ASSERT(sd.getDocOption("NDBandwidth").asString()=="10.0");
	
}

void CNAnalysisMethodChipstreamTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodChipstreamTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodChipstream *cn2=(CNAnalysisMethodChipstream*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("chipstream");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-NDBandwidth");
	CPPUNIT_ASSERT(param1a.GetParameterType()==6);
    CPPUNIT_ASSERT(param1a.GetValueFloat()==10.0f);
	delete cn2;

    //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["NDBandwidth"]="10.111";
	
	CNAnalysisMethod::getParams()->clear();
	SelfCreate *sc=CNAnalysisMethodChipstream::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-NDBandwidth");
	CPPUNIT_ASSERT(param1.GetParameterType()==6);
    CPPUNIT_ASSERT(param1.GetValueFloat()==10.111f);
	
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodChipstream::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="chipstream.NDBandwidth=10.0");
    delete sc;	
    
	//Negative test passed [CNAnalysisMethodChipstream::newObject(params)] Message: SelfCreate::setValue() - '9.99' is not a valid value for parameter: 'NDBandwidth'. The specified range is 10.0 to 20.0
	params["NDBandwidth"]="9.99";
    NEGATIVE_TEST(CNAnalysisMethodChipstream::newObject(params), Except);
	//Negative test passed [CNAnalysisMethodChipstream::newObject(params)] Message: SelfCreate::setValue() - '20.01' is not a valid value for parameter: 'NDBandwidth'. The specified range is 10.0 to 20.0
	params.clear();
	params["NDBandwidth"]="20.01";
	NEGATIVE_TEST(CNAnalysisMethodChipstream::newObject(params),Except);	
}	

void CNAnalysisMethodChipstreamTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodChipstreamTest::runTest****");

    
	CNAnalysisMethodChipstream cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod chipstream is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
