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

#include "copynumber/CNAnalysisMethodCNCyto2Gender.h"
#include "copynumber/CNAnalysisMethodFactory.h"
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
 * @class CNAnalysisMethodCNCyto2GenderTest
 * @brief cppunit class for testing CNAnalysisMethodCNCyto2Gender functions.
 * last change by vliber on 01/28/09
 */

class CNAnalysisMethodCNCyto2GenderTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodCNCyto2GenderTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest);
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodCNCyto2GenderTest );

void CNAnalysisMethodCNCyto2GenderTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodCNCyto2GenderTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodCNCyto2Gender cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="cn-cyto2-gender");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber CNCyto2Gender");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="cn-cyto2-gender");
	
	
	//explainSelf()
	SelfDoc sd=CNAnalysisMethodCNCyto2Gender::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="cn-cyto2-gender.cutoff=0.5");
	CPPUNIT_ASSERT(sd.getDocName()=="cn-cyto2-gender");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber CNCyto2Gender");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==1);
	CPPUNIT_ASSERT(v[0].asString()=="0.5");
	CPPUNIT_ASSERT(sd.getDocOption("cutoff").asString()=="0.5");

	vector<SelfDoc::Opt> sv=CNAnalysisMethodCNCyto2Gender::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==1);
	CPPUNIT_ASSERT(sv[0].name=="cutoff");
	CPPUNIT_ASSERT(sv[0].type==1);
	CPPUNIT_ASSERT(sv[0].value=="0.5");
	CPPUNIT_ASSERT(sv[0].defaultVal=="0.5");
	CPPUNIT_ASSERT(sv[0].minVal=="0");
	CPPUNIT_ASSERT(sv[0].maxVal=="0.5");
	CPPUNIT_ASSERT(sv[0].descript=="Allele Peaks cutoff.");
}
void CNAnalysisMethodCNCyto2GenderTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodCNCyto2GenderTest::newObjectTest****");
     //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodCNCyto2Gender *cn2=(CNAnalysisMethodCNCyto2Gender*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("cn-cyto2-gender");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-cutoff");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1a.GetValueFloat()==0.5f);
	delete cn2;

    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params;
	params["cutoff"]="0.1";
	
	SelfCreate *sc=CNAnalysisMethodCNCyto2Gender::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-cutoff");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1.GetValueFloat()==0.1f); 
		
	 //always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodCNCyto2Gender::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="cn-cyto2-gender.cutoff=0.5");
    delete sc;
    
	//boundary test
	params.clear();
	//SelfCreate::setValue() - '-0.1' is not a valid value for parameter: 'cutoff'. The specified range is 0 to 0.5
	params["cutoff"]="-0.1";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2Gender::newObject(params),Except);
	params.clear();
	params["cutoff"]="0.0";
    POSITIVE_TEST(CNAnalysisMethodCNCyto2Gender::newObject(params));
	params.clear();
	//SelfCreate::setValue() - '0.5001' is not a valid value for parameter: 'cutoff'. The specified range is 0 to 0.5
	params["cutoff"]="0.5001";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2Gender::newObject(params),Except);
	params.clear();
	params["cutoff"]="0.5";
    POSITIVE_TEST(CNAnalysisMethodCNCyto2Gender::newObject(params));
}	


void CNAnalysisMethodCNCyto2GenderTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodCNCyto2GenderTest::runTest****");
    
	CNAnalysisMethodCNCyto2Gender cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod cn-cyto2-gender is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
