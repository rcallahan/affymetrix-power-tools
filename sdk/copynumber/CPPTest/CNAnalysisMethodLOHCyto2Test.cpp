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
#include "copynumber/CNAnalysisMethodLOHCyto2.h"
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
 * @class CNAnalysisMethodLOHCyto2Test
 * @brief cppunit class for testing CNAnalysisMethodLOHCyto2 functions.
 * last change by vliber on 02/09/09
 */

class CNAnalysisMethodLOHCyto2Test : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodLOHCyto2Test);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodLOHCyto2Test );

void CNAnalysisMethodLOHCyto2Test::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodLOHCyto2Test::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodLOHCyto2 cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="loh-cyto2");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="CopyNumber LOH Cyto2");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="loh-cyto2");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodLOHCyto2::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==3);
	CPPUNIT_ASSERT(sv[0].name=="lohCNSegSeparation");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="1000000");
	CPPUNIT_ASSERT(sv[0].defaultVal=="1000000");
	CPPUNIT_ASSERT(sv[0].minVal=="10");
	CPPUNIT_ASSERT(sv[0].maxVal=="10000000");
	CPPUNIT_ASSERT(sv[0].descript=="LOH CN Separation");
	CPPUNIT_ASSERT(sv[1].name=="minInformation");
	CPPUNIT_ASSERT(sv[1].type==2);
	CPPUNIT_ASSERT(sv[1].value=="100");
	CPPUNIT_ASSERT(sv[1].defaultVal=="100");
	CPPUNIT_ASSERT(sv[1].minVal=="NA");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="A control of window size.");
	CPPUNIT_ASSERT(sv[2].name=="lambdaCritical");
	CPPUNIT_ASSERT(sv[2].type==2);
	CPPUNIT_ASSERT(sv[2].value=="8.0");
	CPPUNIT_ASSERT(sv[2].defaultVal=="8.0");
	CPPUNIT_ASSERT(sv[2].minVal=="NA");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="A measure of required likelihood.");

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodLOHCyto2::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0");
    CPPUNIT_ASSERT(sd.getDocName()=="loh-cyto2");
	CPPUNIT_ASSERT(sd.getDocDescription()=="CopyNumber LOH Cyto2");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==3);
	CPPUNIT_ASSERT(v[0].asString()=="1000000");
	CPPUNIT_ASSERT(v[1].asString()=="100");
	CPPUNIT_ASSERT(v[2].asString()=="8.0");
	CPPUNIT_ASSERT(sd.getDocOption("lohCNSegSeparation").asString()=="1000000");
	CPPUNIT_ASSERT(sd.getDocOption("minInformation").asString()=="100");
	CPPUNIT_ASSERT(sd.getDocOption("lambdaCritical").asString()=="8.0");
}

void CNAnalysisMethodLOHCyto2Test::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodLOHCyto2Test::newObjectTest****");

     //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodLOHCyto2 *cn2=(CNAnalysisMethodLOHCyto2*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("loh-cyto2");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-lohCNSegSeparation");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1a.GetValueInt32()==1000000);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-minInformation");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2a.GetValueFloat()==100);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-lambdaCritical");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param3a.GetValueFloat()==8.0f);
	delete cn2;

    //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["lohCNSegSeparation"]="200000";
	params["minInformation"]="1.3";
	params["lambdaCritical"]="2.3";

	CNAnalysisMethod::getParams()->clear();
	SelfCreate *sc=CNAnalysisMethodLOHCyto2::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-lohCNSegSeparation");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1.GetValueInt32()==200000);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-minInformation");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2.GetValueFloat()==1.3f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-lambdaCritical");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param3.GetValueFloat()==2.3f);
		
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodLOHCyto2::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0");
    delete sc;	

	//test parameters
    params["lohCNSegSeparation"]="1000000";
	POSITIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params));
	//FATAL ERROR: SelfCreate::setValue() - '-10' is not a valid value for parameter: 'lohCNSegSeparation'.The specified range is 10 to 10000000
	params["lohCNSegSeparation"]="-10";
	NEGATIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params),Except);
	//FATAL ERROR: SelfCreate::setValue() - '200000000' is not a valid value for parameter: 'lohCNSegSeparation'. The specified range is 10 to 10000000
	params["lohCNSegSeparation"]="200000000";
    NEGATIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params),Except);
	//FATAL ERROR: SelfCreate::setValue() - '0' is not a valid value for parameter: 'lohCNSegSeparation'. The specified range is 10 to 10000000
	params["lohCNSegSeparation"]="0";
	NEGATIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params),Except);

	
	params.clear();
	params["lohCNSegSeparation"]="200000";
	params["minInformation"]="0";
	POSITIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params));
	params["minInformation"]="150.0";
	POSITIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params));
	//FATAL ERROR: Could not convert 'a' to a float.
    params["minInformation"]="a";
	NEGATIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params),Except);

	params.clear();
	params["lohCNSegSeparation"]="200000";
	params["minInformation"]="2.3";
	params["lambdaCritical"]="-11111";
	POSITIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params));
	params["lambdaCritical"]="0";
	POSITIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params));
	//FATAL ERROR: Could not convert 'a' to a float.
    params["lambdaCritical"]="a";
	NEGATIVE_TEST(CNAnalysisMethodLOHCyto2::newObject(params),Except);
	
}	

void CNAnalysisMethodLOHCyto2Test::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodLOHCyto2Test::runTest****");
    
	CNAnalysisMethodLOHCyto2 cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod loh-cyto2 is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
