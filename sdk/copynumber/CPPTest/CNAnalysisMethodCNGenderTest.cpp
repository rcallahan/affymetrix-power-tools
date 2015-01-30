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

#include "copynumber/CNAnalysisMethodCNGender.h"
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
 * @class CNAnalysisMethodCNGenderTest
 * @brief cppunit class for testing CNAnalysisMethodCNGender functions.
 * last change by vliber on 02/10/09
 */

class CNAnalysisMethodCNGenderTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodCNGenderTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodCNGenderTest );

void CNAnalysisMethodCNGenderTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodCNGenderTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodCNGender cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="cn-gender");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber CNGender");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="cn-gender");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodCNGender::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==9);
	CPPUNIT_ASSERT(sv[0].name=="male-chrX-lower-threshold");
	CPPUNIT_ASSERT(sv[0].type==1);
	CPPUNIT_ASSERT(sv[0].value=="0.8");
	CPPUNIT_ASSERT(sv[0].defaultVal=="0.8");
	CPPUNIT_ASSERT(sv[0].minVal=="0");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="Male ChrX Lower Threshold");
	CPPUNIT_ASSERT(sv[1].name=="male-chrX-upper-threshold");
	CPPUNIT_ASSERT(sv[1].type==1);
	CPPUNIT_ASSERT(sv[1].value=="1.3");
	CPPUNIT_ASSERT(sv[1].defaultVal=="1.3");
	CPPUNIT_ASSERT(sv[1].minVal=="0");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="Male ChrX Upper Threshold");
	CPPUNIT_ASSERT(sv[2].name=="male-chrY-lower-threshold");
	CPPUNIT_ASSERT(sv[2].type==1);
	CPPUNIT_ASSERT(sv[2].value=="0.8");
	CPPUNIT_ASSERT(sv[2].defaultVal=="0.8");
	CPPUNIT_ASSERT(sv[2].minVal=="0");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="Male ChrY Lower Threshold");
	CPPUNIT_ASSERT(sv[3].name=="male-chrY-upper-threshold");
	CPPUNIT_ASSERT(sv[3].type==1);
	CPPUNIT_ASSERT(sv[3].value=="1.2");
	CPPUNIT_ASSERT(sv[3].defaultVal=="1.2");
	CPPUNIT_ASSERT(sv[3].minVal=="0");
	CPPUNIT_ASSERT(sv[3].maxVal=="NA");
	CPPUNIT_ASSERT(sv[3].descript=="Male ChrY Upper Threshold");
	CPPUNIT_ASSERT(sv[4].name=="female-chrX-lower-threshold");
	CPPUNIT_ASSERT(sv[4].type==1);
	CPPUNIT_ASSERT(sv[4].value=="1.9");
	CPPUNIT_ASSERT(sv[4].defaultVal=="1.9");
	CPPUNIT_ASSERT(sv[4].minVal=="0");
	CPPUNIT_ASSERT(sv[4].maxVal=="NA");
	CPPUNIT_ASSERT(sv[4].descript=="Female ChrX Lower Threshold");
	CPPUNIT_ASSERT(sv[5].name=="female-chrX-upper-threshold");
	CPPUNIT_ASSERT(sv[5].type==1);
	CPPUNIT_ASSERT(sv[5].value=="2.1");
	CPPUNIT_ASSERT(sv[5].defaultVal=="2.1");
	CPPUNIT_ASSERT(sv[5].minVal=="0");
	CPPUNIT_ASSERT(sv[5].maxVal=="NA");
	CPPUNIT_ASSERT(sv[5].descript=="Female ChrX Upper Threshold");
	CPPUNIT_ASSERT(sv[6].name=="female-chrY-lower-threshold");
	CPPUNIT_ASSERT(sv[6].type==1);
	CPPUNIT_ASSERT(sv[6].value=="0");
	CPPUNIT_ASSERT(sv[6].defaultVal=="0");
	CPPUNIT_ASSERT(sv[6].minVal=="0");
	CPPUNIT_ASSERT(sv[6].maxVal=="NA");
	CPPUNIT_ASSERT(sv[6].descript=="Female ChrY Lower Threshold");
	CPPUNIT_ASSERT(sv[7].name=="female-chrY-upper-threshold");
	CPPUNIT_ASSERT(sv[7].type==1);
	CPPUNIT_ASSERT(sv[7].value=="0.4");
	CPPUNIT_ASSERT(sv[7].defaultVal=="0.4");
	CPPUNIT_ASSERT(sv[7].minVal=="0");
	CPPUNIT_ASSERT(sv[7].maxVal=="NA");
	CPPUNIT_ASSERT(sv[7].descript=="Female ChrY Upper Threshold");
	CPPUNIT_ASSERT(sv[8].name=="mapd-threshold");
	CPPUNIT_ASSERT(sv[8].type==1);
	CPPUNIT_ASSERT(sv[8].value=="0.5");
	CPPUNIT_ASSERT(sv[8].defaultVal=="0.5");
	CPPUNIT_ASSERT(sv[8].minVal=="0");
	CPPUNIT_ASSERT(sv[8].maxVal=="1");
	CPPUNIT_ASSERT(sv[8].descript=="MAPD Threshold");

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodCNGender::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="cn-gender.male-chrX-lower-threshold=0.8.male-chrX-upper-threshold=1.3.male-chrY-lower-threshold=0.8.male-chrY-upper-threshold=1.2.female-chrX-lower-threshold=1.9.female-chrX-upper-threshold=2.1.female-chrY-lower-threshold=0.female-chrY-upper-threshold=0.4.mapd-threshold=0.5");
    CPPUNIT_ASSERT(sd.getDocName()=="cn-gender");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber CNGender");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==9);
	CPPUNIT_ASSERT(v[0].asString()=="0.8");
	CPPUNIT_ASSERT(v[1].asString()=="1.3");
	CPPUNIT_ASSERT(v[2].asString()=="0.8");
	CPPUNIT_ASSERT(v[3].asString()=="1.2");
	CPPUNIT_ASSERT(v[4].asString()=="1.9");
	CPPUNIT_ASSERT(v[5].asString()=="2.1");
	CPPUNIT_ASSERT(v[6].asString()=="0");
	CPPUNIT_ASSERT(v[7].asString()=="0.4");
	CPPUNIT_ASSERT(v[8].asString()=="0.5");
	CPPUNIT_ASSERT(sd.getDocOption("male-chrX-lower-threshold").asString()=="0.8");
	CPPUNIT_ASSERT(sd.getDocOption("male-chrX-upper-threshold").asString()=="1.3");
	CPPUNIT_ASSERT(sd.getDocOption("male-chrY-lower-threshold").asString()=="0.8");
	CPPUNIT_ASSERT(sd.getDocOption("male-chrY-upper-threshold").asString()=="1.2");
	CPPUNIT_ASSERT(sd.getDocOption("female-chrX-lower-threshold").asString()=="1.9");
	CPPUNIT_ASSERT(sd.getDocOption("female-chrX-upper-threshold").asString()=="2.1");
	CPPUNIT_ASSERT(sd.getDocOption("female-chrY-lower-threshold").asString()=="0");
	CPPUNIT_ASSERT(sd.getDocOption("female-chrY-upper-threshold").asString()=="0.4");
	CPPUNIT_ASSERT(sd.getDocOption("mapd-threshold").asString()=="0.5");
	
}

void CNAnalysisMethodCNGenderTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodCNGenderTest::newObjectTest****");

     //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodCNGender *cn2=(CNAnalysisMethodCNGender*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("cn-gender");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==9);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj1->at(0).GetName())=="affymetrix-algorithm-param-male-chrX-lower-threshold");
	CPPUNIT_ASSERT(obj1->at(0).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj1->at(0).GetValueFloat()==0.8f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj1->at(1).GetName())=="affymetrix-algorithm-param-male-chrX-upper-threshold");
	CPPUNIT_ASSERT(obj1->at(1).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj1->at(1).GetValueFloat()==1.3f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj1->at(2).GetName())=="affymetrix-algorithm-param-male-chrY-lower-threshold");
	CPPUNIT_ASSERT(obj1->at(2).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj1->at(2).GetValueFloat()==0.8f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj1->at(3).GetName())=="affymetrix-algorithm-param-male-chrY-upper-threshold");
	CPPUNIT_ASSERT(obj1->at(3).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj1->at(3).GetValueFloat()==1.2f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj1->at(4).GetName())=="affymetrix-algorithm-param-female-chrX-lower-threshold");
	CPPUNIT_ASSERT(obj1->at(4).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj1->at(4).GetValueFloat()==1.9f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj1->at(5).GetName())=="affymetrix-algorithm-param-female-chrX-upper-threshold");
	CPPUNIT_ASSERT(obj1->at(5).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj1->at(5).GetValueFloat()==2.1f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj1->at(6).GetName())=="affymetrix-algorithm-param-female-chrY-lower-threshold");
	CPPUNIT_ASSERT(obj1->at(6).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj1->at(6).GetValueFloat()==0.0f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj1->at(7).GetName())=="affymetrix-algorithm-param-female-chrY-upper-threshold");
	CPPUNIT_ASSERT(obj1->at(7).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj1->at(7).GetValueFloat()==0.4f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj1->at(8).GetName())=="affymetrix-algorithm-param-mapd-threshold");
	CPPUNIT_ASSERT(obj1->at(8).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj1->at(8).GetValueFloat()==0.5f);
	delete cn2;


    //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["male-chrX-lower-threshold"]="0.1";
	params["male-chrX-upper-threshold"]="0.2";
	params["male-chrY-lower-threshold"]="0.3";
	params["male-chrY-upper-threshold"]="0.4";
	params["female-chrX-lower-threshold"]="0.5";
	params["female-chrX-upper-threshold"]="0.6";
	params["female-chrY-lower-threshold"]="0.7";
	params["female-chrY-upper-threshold"]="0.8";
	params["mapd-threshold"]="0.9";
	
	CNAnalysisMethod::getParams()->clear();
	SelfCreate *sc=CNAnalysisMethodCNGender::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==9);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj->at(0).GetName())=="affymetrix-algorithm-param-male-chrX-lower-threshold");
	CPPUNIT_ASSERT(obj->at(0).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj->at(0).GetValueFloat()==0.1f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj->at(1).GetName())=="affymetrix-algorithm-param-male-chrX-upper-threshold");
	CPPUNIT_ASSERT(obj->at(1).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj->at(1).GetValueFloat()==0.2f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj->at(2).GetName())=="affymetrix-algorithm-param-male-chrY-lower-threshold");
	CPPUNIT_ASSERT(obj->at(2).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj->at(2).GetValueFloat()==0.3f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj->at(3).GetName())=="affymetrix-algorithm-param-male-chrY-upper-threshold");
	CPPUNIT_ASSERT(obj->at(3).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj->at(3).GetValueFloat()==0.4f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj->at(4).GetName())=="affymetrix-algorithm-param-female-chrX-lower-threshold");
	CPPUNIT_ASSERT(obj->at(4).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj->at(4).GetValueFloat()==0.5f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj->at(5).GetName())=="affymetrix-algorithm-param-female-chrX-upper-threshold");
	CPPUNIT_ASSERT(obj->at(5).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj->at(5).GetValueFloat()==0.6f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj->at(6).GetName())=="affymetrix-algorithm-param-female-chrY-lower-threshold");
	CPPUNIT_ASSERT(obj->at(6).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj->at(6).GetValueFloat()==0.7f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj->at(7).GetName())=="affymetrix-algorithm-param-female-chrY-upper-threshold");
	CPPUNIT_ASSERT(obj->at(7).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj->at(7).GetValueFloat()==0.8f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(obj->at(8).GetName())=="affymetrix-algorithm-param-mapd-threshold");
	CPPUNIT_ASSERT(obj->at(8).GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(obj->at(8).GetValueFloat()==0.9f);
	
		
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodCNGender::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="cn-gender.male-chrX-lower-threshold=0.8.male-chrX-upper-threshold=1.3.male-chrY-lower-threshold=0.8.male-chrY-upper-threshold=1.2.female-chrX-lower-threshold=1.9.female-chrX-upper-threshold=2.1.female-chrY-lower-threshold=0.female-chrY-upper-threshold=0.4.mapd-threshold=0.5");
	delete sc;

	//Negative parameters tests
   	//FATAL ERROR: SelfCreate::setValue() - '-0.1' is not a valid value for parameter: 'male-chrX-lower-threshold'. The specified range is 0 to NA
	params["male-chrX-lower-threshold"]="-0.1";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
	//FATAL ERROR: Could not convert 'a' to a double.
	params["male-chrX-lower-threshold"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);

	params.clear();
	//FATAL ERROR: SelfCreate::setValue() - '-0.2' is not a valid value for parameter: 'male-chrX-upper-threshold'. The specified range is 0 to NA
	params["male-chrX-upper-threshold"]="-0.2";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
	//FATAL ERROR: Could not convert 'a' to a double.
	params["male-chrX-upper-threshold"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);

	params.clear();
	//FATAL ERROR: SelfCreate::setValue() - '-0.2' is not a valid value for parameter: 'male-chrY-lower-threshold'. The specified range is 0 to NA
	params["male-chrY-lower-threshold"]="-0.2";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
	//FATAL ERROR: Could not convert 'a' to a double.
	params["male-chrY-lower-threshold"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);

	params.clear();
    //FATAL ERROR: SelfCreate::setValue() - '-0.1' is not a valid value for parameter: 'male-chrY-upper-threshold'. The specified range is 0 to NA
	params["male-chrY-upper-threshold"]="-0.1";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
    //FATAL ERROR: Could not convert 'a' to a double.
	params["male-chrY-upper-threshold"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);

    
	params.clear();
	//FATAL ERROR: SelfCreate::setValue() - '-0.1' is not a valid value for parameter: 'female-chrX-lower-threshold'. The specified range is 0 to NA
	params["female-chrX-lower-threshold"]="-0.1";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
	//FATAL ERROR: Could not convert 'a' to a double.
	params["female-chrX-lower-threshold"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);

	params.clear();
	//FATAL ERROR: SelfCreate::setValue() - '-0.1' is not a valid value for parameter: 'female-chrX-upper-threshold'. The specified range is 0 to NA
	params["female-chrX-upper-threshold"]="-0.1";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
	//FATAL ERROR: Could not convert 'a' to a double.
	params["female-chrX-upper-threshold"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);

	params.clear();
	//FATAL ERROR: SelfCreate::setValue() - '-0.1' is not a valid value for parameter: 'female-chrY-lower-threshold'. The specified range is 0 to NA
	params["female-chrY-lower-threshold"]="-0.1";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
	//FATAL ERROR: Could not convert 'a' to a double.
	params["female-chrY-lower-threshold"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);

	params.clear();
	//FATAL ERROR: SelfCreate::setValue() - '-0.1' is not a valid value for parameter: 'female-chrY-upper-threshold'. The specified range is 0 to NA
	params["female-chrY-upper-threshold"]="-0.1";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
	//FATAL ERROR: Could not convert 'a' to a double.
	params["female-chrY-upper-threshold"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);

	params.clear();
	//FATAL ERROR: SelfCreate::setValue() - '-0.1' is not a valid value for parameter: 'mapd-threshold'. The specified range is 0 to 1
	params["mapd-threshold"]="-0.1";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
	//FATAL ERROR: Could not convert 'a' to a double.
	params["mapd-threshold"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNGender::newObject(params),Except);
	
}	

void CNAnalysisMethodCNGenderTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodCNGenderTest::runTest****");

    
	CNAnalysisMethodCNGender cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod cn-gender is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
