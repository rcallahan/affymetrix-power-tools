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
#include "copynumber/CNAnalysisMethodMosaicism.h"
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
 * @class CNAnalysisMethodMosaicismTest
 * @brief cppunit class for testing CNAnalysisMethodMosaicism functions.
 * last change by vliber on 03/17/09
 */

class CNAnalysisMethodMosaicismTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodMosaicismTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodMosaicismTest );

void CNAnalysisMethodMosaicismTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodMosaicismTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodMosaicism cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="mosaicism");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber Mosaicism");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="mosaicism");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodMosaicism::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==5);
	CPPUNIT_ASSERT(sv[0].name=="gains-boundries");
	CPPUNIT_ASSERT(sv[0].type==0);
	CPPUNIT_ASSERT(sv[0].value=="0.08764945,0.15380349,0.21465931,0.27100300");
	CPPUNIT_ASSERT(sv[0].defaultVal=="0.08764945,0.15380349,0.21465931,0.27100300");
	CPPUNIT_ASSERT(sv[0].minVal=="");
	CPPUNIT_ASSERT(sv[0].maxVal=="");
	CPPUNIT_ASSERT(sv[0].descript=="Mosaicism Gains Boundries");
	CPPUNIT_ASSERT(sv[1].name=="losses-boundries");
	CPPUNIT_ASSERT(sv[1].type==0);
	CPPUNIT_ASSERT(sv[1].value=="-0.08293345,-0.17551812,-0.28048196,-0.40165383");
	CPPUNIT_ASSERT(sv[1].defaultVal=="-0.08293345,-0.17551812,-0.28048196,-0.40165383");
	CPPUNIT_ASSERT(sv[1].minVal=="");
	CPPUNIT_ASSERT(sv[1].maxVal=="");
	CPPUNIT_ASSERT(sv[1].descript=="Mosaicism Losses Boundries");
	CPPUNIT_ASSERT(sv[2].name=="marker-bandwidth");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="6000");
	CPPUNIT_ASSERT(sv[2].defaultVal=="6000");
	CPPUNIT_ASSERT(sv[2].minVal=="1");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="Mosaicism Marker Bandwitdth");
	CPPUNIT_ASSERT(sv[3].name=="confidence-window");
	CPPUNIT_ASSERT(sv[3].type==3);
	CPPUNIT_ASSERT(sv[3].value=="251");
	CPPUNIT_ASSERT(sv[3].defaultVal=="251");
	CPPUNIT_ASSERT(sv[3].minVal=="31");
	CPPUNIT_ASSERT(sv[3].maxVal=="NA");
	CPPUNIT_ASSERT(sv[3].descript=="Mosaicism confidence running median window size");
	CPPUNIT_ASSERT(sv[4].name=="run-y-chromosome");
	CPPUNIT_ASSERT(sv[4].type==4);
	CPPUNIT_ASSERT(sv[4].value=="true");
	CPPUNIT_ASSERT(sv[4].defaultVal=="true");
	CPPUNIT_ASSERT(sv[4].minVal=="NA");
	CPPUNIT_ASSERT(sv[4].maxVal=="NA");
	CPPUNIT_ASSERT(sv[4].descript=="Run Mosaicism analysis on Y Chromosome");

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodMosaicism::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="mosaicism.gains-boundries=0.08764945,0.15380349,0.21465931,0.27100300.losses-boundries=-0.08293345,-0.17551812,-0.28048196,-0.40165383.marker-bandwidth=6000.confidence-window=251.run-y-chromosome=true");
    CPPUNIT_ASSERT(sd.getDocName()=="mosaicism");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber Mosaicism");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==5);
	CPPUNIT_ASSERT(v[0].asString()=="0.08764945,0.15380349,0.21465931,0.27100300");
	CPPUNIT_ASSERT(v[1].asString()=="-0.08293345,-0.17551812,-0.28048196,-0.40165383");
	CPPUNIT_ASSERT(v[2].asString()=="6000");
	CPPUNIT_ASSERT(v[3].asString()=="251");
	CPPUNIT_ASSERT(v[4].asString()=="true");
	CPPUNIT_ASSERT(sd.getDocOption("gains-boundries").asString()=="0.08764945,0.15380349,0.21465931,0.27100300");
	CPPUNIT_ASSERT(sd.getDocOption("losses-boundries").asString()=="-0.08293345,-0.17551812,-0.28048196,-0.40165383");
	CPPUNIT_ASSERT(sd.getDocOption("marker-bandwidth").asString()=="6000");
	CPPUNIT_ASSERT(sd.getDocOption("confidence-window").asString()=="251");
	CPPUNIT_ASSERT(sd.getDocOption("run-y-chromosome").asString()=="true");
	
}

void CNAnalysisMethodMosaicismTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodMosaicismTest::newObjectTest****");

     //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodMosaicism *cn2=(CNAnalysisMethodMosaicism*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("mosaicism");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==5);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5a = obj1->at(4);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-gains-boundries");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param1a.GetValueAscii()=="0.08764945,0.15380349,0.21465931,0.27100300");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-losses-boundries");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param2a.GetValueAscii()=="-0.08293345,-0.17551812,-0.28048196,-0.40165383");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-marker-bandwidth");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==6000);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-confidence-window");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param4a.GetValueInt32()==251);
	   
	delete cn2;

    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params;
	params["gains-boundries"]="0.1,0.2,0.3,0.4";
	params["losses-boundries"]="0.2,-2.2,-3.3,-4.4";
	params["marker-bandwidth"]="2";
	params["confidence-window"]="41";
	
	SelfCreate *sc=CNAnalysisMethodMosaicism::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==5);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(3);

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-gains-boundries");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param1.GetValueAscii()=="0.1,0.2,0.3,0.4"); 
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-losses-boundries");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param2.GetValueAscii()=="0.2,-2.2,-3.3,-4.4");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-marker-bandwidth");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3.GetValueInt32()==2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-confidence-window");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param4.GetValueInt32()==41);

	//Exceptions 
    //losses-boundries must conform to gains-boundries.
    std::map<std::string,std::string> params1;
	params1["gains-boundries"]="0,1,7";
	params1["losses-boundries"]="2.5,-0.62";
	NEGATIVE_TEST(CNAnalysisMethodMosaicism::newObject(params1),Except);

	//gains-boundries and losses-boundries must contain four floating point numbers.
    std::map<std::string,std::string> params2;
	params2["gains-boundries"]="0,1,7";
	params2["losses-boundries"]="0.3,-0.2,-0.3";
	NEGATIVE_TEST(CNAnalysisMethodMosaicism::newObject(params2),Except);

	//gains-boundries parameters must ordered increasing.
    std::map<std::string,std::string> params3;
	params3["gains-boundries"]="0.3,-0.2,0.3,0.2";
	NEGATIVE_TEST(CNAnalysisMethodMosaicism::newObject(params3),Except);

	//losses-boundries parameters must ordered decreasing.
    std::map<std::string,std::string> params4;
	params4["losses-boundries"]="-2.5,-0.62,0.0,2.49";
	NEGATIVE_TEST(CNAnalysisMethodMosaicism::newObject(params4),Except); 
			
	 //always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodMosaicism::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="mosaicism.gains-boundries=0.08764945,0.15380349,0.21465931,0.27100300.losses-boundries=-0.08293345,-0.17551812,-0.28048196,-0.40165383.marker-bandwidth=6000.confidence-window=251.run-y-chromosome=true");
    delete sc;

	params.clear();
	params["gains-boundries"]="0.1,0.2,0.3,0.4";
	params["losses-boundries"]="0.2,-2.2,-3.3,-4.4";
	params["marker-bandwidth"]="2";
	params["confidence-window"]="41";
	POSITIVE_TEST(CNAnalysisMethodMosaicism::newObject(params));
	//SelfCreate::setValue() - '0' is not a valid value for parameter: 'marker-bandwidth'. The specified range is 1 to NA
    params["marker-bandwidth"]="0";
	NEGATIVE_TEST(CNAnalysisMethodMosaicism::newObject(params),Except);
	params["marker-bandwidth"]="100000";
	POSITIVE_TEST(CNAnalysisMethodMosaicism::newObject(params));
	//SelfCreate::setValue() - '30' is not a valid value for parameter: 'confidence-window'. The specified range is 31 to NA
	params["confidence-window"]="30";
	NEGATIVE_TEST(CNAnalysisMethodMosaicism::newObject(params),Except);
	params["confidence-window"]="100000";
	POSITIVE_TEST(CNAnalysisMethodMosaicism::newObject(params));
	

	

}	

void CNAnalysisMethodMosaicismTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodMosaicismTest::runTest****");
    
	CNAnalysisMethodMosaicism cnrfCNam;
    //CNAnalysisMethod mosaicism is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
