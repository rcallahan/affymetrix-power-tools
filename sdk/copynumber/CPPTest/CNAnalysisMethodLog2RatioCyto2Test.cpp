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
#include "copynumber/CNAnalysisMethodLog2RatioCyto2.h"
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
 * @class CNAnalysisMethodLog2RatioCyto2Test
 * @brief cppunit class for testing CNAnalysisMethodLog2RatioCyto2 functions.
 * last change by vliber on 03/24/09
 */

class CNAnalysisMethodLog2RatioCyto2Test : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodLog2RatioCyto2Test);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest(); 
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodLog2RatioCyto2Test );

void CNAnalysisMethodLog2RatioCyto2Test::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodLog2RatioCyto2Test::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodLog2RatioCyto2 cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="log2-ratio-cyto2");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber Log2RatioCyto2");
		
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodLog2RatioCyto2::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==5);
	CPPUNIT_ASSERT(sv[0].name=="gc-correction");
	CPPUNIT_ASSERT(sv[0].type==4);
	CPPUNIT_ASSERT(sv[0].value=="true");
	CPPUNIT_ASSERT(sv[0].defaultVal=="true");
	CPPUNIT_ASSERT(sv[0].minVal=="NA");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="Log2RatioCyto2 GC Correction");
	CPPUNIT_ASSERT(sv[1].name=="median-autosome-median-normalization");
	CPPUNIT_ASSERT(sv[1].type==4);
	CPPUNIT_ASSERT(sv[1].value=="true");
	CPPUNIT_ASSERT(sv[1].defaultVal=="true");
	CPPUNIT_ASSERT(sv[1].minVal=="NA");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="Log2RatioCyto2 Median Autosmome Median Normalization");
	CPPUNIT_ASSERT(sv[2].name=="median-smooth-marker-count");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="5");
	CPPUNIT_ASSERT(sv[2].defaultVal=="5");
	CPPUNIT_ASSERT(sv[2].minVal=="3");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="Log2Ratio Median Smooth marker count.");
	CPPUNIT_ASSERT(sv[3].name=="trim-high");
	CPPUNIT_ASSERT(sv[3].type==1);
	CPPUNIT_ASSERT(sv[3].value=="2.0");
	CPPUNIT_ASSERT(sv[3].defaultVal=="2.0");
	CPPUNIT_ASSERT(sv[3].minVal=="0");
	CPPUNIT_ASSERT(sv[3].maxVal=="NA");
	CPPUNIT_ASSERT(sv[3].descript=="High trim value for Log2Ratios.");
	CPPUNIT_ASSERT(sv[4].name=="trim-low");
	CPPUNIT_ASSERT(sv[4].type==1);
	CPPUNIT_ASSERT(sv[4].value=="-2.5");
	CPPUNIT_ASSERT(sv[4].defaultVal=="-2.5");
	CPPUNIT_ASSERT(sv[4].minVal=="NA");
	CPPUNIT_ASSERT(sv[4].maxVal=="0");
	CPPUNIT_ASSERT(sv[4].descript=="Low trim value for Log2Ratios.");

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodLog2RatioCyto2::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5");
    CPPUNIT_ASSERT(sd.getDocName()=="log2-ratio-cyto2");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber Log2RatioCyto2");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==5);
	CPPUNIT_ASSERT(v[0].asString()=="true");
	CPPUNIT_ASSERT(v[1].asString()=="true");
	CPPUNIT_ASSERT(v[2].asString()=="5");
	CPPUNIT_ASSERT(v[3].asString()=="2.0");
	CPPUNIT_ASSERT(v[4].asString()=="-2.5");
	CPPUNIT_ASSERT(sd.getDocOption("gc-correction").asString()=="true");
	CPPUNIT_ASSERT(sd.getDocOption("median-autosome-median-normalization").asString()=="true");	
	CPPUNIT_ASSERT(sd.getDocOption("median-smooth-marker-count").asString()=="5");
	CPPUNIT_ASSERT(sd.getDocOption("trim-high").asString()=="2.0");
	CPPUNIT_ASSERT(sd.getDocOption("trim-low").asString()=="-2.5");
}

void CNAnalysisMethodLog2RatioCyto2Test::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodLog2RatioCyto2Test::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodLog2RatioCyto2 *cn2=(CNAnalysisMethodLog2RatioCyto2*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("log2-ratio-cyto2");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==5);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param5a = obj1->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param6a = obj1->at(4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-gc-correction");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param2a.GetValueInt8()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-median-autosome-median-normalization");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param3a.GetValueInt8()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-median-smooth-marker-count");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param4a.GetValueInt32()==5);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5a.GetName())=="affymetrix-algorithm-param-trim-high");
	CPPUNIT_ASSERT(param5a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5a.GetValueFloat()==2.0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6a.GetName())=="affymetrix-algorithm-param-trim-low");
	CPPUNIT_ASSERT(param6a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param6a.GetValueFloat()==-2.5);
	delete cn2;


    //create newObject with different parameters and test them 
	std::map<std::string,std::string> params;
	params["gc-correction"]="false";
	params["median-autosome-median-normalization"]="true";
	params["median-smooth-marker-count"]="5";
	params["trim-high"]="15";
	params["trim-low"]="-5.0";
	CNAnalysisMethod::getParams()->clear();
	SelfCreate *sc=CNAnalysisMethodLog2RatioCyto2::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==5);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param5 = obj->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param6 = obj->at(4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-gc-correction");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param2.GetValueInt8()==0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-median-autosome-median-normalization");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param3.GetValueInt8()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-median-smooth-marker-count");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param4.GetValueInt32()==5);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5.GetName())=="affymetrix-algorithm-param-trim-high");
	CPPUNIT_ASSERT(param5.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5.GetValueFloat()==15);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6.GetName())=="affymetrix-algorithm-param-trim-low");
	CPPUNIT_ASSERT(param6.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param6.GetValueFloat()==-5.0);
	delete sc;

	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodLog2RatioCyto2::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5");
    
    //test parameters
	//positive
	params["gc-correction"]="true";
	params["median-autosome-median-normalization"]="true";
	POSITIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params));
	params["gc-correction"]="false";
	params["median-autosome-median-normalization"]="false";
	POSITIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params));
	//negative
	//FATAL ERROR: SelfCreate::setValue() - '1' is not a valid value for parameter: 'gc-correction'. The specified range is NA to NA
    params["gc-correction"]="1";
	NEGATIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params),Except);
	//FATAL ERROR: SelfCreate::setValue() - 'a' is not a valid value for parameter: 'gc-correction'. The specified range is NA to NA
    params["gc-correction"]="a";
	NEGATIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params),Except);
	//FATAL ERROR: SelfCreate::setValue() - '0' is not a valid value for parameter: 'median-autosome-median-normalization'. The specified range is NA to NA
	params["gc-correction"]="true";
	params["median-autosome-median-normalization"]="0";
	NEGATIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params),Except);
	//FATAL ERROR: SelfCreate::setValue() - 'a' is not a valid value for parameter: 'median-autosome-median-normalization'. The specified range is NA to NA
	params["median-autosome-median-normalization"]="a";
	NEGATIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params),Except);
	//FATAL ERROR:Message: SelfCreate::setValue() - '2' is not a valid value for parameter: 'median-smooth-marker-count'. The specified range is 3 to NA
	params["median-autosome-median-normalization"]="true";
	params["median-smooth-marker-count"]="2";
	NEGATIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params),Except);
	params["median-smooth-marker-count"]="3";
	POSITIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params));
	params["median-smooth-marker-count"]="33333";
	POSITIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params));
    //FATAL ERROR:Message: SelfCreate::setValue() - '-15' is not a valid value for parameter: 'trim-high'. The specified range is 0 to NA
	params["trim-high"]="-15";
    NEGATIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params),Except);
	params["trim-high"]="0";
    POSITIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params));
	params["trim-high"]="111110";
    POSITIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params));
	//FATAL ERROR: Message: SelfCreate::setValue() - '0.0001' is not a valid value for parameter: 'trim-low'. The specified range is NA to 0
	params["trim-low"]="0.0001";
	NEGATIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params),Except);
	params["trim-low"]="0.0";
	POSITIVE_TEST(CNAnalysisMethodLog2RatioCyto2::newObject(params));
}	

