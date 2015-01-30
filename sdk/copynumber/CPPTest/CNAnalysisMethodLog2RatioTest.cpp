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
#include "copynumber/CNAnalysisMethodLog2Ratio.h"
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
 * @class CNAnalysisMethodLog2RatioTest
 * @brief cppunit class for testing CNAnalysisMethodLog2Ratio functions.
 * last change by vliber on 03/23/09
 */

class CNAnalysisMethodLog2RatioTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodLog2RatioTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest);
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodLog2RatioTest );

void CNAnalysisMethodLog2RatioTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodLog2RatioTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodLog2Ratio cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="log2-ratio");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber Log2Ratio");
		
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodLog2Ratio::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==3);
	CPPUNIT_ASSERT(sv[0].name=="gc-correction");
	CPPUNIT_ASSERT(sv[0].type==4);
	CPPUNIT_ASSERT(sv[0].value=="true");
	CPPUNIT_ASSERT(sv[0].defaultVal=="true");
	CPPUNIT_ASSERT(sv[0].minVal=="NA");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="Log2Ratio GC Correction");
	CPPUNIT_ASSERT(sv[1].name=="median-autosome-median-normalization");
	CPPUNIT_ASSERT(sv[1].type==4);
	CPPUNIT_ASSERT(sv[1].value=="true");
	CPPUNIT_ASSERT(sv[1].defaultVal=="true");
	CPPUNIT_ASSERT(sv[1].minVal=="NA");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="Log2Ratio Median Autosmome Median Normalization");
	CPPUNIT_ASSERT(sv[2].name=="median-smooth-marker-count");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="5");
	CPPUNIT_ASSERT(sv[2].defaultVal=="5");
	CPPUNIT_ASSERT(sv[2].minVal=="3");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="Log2Ratio Median Smooth marker count.");
	
	//explainSelf()
	SelfDoc sd=CNAnalysisMethodLog2Ratio::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="log2-ratio.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5");
    CPPUNIT_ASSERT(sd.getDocName()=="log2-ratio");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber Log2Ratio");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==3);
	CPPUNIT_ASSERT(v[0].asString()=="true");
	CPPUNIT_ASSERT(v[1].asString()=="true");
	CPPUNIT_ASSERT(v[2].asString()=="5");
	CPPUNIT_ASSERT(sd.getDocOption("gc-correction").asString()=="true");
	CPPUNIT_ASSERT(sd.getDocOption("median-autosome-median-normalization").asString()=="true");	
	CPPUNIT_ASSERT(sd.getDocOption("median-smooth-marker-count").asString()=="5");	

    
}

void CNAnalysisMethodLog2RatioTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodLog2RatioTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodLog2Ratio *cn2=(CNAnalysisMethodLog2Ratio*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("log2-ratio");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-gc-correction");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.ToString())=="1");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-median-autosome-median-normalization");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.ToString())=="1");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-median-smooth-marker-count");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.ToString())=="5");
	
    delete cn2;

    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params;
	params["gc-correction"]="false";
	params["median-autosome-median-normalization"]="true";
	params["median-smooth-marker-count"]="14";
	SelfCreate *sc=CNAnalysisMethodLog2Ratio::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-gc-correction");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.ToString())=="0");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-median-autosome-median-normalization");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.ToString())=="1");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-median-smooth-marker-count");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.ToString())=="14");
    
	//positive tests with different option values
    POSITIVE_TEST(CNAnalysisMethodLog2Ratio::newObject(params));
    params["gc-correction"]="true";
    POSITIVE_TEST(CNAnalysisMethodLog2Ratio::newObject(params));
	params["gc-correction"]="false";
    POSITIVE_TEST(CNAnalysisMethodLog2Ratio::newObject(params));
	params["median-autosome-median-normalization"]="true";
    POSITIVE_TEST(CNAnalysisMethodLog2Ratio::newObject(params));
	params["median-autosome-median-normalization"]="false";
    POSITIVE_TEST(CNAnalysisMethodLog2Ratio::newObject(params));
    
	//Negative tests to check option value

	//FATAL ERROR: SelfCreate::setValue() - 'a' is not a valid value for parameter: 'gc-correction'. The specified range is NA to NA
	params["gc-correction"]="a";
	NEGATIVE_TEST(CNAnalysisMethodLog2Ratio::newObject(params),Except);

	//FATAL ERROR: SelfCreate::setValue() - 'a' is not a valid value for parameter: 'median-autosome-median-normalization'. The specified range is NA to NA
	params["gc-correction"]="false";
	params["median-autosome-median-normalization"]="a";
	NEGATIVE_TEST(CNAnalysisMethodLog2Ratio::newObject(params),Except);
 
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodLog2Ratio::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="log2-ratio.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5");
    delete sc;	
}

void CNAnalysisMethodLog2RatioTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodLog2RatioTest::runTest****");
    
	CNAnalysisMethodLog2Ratio cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod log2-ratio is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}

