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

#include "copynumber/CNAnalysisMethodCN.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNLog2RatioAdjustmentMethodWaveCorrection.h"
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
 * @class CNLog2RatioAdjustmentMethodWaveCorrectionTest
 * @brief cppunit class for testing CNLog2RatioAdjustmentMethodWaveCorrection functions.
 * last change by vliber on 03/31/09
 */

class CNLog2RatioAdjustmentMethodWaveCorrectionTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNLog2RatioAdjustmentMethodWaveCorrectionTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest);
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNLog2RatioAdjustmentMethodWaveCorrectionTest);

void CNLog2RatioAdjustmentMethodWaveCorrectionTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNLog2RatioAdjustmentMethodWaveCorrectionTest::functionsTest****");
	// constructor && m_vParams()
	CNLog2RatioAdjustmentMethodWaveCorrection cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="wave-correction-log2ratio-adjustment-method");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber Wave Correction Log2Ratio Adjustment");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="wave-correction-log2ratio-adjustment-method");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNLog2RatioAdjustmentMethodWaveCorrection::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==4);
	CPPUNIT_ASSERT(sv[0].name=="bandwidth");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="101");
	CPPUNIT_ASSERT(sv[0].defaultVal=="101");
	CPPUNIT_ASSERT(sv[0].minVal=="NA");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="Wave Correction bandwidth.");
	CPPUNIT_ASSERT(sv[1].name=="bin-count");
	CPPUNIT_ASSERT(sv[1].type==3);
	CPPUNIT_ASSERT(sv[1].value=="25");
	CPPUNIT_ASSERT(sv[1].defaultVal=="25");
	CPPUNIT_ASSERT(sv[1].minVal=="NA");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="Wave Correction bin count.");
	CPPUNIT_ASSERT(sv[2].name=="wave-count");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="-1");
	CPPUNIT_ASSERT(sv[2].defaultVal=="-1");
	CPPUNIT_ASSERT(sv[2].minVal=="-1");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="Wave Correction count.");
	CPPUNIT_ASSERT(sv[3].name=="wave-smooth");
	CPPUNIT_ASSERT(sv[3].type==4);
	CPPUNIT_ASSERT(sv[3].value=="true");
	CPPUNIT_ASSERT(sv[3].defaultVal=="true");
	CPPUNIT_ASSERT(sv[3].minVal=="NA");
	CPPUNIT_ASSERT(sv[3].maxVal=="NA");
	CPPUNIT_ASSERT(sv[3].descript=="Turn the non-parametric smoothing on or off.");

	

	//explainSelf()
	SelfDoc sd=CNLog2RatioAdjustmentMethodWaveCorrection::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=-1.wave-smooth=true");
    CPPUNIT_ASSERT(sd.getDocName()=="wave-correction-log2ratio-adjustment-method");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber Wave Correction Log2Ratio Adjustment");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==4);
	CPPUNIT_ASSERT(v[0].asString()=="101");
	CPPUNIT_ASSERT(v[1].asString()=="25");
	CPPUNIT_ASSERT(v[2].asString()=="-1");
    CPPUNIT_ASSERT(v[3].asString()=="true");
	CPPUNIT_ASSERT(sd.getDocOption("bandwidth").asString()=="101");
	CPPUNIT_ASSERT(sd.getDocOption("bin-count").asString()=="25");
	CPPUNIT_ASSERT(sd.getDocOption("wave-count").asString()=="-1");
    CPPUNIT_ASSERT(sd.getDocOption("wave-smooth").asString()=="true");
}

void CNLog2RatioAdjustmentMethodWaveCorrectionTest::newObjectTest()
{
	
	Verbose::out(1, "****CNLog2RatioAdjustmentMethodWaveCorrectionTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNLog2RatioAdjustmentMethodWaveCorrection *cn2=(CNLog2RatioAdjustmentMethodWaveCorrection*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("wave-correction-log2ratio-adjustment-method");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==4);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);	
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
    affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(3);	
		
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-bandwidth");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1a.GetValueInt32()==101);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-bin-count");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param2a.GetValueInt32()==25);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-wave-count");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==-1);

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-wave-smooth");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param4a.GetValueInt8()==1);

	delete cn2;

    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params1;
	params1["bandwidth"]="200";
	params1["bin-count"]="20";
	params1["wave-count"]="2";

	SelfCreate *sc=CNLog2RatioAdjustmentMethodWaveCorrection::newObject(params1);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==4);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
    affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(3);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-bandwidth");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1.GetValueInt32()==200);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-bin-count");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param2.GetValueInt32()==20);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-wave-count");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3.GetValueInt32()==2);

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-wave-smooth");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param4.GetValueInt8()==1);
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNLog2RatioAdjustmentMethodWaveCorrection::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=-1.wave-smooth=true");
    delete sc;
	
}	

void CNLog2RatioAdjustmentMethodWaveCorrectionTest::runTest()
{
    Verbose::out(1, "****CNLog2RatioAdjustmentMethodWaveCorrectionTest::runTest****");
    
	CNLog2RatioAdjustmentMethodWaveCorrection cnrfCNam;
    //CNAnalysisMethod wave-correction-log2ratio-adjustment-method is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
