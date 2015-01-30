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
#include "copynumber/CNIntensityAdjustmentMethodHighPassFilter.h"
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
 * @class CNIntensityAdjustmentMethodHighPassFilterTest
 * @brief cppunit class for testing CNIntensityAdjustmentMethodHighPassFilter functions.
 * last change by vliber on 03/24/09
 */

class CNIntensityAdjustmentMethodHighPassFilterTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNIntensityAdjustmentMethodHighPassFilterTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNIntensityAdjustmentMethodHighPassFilterTest);

void CNIntensityAdjustmentMethodHighPassFilterTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNIntensityAdjustmentMethodHighPassFilterTest::functionsTest****");
	// constructor && m_vParams()
	CNIntensityAdjustmentMethodHighPassFilter cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="high-pass-filter-intensity-adjustment-method");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="CopyNumber HighPassFilter");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="high-pass-filter-intensity-adjustment-method");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNIntensityAdjustmentMethodHighPassFilter::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==8);
	CPPUNIT_ASSERT(sv[0].name=="data-block-rows");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="320");
	CPPUNIT_ASSERT(sv[0].defaultVal=="320");
	CPPUNIT_ASSERT(sv[0].minVal=="NA");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="Data block rows.");
	CPPUNIT_ASSERT(sv[1].name=="data-block-cols");
	CPPUNIT_ASSERT(sv[1].type==3);
	CPPUNIT_ASSERT(sv[1].value=="2015");
	CPPUNIT_ASSERT(sv[1].defaultVal=="2015");
	CPPUNIT_ASSERT(sv[1].minVal=="NA");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="Data block columns.");
	CPPUNIT_ASSERT(sv[2].name=="mini-block-rows");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="8");
	CPPUNIT_ASSERT(sv[2].defaultVal=="8");
	CPPUNIT_ASSERT(sv[2].minVal=="NA");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="Mini block rows.");
	CPPUNIT_ASSERT(sv[3].name=="mini-block-cols");
	CPPUNIT_ASSERT(sv[3].type==3);
	CPPUNIT_ASSERT(sv[3].value=="8");
	CPPUNIT_ASSERT(sv[3].defaultVal=="8");
	CPPUNIT_ASSERT(sv[3].minVal=="NA");
	CPPUNIT_ASSERT(sv[3].maxVal=="NA");
	CPPUNIT_ASSERT(sv[3].descript=="Mini block columns.");
	CPPUNIT_ASSERT(sv[4].name=="global-smooth-weight");
	CPPUNIT_ASSERT(sv[4].type==1);
	CPPUNIT_ASSERT(sv[4].value=="256.0");
	CPPUNIT_ASSERT(sv[4].defaultVal=="256.0");
	CPPUNIT_ASSERT(sv[4].minVal=="NA");
	CPPUNIT_ASSERT(sv[4].maxVal=="NA");
	CPPUNIT_ASSERT(sv[4].descript=="Global smooth weight.");
	CPPUNIT_ASSERT(sv[5].name=="local-smooth-weight");
	CPPUNIT_ASSERT(sv[5].type==1);
	CPPUNIT_ASSERT(sv[5].value=="64.0");
	CPPUNIT_ASSERT(sv[5].defaultVal=="64.0");
	CPPUNIT_ASSERT(sv[5].minVal=="NA");
	CPPUNIT_ASSERT(sv[5].maxVal=="NA");
	CPPUNIT_ASSERT(sv[5].descript=="Local smooth weight.");
	CPPUNIT_ASSERT(sv[6].name=="converged");
	CPPUNIT_ASSERT(sv[6].type==1);
	CPPUNIT_ASSERT(sv[6].value=="0.0001");
	CPPUNIT_ASSERT(sv[6].defaultVal=="0.0001");
	CPPUNIT_ASSERT(sv[6].minVal=="NA");
	CPPUNIT_ASSERT(sv[6].maxVal=="NA");
	CPPUNIT_ASSERT(sv[6].descript=="Converged.");
		

	//explainSelf()
	SelfDoc sd=CNIntensityAdjustmentMethodHighPassFilter::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.0.local-smooth-weight=64.0.converged=0.0001.use-single-block=false");
    CPPUNIT_ASSERT(sd.getDocName()=="high-pass-filter-intensity-adjustment-method");
	CPPUNIT_ASSERT(sd.getDocDescription()=="CopyNumber HighPassFilter");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==8);
	CPPUNIT_ASSERT(v[0].asString()=="320");
	CPPUNIT_ASSERT(v[1].asString()=="2015");
	CPPUNIT_ASSERT(v[2].asString()=="8");
	CPPUNIT_ASSERT(v[3].asString()=="8");
	CPPUNIT_ASSERT(v[4].asString()=="256.0");
	CPPUNIT_ASSERT(v[5].asString()=="64.0");
	CPPUNIT_ASSERT(v[6].asString()=="0.0001");
	CPPUNIT_ASSERT(sd.getDocOption("data-block-rows").asString()=="320");
	CPPUNIT_ASSERT(sd.getDocOption("data-block-cols").asString()=="2015");
	CPPUNIT_ASSERT(sd.getDocOption("mini-block-rows").asString()=="8");
	CPPUNIT_ASSERT(sd.getDocOption("mini-block-cols").asString()=="8");
	CPPUNIT_ASSERT(sd.getDocOption("global-smooth-weight").asString()=="256.0");
	CPPUNIT_ASSERT(sd.getDocOption("local-smooth-weight").asString()=="64.0");
	CPPUNIT_ASSERT(sd.getDocOption("converged").asString()=="0.0001");
}
void CNIntensityAdjustmentMethodHighPassFilterTest::newObjectTest()
{
	
	Verbose::out(1, "****CNIntensityAdjustmentMethodHighPassFilterTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNIntensityAdjustmentMethodHighPassFilter *cn2=(CNIntensityAdjustmentMethodHighPassFilter*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("high-pass-filter-intensity-adjustment-method");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==8);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5a = obj1->at(4);
	affymetrix_calvin_parameter::ParameterNameValueType param6a = obj1->at(5);
	affymetrix_calvin_parameter::ParameterNameValueType param7a = obj1->at(6);
			
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-data-block-rows");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1a.GetValueInt32()==320);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-data-block-cols");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param2a.GetValueInt32()==2015);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-mini-block-rows");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==8);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-mini-block-cols");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param4a.GetValueInt32()==8);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5a.GetName())=="affymetrix-algorithm-param-global-smooth-weight");
	CPPUNIT_ASSERT(param5a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5a.GetValueFloat()==256.0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6a.GetName())=="affymetrix-algorithm-param-local-smooth-weight");
	CPPUNIT_ASSERT(param6a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param6a.GetValueFloat()==64.0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7a.GetName())=="affymetrix-algorithm-param-converged");
	CPPUNIT_ASSERT(param7a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param7a.GetValueFloat()==0.0001f);
	delete cn2;
   
    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params1;
	params1["data-block-rows"]="350";
	params1["data-block-cols"]="450";
	params1["mini-block-rows"]="6";
	params1["mini-block-cols"]="4";
	params1["global-smooth-weight"]="10.50";
	params1["local-smooth-weight"]="12.50";
	params1["converged"]="0.18";
	SelfCreate *sc=CNIntensityAdjustmentMethodHighPassFilter::newObject(params1);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==8);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5 = obj->at(4);
	affymetrix_calvin_parameter::ParameterNameValueType param6 = obj->at(5);
	affymetrix_calvin_parameter::ParameterNameValueType param7 = obj->at(6);
			
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-data-block-rows");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1.GetValueInt32()==350);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-data-block-cols");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param2.GetValueInt32()==450);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-mini-block-rows");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3.GetValueInt32()==6);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-mini-block-cols");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param4.GetValueInt32()==4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5.GetName())=="affymetrix-algorithm-param-global-smooth-weight");
	CPPUNIT_ASSERT(param5.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5.GetValueFloat()==10.50);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6.GetName())=="affymetrix-algorithm-param-local-smooth-weight");
	CPPUNIT_ASSERT(param6.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param6.GetValueFloat()==12.50);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7.GetName())=="affymetrix-algorithm-param-converged");
	CPPUNIT_ASSERT(param7.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param7.GetValueFloat()==0.18f);

	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNIntensityAdjustmentMethodHighPassFilter::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.0.local-smooth-weight=64.0.converged=0.0001.use-single-block=false");
    delete sc;	
}

