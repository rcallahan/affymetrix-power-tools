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
#include "copynumber/CNIntensityAdjustmentMethodPDNN.h"
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
 * @class CNIntensityAdjustmentMethodPDNNTest
 * @brief cppunit class for testing CNIntensityAdjustmentMethodPDNN functions.
 * last change by vliber on 03/24/09
 */

class CNIntensityAdjustmentMethodPDNNTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNIntensityAdjustmentMethodPDNNTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest);
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNIntensityAdjustmentMethodPDNNTest);

void CNIntensityAdjustmentMethodPDNNTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNIntensityAdjustmentMethodPDNNTest::functionsTest****");
	// constructor && m_vParams()
	CNIntensityAdjustmentMethodPDNN cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="pdnn-intensity-adjustment-method");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber PDNN Intensity Adjustment");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="pdnn-intensity-adjustment-method");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNIntensityAdjustmentMethodPDNN::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==3);
	CPPUNIT_ASSERT(sv[0].name=="predicted-intensity-bin-count");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="20");
	CPPUNIT_ASSERT(sv[0].defaultVal=="20");
	CPPUNIT_ASSERT(sv[0].minVal=="1");
	CPPUNIT_ASSERT(sv[0].maxVal=="126");
	CPPUNIT_ASSERT(sv[0].descript=="Predicted Intensity bin count.");
	CPPUNIT_ASSERT(sv[1].name=="gc-bin-count");
	CPPUNIT_ASSERT(sv[1].type==3);
	CPPUNIT_ASSERT(sv[1].value=="20");
	CPPUNIT_ASSERT(sv[1].defaultVal=="20");
	CPPUNIT_ASSERT(sv[1].minVal=="1");
	CPPUNIT_ASSERT(sv[1].maxVal=="126");
	CPPUNIT_ASSERT(sv[1].descript=="Predicted Intensity bin count.");
	CPPUNIT_ASSERT(sv[2].name=="residual-trim");
	CPPUNIT_ASSERT(sv[2].type==1);
	CPPUNIT_ASSERT(sv[2].value=="2.0");
	CPPUNIT_ASSERT(sv[2].defaultVal=="2.0");
	CPPUNIT_ASSERT(sv[2].minVal=="NA");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="Residual Trim value.");
		

	//explainSelf()
	SelfDoc sd=CNIntensityAdjustmentMethodPDNN::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0");
    CPPUNIT_ASSERT(sd.getDocName()=="pdnn-intensity-adjustment-method");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber PDNN Intensity Adjustment");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==3);
	CPPUNIT_ASSERT(v[0].asString()=="20");
	CPPUNIT_ASSERT(v[1].asString()=="20");
	CPPUNIT_ASSERT(v[2].asString()=="2.0");
	CPPUNIT_ASSERT(sd.getDocOption("predicted-intensity-bin-count").asString()=="20");
	CPPUNIT_ASSERT(sd.getDocOption("gc-bin-count").asString()=="20");
	CPPUNIT_ASSERT(sd.getDocOption("residual-trim").asString()=="2.0");
}
//vliber inconsistency in data type (double-int)
void CNIntensityAdjustmentMethodPDNNTest::newObjectTest()
{
	
	Verbose::out(1, "****CNIntensityAdjustmentMethodPDNNTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNIntensityAdjustmentMethodPDNN *cn2=(CNIntensityAdjustmentMethodPDNN*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("pdnn-intensity-adjustment-method");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
	
		
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-predicted-intensity-bin-count");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1a.GetValueInt32()==20);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-gc-bin-count");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param2a.GetValueInt32()==20);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-residual-trim");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param3a.GetValueFloat()==2.0);
	delete cn2;
   
    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params1;
	params1["predicted-intensity-bin-count"]="50";
	params1["gc-bin-count"]="45";
	params1["residual-trim"]="1.50";
	SelfCreate *sc=CNIntensityAdjustmentMethodPDNN::newObject(params1);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-predicted-intensity-bin-count");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1.GetValueInt32()==50);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-gc-bin-count");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param2.GetValueInt32()==45);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-residual-trim");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param3.GetValueFloat()==1.50);
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNIntensityAdjustmentMethodPDNN::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0");
    delete sc;	
}
void CNIntensityAdjustmentMethodPDNNTest::runTest()
{
   Verbose::out(1, "****CNIntensityAdjustmentMethodPDNNTest::runTest****");
    
	CNIntensityAdjustmentMethodPDNN cnrfCNam;
    //FATAL ERROR: CNIntensityAdjustmentMethodPDNN cn-state is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
