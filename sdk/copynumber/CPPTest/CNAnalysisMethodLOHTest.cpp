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
#include "copynumber/CNAnalysisMethodLOH.h"
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
 * @class CNAnalysisMethodLOHTest
 * @brief cppunit class for testing CNAnalysisMethodLOH functions.
 * last change by vliber on 01/28/09
 */

class CNAnalysisMethodLOHTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodLOHTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodLOHTest );

void CNAnalysisMethodLOHTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodLOHTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodLOH cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="loh");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="CopyNumber LOH");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="loh");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodLOH::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==7);
	CPPUNIT_ASSERT(sv[0].name=="lohCN_errorrate");
	CPPUNIT_ASSERT(sv[0].type==1);
	CPPUNIT_ASSERT(sv[0].value=="0.05");
	CPPUNIT_ASSERT(sv[0].defaultVal=="0.05");
	CPPUNIT_ASSERT(sv[0].minVal=="NA");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="LOH CN Error Rate");
	CPPUNIT_ASSERT(sv[1].name=="lohCN_beta");
	CPPUNIT_ASSERT(sv[1].type==1);
	CPPUNIT_ASSERT(sv[1].value=="0.001");
	CPPUNIT_ASSERT(sv[1].defaultVal=="0.001");
	CPPUNIT_ASSERT(sv[1].minVal=="NA");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="LOH CN Beta");
	CPPUNIT_ASSERT(sv[2].name=="lohCN_alpha");
	CPPUNIT_ASSERT(sv[2].type==1);
	CPPUNIT_ASSERT(sv[2].value=="0.01");
	CPPUNIT_ASSERT(sv[2].defaultVal=="0.01");
	CPPUNIT_ASSERT(sv[2].minVal=="NA");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="LOH CN Alpha");
	CPPUNIT_ASSERT(sv[3].name=="lohCN_separation");
	CPPUNIT_ASSERT(sv[3].type==3);
	CPPUNIT_ASSERT(sv[3].value=="1000000");
	CPPUNIT_ASSERT(sv[3].defaultVal=="1000000");
	CPPUNIT_ASSERT(sv[3].minVal=="NA");
	CPPUNIT_ASSERT(sv[3].maxVal=="NA");
	CPPUNIT_ASSERT(sv[3].descript=="LOH CN Separation");
	CPPUNIT_ASSERT(sv[4].name=="lohCN_nMinMarkers");
	CPPUNIT_ASSERT(sv[4].type==3);
	CPPUNIT_ASSERT(sv[4].value=="10");
	CPPUNIT_ASSERT(sv[4].defaultVal=="10");
	CPPUNIT_ASSERT(sv[4].minVal=="NA");
	CPPUNIT_ASSERT(sv[4].maxVal=="NA");
	CPPUNIT_ASSERT(sv[4].descript=="LOH CN Minimum Marker Count");
	CPPUNIT_ASSERT(sv[5].name=="lohCN_NoCallThreshold");
	CPPUNIT_ASSERT(sv[5].type==1);
	CPPUNIT_ASSERT(sv[5].value=="0.05");
	CPPUNIT_ASSERT(sv[5].defaultVal=="0.05");
	CPPUNIT_ASSERT(sv[5].minVal=="NA");
	CPPUNIT_ASSERT(sv[5].maxVal=="NA");
	CPPUNIT_ASSERT(sv[5].descript=="LOH CN No Call Threshold");
	CPPUNIT_ASSERT(sv[6].name=="lohCN_minGenomicSpan");
	CPPUNIT_ASSERT(sv[6].type==3);
	CPPUNIT_ASSERT(sv[6].value=="1000000");
	CPPUNIT_ASSERT(sv[6].defaultVal=="1000000");
	CPPUNIT_ASSERT(sv[6].minVal=="NA");
	CPPUNIT_ASSERT(sv[6].maxVal=="NA");
	CPPUNIT_ASSERT(sv[6].descript=="LOH CN Minimum Genomic Span");
	//explainSelf()
	SelfDoc sd=CNAnalysisMethodLOH::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="loh.lohCN_errorrate=0.05.lohCN_beta=0.001.lohCN_alpha=0.01.lohCN_separation=1000000.lohCN_nMinMarkers=10.lohCN_NoCallThreshold=0.05.lohCN_minGenomicSpan=1000000");
    CPPUNIT_ASSERT(sd.getDocName()=="loh");
	CPPUNIT_ASSERT(sd.getDocDescription()=="CopyNumber LOH");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==7);
	CPPUNIT_ASSERT(v[0].asString()=="0.05");
	CPPUNIT_ASSERT(v[1].asString()=="0.001");
	CPPUNIT_ASSERT(v[2].asString()=="0.01");
	CPPUNIT_ASSERT(v[3].asString()=="1000000");
	CPPUNIT_ASSERT(v[4].asString()=="10");
	CPPUNIT_ASSERT(v[5].asString()=="0.05");
	CPPUNIT_ASSERT(v[6].asString()=="1000000");
	CPPUNIT_ASSERT(sd.getDocOption("lohCN_errorrate").asString()=="0.05");
	CPPUNIT_ASSERT(sd.getDocOption("lohCN_beta").asString()=="0.001");
	CPPUNIT_ASSERT(sd.getDocOption("lohCN_alpha").asString()=="0.01");
	CPPUNIT_ASSERT(sd.getDocOption("lohCN_separation").asString()=="1000000");
	CPPUNIT_ASSERT(sd.getDocOption("lohCN_nMinMarkers").asString()=="10");
	CPPUNIT_ASSERT(sd.getDocOption("lohCN_NoCallThreshold").asString()=="0.05");
	CPPUNIT_ASSERT(sd.getDocOption("lohCN_minGenomicSpan").asString()=="1000000");
}

void CNAnalysisMethodLOHTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodLOHTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodLOH *cn2=(CNAnalysisMethodLOH*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("loh");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==7);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5a = obj1->at(4);
	affymetrix_calvin_parameter::ParameterNameValueType param6a = obj1->at(5);
	affymetrix_calvin_parameter::ParameterNameValueType param7a = obj1->at(6);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-lohCN_errorrate");
	CPPUNIT_ASSERT(param1a.GetParameterType()==6);
    CPPUNIT_ASSERT(param1a.GetValueFloat()==0.05f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-lohCN_beta");
	CPPUNIT_ASSERT(param2a.GetParameterType()==6);
	CPPUNIT_ASSERT(param2a.GetValueFloat()==0.001f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-lohCN_alpha");
	CPPUNIT_ASSERT(param3a.GetParameterType()==6);
	CPPUNIT_ASSERT(param3a.GetValueFloat()==0.01f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-lohCN_separation");
	CPPUNIT_ASSERT(param4a.GetParameterType()==4);
	CPPUNIT_ASSERT(param4a.GetValueInt32()==1000000);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5a.GetName())=="affymetrix-algorithm-param-lohCN_nMinMarkers");
	CPPUNIT_ASSERT(param5a.GetParameterType()==4);
	CPPUNIT_ASSERT(param5a.GetValueInt32()==10);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6a.GetName())=="affymetrix-algorithm-param-lohCN_NoCallThreshold");
	CPPUNIT_ASSERT(param6a.GetParameterType()==6);
	CPPUNIT_ASSERT(param6a.GetValueFloat()==0.05f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7a.GetName())=="affymetrix-algorithm-param-lohCN_minGenomicSpan");
	CPPUNIT_ASSERT(param7a.GetParameterType()==4);
	CPPUNIT_ASSERT(param7a.GetValueInt32()==1000000);
	delete cn2;


   //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["lohCN_errorrate"]="0.111";
	params["lohCN_beta"]="0.222";
	params["lohCN_alpha"]="0.100";
	params["lohCN_separation"]="200";
	params["lohCN_nMinMarkers"]="100";
	params["lohCN_NoCallThreshold"]="0.100";
	params["lohCN_minGenomicSpan"]="1000";
	CNAnalysisMethod::getParams()->clear();
	SelfCreate *sc=CNAnalysisMethodLOH::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==7);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5 = obj->at(4);
	affymetrix_calvin_parameter::ParameterNameValueType param6 = obj->at(5);
	affymetrix_calvin_parameter::ParameterNameValueType param7 = obj->at(6);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-lohCN_errorrate");
	CPPUNIT_ASSERT(param1.GetParameterType()==6);
    CPPUNIT_ASSERT(param1.GetValueFloat()==0.111f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-lohCN_beta");
	CPPUNIT_ASSERT(param2.GetParameterType()==6);
	CPPUNIT_ASSERT(param2.GetValueFloat()==0.222f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-lohCN_alpha");
	CPPUNIT_ASSERT(param3.GetParameterType()==6);
	CPPUNIT_ASSERT(param3.GetValueFloat()==0.100f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-lohCN_separation");
	CPPUNIT_ASSERT(param4.GetParameterType()==4);
	CPPUNIT_ASSERT(param4.GetValueInt32()==200);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5.GetName())=="affymetrix-algorithm-param-lohCN_nMinMarkers");
	CPPUNIT_ASSERT(param5.GetParameterType()==4);
	CPPUNIT_ASSERT(param5.GetValueInt32()==100);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6.GetName())=="affymetrix-algorithm-param-lohCN_NoCallThreshold");
	CPPUNIT_ASSERT(param6.GetParameterType()==6);
	CPPUNIT_ASSERT(param6.GetValueFloat()==0.100f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7.GetName())=="affymetrix-algorithm-param-lohCN_minGenomicSpan");
	CPPUNIT_ASSERT(param7.GetParameterType()==4);
	CPPUNIT_ASSERT(param7.GetValueInt32()==1000);
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodLOH::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="loh.lohCN_errorrate=0.05.lohCN_beta=0.001.lohCN_alpha=0.01.lohCN_separation=1000000.lohCN_nMinMarkers=10.lohCN_NoCallThreshold=0.05.lohCN_minGenomicSpan=1000000");
    delete sc;	
}	

void CNAnalysisMethodLOHTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodLOHTest::runTest****");

    
	CNAnalysisMethodLOH cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod loh is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
