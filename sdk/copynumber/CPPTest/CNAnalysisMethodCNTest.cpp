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
 * @class CNAnalysisMethodCNTest
 * @brief cppunit class for testing CNAnalysisMethodCN functions.
 * last change by vliber on 01/29/09
 */

class CNAnalysisMethodCNTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodCNTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodCNTest );

void CNAnalysisMethodCNTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodCNTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodCN cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="cn-state");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="CopyNumber CNState");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="cn-state");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodCN::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==14);
	CPPUNIT_ASSERT(sv[0].name=="hmmCN_state");
	CPPUNIT_ASSERT(sv[0].type==0);
	CPPUNIT_ASSERT(sv[0].value=="0,1,2,3,4");
	CPPUNIT_ASSERT(sv[0].defaultVal=="0,1,2,3,4");
	CPPUNIT_ASSERT(sv[0].minVal=="");
	CPPUNIT_ASSERT(sv[0].maxVal=="");
	CPPUNIT_ASSERT(sv[0].descript=="CN HMM State");
	CPPUNIT_ASSERT(sv[1].name=="hmmCN_prior_prob");
	CPPUNIT_ASSERT(sv[1].type==0);
	CPPUNIT_ASSERT(sv[1].value=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(sv[1].defaultVal=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(sv[1].minVal=="");
	CPPUNIT_ASSERT(sv[1].maxVal=="");
	CPPUNIT_ASSERT(sv[1].descript=="CN HMM Prior Probability");
	CPPUNIT_ASSERT(sv[2].name=="hmmCN_mu");
	CPPUNIT_ASSERT(sv[2].type==0);
	CPPUNIT_ASSERT(sv[2].value=="-2,-0.533,0,0.363,0.567");
	CPPUNIT_ASSERT(sv[2].defaultVal=="-2,-0.533,0,0.363,0.567");
	CPPUNIT_ASSERT(sv[2].minVal=="");
	CPPUNIT_ASSERT(sv[2].maxVal=="");
	CPPUNIT_ASSERT(sv[2].descript=="CN HMM Mu");
	CPPUNIT_ASSERT(sv[3].name=="hmmCN_sigma");
	CPPUNIT_ASSERT(sv[3].type==0);
	CPPUNIT_ASSERT(sv[3].value=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(sv[3].defaultVal=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(sv[3].minVal=="");
	CPPUNIT_ASSERT(sv[3].maxVal=="");
	CPPUNIT_ASSERT(sv[3].descript=="CN HMM Sigma");
	CPPUNIT_ASSERT(sv[4].name=="hmmCN_TransitionDecay");
	CPPUNIT_ASSERT(sv[4].type==1);
	CPPUNIT_ASSERT(sv[4].value=="1e+9");
	CPPUNIT_ASSERT(sv[4].defaultVal=="1e+9");
	CPPUNIT_ASSERT(sv[4].minVal=="NA");
	CPPUNIT_ASSERT(sv[4].maxVal=="NA");
	CPPUNIT_ASSERT(sv[4].descript=="CN HMM Transition Decay");
	CPPUNIT_ASSERT(sv[5].name=="hmmCN_StateEstimationMethod");
	CPPUNIT_ASSERT(sv[5].type==0);
	CPPUNIT_ASSERT(sv[5].value=="EM");
	CPPUNIT_ASSERT(sv[5].defaultVal=="EM");
	CPPUNIT_ASSERT(sv[5].minVal=="");
	CPPUNIT_ASSERT(sv[5].maxVal=="");
	CPPUNIT_ASSERT(sv[5].descript=="CN HMM StateEstimation Method");
	CPPUNIT_ASSERT(sv[6].name=="hmmCN_EMIterations");
	CPPUNIT_ASSERT(sv[6].type==3);
	CPPUNIT_ASSERT(sv[6].value=="1");
	CPPUNIT_ASSERT(sv[6].defaultVal=="1");
	CPPUNIT_ASSERT(sv[6].minVal=="NA");
	CPPUNIT_ASSERT(sv[6].maxVal=="NA");
	CPPUNIT_ASSERT(sv[6].descript=="CN HMM Estimation Method Iterations");
	CPPUNIT_ASSERT(sv[7].name=="hmmCN_EMConvergenceThreshold");
	CPPUNIT_ASSERT(sv[7].type==1);
	CPPUNIT_ASSERT(sv[7].value=="0.0001");
	CPPUNIT_ASSERT(sv[7].defaultVal=="0.0001");
	CPPUNIT_ASSERT(sv[7].minVal=="NA");
	CPPUNIT_ASSERT(sv[7].maxVal=="NA");
	CPPUNIT_ASSERT(sv[7].descript=="CN HMM Estimation Method Convergence Threshold");
	CPPUNIT_ASSERT(sv[8].name=="hmmCN_NormalState");
	CPPUNIT_ASSERT(sv[8].type==3);
	CPPUNIT_ASSERT(sv[8].value=="2");
	CPPUNIT_ASSERT(sv[8].defaultVal=="2");
	CPPUNIT_ASSERT(sv[8].minVal=="NA");
	CPPUNIT_ASSERT(sv[8].maxVal=="NA");
	CPPUNIT_ASSERT(sv[8].descript=="CN HMM Normal State");
	CPPUNIT_ASSERT(sv[9].name=="hmmCN_ForwardOnly");
	CPPUNIT_ASSERT(sv[9].type==3);
	CPPUNIT_ASSERT(sv[9].value=="0");
	CPPUNIT_ASSERT(sv[9].defaultVal=="0");
	CPPUNIT_ASSERT(sv[9].minVal=="NA");
	CPPUNIT_ASSERT(sv[9].maxVal=="NA");
	CPPUNIT_ASSERT(sv[9].descript=="CN HMM Forward Only");
	CPPUNIT_ASSERT(sv[10].name=="hmmCN_NormalStateMinObservations");
	CPPUNIT_ASSERT(sv[10].type==3);
	CPPUNIT_ASSERT(sv[10].value=="2");
	CPPUNIT_ASSERT(sv[10].defaultVal=="2");
	CPPUNIT_ASSERT(sv[10].minVal=="NA");
	CPPUNIT_ASSERT(sv[10].maxVal=="NA");
	CPPUNIT_ASSERT(sv[10].descript=="CN HMM Normal State Min Observations");
    CPPUNIT_ASSERT(sv[11].name=="hmmCN_SmoothOutliers");
	CPPUNIT_ASSERT(sv[11].type==3);
	CPPUNIT_ASSERT(sv[11].value=="1");
	CPPUNIT_ASSERT(sv[11].defaultVal=="1");
	CPPUNIT_ASSERT(sv[11].minVal=="NA");
	CPPUNIT_ASSERT(sv[11].maxVal=="NA");
	CPPUNIT_ASSERT(sv[11].descript=="CN HMM Smooth Outliers");
    CPPUNIT_ASSERT(sv[12].name=="hmmCN_TransTypeStat");
	CPPUNIT_ASSERT(sv[12].type==3);
	CPPUNIT_ASSERT(sv[12].value=="0");
	CPPUNIT_ASSERT(sv[12].defaultVal=="0");
	CPPUNIT_ASSERT(sv[12].minVal=="NA");
	CPPUNIT_ASSERT(sv[12].maxVal=="NA");
	CPPUNIT_ASSERT(sv[12].descript=="CN HMM Trans Type Stat");
    CPPUNIT_ASSERT(sv[13].name=="PostCNFitMaxOutlierRemoveRunSize");
	CPPUNIT_ASSERT(sv[13].type==3);
	CPPUNIT_ASSERT(sv[13].value=="1");
	CPPUNIT_ASSERT(sv[13].defaultVal=="1");
	CPPUNIT_ASSERT(sv[13].minVal=="NA");
	CPPUNIT_ASSERT(sv[13].maxVal=="NA");
	CPPUNIT_ASSERT(sv[13].descript=="Post CN Fit Maximum Outlier Remove Run Size");

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodCN::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="cn-state.hmmCN_state=0,1,2,3,4.hmmCN_prior_prob=0.2,0.2,0.2,0.2,0.2.hmmCN_mu=-2,-0.533,0,0.363,0.567.hmmCN_sigma=0.2,0.2,0.2,0.2,0.2.hmmCN_TransitionDecay=1e+9.hmmCN_StateEstimationMethod=EM.hmmCN_EMIterations=1.hmmCN_EMConvergenceThreshold=0.0001.hmmCN_NormalState=2.hmmCN_ForwardOnly=0.hmmCN_NormalStateMinObservations=2.hmmCN_SmoothOutliers=1.hmmCN_TransTypeStat=0.PostCNFitMaxOutlierRemoveRunSize=1");
    CPPUNIT_ASSERT(sd.getDocName()=="cn-state");
	CPPUNIT_ASSERT(sd.getDocDescription()=="CopyNumber CNState");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==14);
	CPPUNIT_ASSERT(v[0].asString()=="0,1,2,3,4");
	CPPUNIT_ASSERT(v[1].asString()=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(v[2].asString()=="-2,-0.533,0,0.363,0.567");
	CPPUNIT_ASSERT(v[3].asString()=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(v[4].asString()=="1e+9");
	CPPUNIT_ASSERT(v[5].asString()=="EM");
	CPPUNIT_ASSERT(v[6].asString()=="1");
	CPPUNIT_ASSERT(v[7].asString()=="0.0001");
	CPPUNIT_ASSERT(v[8].asString()=="2");
	CPPUNIT_ASSERT(v[9].asString()=="0");
	CPPUNIT_ASSERT(v[10].asString()=="2");
	CPPUNIT_ASSERT(v[11].asString()=="1");
	CPPUNIT_ASSERT(v[12].asString()=="0");
	CPPUNIT_ASSERT(v[13].asString()=="1");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_state").asString()=="0,1,2,3,4");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_prior_prob").asString()=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_mu").asString()=="-2,-0.533,0,0.363,0.567");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_sigma").asString()=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_TransitionDecay").asString()=="1e+9");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_StateEstimationMethod").asString()=="EM");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_EMIterations").asString()=="1");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_EMConvergenceThreshold").asString()=="0.0001");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_NormalState").asString()=="2");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_ForwardOnly").asString()=="0");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_NormalStateMinObservations").asString()=="2");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_SmoothOutliers").asString()=="1");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_TransTypeStat").asString()=="0");
	CPPUNIT_ASSERT(sd.getDocOption("PostCNFitMaxOutlierRemoveRunSize").asString()=="1");
}

void CNAnalysisMethodCNTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodCNTest::newObjectTest****");

     //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodCN *cn2=(CNAnalysisMethodCN*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("cn-state");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==14);
	
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5a = obj1->at(4);
	affymetrix_calvin_parameter::ParameterNameValueType param6a = obj1->at(5);
	affymetrix_calvin_parameter::ParameterNameValueType param7a = obj1->at(6);
	affymetrix_calvin_parameter::ParameterNameValueType param8a = obj1->at(7);
	affymetrix_calvin_parameter::ParameterNameValueType param9a = obj1->at(8);
	affymetrix_calvin_parameter::ParameterNameValueType param10a = obj1->at(9);
	affymetrix_calvin_parameter::ParameterNameValueType param11a = obj1->at(10);
	affymetrix_calvin_parameter::ParameterNameValueType param12a = obj1->at(11);
	affymetrix_calvin_parameter::ParameterNameValueType param13a = obj1->at(12);
	affymetrix_calvin_parameter::ParameterNameValueType param14a = obj1->at(13);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-hmmCN_state");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param1a.GetValueAscii()=="0,1,2,3,4");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-hmmCN_prior_prob");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param2a.GetValueAscii()=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-hmmCN_mu");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param3a.GetValueAscii()=="-2,-0.533,0,0.363,0.567");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-hmmCN_sigma");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param4a.GetValueAscii()=="0.2,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5a.GetName())=="affymetrix-algorithm-param-hmmCN_TransitionDecay");
	CPPUNIT_ASSERT(param5a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5a.GetValueFloat()==1e+9);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6a.GetName())=="affymetrix-algorithm-param-hmmCN_StateEstimationMethod");
	CPPUNIT_ASSERT(param6a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param6a.GetValueAscii()=="EM");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7a.GetName())=="affymetrix-algorithm-param-hmmCN_ForwardOnly");
	CPPUNIT_ASSERT(param7a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param7a.GetValueInt32()==0);
    CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param8a.GetName())=="affymetrix-algorithm-param-hmmCN_NormalStateMinObservations");
	CPPUNIT_ASSERT(param8a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param8a.GetValueInt32()==2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param9a.GetName())=="affymetrix-algorithm-param-hmmCN_SmoothOutliers");
	CPPUNIT_ASSERT(param9a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param9a.GetValueInt32()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param10a.GetName())=="affymetrix-algorithm-param-hmmCN_TransTypeStat");
	CPPUNIT_ASSERT(param10a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param10a.GetValueInt32()==0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param11a.GetName())=="affymetrix-algorithm-param-hmmCN_EMIterations");
	CPPUNIT_ASSERT(param11a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param11a.GetValueInt32()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param12a.GetName())=="affymetrix-algorithm-param-hmmCN_EMConvergenceThreshold");
	CPPUNIT_ASSERT(param12a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param12a.GetValueFloat()==0.0001f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param13a.GetName())=="affymetrix-algorithm-param-hmmCN_NormalState");
	CPPUNIT_ASSERT(param13a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param13a.GetValueInt32()==2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param14a.GetName())=="affymetrix-algorithm-param-PostCNFitMaxOutlierRemoveRunSize");
	CPPUNIT_ASSERT(param14a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param14a.GetValueInt32()==1);
	
	delete cn2;

    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params1;
	params1["hmmCN_state"]="0,1,3,3,4";
	params1["hmmCN_prior_prob"]="0.3,0.2,0.2,0.2,0.2";
	params1["hmmCN_mu"]="-2,-0.333,0,0.363,0.567";
	params1["hmmCN_sigma"]="0.2,0.3,0.2,0.2,0.2";
	params1["hmmCN_TransitionDecay"]="1e+8";
	params1["hmmCN_StateEstimationMethod"]="EM1";
	params1["hmmCN_ForwardOnly"]="10";
	params1["hmmCN_NormalStateMinObservations"]="21";
	params1["hmmCN_SmoothOutliers"]="11";
    params1["hmmCN_TransTypeStat"]="10"; 
	params1["hmmCN_EMIterations"]="10";
	params1["hmmCN_EMConvergenceThreshold"]="0.0002";
	params1["hmmCN_NormalState"]="3";	
	params1["PostCNFitMaxOutlierRemoveRunSize"]="21";
	SelfCreate *sc=CNAnalysisMethodCN::newObject(params1);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==14);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5 = obj->at(4);
	affymetrix_calvin_parameter::ParameterNameValueType param6 = obj->at(5);
	affymetrix_calvin_parameter::ParameterNameValueType param7 = obj->at(6);
	affymetrix_calvin_parameter::ParameterNameValueType param8 = obj->at(7);
	affymetrix_calvin_parameter::ParameterNameValueType param9 = obj->at(8);
	affymetrix_calvin_parameter::ParameterNameValueType param10 = obj->at(9);
	affymetrix_calvin_parameter::ParameterNameValueType param11 = obj->at(10);
	affymetrix_calvin_parameter::ParameterNameValueType param12 = obj->at(11);
	affymetrix_calvin_parameter::ParameterNameValueType param13 = obj->at(12);
	affymetrix_calvin_parameter::ParameterNameValueType param14 = obj->at(13);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-hmmCN_state");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param1.GetValueAscii()=="0,1,3,3,4");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-hmmCN_prior_prob");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param2.GetValueAscii()=="0.3,0.2,0.2,0.2,0.2");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-hmmCN_mu");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param3.GetValueAscii()=="-2,-0.333,0,0.363,0.567");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-hmmCN_sigma");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param4.GetValueAscii()=="0.2,0.3,0.2,0.2,0.2");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5.GetName())=="affymetrix-algorithm-param-hmmCN_TransitionDecay");
	CPPUNIT_ASSERT(param5.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5.GetValueFloat()==1e+8);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6.GetName())=="affymetrix-algorithm-param-hmmCN_StateEstimationMethod");
	CPPUNIT_ASSERT(param6.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param6.GetValueAscii()=="EM1"); 
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7.GetName())=="affymetrix-algorithm-param-hmmCN_ForwardOnly");
	CPPUNIT_ASSERT(param7.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param7.GetValueInt32()==10);
    CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param8.GetName())=="affymetrix-algorithm-param-hmmCN_NormalStateMinObservations");
	CPPUNIT_ASSERT(param8.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param8.GetValueInt32()==21);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param9.GetName())=="affymetrix-algorithm-param-hmmCN_SmoothOutliers");
	CPPUNIT_ASSERT(param9.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param9.GetValueInt32()==11);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param10.GetName())=="affymetrix-algorithm-param-hmmCN_TransTypeStat");
	CPPUNIT_ASSERT(param10.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param10.GetValueInt32()==10);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param11.GetName())=="affymetrix-algorithm-param-hmmCN_EMIterations");
	CPPUNIT_ASSERT(param11.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param11.GetValueInt32()==10);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param12.GetName())=="affymetrix-algorithm-param-hmmCN_EMConvergenceThreshold");
	CPPUNIT_ASSERT(param12.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param12.GetValueFloat()==0.0002f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param13.GetName())=="affymetrix-algorithm-param-hmmCN_NormalState");
	CPPUNIT_ASSERT(param13.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param13.GetValueInt32()==3);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param14.GetName())=="affymetrix-algorithm-param-PostCNFitMaxOutlierRemoveRunSize");
	CPPUNIT_ASSERT(param14.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param14.GetValueInt32()==21);
	
    //Positive test. Application cares only size not content
	std::map<std::string,std::string> params2;
	params2["hmmCN_state"]="t,r,u,e,t";
	POSITIVE_TEST(CNAnalysisMethodCN::newObject(params2));
	params2.clear();
	//FATAL ERROR: Could not convert 'test' to a double.
    params2["hmmCN_TransitionDecay"]="test";
	NEGATIVE_TEST(CNAnalysisMethodCN::newObject(params2),Except);
	//positive with 0 value
	params2.clear();
    params2["hmmCN_StateEstimationMethod"]="0";
	POSITIVE_TEST(CNAnalysisMethodCN::newObject(params2));
	

	//Exceptions
    //PriorProb parameter must be the same size as the CNState parameter
	std::map<std::string,std::string> params3;
	params3["hmmCN_state"]="0,1,3,3,4";
	params3["hmmCN_prior_prob"]="0.3,0.2,0.2,0.2";
	NEGATIVE_TEST(CNAnalysisMethodCN::newObject(params3),Except);
  
	//Mu parameter must be the same size as the CNState parameter
    std::map<std::string,std::string> params4;
	params4["hmmCN_state"]="0,1,3,3,4";
	params4["hmmCN_mu"]="0.3,0.2,0.3,0.2";
	NEGATIVE_TEST(CNAnalysisMethodCN::newObject(params4),Except);

	//Sigma parameter must be the same size as the CNState parameter
    std::map<std::string,std::string> params5;
	params5["hmmCN_state"]="0,1,3,3,4";
	params5["hmmCN_sigma"]="0.3,0.2,0.3,0.2";
	NEGATIVE_TEST(CNAnalysisMethodCN::newObject(params5),Except);

	//PriorProb parameters must be greater than or equal to zero
    std::map<std::string,std::string> params6;
	params6["hmmCN_prior_prob"]="0.3,-0.2,0.2,0.2,0.2";
	NEGATIVE_TEST(CNAnalysisMethodCN::newObject(params6),Except);

	//Sigma parameters must be greater than or equal to zero
    std::map<std::string,std::string> params7;
	params7["hmmCN_sigma"]="0.2,0.3,0.2,0.2,-0.2";
	NEGATIVE_TEST(CNAnalysisMethodCN::newObject(params7),Except);

	//Mu parameters must be sorted and increasing in order.
    std::map<std::string,std::string> params8;
	params8["hmmCN_mu"]="-2,-0.333,1,0.363,0.567";
	NEGATIVE_TEST(CNAnalysisMethodCN::newObject(params8),Except);

	//Normal State parameter must be less than the number of CN States specified.
    std::map<std::string,std::string> params9;
	params9["hmmCN_NormalState"]="13";
	NEGATIVE_TEST(CNAnalysisMethodCN::newObject(params9),Except);

	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodCN::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="cn-state.hmmCN_state=0,1,2,3,4.hmmCN_prior_prob=0.2,0.2,0.2,0.2,0.2.hmmCN_mu=-2,-0.533,0,0.363,0.567.hmmCN_sigma=0.2,0.2,0.2,0.2,0.2.hmmCN_TransitionDecay=1e+9.hmmCN_StateEstimationMethod=EM.hmmCN_EMIterations=1.hmmCN_EMConvergenceThreshold=0.0001.hmmCN_NormalState=2.hmmCN_ForwardOnly=0.hmmCN_NormalStateMinObservations=2.hmmCN_SmoothOutliers=1.hmmCN_TransTypeStat=0.PostCNFitMaxOutlierRemoveRunSize=1");
    delete sc;
	
}	

void CNAnalysisMethodCNTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodCNTest::runTest****");
    
	CNAnalysisMethodCN cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod cn-state is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
