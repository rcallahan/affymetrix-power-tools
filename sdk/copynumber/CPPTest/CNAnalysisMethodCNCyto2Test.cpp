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

#include "copynumber/CNAnalysisMethodCNCyto2.h"
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
 * @class CNAnalysisMethodCNCyto2Test
 * @brief cppunit class for testing CNAnalysisMethodCNCyto2 functions.
 * last change by vliber on 04/29/09
 */

class CNAnalysisMethodCNCyto2Test : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodCNCyto2Test);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodCNCyto2Test );

void CNAnalysisMethodCNCyto2Test::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodCNCyto2Test::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodCNCyto2 cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="cn-cyto2");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="CopyNumber CNCyto2");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="cn-cyto2");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodCNCyto2::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==27);
	CPPUNIT_ASSERT(sv[0].name=="hmmCN_state");
	CPPUNIT_ASSERT(sv[0].type==0);
	CPPUNIT_ASSERT(sv[0].value=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(sv[0].defaultVal=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(sv[0].minVal=="");
	CPPUNIT_ASSERT(sv[0].maxVal=="");
	CPPUNIT_ASSERT(sv[0].descript=="CN HMM State");
	CPPUNIT_ASSERT(sv[1].name=="hmmCN_mu");
	CPPUNIT_ASSERT(sv[1].type==0);
	CPPUNIT_ASSERT(sv[1].value=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(sv[1].defaultVal=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(sv[1].minVal=="");
	CPPUNIT_ASSERT(sv[1].maxVal=="");
	CPPUNIT_ASSERT(sv[1].descript=="CN HMM Mu");
	CPPUNIT_ASSERT(sv[2].name=="hmmCN_sigma");
	CPPUNIT_ASSERT(sv[2].type==0);
	CPPUNIT_ASSERT(sv[2].value=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(sv[2].defaultVal=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(sv[2].minVal=="");
	CPPUNIT_ASSERT(sv[2].maxVal=="");
	CPPUNIT_ASSERT(sv[2].descript=="CN HMM Sigma");
	CPPUNIT_ASSERT(sv[3].name=="diagonal-weight");
	CPPUNIT_ASSERT(sv[3].type==1);
	CPPUNIT_ASSERT(sv[3].value=="0.995");
	CPPUNIT_ASSERT(sv[3].defaultVal=="0.995");
	CPPUNIT_ASSERT(sv[3].minVal=="0.0");
	CPPUNIT_ASSERT(sv[3].maxVal=="1.0");
	CPPUNIT_ASSERT(sv[3].descript=="CN HMM Diagonal Weight");
	CPPUNIT_ASSERT(sv[4].name=="mapd-weight");
	CPPUNIT_ASSERT(sv[4].type==1);
	CPPUNIT_ASSERT(sv[4].value=="0.22");
	CPPUNIT_ASSERT(sv[4].defaultVal=="0.22");
	CPPUNIT_ASSERT(sv[4].minVal=="0.0");
	CPPUNIT_ASSERT(sv[4].maxVal=="10000.0");
	CPPUNIT_ASSERT(sv[4].descript=="CN MAPD Weight");
	CPPUNIT_ASSERT(sv[5].name=="min-segment-size");
	CPPUNIT_ASSERT(sv[5].type==3);
	CPPUNIT_ASSERT(sv[5].value=="5");
	CPPUNIT_ASSERT(sv[5].defaultVal=="5");
	CPPUNIT_ASSERT(sv[5].minVal=="1");
	CPPUNIT_ASSERT(sv[5].maxVal=="NA");
	CPPUNIT_ASSERT(sv[5].descript=="CN Min Segment Size");
	CPPUNIT_ASSERT(sv[6].name=="hmm-confidence-weight");
	CPPUNIT_ASSERT(sv[6].type==1);
	CPPUNIT_ASSERT(sv[6].value=="0.6");
	CPPUNIT_ASSERT(sv[6].defaultVal=="0.6");
	CPPUNIT_ASSERT(sv[6].minVal=="0.0");
	CPPUNIT_ASSERT(sv[6].maxVal=="1.0");
	CPPUNIT_ASSERT(sv[6].descript=="HMM influence on confidence");

	CPPUNIT_ASSERT(sv[7].name=="hmmCN_state-X");
	CPPUNIT_ASSERT(sv[7].type==0);
	CPPUNIT_ASSERT(sv[7].value=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(sv[7].defaultVal=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(sv[7].minVal=="");
	CPPUNIT_ASSERT(sv[7].maxVal=="");
	CPPUNIT_ASSERT(sv[7].descript=="CN HMM State for X");
	CPPUNIT_ASSERT(sv[8].name=="hmmCN_mu-X");
	CPPUNIT_ASSERT(sv[8].type==0);
	CPPUNIT_ASSERT(sv[8].value=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(sv[8].defaultVal=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(sv[8].minVal=="");
	CPPUNIT_ASSERT(sv[8].maxVal=="");
	CPPUNIT_ASSERT(sv[8].descript=="CN HMM Mu for X");
	CPPUNIT_ASSERT(sv[9].name=="hmmCN_sigma-X");
	CPPUNIT_ASSERT(sv[9].type==0);
	CPPUNIT_ASSERT(sv[9].value=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(sv[9].defaultVal=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(sv[9].minVal=="");
	CPPUNIT_ASSERT(sv[9].maxVal=="");
	CPPUNIT_ASSERT(sv[9].descript=="CN HMM Sigma for X");
	CPPUNIT_ASSERT(sv[10].name=="diagonal-weight-X");
	CPPUNIT_ASSERT(sv[10].type==1);
	CPPUNIT_ASSERT(sv[10].value=="0.995");
	CPPUNIT_ASSERT(sv[10].defaultVal=="0.995");
	CPPUNIT_ASSERT(sv[10].minVal=="0.0");
	CPPUNIT_ASSERT(sv[10].maxVal=="1.0");
	CPPUNIT_ASSERT(sv[10].descript=="CN HMM Diagonal Weight for X");
	CPPUNIT_ASSERT(sv[11].name=="mapd-weight-X");
	CPPUNIT_ASSERT(sv[11].type==1);
	CPPUNIT_ASSERT(sv[11].value=="0.22");
	CPPUNIT_ASSERT(sv[11].defaultVal=="0.22");
	CPPUNIT_ASSERT(sv[11].minVal=="0.0");
	CPPUNIT_ASSERT(sv[11].maxVal=="10000.0");
	CPPUNIT_ASSERT(sv[11].descript=="CN MAPD Weight for X");
	CPPUNIT_ASSERT(sv[12].name=="min-segment-size-X");
	CPPUNIT_ASSERT(sv[12].type==3);
	CPPUNIT_ASSERT(sv[12].value=="5");
	CPPUNIT_ASSERT(sv[12].defaultVal=="5");
	CPPUNIT_ASSERT(sv[12].minVal=="1");
	CPPUNIT_ASSERT(sv[12].maxVal=="NA");
	CPPUNIT_ASSERT(sv[12].descript=="CN Min Segment Size for X");
	CPPUNIT_ASSERT(sv[13].name=="hmm-confidence-weight-X");
	CPPUNIT_ASSERT(sv[13].type==1);
	CPPUNIT_ASSERT(sv[13].value=="0.6");
	CPPUNIT_ASSERT(sv[13].defaultVal=="0.6");
	CPPUNIT_ASSERT(sv[13].minVal=="0.0");
	CPPUNIT_ASSERT(sv[13].maxVal=="1.0");
	CPPUNIT_ASSERT(sv[13].descript=="HMM influence on confidence for X");

	CPPUNIT_ASSERT(sv[14].name=="hmmCN_state-Y");
	CPPUNIT_ASSERT(sv[14].type==0);
	CPPUNIT_ASSERT(sv[14].value=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(sv[14].defaultVal=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(sv[14].minVal=="");
	CPPUNIT_ASSERT(sv[14].maxVal=="");
	CPPUNIT_ASSERT(sv[14].descript=="CN HMM State for Y");
	CPPUNIT_ASSERT(sv[15].name=="hmmCN_mu-Y");
	CPPUNIT_ASSERT(sv[15].type==0);
	CPPUNIT_ASSERT(sv[15].value=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(sv[15].defaultVal=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(sv[15].minVal=="");
	CPPUNIT_ASSERT(sv[15].maxVal=="");
	CPPUNIT_ASSERT(sv[15].descript=="CN HMM Mu for Y");
	CPPUNIT_ASSERT(sv[16].name=="hmmCN_sigma-Y");
	CPPUNIT_ASSERT(sv[16].type==0);
	CPPUNIT_ASSERT(sv[16].value=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(sv[16].defaultVal=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(sv[16].minVal=="");
	CPPUNIT_ASSERT(sv[16].maxVal=="");
	CPPUNIT_ASSERT(sv[16].descript=="CN HMM Sigma for Y");
	CPPUNIT_ASSERT(sv[17].name=="diagonal-weight-Y");
	CPPUNIT_ASSERT(sv[17].type==1);
	CPPUNIT_ASSERT(sv[17].value=="0.995");
	CPPUNIT_ASSERT(sv[17].defaultVal=="0.995");
	CPPUNIT_ASSERT(sv[17].minVal=="0.0");
	CPPUNIT_ASSERT(sv[17].maxVal=="1.0");
	CPPUNIT_ASSERT(sv[17].descript=="CN HMM Diagonal Weight for Y");
	CPPUNIT_ASSERT(sv[18].name=="mapd-weight-Y");
	CPPUNIT_ASSERT(sv[18].type==1);
	CPPUNIT_ASSERT(sv[18].value=="0.22");
	CPPUNIT_ASSERT(sv[18].defaultVal=="0.22");
	CPPUNIT_ASSERT(sv[18].minVal=="0.0");
	CPPUNIT_ASSERT(sv[18].maxVal=="10000.0");
	CPPUNIT_ASSERT(sv[18].descript=="CN MAPD Weight for Y");
	CPPUNIT_ASSERT(sv[19].name=="min-segment-size-Y");
	CPPUNIT_ASSERT(sv[19].type==3);
	CPPUNIT_ASSERT(sv[19].value=="5");
	CPPUNIT_ASSERT(sv[19].defaultVal=="5");
	CPPUNIT_ASSERT(sv[19].minVal=="1");
	CPPUNIT_ASSERT(sv[19].maxVal=="NA");
	CPPUNIT_ASSERT(sv[19].descript=="CN Min Segment Size for Y");
	CPPUNIT_ASSERT(sv[20].name=="hmm-confidence-weight-Y");
	CPPUNIT_ASSERT(sv[20].type==1);
	CPPUNIT_ASSERT(sv[20].value=="0.6");
	CPPUNIT_ASSERT(sv[20].defaultVal=="0.6");
	CPPUNIT_ASSERT(sv[20].minVal=="0.0");
	CPPUNIT_ASSERT(sv[20].maxVal=="1.0");
	CPPUNIT_ASSERT(sv[20].descript=="HMM influence on confidence for Y");

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodCNCyto2::explainSelf();
	CPPUNIT_ASSERT(sd.getState()==string("") +
        "cn-cyto2.hmmCN_state=0,1,2,3,4,5.hmmCN_mu=-1.63,-0.58,0,0.45,0.72,0.93.hmmCN_sigma=0.3,0.3,0.3,0.3,0.3,0.3.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6" +
        ".hmmCN_state-X=0,1,2,3,4,5.hmmCN_mu-X=-1.63,-0.58,0,0.45,0.72,0.93.hmmCN_sigma-X=0.3,0.3,0.3,0.3,0.3,0.3.diagonal-weight-X=0.995.mapd-weight-X=0.22.min-segment-size-X=5.hmm-confidence-weight-X=0.6" +
        ".hmmCN_state-Y=0,1,2,3,4,5.hmmCN_mu-Y=-1.63,-0.58,0,0.45,0.72,0.93.hmmCN_sigma-Y=0.3,0.3,0.3,0.3,0.3,0.3.diagonal-weight-Y=0.995.mapd-weight-Y=0.22.min-segment-size-Y=5.hmm-confidence-weight-Y=0.6" +
        ".shrink=false.shrink-lprec=0.5,4.0.shrink-converge=0.000001.shrink-downweight-outlier=true.shrink-downweight-df=6.5.shrink-downweight-maxiter=15"
        );
    CPPUNIT_ASSERT(sd.getDocName()=="cn-cyto2");
	CPPUNIT_ASSERT(sd.getDocDescription()=="CopyNumber CNCyto2");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==27);
	CPPUNIT_ASSERT(v[0].asString()=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(v[1].asString()=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(v[2].asString()=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(v[3].asString()=="0.995");
	CPPUNIT_ASSERT(v[4].asString()=="0.22");
	CPPUNIT_ASSERT(v[5].asString()=="5");
	CPPUNIT_ASSERT(v[6].asString()=="0.6");

	CPPUNIT_ASSERT(v[7].asString()=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(v[8].asString()=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(v[9].asString()=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(v[10].asString()=="0.995");
	CPPUNIT_ASSERT(v[11].asString()=="0.22");
	CPPUNIT_ASSERT(v[12].asString()=="5");
	CPPUNIT_ASSERT(v[13].asString()=="0.6");

	CPPUNIT_ASSERT(v[14].asString()=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(v[15].asString()=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(v[16].asString()=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(v[17].asString()=="0.995");
	CPPUNIT_ASSERT(v[18].asString()=="0.22");
	CPPUNIT_ASSERT(v[19].asString()=="5");
	CPPUNIT_ASSERT(v[20].asString()=="0.6");

	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_state").asString()=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_mu").asString()=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_sigma").asString()=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(sd.getDocOption("diagonal-weight").asString()=="0.995");
	CPPUNIT_ASSERT(sd.getDocOption("mapd-weight").asString()=="0.22");
	CPPUNIT_ASSERT(sd.getDocOption("min-segment-size").asString()=="5");
	CPPUNIT_ASSERT(sd.getDocOption("hmm-confidence-weight").asString()=="0.6");

	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_state-X").asString()=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_mu-X").asString()=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_sigma-X").asString()=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(sd.getDocOption("diagonal-weight-X").asString()=="0.995");
	CPPUNIT_ASSERT(sd.getDocOption("mapd-weight-X").asString()=="0.22");
	CPPUNIT_ASSERT(sd.getDocOption("min-segment-size-X").asString()=="5");
	CPPUNIT_ASSERT(sd.getDocOption("hmm-confidence-weight-X").asString()=="0.6");

	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_state-Y").asString()=="0,1,2,3,4,5");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_mu-Y").asString()=="-1.63,-0.58,0,0.45,0.72,0.93");
	CPPUNIT_ASSERT(sd.getDocOption("hmmCN_sigma-Y").asString()=="0.3,0.3,0.3,0.3,0.3,0.3");
	CPPUNIT_ASSERT(sd.getDocOption("diagonal-weight-Y").asString()=="0.995");
	CPPUNIT_ASSERT(sd.getDocOption("mapd-weight-Y").asString()=="0.22");
	CPPUNIT_ASSERT(sd.getDocOption("min-segment-size-Y").asString()=="5");
	CPPUNIT_ASSERT(sd.getDocOption("hmm-confidence-weight-Y").asString()=="0.6");
}

void CNAnalysisMethodCNCyto2Test::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodCNCyto2Test::newObjectTest****");

     //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodCNCyto2 *cn2=(CNAnalysisMethodCNCyto2*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("cn-cyto2");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==27);
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

	affymetrix_calvin_parameter::ParameterNameValueType param15a = obj1->at(14);
	affymetrix_calvin_parameter::ParameterNameValueType param16a = obj1->at(15);
	affymetrix_calvin_parameter::ParameterNameValueType param17a = obj1->at(16);
	affymetrix_calvin_parameter::ParameterNameValueType param18a = obj1->at(17);
	affymetrix_calvin_parameter::ParameterNameValueType param19a = obj1->at(18);
	affymetrix_calvin_parameter::ParameterNameValueType param20a = obj1->at(19);
	affymetrix_calvin_parameter::ParameterNameValueType param21a = obj1->at(20);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-diagonal-weight");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1a.GetValueFloat()==0.995f);

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-mapd-weight");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2a.GetValueFloat()==0.22f);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-min-segment-size");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==5);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-hmmCN_state");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param4a.GetValueAscii()=="0,1,2,3,4,5");
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5a.GetName())=="affymetrix-algorithm-param-hmmCN_mu");
	CPPUNIT_ASSERT(param5a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param5a.GetValueAscii()=="-1.63,-0.58,0,0.45,0.72,0.93");
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6a.GetName())=="affymetrix-algorithm-param-hmmCN_sigma");
	CPPUNIT_ASSERT(param6a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param6a.GetValueAscii()=="0.3,0.3,0.3,0.3,0.3,0.3");

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7a.GetName())=="affymetrix-algorithm-param-hmm-confidence-weight");
	CPPUNIT_ASSERT(param7a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param7a.GetValueFloat()==0.6f);

////////////////
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param8a.GetName())=="affymetrix-algorithm-param-diagonal-weight-X");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1a.GetValueFloat()==0.995f);

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param9a.GetName())=="affymetrix-algorithm-param-mapd-weight-X");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2a.GetValueFloat()==0.22f);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param10a.GetName())=="affymetrix-algorithm-param-min-segment-size-X");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==5);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param11a.GetName())=="affymetrix-algorithm-param-hmmCN_state-X");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param4a.GetValueAscii()=="0,1,2,3,4,5");
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param12a.GetName())=="affymetrix-algorithm-param-hmmCN_mu-X");
	CPPUNIT_ASSERT(param5a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param5a.GetValueAscii()=="-1.63,-0.58,0,0.45,0.72,0.93");
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param13a.GetName())=="affymetrix-algorithm-param-hmmCN_sigma-X");
	CPPUNIT_ASSERT(param6a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param6a.GetValueAscii()=="0.3,0.3,0.3,0.3,0.3,0.3");

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param14a.GetName())=="affymetrix-algorithm-param-hmm-confidence-weight-X");
	CPPUNIT_ASSERT(param7a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param7a.GetValueFloat()==0.6f);

////////////////
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param15a.GetName())=="affymetrix-algorithm-param-diagonal-weight-Y");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1a.GetValueFloat()==0.995f);

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param16a.GetName())=="affymetrix-algorithm-param-mapd-weight-Y");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2a.GetValueFloat()==0.22f);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param17a.GetName())=="affymetrix-algorithm-param-min-segment-size-Y");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==5);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param18a.GetName())=="affymetrix-algorithm-param-hmmCN_state-Y");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param4a.GetValueAscii()=="0,1,2,3,4,5");
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param19a.GetName())=="affymetrix-algorithm-param-hmmCN_mu-Y");
	CPPUNIT_ASSERT(param5a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param5a.GetValueAscii()=="-1.63,-0.58,0,0.45,0.72,0.93");
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param20a.GetName())=="affymetrix-algorithm-param-hmmCN_sigma-Y");
	CPPUNIT_ASSERT(param6a.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param6a.GetValueAscii()=="0.3,0.3,0.3,0.3,0.3,0.3");

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param21a.GetName())=="affymetrix-algorithm-param-hmm-confidence-weight-Y");
	CPPUNIT_ASSERT(param7a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param7a.GetValueFloat()==0.6f);
	
   
	delete cn2;

    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params;
	params["diagonal-weight"]="0.1";
	params["mapd-weight"]="0.2";
	params["min-segment-size"]="2";
	params["hmmCN_state"]="0,1,7,3,9,5";
	params["hmmCN_mu"]="-2.5,-0.62,0.0,0.49,0.69,1.94";
	params["hmmCN_sigma"]="0.3,0.2,0.3,0.2,0.3,0.2";
	params["hmm-confidence-weight"]="0.9";
	SelfCreate *sc=CNAnalysisMethodCNCyto2::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==27);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5 = obj->at(4);
	affymetrix_calvin_parameter::ParameterNameValueType param6 = obj->at(5);
	affymetrix_calvin_parameter::ParameterNameValueType param7 = obj->at(6);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-diagonal-weight");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1.GetValueFloat()==0.1f); 
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-mapd-weight");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2.GetValueFloat()==0.2f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-min-segment-size");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3.GetValueInt32()==2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-hmmCN_state");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param4.GetValueAscii()=="0,1,7,3,9,5");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5.GetName())=="affymetrix-algorithm-param-hmmCN_mu");
	CPPUNIT_ASSERT(param5.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param5.GetValueAscii()=="-2.5,-0.62,0.0,0.49,0.69,1.94");
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6.GetName())=="affymetrix-algorithm-param-hmmCN_sigma");
	CPPUNIT_ASSERT(param6.GetParameterType()==ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT(param6.GetValueAscii()=="0.3,0.2,0.3,0.2,0.3,0.2");
    CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7.GetName())=="affymetrix-algorithm-param-hmm-confidence-weight");
	CPPUNIT_ASSERT(param7.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param7.GetValueFloat()==0.9f);
	//Exceptions 
    //Mu must conform to CNState.
    std::map<std::string,std::string> params1;
	params1["hmmCN_state"]="0,1,7,3,9,5,7";
	params1["hmmCN_mu"]="-2.5,-0.62,0.0,0.49,0.69,1.94";
	NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params1),Except);

	//Sigma must conform to CNState.
    std::map<std::string,std::string> params2;
	params2["hmmCN_state"]="0,1,7,3,9,5";
	params2["hmmCN_sigma"]="0.3,0.2,0.3,0.2,0.3,0.2,0.3";
	NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params2),Except);

	//Sigma parameters must be > 0.
    std::map<std::string,std::string> params3;
	params3["hmmCN_sigma"]="0.3,-0.2,0.3,0.2,0.3,0.2";
	NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params3),Except);

	//Mu parameters must ordered increasing
    std::map<std::string,std::string> params4;
	params4["hmmCN_mu"]="-2.5,-0.62,0.0,2.49,0.69,1.94";
	NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params4),Except);
    
	
	std::map<std::string,std::string> params5;
	//FATAL ERROR: SelfCreate::setValue() - '111.11' is not a valid value for parameter: 'diagonal-weight'. The specified range is 0.0 to 1.0
	params5["diagonal-weight"]="111.11";
	NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);

	//FATAL ERROR: SelfCreate::setValue() - '-111.11' is not a valid value for parameter: 'diagonal-weight'. The specified range is 0.0 to 1.0
	params5["diagonal-weight"]="-111.11";
	NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);

	//FATAL ERROR: Could not convert 'a' to a double.	
    params5["diagonal-weight"]="a";
	NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);

    
	params5["diagonal-weight"]="0.1";
	//FATAL ERROR: SelfCreate::setValue() - '1.8' is not a valid value for parameter: 'mapd-weight'. The specified range is 0.0 to 10000.0
    params5["mapd-weight"]="-1.8";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);
	//FATAL ERROR: SelfCreate::setValue() - '-1111111.1' is not a valid value for parameter: 'mapd-weight'. The specified range is 0.0 to 10000.0
	
	params5["mapd-weight"]="-1111111.1";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);
    
	//FATAL ERROR: Could not convert 'a' to a double.
	params5["mapd-weight"]="a";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);

	params5["mapd-weight"]="0.1";
	//FATAL ERROR: SelfCreate::setValue() - '0' is not a valid value for parameter: 'min-segment-size'. The specified range is 1 to NA
    params5["min-segment-size"]="0";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);
	//FATAL ERROR: SelfCreate::setValue() - '-10' is not a valid value for parameter: 'min-segment-size'. The specified range is 1 to NA
	params5["min-segment-size"]="-10";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);
	//FATAL ERROR: Could not convert 'a' to a int.
	params5["min-segment-size"]="a";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);

	//HMM Confidence weight must be in [0,1]
	params5["min-segment-size"]="5";
	params5["hmm-confidence-weight"]="0.0";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);
	params5["hmm-confidence-weight"]="1.0";
    NEGATIVE_TEST(CNAnalysisMethodCNCyto2::newObject(params5),Except);

		
	 //always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodCNCyto2::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()==string("") +
        "cn-cyto2.hmmCN_state=0,1,2,3,4,5.hmmCN_mu=-1.63,-0.58,0,0.45,0.72,0.93.hmmCN_sigma=0.3,0.3,0.3,0.3,0.3,0.3.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6" +
        ".hmmCN_state-X=0,1,2,3,4,5.hmmCN_mu-X=-1.63,-0.58,0,0.45,0.72,0.93.hmmCN_sigma-X=0.3,0.3,0.3,0.3,0.3,0.3.diagonal-weight-X=0.995.mapd-weight-X=0.22.min-segment-size-X=5.hmm-confidence-weight-X=0.6" +
        ".hmmCN_state-Y=0,1,2,3,4,5.hmmCN_mu-Y=-1.63,-0.58,0,0.45,0.72,0.93.hmmCN_sigma-Y=0.3,0.3,0.3,0.3,0.3,0.3.diagonal-weight-Y=0.995.mapd-weight-Y=0.22.min-segment-size-Y=5.hmm-confidence-weight-Y=0.6" +
        ".shrink=false.shrink-lprec=0.5,4.0.shrink-converge=0.000001.shrink-downweight-outlier=true.shrink-downweight-df=6.5.shrink-downweight-maxiter=15"
        );
 
    delete sc;
	

	

}	

void CNAnalysisMethodCNCyto2Test::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodCNCyto2Test::runTest****");
    
	CNAnalysisMethodCNCyto2 cnrfCNam;
    //CNAnalysisMethod cn-cyto2 is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
