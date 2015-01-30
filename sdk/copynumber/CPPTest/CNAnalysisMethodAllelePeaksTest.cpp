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

#include "copynumber/CNAnalysisMethodAllelePeaks.h"
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
 * @class CNAnalysisMethodAllelePeaksTest
 * @brief cppunit class for testing CNAnalysisMethodAllelePeaks functions.
 * last change by vliber on 03/25/09
 */

class CNAnalysisMethodAllelePeaksTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodAllelePeaksTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest);
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodAllelePeaksTest );

void CNAnalysisMethodAllelePeaksTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodAllelePeaksTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodAllelePeaks cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="allele-peaks");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber AllelePeaks");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="allele-peaks");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodAllelePeaks::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==7);
	CPPUNIT_ASSERT(sv[0].name=="step");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="50");
	CPPUNIT_ASSERT(sv[0].defaultVal=="50");
	CPPUNIT_ASSERT(sv[0].minVal=="1");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="Allele Peaks step size.");
	CPPUNIT_ASSERT(sv[1].name=="window");
	CPPUNIT_ASSERT(sv[1].type==3);
	CPPUNIT_ASSERT(sv[1].value=="250");
	CPPUNIT_ASSERT(sv[1].defaultVal=="250");
	CPPUNIT_ASSERT(sv[1].minVal=="30");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="Allele Peaks window.");
	CPPUNIT_ASSERT(sv[2].name=="point-count");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="128");
	CPPUNIT_ASSERT(sv[2].defaultVal=="128");
	CPPUNIT_ASSERT(sv[2].minVal=="30");
	CPPUNIT_ASSERT(sv[2].maxVal=="1024");
	CPPUNIT_ASSERT(sv[2].descript=="Allele Peaks number of points.");
	CPPUNIT_ASSERT(sv[3].name=="bandwidth");
	CPPUNIT_ASSERT(sv[3].type==1);
	CPPUNIT_ASSERT(sv[3].value=="0.45");
	CPPUNIT_ASSERT(sv[3].defaultVal=="0.45");
	CPPUNIT_ASSERT(sv[3].minVal=="0");
	CPPUNIT_ASSERT(sv[3].maxVal=="1");
	CPPUNIT_ASSERT(sv[3].descript=="Allele Peaks bandwidth.");
	CPPUNIT_ASSERT(sv[4].name=="cutoff");
	CPPUNIT_ASSERT(sv[4].type==1);
	CPPUNIT_ASSERT(sv[4].value=="0.05");
	CPPUNIT_ASSERT(sv[4].defaultVal=="0.05");
	CPPUNIT_ASSERT(sv[4].minVal=="0");
	CPPUNIT_ASSERT(sv[4].maxVal=="0.5");
	CPPUNIT_ASSERT(sv[4].descript=="Allele Peaks cutoff.");
	

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodAllelePeaks::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="allele-peaks.step=50.window=250.point-count=128.bandwidth=0.45.cutoff=0.05.clean-threshold=0.25.symmetry=true");
    CPPUNIT_ASSERT(sd.getDocName()=="allele-peaks");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber AllelePeaks");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==7);
	CPPUNIT_ASSERT(v[0].asString()=="50");
	CPPUNIT_ASSERT(v[1].asString()=="250");
	CPPUNIT_ASSERT(v[2].asString()=="128");
	CPPUNIT_ASSERT(v[3].asString()=="0.45");
	CPPUNIT_ASSERT(v[4].asString()=="0.05");
	CPPUNIT_ASSERT(sd.getDocOption("step").asString()=="50");
	CPPUNIT_ASSERT(sd.getDocOption("window").asString()=="250");
	CPPUNIT_ASSERT(sd.getDocOption("point-count").asString()=="128");
	CPPUNIT_ASSERT(sd.getDocOption("bandwidth").asString()=="0.45");
	CPPUNIT_ASSERT(sd.getDocOption("cutoff").asString()=="0.05");
}

void CNAnalysisMethodAllelePeaksTest::newObjectTest()
{	
	Verbose::out(1, "****CNAnalysisMethodAllelePeaksTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodAllelePeaks *cn2=(CNAnalysisMethodAllelePeaks*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("allele-peaks");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==7);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5a = obj1->at(4);
	
		
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-step");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1a.GetValueInt32()==50);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-window");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param2a.GetValueInt32()==250);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-point-count");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==128);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-bandwidth");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param4a.GetValueFloat()==0.45f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5a.GetName())=="affymetrix-algorithm-param-cutoff");
	CPPUNIT_ASSERT(param5a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5a.GetValueFloat()==0.05f);
	delete cn2;

    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params1;
	params1["step"]="50";
	params1["window"]="100";
	params1["point-count"]="150";
	params1["bandwidth"]="0.2";
	params1["cutoff"]="0.003";
	SelfCreate *sc=CNAnalysisMethodAllelePeaks::newObject(params1);
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==7);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5 = obj->at(4);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-step");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param1.GetValueInt32()==50);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-window");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param2.GetValueInt32()==100);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-point-count");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3.GetValueInt32()==150);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-bandwidth");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param4.GetValueFloat()==0.2f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5.GetName())=="affymetrix-algorithm-param-cutoff");
	CPPUNIT_ASSERT(param5.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5.GetValueFloat()==0.003f);
	   
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodAllelePeaks::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="allele-peaks.step=50.window=250.point-count=128.bandwidth=0.45.cutoff=0.05.clean-threshold=0.25.symmetry=true");
    delete sc;
    
	//parameter tests 
    //step 1-NA
    params1.clear();
	params1["step"]="0";
	params1["window"]="100";
	params1["point-count"]="150";
	params1["bandwidth"]="0.2";
	params1["cutoff"]="0.003";
	NEGATIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1),Except);
	params1.clear();
	params1["step"]="1";
	params1["window"]="100";
	params1["point-count"]="150";
	params1["bandwidth"]="0.2";
	params1["cutoff"]="0.003";
	POSITIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1));
    //window 30-NA
	params1.clear();
	params1["step"]="10";
	params1["window"]="29";
	params1["point-count"]="31";
	params1["bandwidth"]="0.2";
	params1["cutoff"]="0.003";
	NEGATIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1),Except);
	params1.clear();
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="31";
	params1["bandwidth"]="0.2";
	params1["cutoff"]="0.003";
	POSITIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1));
	//point-count 30-1024
	params1.clear();
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="29";
	params1["bandwidth"]="0.2";
	params1["cutoff"]="0.003";
	NEGATIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1),Except);
	params1.clear();
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="30";
	params1["bandwidth"]="0.2";
	params1["cutoff"]="0.003";
	POSITIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1));
	params1.clear();
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1025";
	params1["bandwidth"]="0.2";
	params1["cutoff"]="0.003";
	NEGATIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1),Except);
	params1.clear();
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1024";
	params1["bandwidth"]="0.2";
	params1["cutoff"]="0.003";
	POSITIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1));
    //bandwidth 0-1
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1024";
	params1["bandwidth"]="-0.2";
	params1["cutoff"]="0.003";
	NEGATIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1),Except);
	params1.clear();
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1024";
	params1["bandwidth"]="0.0";
	params1["cutoff"]="0.003";
	POSITIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1));
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1024";
	params1["bandwidth"]="1.001";
	params1["cutoff"]="0.003";
	NEGATIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1),Except);
	params1.clear();
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1024";
	params1["bandwidth"]="1.0";
	params1["cutoff"]="0.003";
	POSITIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1));
	//cutoff 0-0.5
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1024";
	params1["bandwidth"]="0.5";
	params1["cutoff"]="-0.003";
	NEGATIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1),Except);
	params1.clear();
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1024";
	params1["bandwidth"]="0.5";
	params1["cutoff"]="0.0";
	POSITIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1));
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1024";
	params1["bandwidth"]="1.0";
	params1["cutoff"]="0.503";
	NEGATIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1),Except);
	params1.clear();
	params1["step"]="10";
	params1["window"]="30";
	params1["point-count"]="1024";
	params1["bandwidth"]="1.0";
	params1["cutoff"]="0.5";
	POSITIVE_TEST(CNAnalysisMethodAllelePeaks::newObject(params1));
		
}
void CNAnalysisMethodAllelePeaksTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodAllelePeaksTest::runTest****");
    
	CNAnalysisMethodAllelePeaks cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod allele-peaks is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
