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

//#include "chipstream/QuantLabelZ.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNAnalysisMethodGaussianSmooth.h"
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
 * @class CNAnalysisMethodGaussianSmoothTest
 * @brief cppunit class for testing CNAnalysisMethodGaussianSmooth functions.
 * last change by vliber on 01/28/09
 */

class CNAnalysisMethodGaussianSmoothTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodGaussianSmoothTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodGaussianSmoothTest );

void CNAnalysisMethodGaussianSmoothTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodGaussianSmoothTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodGaussianSmooth cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="gaussian-smooth");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="CopyNumber GaussianSmooth");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="gaussian-smooth");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodGaussianSmooth::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==3);
	CPPUNIT_ASSERT(sv[0].name=="expSmoothSignal");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="1");
	CPPUNIT_ASSERT(sv[0].defaultVal=="1");
	CPPUNIT_ASSERT(sv[0].minVal=="NA");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="Exp Smooth Signal");
	CPPUNIT_ASSERT(sv[1].name=="smooth_sigma_multiplier");
	CPPUNIT_ASSERT(sv[1].type==1);
	CPPUNIT_ASSERT(sv[1].value=="2");
	CPPUNIT_ASSERT(sv[1].defaultVal=="2");
	CPPUNIT_ASSERT(sv[1].minVal=="0.0");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="Smooth Sigma Multiplier");
	CPPUNIT_ASSERT(sv[2].name=="smooth_bw");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="50000");
	CPPUNIT_ASSERT(sv[2].defaultVal=="50000");
	CPPUNIT_ASSERT(sv[2].minVal=="1");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="Smooth Gaussian BW");
	//explainSelf()
	SelfDoc sd=CNAnalysisMethodGaussianSmooth::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="gaussian-smooth.expSmoothSignal=1.smooth_sigma_multiplier=2.smooth_bw=50000");
    CPPUNIT_ASSERT(sd.getDocName()=="gaussian-smooth");
	CPPUNIT_ASSERT(sd.getDocDescription()=="CopyNumber GaussianSmooth");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==3);
	CPPUNIT_ASSERT(v[0].asString()=="1");
	CPPUNIT_ASSERT(v[1].asString()=="2");
	CPPUNIT_ASSERT(v[2].asString()=="50000");
	CPPUNIT_ASSERT(sd.getDocOption("expSmoothSignal").asString()=="1");
	CPPUNIT_ASSERT(sd.getDocOption("smooth_sigma_multiplier").asString()=="2");
	CPPUNIT_ASSERT(sd.getDocOption("smooth_bw").asString()=="50000");
}

void CNAnalysisMethodGaussianSmoothTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodGaussianSmoothTest::newObjectTest****");

      //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodGaussianSmooth *cn2=(CNAnalysisMethodGaussianSmooth*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("gaussian-smooth");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-expSmoothSignal");
	CPPUNIT_ASSERT(param1a.GetParameterType()==4);
    CPPUNIT_ASSERT(param1a.GetValueInt32()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-smooth_sigma_multiplier");
	CPPUNIT_ASSERT(param2a.GetParameterType()==6);
	CPPUNIT_ASSERT(param2a.GetValueFloat()==2.0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-smooth_bw");
	CPPUNIT_ASSERT(param3a.GetParameterType()==4);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==50000);
	delete cn2;


     //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["expSmoothSignal"]="111";
	params["smooth_sigma_multiplier"]="222.0";
	params["smooth_bw"]="100";
	CNAnalysisMethod::getParams()->clear();
	SelfCreate *sc=CNAnalysisMethodGaussianSmooth::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-expSmoothSignal");
	CPPUNIT_ASSERT(param1.GetParameterType()==4);
    CPPUNIT_ASSERT(param1.GetValueInt32()==111);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-smooth_sigma_multiplier");
	CPPUNIT_ASSERT(param2.GetParameterType()==6);
	CPPUNIT_ASSERT(param2.GetValueFloat()==222.0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-smooth_bw");
	CPPUNIT_ASSERT(param3.GetParameterType()==4);
	CPPUNIT_ASSERT(param3.GetValueInt32()==100);
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodGaussianSmooth::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="gaussian-smooth.expSmoothSignal=1.smooth_sigma_multiplier=2.smooth_bw=50000");
    
    // Negative test passed [CNAnalysisMethodGaussianSmooth::newObject(params)] Message: SelfCreate::setValue() - '-222.0' is not a valid value for parameter: 'smooth_sigma_multiplier'. The specified range is 0.0 to NA
	params["smooth_sigma_multiplier"]="-222.0";
	NEGATIVE_TEST(CNAnalysisMethodGaussianSmooth::newObject(params),Except);
	params["smooth_sigma_multiplier"]="0.0";
	POSITIVE_TEST(CNAnalysisMethodGaussianSmooth::newObject(params));
	delete sc;
}	

void CNAnalysisMethodGaussianSmoothTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodGaussianSmoothTest::runTest****");

    
	CNAnalysisMethodGaussianSmooth cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod gaussian-smooth is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
