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
#include "copynumber/CNAnalysisMethodKernelSmooth.h"
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
 * @class CNAnalysisMethodKernelSmoothTest
 * @brief cppunit class for testing CNAnalysisMethodKernelSmooth functions.
 * last change by vliber on 01/28/09
 */

class CNAnalysisMethodKernelSmoothTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodKernelSmoothTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodKernelSmoothTest );

void CNAnalysisMethodKernelSmoothTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodKernelSmoothTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodKernelSmooth cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="kernel-smooth");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="CopyNumber KernelSmooth");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="kernel-smooth");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodKernelSmooth::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==1);
	CPPUNIT_ASSERT(sv[0].name=="sigma_span");
	CPPUNIT_ASSERT(sv[0].type==1);
	CPPUNIT_ASSERT(sv[0].value=="50.0");
	CPPUNIT_ASSERT(sv[0].defaultVal=="50.0");
	CPPUNIT_ASSERT(sv[0].minVal=="0.0");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="Probes spanned by one stdev");

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodKernelSmooth::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="kernel-smooth.sigma_span=50.0");
    CPPUNIT_ASSERT(sd.getDocName()=="kernel-smooth");
	CPPUNIT_ASSERT(sd.getDocDescription()=="CopyNumber KernelSmooth");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==1);
	CPPUNIT_ASSERT(v[0].asString()=="50.0");
	CPPUNIT_ASSERT(sd.getDocOption("sigma_span").asString()=="50.0");
	
}

void CNAnalysisMethodKernelSmoothTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodKernelSmoothTest::newObjectTest****");

     //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodKernelSmooth *cn2=(CNAnalysisMethodKernelSmooth*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("kernel-smooth");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-sigma_span");
	CPPUNIT_ASSERT(param1a.GetParameterType()==6);
    CPPUNIT_ASSERT(param1a.GetValueFloat()==50.0f);
	delete cn2;


    //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["sigma_span"]="0.111";
	
	CNAnalysisMethod::getParams()->clear();
	SelfCreate *sc=CNAnalysisMethodKernelSmooth::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-sigma_span");
	CPPUNIT_ASSERT(param1.GetParameterType()==6);
    CPPUNIT_ASSERT(param1.GetValueFloat()==0.111f);
	
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodKernelSmooth::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="kernel-smooth.sigma_span=50.0");
    delete sc;	
    
	//FATAL ERROR: SelfCreate::setValue() - '-0.111' is not a valid value for parameter: 'sigma_span'. The specified range is 0.0 to NA
	params["sigma_span"]="-0.111";
    NEGATIVE_TEST(CNAnalysisMethodKernelSmooth::newObject(params), Except);
	//FATAL ERROR: Could not convert 'a' to a double.
	params["sigma_span"]="a";
	NEGATIVE_TEST(CNAnalysisMethodKernelSmooth::newObject(params),Except);
	params["sigma_span"]="0";
    POSITIVE_TEST(CNAnalysisMethodKernelSmooth::newObject(params));

}	

void CNAnalysisMethodKernelSmoothTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodKernelSmoothTest::runTest****");

    
	CNAnalysisMethodKernelSmooth cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod kernel-smooth is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
