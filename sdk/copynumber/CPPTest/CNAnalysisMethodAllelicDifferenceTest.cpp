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

#include "copynumber/CNAnalysisMethodAllelicDifference.h"
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
 * @class CNAnalysisMethodAllelicDifferenceTest
 * @brief cppunit class for testing CNAnalysisMethodAllelicDifference functions.
 * last change by vliber on 01/28/09
 */

class CNAnalysisMethodAllelicDifferenceTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodAllelicDifferenceTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodAllelicDifferenceTest );

void CNAnalysisMethodAllelicDifferenceTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodAllelicDifferenceTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodAllelicDifference cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="allelic-difference");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber AllelicDifference");
	
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodAllelicDifference::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==1);
	CPPUNIT_ASSERT(sv[0].name=="outlier-trim");
	CPPUNIT_ASSERT(sv[0].type==1);
	CPPUNIT_ASSERT(sv[0].value=="3.0");
	CPPUNIT_ASSERT(sv[0].defaultVal=="3.0");
	CPPUNIT_ASSERT(sv[0].minVal=="NA");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="AllelicDifference Outlier Trim");
#if 0
	CPPUNIT_ASSERT(sv[1].name=="step");
	CPPUNIT_ASSERT(sv[1].type==3);
	CPPUNIT_ASSERT(sv[1].value=="20");
	CPPUNIT_ASSERT(sv[1].defaultVal=="20");
	CPPUNIT_ASSERT(sv[1].minVal=="1");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="AllelicDifference step size");

	CPPUNIT_ASSERT(sv[2].name=="window");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="100");
	CPPUNIT_ASSERT(sv[2].defaultVal=="100");
	CPPUNIT_ASSERT(sv[2].minVal=="30");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="AllelicDifference window");

	CPPUNIT_ASSERT(sv[3].name=="point-count");
	CPPUNIT_ASSERT(sv[3].type==3);
	CPPUNIT_ASSERT(sv[3].value=="128");
	CPPUNIT_ASSERT(sv[3].defaultVal=="128");
	CPPUNIT_ASSERT(sv[3].minVal=="30");
	CPPUNIT_ASSERT(sv[3].maxVal=="1024");
	CPPUNIT_ASSERT(sv[3].descript=="AllelicDifference number of points");

	CPPUNIT_ASSERT(sv[4].name=="bandwidth");
	CPPUNIT_ASSERT(sv[4].type==1);
	CPPUNIT_ASSERT(sv[4].value=="0.25");
	CPPUNIT_ASSERT(sv[4].defaultVal=="0.25");
	CPPUNIT_ASSERT(sv[4].minVal=="0");
	CPPUNIT_ASSERT(sv[4].maxVal=="1");
	CPPUNIT_ASSERT(sv[4].descript=="AllelicDifference bandwidth");

	CPPUNIT_ASSERT(sv[5].name=="cutoff");
	CPPUNIT_ASSERT(sv[5].type==1);
	CPPUNIT_ASSERT(sv[5].value=="0.05");
	CPPUNIT_ASSERT(sv[5].defaultVal=="0.05");
	CPPUNIT_ASSERT(sv[5].minVal=="0");
	CPPUNIT_ASSERT(sv[5].maxVal=="0.5");
	CPPUNIT_ASSERT(sv[5].descript=="AllelicDifference cutoff");

	CPPUNIT_ASSERT(sv[6].name=="clean-threshold");
	CPPUNIT_ASSERT(sv[6].type==1);
	CPPUNIT_ASSERT(sv[6].value=="0.35");
	CPPUNIT_ASSERT(sv[6].defaultVal=="0.35");
	CPPUNIT_ASSERT(sv[6].minVal=="0");
	CPPUNIT_ASSERT(sv[6].maxVal=="1.0");
	CPPUNIT_ASSERT(sv[6].descript=="AllelicDifference clean threshold");

	CPPUNIT_ASSERT(sv[7].name=="symmetry");
	CPPUNIT_ASSERT(sv[7].type==4);
	CPPUNIT_ASSERT(sv[7].value=="true");
	CPPUNIT_ASSERT(sv[7].defaultVal=="true");
	CPPUNIT_ASSERT(sv[7].minVal=="NA");
	CPPUNIT_ASSERT(sv[7].maxVal=="NA");
	CPPUNIT_ASSERT(sv[7].descript=="AllelicDifference SCAR mirror flag");
#endif

	//explainSelf()
	SelfDoc sd=CNAnalysisMethodAllelicDifference::explainSelf();
	//CPPUNIT_ASSERT(sd.getState()=="allelic-difference.outlier-trim=3.0.step=20.window=100.point-count=128.bandwidth=0.25.cutoff=0.05.clean-threshold=0.35.symmetry=true");
	CPPUNIT_ASSERT(sd.getState()=="allelic-difference.outlier-trim=3.0");
    CPPUNIT_ASSERT(sd.getDocName()=="allelic-difference");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Copynumber AllelicDifference");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==1);
	CPPUNIT_ASSERT(v[0].asString()=="3.0");
#if 0
    CPPUNIT_ASSERT(v[1].asString()=="20");
    CPPUNIT_ASSERT(v[2].asString()=="100");
    CPPUNIT_ASSERT(v[3].asString()=="128");
    CPPUNIT_ASSERT(v[4].asString()=="0.25");
    CPPUNIT_ASSERT(v[5].asString()=="0.05");
    CPPUNIT_ASSERT(v[6].asString()=="0.35");
    CPPUNIT_ASSERT(v[7].asString()=="true");
#endif
	CPPUNIT_ASSERT(sd.getDocOption("outlier-trim").asString()=="3.0");
#if 0
    CPPUNIT_ASSERT(sd.getDocOption("step").asString()=="20");
    CPPUNIT_ASSERT(sd.getDocOption("window").asString()=="100");
    CPPUNIT_ASSERT(sd.getDocOption("point-count").asString()=="128");
    CPPUNIT_ASSERT(sd.getDocOption("bandwidth").asString()=="0.25");
    CPPUNIT_ASSERT(sd.getDocOption("cutoff").asString()=="0.05");
    CPPUNIT_ASSERT(sd.getDocOption("clean-threshold").asString()=="0.35");
    CPPUNIT_ASSERT(sd.getDocOption("symmetry").asString()=="true");
#endif
}

void CNAnalysisMethodAllelicDifferenceTest::newObjectTest()
{
	
	Verbose::out(1, "****CNAnalysisMethodAllelicDifferenceTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNAnalysisMethodAllelicDifference *cn2=(CNAnalysisMethodAllelicDifference*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("allelic-difference");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param = obj->at(0);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param.GetName())=="affymetrix-algorithm-param-outlier-trim");
	CPPUNIT_ASSERT(param.GetParameterType()==6);
    CPPUNIT_ASSERT(param.GetValueFloat()==3.0f);
	delete cn2;


   //create newObject with different parameters and test them 
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params;
	params["outlier-trim"]="0.111";	
	SelfCreate *sc=CNAnalysisMethodAllelicDifference::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj1->at(0);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-outlier-trim");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT(param1.GetValueFloat()==0.111f);
	
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNAnalysisMethodAllelicDifference::explainSelf();
	//CPPUNIT_ASSERT(sd1.getState()=="allelic-difference.outlier-trim=3.0.step=20.window=100.point-count=128.bandwidth=0.25.cutoff=0.05.clean-threshold=0.35.symmetry=true");
	CPPUNIT_ASSERT(sd1.getState()=="allelic-difference.outlier-trim=3.0");
    delete sc;	
}	

void CNAnalysisMethodAllelicDifferenceTest::runTest()
{
    Verbose::out(1, "****CNAnalysisMethodAllelicDifferenceTest::runTest****");

    
	CNAnalysisMethodAllelicDifference cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod allelic-difference is not setup properly.
	NEGATIVE_TEST(cnrfCNam.run(),Except);
}
