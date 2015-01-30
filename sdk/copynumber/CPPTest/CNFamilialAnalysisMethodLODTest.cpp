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

#include "copynumber/CNFamilialAnalysisMethodFactory.h"
#include "copynumber/CNFamilialAnalysisMethodLOD.h"
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
 * @class CNFamilialAnalysisMethodLODTest
 * @brief cppunit class for testing CNFamilialAnalysisMethodLOD functions.
 * last change by vliber on 11/03/09
  */

class CNFamilialAnalysisMethodLODTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNFamilialAnalysisMethodLODTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(alleleFrequencyFileTest); 
  CPPUNIT_TEST(runTest); 
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void alleleFrequencyFileTest();
  void runTest();
 
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNFamilialAnalysisMethodLODTest );

void CNFamilialAnalysisMethodLODTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNFamilialAnalysisMethodLODTest::functionsTest****");
	// constructor && m_vParams()
	CNFamilialAnalysisMethodLOD cnrfCNam;
	CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="paternity");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Familial LOD");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="paternity");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNFamilialAnalysisMethodLOD::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==2);
	CPPUNIT_ASSERT(sv[0].name=="role-validity-threshold");
	CPPUNIT_ASSERT(sv[0].type==2);
	CPPUNIT_ASSERT(sv[0].value=="1000.0");
	CPPUNIT_ASSERT(sv[0].defaultVal=="1000.0");
	CPPUNIT_ASSERT(sv[0].minVal=="0");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript==" Paternity Role Validity Threshold");
	CPPUNIT_ASSERT(sv[1].name=="error-rate");
	CPPUNIT_ASSERT(sv[1].type==2);
	CPPUNIT_ASSERT(sv[1].value=="0.01");
	CPPUNIT_ASSERT(sv[1].defaultVal=="0.01");
	CPPUNIT_ASSERT(sv[1].minVal=="0.0");
	CPPUNIT_ASSERT(sv[1].maxVal=="1.0");
	CPPUNIT_ASSERT(sv[1].descript==" Paternity Error Rate");
		
	//explainSelf()
	SelfDoc sd=CNFamilialAnalysisMethodLOD::explainSelf();
	//std::cout<<sd.getState()<<std::endl;
	CPPUNIT_ASSERT(sd.getState()=="paternity.role-validity-threshold=1000.0.error-rate=0.01");
	CPPUNIT_ASSERT(sd.getDocName()=="paternity");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Familial LOD");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==2);
	CPPUNIT_ASSERT(v[0].asString()=="1000.0");
	CPPUNIT_ASSERT(sd.getDocOption("role-validity-threshold").asString()=="1000.0");
	CPPUNIT_ASSERT(v[1].asString()=="0.01");
	CPPUNIT_ASSERT(sd.getDocOption("error-rate").asString()=="0.01");
}
void CNFamilialAnalysisMethodLODTest::newObjectTest()
{
	
	Verbose::out(1, "****CNFamilialAnalysisMethodLODTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNFamilialAnalysisMethod::getParams().clear();
	CNFamilialAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	POSITIVE_TEST((CNFamilialAnalysisMethodLOD*)obj_CNAnalysisMethodFactory.CNFamilialAnalysisMethodForString("paternity"));
	vector<affymetrix_calvin_parameter::ParameterNameValueType> obj1=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1.size()==2);
	affymetrix_calvin_parameter::ParameterNameValueType param = obj1.at(0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param.GetName())=="affymetrix-algorithm-param-role-validity-threshold");
	CPPUNIT_ASSERT(param.GetParameterType()==6);
  CPPUNIT_ASSERT(param.GetValueFloat()==1000.0);
	param = obj1.at(1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param.GetName())=="affymetrix-algorithm-param-error-rate");
	CPPUNIT_ASSERT(param.GetParameterType()==6);
    //CPPUNIT_ASSERT(param.GetValueFloat()==0.01);

   //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["role-validity-threshold"]="31.0";
	params["error-rate"]="0.7";
	CNFamilialAnalysisMethod::getParams().clear();
	SelfCreate *sc=CNFamilialAnalysisMethodLOD::newObject(params);
    vector<affymetrix_calvin_parameter::ParameterNameValueType> obj=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj.size()==2);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj.at(0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-role-validity-threshold");
	CPPUNIT_ASSERT(param1.GetParameterType()==6);
    CPPUNIT_ASSERT(param1.GetValueFloat()==31.0);
	param1 = obj.at(1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-error-rate");
	CPPUNIT_ASSERT(param1.GetParameterType()==6);
    //CPPUNIT_ASSERT(param1.GetValueFloat()==0.7);
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNFamilialAnalysisMethodLOD::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="paternity.role-validity-threshold=1000.0.error-rate=0.01");
    delete sc;
}


void CNFamilialAnalysisMethodLODTest::alleleFrequencyFileTest()
{
	Verbose::out(1, "****CNFamilialAnalysisMethodLODTest::alleleFrequencyFileTest****");
    CNFamilialAnalysisMethodLOD af;
    
	CPPUNIT_ASSERT(!af.readFile(INPUT+"/Cyto/alleleFrequencyFile_test_1.txt"));
	CPPUNIT_ASSERT(af.readFile(INPUT+"/Cyto/alleleFrequencyFile_test.txt"));
    CPPUNIT_ASSERT(af.getAlleleFreqSize()==3);
	CPPUNIT_ASSERT(af.getAlleleFrequency(0)->SNPID=="11111");
	CPPUNIT_ASSERT(af.getAlleleFrequency(1)->SNPID=="22222");
	CPPUNIT_ASSERT(af.getAlleleFrequency(2)->SNPID=="33333");
	CPPUNIT_ASSERT(af.getAlleleFrequency(0)->AAlleleFrequency==0.2);
	CPPUNIT_ASSERT(af.getAlleleFrequency(1)->AAlleleFrequency==0.3);
	CPPUNIT_ASSERT(af.getAlleleFrequency(2)->AAlleleFrequency==0.1);
    CPPUNIT_ASSERT(af.getAlleleFrequency(0)->BAlleleFrequency==0.5);
	CPPUNIT_ASSERT(af.getAlleleFrequency(1)->BAlleleFrequency==0.4);
	CPPUNIT_ASSERT(af.getAlleleFrequency(2)->BAlleleFrequency==0.2);
}
void CNFamilialAnalysisMethodLODTest::runTest()
{
	Verbose::out(1, "****CNFamilialAnalysisMethodLODTest::runTest****");
	CNFamilialAnalysisMethodLOD m_obj;
	//Message: Cannot open allele frequency file:
    NEGATIVE_TEST(m_obj.run(),Except);
}
