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
#include "copynumber/CNFamilialAnalysisMethodGenotypeConcordance.h"
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
 * @class CNFamilialAnalysisMethodGenotypeConcordanceTest
 * @brief cppunit class for testing CNFamilialAnalysisMethodGenotypeConcordance functions.
 * last change by vliber on 10/08/09
  */

class CNFamilialAnalysisMethodGenotypeConcordanceTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNFamilialAnalysisMethodGenotypeConcordanceTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
 
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNFamilialAnalysisMethodGenotypeConcordanceTest );

void CNFamilialAnalysisMethodGenotypeConcordanceTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNFamilialAnalysisMethodGenotypeConcordanceTest::functionsTest****");
	// constructor && m_vParams()
	CNFamilialAnalysisMethodGenotypeConcordance cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="genotype-concordance");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Familial GenotypeConcordance");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="genotype-concordance");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNFamilialAnalysisMethodGenotypeConcordance::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==3);
	CPPUNIT_ASSERT(sv[0].name=="genotype-concordance-marker-count-cutoff");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="21");
	CPPUNIT_ASSERT(sv[0].defaultVal=="21");
	CPPUNIT_ASSERT(sv[0].minVal=="1");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="GenotypeConcordance Marker Count Cutoff");
	CPPUNIT_ASSERT(sv[1].name=="genotype-concordance-call-cutoff");
	CPPUNIT_ASSERT(sv[1].type==3);
	CPPUNIT_ASSERT(sv[1].value=="4");
	CPPUNIT_ASSERT(sv[1].defaultVal=="4");
	CPPUNIT_ASSERT(sv[1].minVal=="0");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="GenotypeConcordance Call Cutoff");
	CPPUNIT_ASSERT(sv[2].name=="genotype-concordance-min-genomic-span");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="1000000");
	CPPUNIT_ASSERT(sv[2].defaultVal=="1000000");
	CPPUNIT_ASSERT(sv[2].minVal=="1");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="GenotypeConcordance Minimum Genomic Span");
	
	
	//explainSelf()
	SelfDoc sd=CNFamilialAnalysisMethodGenotypeConcordance::explainSelf();
	//std::cout<<sd.getState()<<std::endl;
	CPPUNIT_ASSERT(sd.getState()=="genotype-concordance.genotype-concordance-marker-count-cutoff=21.genotype-concordance-call-cutoff=4.genotype-concordance-min-genomic-span=1000000");
	CPPUNIT_ASSERT(sd.getDocName()=="genotype-concordance");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Familial GenotypeConcordance");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==3);
	CPPUNIT_ASSERT(v[0].asString()=="21");
	CPPUNIT_ASSERT(v[1].asString()=="4");
	CPPUNIT_ASSERT(v[2].asString()=="1000000");
	CPPUNIT_ASSERT(sd.getDocOption("genotype-concordance-marker-count-cutoff").asString()=="21");
	CPPUNIT_ASSERT(sd.getDocOption("genotype-concordance-call-cutoff").asString()=="4");
	CPPUNIT_ASSERT(sd.getDocOption("genotype-concordance-min-genomic-span").asString()=="1000000");
}
void CNFamilialAnalysisMethodGenotypeConcordanceTest::newObjectTest()
{
	
	Verbose::out(1, "****CNFamilialAnalysisMethodGenotypeConcordanceTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNFamilialAnalysisMethod::getParams().clear();
	CNFamilialAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNFamilialAnalysisMethodGenotypeConcordance *cn2=(CNFamilialAnalysisMethodGenotypeConcordance*)obj_CNAnalysisMethodFactory.CNFamilialAnalysisMethodForString("genotype-concordance");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> obj1=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1.size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1.at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1.at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1.at(2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-genotype-concordance-marker-count-cutoff");
	CPPUNIT_ASSERT(param1a.GetParameterType()==4);
    CPPUNIT_ASSERT(param1a.GetValueInt32()==21);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-genotype-concordance-call-cutoff");
	CPPUNIT_ASSERT(param2a.GetParameterType()==4);
	CPPUNIT_ASSERT(param2a.GetValueInt32()==4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-genotype-concordance-min-genomic-span");
	CPPUNIT_ASSERT(param3a.GetParameterType()==4);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==1000000);
	delete cn2;


   //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["genotype-concordance-marker-count-cutoff"]="31";
	params["genotype-concordance-call-cutoff"]="14";
	params["genotype-concordance-min-genomic-span"]="100";
	CNFamilialAnalysisMethod::getParams().clear();
	SelfCreate *sc=CNFamilialAnalysisMethodGenotypeConcordance::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> obj=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj.size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj.at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj.at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj.at(2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-genotype-concordance-marker-count-cutoff");
	CPPUNIT_ASSERT(param1.GetParameterType()==4);
    CPPUNIT_ASSERT(param1.GetValueInt32()==31);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-genotype-concordance-call-cutoff");
	CPPUNIT_ASSERT(param2.GetParameterType()==4);
	CPPUNIT_ASSERT(param2.GetValueInt32()==14);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-genotype-concordance-min-genomic-span");
	CPPUNIT_ASSERT(param3.GetParameterType()==4);
	CPPUNIT_ASSERT(param3.GetValueInt32()==100);
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNFamilialAnalysisMethodGenotypeConcordance::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="genotype-concordance.genotype-concordance-marker-count-cutoff=21.genotype-concordance-call-cutoff=4.genotype-concordance-min-genomic-span=1000000");
	delete sc;	
}

