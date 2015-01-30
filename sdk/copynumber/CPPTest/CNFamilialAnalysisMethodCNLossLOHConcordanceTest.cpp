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
#include "copynumber/CNFamilialAnalysisMethodCNLossLOHConcordance.h"
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
 * @class CNFamilialAnalysisMethodCNLossLOHConcordanceTest
 * @brief cppunit class for testing CNFamilialAnalysisMethodCNLossLOHConcordance functions.
 * last change by vliber on 10/07/09
  */

class CNFamilialAnalysisMethodCNLossLOHConcordanceTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNFamilialAnalysisMethodCNLossLOHConcordanceTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST_SUITE_END();
  
public:  
  void functionsTest();
  void newObjectTest();
 
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNFamilialAnalysisMethodCNLossLOHConcordanceTest );

void CNFamilialAnalysisMethodCNLossLOHConcordanceTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNFamilialAnalysisMethodCNLossLOHConcordanceTest::functionsTest****");
	// constructor && m_vParams()
	CNFamilialAnalysisMethodCNLossLOHConcordance cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="cn-loss-loh-concordance");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Familial CNLossLOHConcordance");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="cn-loss-loh-concordance");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNFamilialAnalysisMethodCNLossLOHConcordance::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==5);
	CPPUNIT_ASSERT(sv[0].name=="cn-loss-loh-concordance-marker-count-cutoff");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="21");
	CPPUNIT_ASSERT(sv[0].defaultVal=="21");
	CPPUNIT_ASSERT(sv[0].minVal=="1");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="CNLossLOHConcordance Marker Count Cutoff");
	CPPUNIT_ASSERT(sv[1].name=="cn-loss-loh-concordance-call-cutoff");
	CPPUNIT_ASSERT(sv[1].type==3);
	CPPUNIT_ASSERT(sv[1].value=="4");
	CPPUNIT_ASSERT(sv[1].defaultVal=="4");
	CPPUNIT_ASSERT(sv[1].minVal=="0");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="CNLossLOHConcordance Call Cutoff");
	CPPUNIT_ASSERT(sv[2].name=="cn-loss-loh-concordance-min-genomic-span");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="1000000");
	CPPUNIT_ASSERT(sv[2].defaultVal=="1000000");
	CPPUNIT_ASSERT(sv[2].minVal=="1");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="CNLossLOHConcordance Minimum Genomic Span");
	CPPUNIT_ASSERT(sv[3].name=="cn-loss-loh-concordance-xChromosome");
	CPPUNIT_ASSERT(sv[3].type==3);
	CPPUNIT_ASSERT(sv[3].value=="24");
	CPPUNIT_ASSERT(sv[3].defaultVal=="24");
	CPPUNIT_ASSERT(sv[3].minVal=="1");
	CPPUNIT_ASSERT(sv[3].maxVal=="NA");
	CPPUNIT_ASSERT(sv[3].descript=="CNLossLOHConcordance X Chromosome");
	CPPUNIT_ASSERT(sv[4].name=="cn-loss-loh-concordance-yChromosome");
	CPPUNIT_ASSERT(sv[4].type==3);
	CPPUNIT_ASSERT(sv[4].value=="25");
	CPPUNIT_ASSERT(sv[4].defaultVal=="25");
	CPPUNIT_ASSERT(sv[4].minVal=="1");
	CPPUNIT_ASSERT(sv[4].maxVal=="NA");
	CPPUNIT_ASSERT(sv[4].descript=="CNLossLOHConcordance Y Chromosome");
	
	//explainSelf()
	SelfDoc sd=CNFamilialAnalysisMethodCNLossLOHConcordance::explainSelf();
	//std::cout<<sd.getState()<<std::endl;
	CPPUNIT_ASSERT(sd.getState()=="cn-loss-loh-concordance.cn-loss-loh-concordance-marker-count-cutoff=21.cn-loss-loh-concordance-call-cutoff=4.cn-loss-loh-concordance-min-genomic-span=1000000.cn-loss-loh-concordance-xChromosome=24.cn-loss-loh-concordance-yChromosome=25");
	CPPUNIT_ASSERT(sd.getDocName()=="cn-loss-loh-concordance");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Familial CNLossLOHConcordance");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==5);
	CPPUNIT_ASSERT(v[0].asString()=="21");
	CPPUNIT_ASSERT(v[1].asString()=="4");
	CPPUNIT_ASSERT(v[2].asString()=="1000000");
	CPPUNIT_ASSERT(v[3].asString()=="24");
	CPPUNIT_ASSERT(v[4].asString()=="25");
	CPPUNIT_ASSERT(sd.getDocOption("cn-loss-loh-concordance-marker-count-cutoff").asString()=="21");
	CPPUNIT_ASSERT(sd.getDocOption("cn-loss-loh-concordance-call-cutoff").asString()=="4");
	CPPUNIT_ASSERT(sd.getDocOption("cn-loss-loh-concordance-min-genomic-span").asString()=="1000000");
	CPPUNIT_ASSERT(sd.getDocOption("cn-loss-loh-concordance-xChromosome").asString()=="24");
	CPPUNIT_ASSERT(sd.getDocOption("cn-loss-loh-concordance-yChromosome").asString()=="25");
	
}
void CNFamilialAnalysisMethodCNLossLOHConcordanceTest::newObjectTest()
{
	
	Verbose::out(1, "****CNFamilialAnalysisMethodCNLossLOHConcordanceTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNFamilialAnalysisMethod::getParams().clear();
	CNFamilialAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNFamilialAnalysisMethodCNLossLOHConcordance *cn2=(CNFamilialAnalysisMethodCNLossLOHConcordance*)obj_CNAnalysisMethodFactory.CNFamilialAnalysisMethodForString("cn-loss-loh-concordance");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> obj1=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1.size()==5);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1.at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1.at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1.at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1.at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5a = obj1.at(4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-marker-count-cutoff");
	CPPUNIT_ASSERT(param1a.GetParameterType()==4);
    CPPUNIT_ASSERT(param1a.GetValueInt32()==21);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-call-cutoff");
	CPPUNIT_ASSERT(param2a.GetParameterType()==4);
	CPPUNIT_ASSERT(param2a.GetValueInt32()==4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-min-genomic-span");
	CPPUNIT_ASSERT(param3a.GetParameterType()==4);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==1000000);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-xChromosome");
	CPPUNIT_ASSERT(param4a.GetParameterType()==4);
	CPPUNIT_ASSERT(param4a.GetValueInt32()==24);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5a.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-yChromosome");
	CPPUNIT_ASSERT(param5a.GetParameterType()==4);
	CPPUNIT_ASSERT(param5a.GetValueInt32()==25);
	delete cn2;


   //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["cn-loss-loh-concordance-marker-count-cutoff"]="31";
	params["cn-loss-loh-concordance-call-cutoff"]="14";
	params["cn-loss-loh-concordance-min-genomic-span"]="100";
	params["cn-loss-loh-concordance-xChromosome"]="12";
	params["cn-loss-loh-concordance-yChromosome"]="13";
	CNFamilialAnalysisMethod::getParams().clear();
	SelfCreate *sc=CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> obj=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj.size()==5);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj.at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj.at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj.at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj.at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5 = obj.at(4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-marker-count-cutoff");
	CPPUNIT_ASSERT(param1.GetParameterType()==4);
    CPPUNIT_ASSERT(param1.GetValueInt32()==31);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-call-cutoff");
	CPPUNIT_ASSERT(param2.GetParameterType()==4);
	CPPUNIT_ASSERT(param2.GetValueInt32()==14);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-min-genomic-span");
	CPPUNIT_ASSERT(param3.GetParameterType()==4);
	CPPUNIT_ASSERT(param3.GetValueInt32()==100);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-xChromosome");
	CPPUNIT_ASSERT(param4.GetParameterType()==4);
	CPPUNIT_ASSERT(param4.GetValueInt32()==12);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5.GetName())=="affymetrix-algorithm-param-cn-loss-loh-concordance-yChromosome");
	CPPUNIT_ASSERT(param5.GetParameterType()==4);
	CPPUNIT_ASSERT(param5.GetValueInt32()==13);
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNFamilialAnalysisMethodCNLossLOHConcordance::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="cn-loss-loh-concordance.cn-loss-loh-concordance-marker-count-cutoff=21.cn-loss-loh-concordance-call-cutoff=4.cn-loss-loh-concordance-min-genomic-span=1000000.cn-loss-loh-concordance-xChromosome=24.cn-loss-loh-concordance-yChromosome=25");
    delete sc;	
    //parameter tests
	params.clear();
	params["cn-loss-loh-concordance-marker-count-cutoff"]="31";
	params["cn-loss-loh-concordance-call-cutoff"]="14";
	params["cn-loss-loh-concordance-min-genomic-span"]="100";
	params["cn-loss-loh-concordance-xChromosome"]="12";
	params["cn-loss-loh-concordance-yChromosome"]="13";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
	//SelfCreate::setValue() - '0' is not a valid value for parameter: 'cn-loss-loh-concordance-marker-count-cutoff'. The specified range is 1 to NA
	params["cn-loss-loh-concordance-marker-count-cutoff"]="0";
	NEGATIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params),Except);
	params["cn-loss-loh-concordance-marker-count-cutoff"]="1";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
	params["cn-loss-loh-concordance-marker-count-cutoff"]="11111";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
	//SelfCreate::setValue() - '-14' is not a valid value for parameter: 'cn-loss-loh-concordance-call-cutoff'. The specified range is 0 to NA
	params["cn-loss-loh-concordance-call-cutoff"]="-14";
	NEGATIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params),Except);
	params["cn-loss-loh-concordance-call-cutoff"]="0";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
	params["cn-loss-loh-concordance-call-cutoff"]="11110";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
    //SelfCreate::setValue() - '0' is not a valid value for parameter: 'cn-loss-loh-concordance-min-genomic-span'. The specified range is 1 to NA
	params["cn-loss-loh-concordance-min-genomic-span"]="0";
	NEGATIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params),Except);
	params["cn-loss-loh-concordance-min-genomic-span"]="1";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
	params["cn-loss-loh-concordance-min-genomic-span"]="11110";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
	//SelfCreate::setValue() - '0' is not a valid value for parameter: 'cn-loss-loh-concordance-xChromosome'. The specified range is 1 to NA
	params["cn-loss-loh-concordance-xChromosome"]="0";
	NEGATIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params),Except);
	params["cn-loss-loh-concordance-xChromosome"]="1";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
	params["cn-loss-loh-concordance-xChromosome"]="11110";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
	//SelfCreate::setValue() - '0' is not a valid value for parameter: 'cn-loss-loh-concordance-yChromosome'. The specified range is 1 to NA
	params["cn-loss-loh-concordance-yChromosome"]="0";
	NEGATIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params),Except);
	params["cn-loss-loh-concordance-yChromosome"]="1";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));
	params["cn-loss-loh-concordance-yChromosome"]="11110";
	POSITIVE_TEST(CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(params));

}

