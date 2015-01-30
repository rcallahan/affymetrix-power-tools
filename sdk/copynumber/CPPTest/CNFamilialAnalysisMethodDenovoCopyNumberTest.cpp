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
#include "copynumber/CNFamilialAnalysisMethodDenovoCopyNumber.h"
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
 * @class CNFamilialAnalysisMethodDenovoCopyNumberTest
 * @brief cppunit class for testing CNFamilialAnalysisMethodDenovoCopyNumber functions.
 * last change by vliber on 10/08/09
  */

class CNFamilialAnalysisMethodDenovoCopyNumberTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNFamilialAnalysisMethodDenovoCopyNumberTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
 
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNFamilialAnalysisMethodDenovoCopyNumberTest );

void CNFamilialAnalysisMethodDenovoCopyNumberTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNFamilialAnalysisMethodDenovoCopyNumberTest::functionsTest****");
	// constructor && m_vParams()
	CNFamilialAnalysisMethodDenovoCopyNumber cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="denovo-cn");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Familial DenovoCopyNumber");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="denovo-cn");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNFamilialAnalysisMethodDenovoCopyNumber::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==3);
	CPPUNIT_ASSERT(sv[0].name=="denovo-cn-marker-count-cutoff");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="21");
	CPPUNIT_ASSERT(sv[0].defaultVal=="21");
	CPPUNIT_ASSERT(sv[0].minVal=="1");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="DenovoCopyNumber Marker Count Cutoff");
	CPPUNIT_ASSERT(sv[1].name=="denovo-cn-call-cutoff");
	CPPUNIT_ASSERT(sv[1].type==3);
	CPPUNIT_ASSERT(sv[1].value=="4");
	CPPUNIT_ASSERT(sv[1].defaultVal=="4");
	CPPUNIT_ASSERT(sv[1].minVal=="0");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="DenovoCopyNumber Call Cutoff");
	CPPUNIT_ASSERT(sv[2].name=="denovo-cn-min-genomic-span");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="1000000");
	CPPUNIT_ASSERT(sv[2].defaultVal=="1000000");
	CPPUNIT_ASSERT(sv[2].minVal=="1");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="DenovoCopyNumber Minimum Genomic Span");
	
	
	//explainSelf()
	SelfDoc sd=CNFamilialAnalysisMethodDenovoCopyNumber::explainSelf();
	//std::cout<<sd.getState()<<std::endl;
	CPPUNIT_ASSERT(sd.getState()=="denovo-cn.denovo-cn-marker-count-cutoff=21.denovo-cn-call-cutoff=4.denovo-cn-min-genomic-span=1000000");
	CPPUNIT_ASSERT(sd.getDocName()=="denovo-cn");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Familial DenovoCopyNumber");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==3);
	CPPUNIT_ASSERT(v[0].asString()=="21");
	CPPUNIT_ASSERT(v[1].asString()=="4");
	CPPUNIT_ASSERT(v[2].asString()=="1000000");
	CPPUNIT_ASSERT(sd.getDocOption("denovo-cn-marker-count-cutoff").asString()=="21");
	CPPUNIT_ASSERT(sd.getDocOption("denovo-cn-call-cutoff").asString()=="4");
	CPPUNIT_ASSERT(sd.getDocOption("denovo-cn-min-genomic-span").asString()=="1000000");
}
void CNFamilialAnalysisMethodDenovoCopyNumberTest::newObjectTest()
{
	
	Verbose::out(1, "****CNFamilialAnalysisMethodDenovoCopyNumberTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNFamilialAnalysisMethod::getParams().clear();
	CNFamilialAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNFamilialAnalysisMethodDenovoCopyNumber *cn2=(CNFamilialAnalysisMethodDenovoCopyNumber*)obj_CNAnalysisMethodFactory.CNFamilialAnalysisMethodForString("denovo-cn");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> obj1=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1.size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1.at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1.at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1.at(2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-denovo-cn-marker-count-cutoff");
	CPPUNIT_ASSERT(param1a.GetParameterType()==4);
    CPPUNIT_ASSERT(param1a.GetValueInt32()==21);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-denovo-cn-call-cutoff");
	CPPUNIT_ASSERT(param2a.GetParameterType()==4);
	CPPUNIT_ASSERT(param2a.GetValueInt32()==4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-denovo-cn-min-genomic-span");
	CPPUNIT_ASSERT(param3a.GetParameterType()==4);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==1000000);
	delete cn2;


   //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["denovo-cn-marker-count-cutoff"]="11";
	params["denovo-cn-call-cutoff"]="44";
	params["denovo-cn-min-genomic-span"]="500";
	CNFamilialAnalysisMethod::getParams().clear();
	SelfCreate *sc=CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> obj=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj.size()==3);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj.at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj.at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj.at(2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-denovo-cn-marker-count-cutoff");
	CPPUNIT_ASSERT(param1.GetParameterType()==4);
    CPPUNIT_ASSERT(param1.GetValueInt32()==11);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-denovo-cn-call-cutoff");
	CPPUNIT_ASSERT(param2.GetParameterType()==4);
	CPPUNIT_ASSERT(param2.GetValueInt32()==44);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-denovo-cn-min-genomic-span");
	CPPUNIT_ASSERT(param3.GetParameterType()==4);
	CPPUNIT_ASSERT(param3.GetValueInt32()==500);
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNFamilialAnalysisMethodDenovoCopyNumber::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="denovo-cn.denovo-cn-marker-count-cutoff=21.denovo-cn-call-cutoff=4.denovo-cn-min-genomic-span=1000000");
	delete sc;	
    //parameter tests
	params.clear();
	params["denovo-cn-marker-count-cutoff"]="11";
	params["denovo-cn-call-cutoff"]="44";
	params["denovo-cn-min-genomic-span"]="500";
	POSITIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params));
	//SelfCreate::setValue() - '0' is not a valid value for parameter: 'denovo-cn-marker-count-cutoff'. The specified range is 1 to NA
	params["denovo-cn-marker-count-cutoff"]="0";
	NEGATIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params),Except);
	params["denovo-cn-marker-count-cutoff"]="1";
	POSITIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params));
	params["denovo-cn-marker-count-cutoff"]="11111";
	POSITIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params));
	//SelfCreate::setValue() - '-1' is not a valid value for parameter: 'denovo-cn-call-cutoff'. The specified range is 0 to NA
	params["denovo-cn-call-cutoff"]="-1";
	NEGATIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params),Except);
	params["denovo-cn-call-cutoff"]="0";
	POSITIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params));
	params["denovo-cn-call-cutoff"]="11111";
	POSITIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params));
	//SelfCreate::setValue() - '0' is not a valid value for parameter: 'denovo-cn-min-genomic-span'. The specified range is 1 to NA
	params["denovo-cn-min-genomic-span"]="0";
	NEGATIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params),Except);
	params["denovo-cn-min-genomic-span"]="1";
	POSITIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params));
	params["denovo-cn-min-genomic-span"]="11111";
	POSITIVE_TEST(CNFamilialAnalysisMethodDenovoCopyNumber::newObject(params));
}

