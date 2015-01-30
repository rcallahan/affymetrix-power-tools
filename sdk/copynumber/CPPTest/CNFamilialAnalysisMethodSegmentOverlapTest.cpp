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
#include "copynumber/CNFamilialAnalysisMethodSegmentOverlap.h"
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
 * @class CNFamilialAnalysisMethodSegmentOverlapTest
 * @brief cppunit class for testing CNFamilialAnalysisMethodSegmentOverlap functions.
 * last change by vliber on 10/08/09
  */

class CNFamilialAnalysisMethodSegmentOverlapTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNFamilialAnalysisMethodSegmentOverlapTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
 
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNFamilialAnalysisMethodSegmentOverlapTest );

void CNFamilialAnalysisMethodSegmentOverlapTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNFamilialAnalysisMethodSegmentOverlapTest::functionsTest****");
	// constructor && m_vParams()
	CNFamilialAnalysisMethodSegmentOverlap cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="segment-overlap");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Familial SegmentOverlap");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="segment-overlap");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNFamilialAnalysisMethodSegmentOverlap::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==1);
	CPPUNIT_ASSERT(sv[0].name=="segment-overlap-threshold");
	CPPUNIT_ASSERT(sv[0].type==2);
	CPPUNIT_ASSERT(sv[0].value=="0.75");
	CPPUNIT_ASSERT(sv[0].defaultVal=="0.75");
	CPPUNIT_ASSERT(sv[0].minVal=="0");
	CPPUNIT_ASSERT(sv[0].maxVal=="1");
	CPPUNIT_ASSERT(sv[0].descript=="SegmentOverlap threshold");
	
	//explainSelf()
	SelfDoc sd=CNFamilialAnalysisMethodSegmentOverlap::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="segment-overlap.segment-overlap-threshold=0.75");
	CPPUNIT_ASSERT(sd.getDocName()=="segment-overlap");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Familial SegmentOverlap");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==1);
	CPPUNIT_ASSERT(v[0].asString()=="0.75");
	CPPUNIT_ASSERT(sd.getDocOption("segment-overlap-threshold").asString()=="0.75");
	
	
}
void CNFamilialAnalysisMethodSegmentOverlapTest::newObjectTest()
{
	
	Verbose::out(1, "****CNFamilialAnalysisMethodSegmentOverlapTest::newObjectTest****");
    
    //create default object from factory by string and verify default parameters
	CNFamilialAnalysisMethod::getParams().clear();
	CNFamilialAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNFamilialAnalysisMethodSegmentOverlap *cn2=(CNFamilialAnalysisMethodSegmentOverlap*)obj_CNAnalysisMethodFactory.CNFamilialAnalysisMethodForString("segment-overlap");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> obj1=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1.size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1.at(0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-segment-overlap-threshold");
	CPPUNIT_ASSERT(param1a.GetParameterType()==6);
    CPPUNIT_ASSERT(param1a.GetValueFloat()==0.75f);
	delete cn2;


   //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["segment-overlap-threshold"]="0.55";
	CNFamilialAnalysisMethod::getParams().clear();
	SelfCreate *sc=CNFamilialAnalysisMethodSegmentOverlap::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> obj=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj.size()==1);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj.at(0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-segment-overlap-threshold");
	CPPUNIT_ASSERT(param1.GetParameterType()==6);
    CPPUNIT_ASSERT(param1.GetValueFloat()==0.55f);
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNFamilialAnalysisMethodSegmentOverlap::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="segment-overlap.segment-overlap-threshold=0.75");
    delete sc;	
	
}

