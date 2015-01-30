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

#include "copynumber/CNAnalysisMethodCN.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNReferenceMethodWaveCorrection.h"
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
 * @class CNReferenceMethodWaveCorrectionTest
 * @brief cppunit class for testing CNReferenceMethodWaveCorrection functions.
 * last change by vliber on 07/01/09
 */

class CNReferenceMethodWaveCorrectionTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNReferenceMethodWaveCorrectionTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNReferenceMethodWaveCorrectionTest);

void CNReferenceMethodWaveCorrectionTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNReferenceMethodWaveCorrectionTest::functionsTest****");
	// constructor && m_vParams()
	CNReferenceMethodWaveCorrection cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="wave-correction-reference-method");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="CopyNumber WaveCorrection");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="wave-correction-reference-method");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNReferenceMethodWaveCorrection::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==4);
	CPPUNIT_ASSERT(sv[0].name=="trim");
	CPPUNIT_ASSERT(sv[0].type==1);
	CPPUNIT_ASSERT(sv[0].value=="2.0");
	CPPUNIT_ASSERT(sv[0].defaultVal=="2.0");
	CPPUNIT_ASSERT(sv[0].minVal=="NA");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="Log2Ratio Trim value.");
	CPPUNIT_ASSERT(sv[1].name=="percentile");
	CPPUNIT_ASSERT(sv[1].type==1);
	CPPUNIT_ASSERT(sv[1].value=="0.75");
	CPPUNIT_ASSERT(sv[1].defaultVal=="0.75");
	CPPUNIT_ASSERT(sv[1].minVal=="NA");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="High Percentile value.");
	CPPUNIT_ASSERT(sv[2].name=="wave-count");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="6");
	CPPUNIT_ASSERT(sv[2].defaultVal=="6");
	CPPUNIT_ASSERT(sv[2].minVal=="0");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="Number of waves to add to the reference.");
	CPPUNIT_ASSERT(sv[3].name=="demean");
	CPPUNIT_ASSERT(sv[3].type==4);
	CPPUNIT_ASSERT(sv[3].value=="false");
	CPPUNIT_ASSERT(sv[3].defaultVal=="false");
	CPPUNIT_ASSERT(sv[3].minVal=="NA");
	CPPUNIT_ASSERT(sv[3].maxVal=="NA");
	CPPUNIT_ASSERT(sv[3].descript=="Demean the input to the SVD.");	

	//explainSelf()
	SelfDoc sd=CNReferenceMethodWaveCorrection::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="wave-correction-reference-method.trim=2.0.percentile=0.75.wave-count=6.demean=false");
    CPPUNIT_ASSERT(sd.getDocName()=="wave-correction-reference-method");
	CPPUNIT_ASSERT(sd.getDocDescription()=="CopyNumber WaveCorrection");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==4);
	CPPUNIT_ASSERT(v[0].asString()=="2.0");
	CPPUNIT_ASSERT(v[1].asString()=="0.75");
	CPPUNIT_ASSERT(v[2].asString()=="6");
	CPPUNIT_ASSERT(v[3].asString()=="false");
	CPPUNIT_ASSERT(sd.getDocOption("trim").asString()=="2.0");
	CPPUNIT_ASSERT(sd.getDocOption("percentile").asString()=="0.75");
	CPPUNIT_ASSERT(sd.getDocOption("wave-count").asString()=="6");
	CPPUNIT_ASSERT(sd.getDocOption("demean").asString()=="false");
}

void CNReferenceMethodWaveCorrectionTest::newObjectTest()
{
	
	Verbose::out(1, "****CNReferenceMethodWaveCorrectionTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNReferenceMethodWaveCorrection *cn2=(CNReferenceMethodWaveCorrection*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("wave-correction-reference-method");
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==4);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(3);
	
		
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-trim");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1a.GetValueFloat()==2.0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-percentile");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2a.GetValueFloat()==0.75);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-wave-count");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==6);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-demean");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param4a.GetValueInt8()==0);
	delete cn2;

    //create newObject with different parameters and test them
	CNAnalysisMethod::getParams()->clear();
	std::map<std::string,std::string> params1;
	params1["trim"]="5.0";
	params1["percentile"]="0.100";
	params1["wave-count"]="1";
	params1["demean"]="true";
	SelfCreate *sc=CNReferenceMethodWaveCorrection::newObject(params1);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==4);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(3);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-trim");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1.GetValueFloat()==5.0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-percentile");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2.GetValueFloat()==0.100f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-wave-count");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3.GetValueInt32()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-demean");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param4.GetValueInt8()==1);
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNReferenceMethodWaveCorrection::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="wave-correction-reference-method.trim=2.0.percentile=0.75.wave-count=6.demean=false");
    delete sc;
	
}	
