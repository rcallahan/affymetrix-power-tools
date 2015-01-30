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
#include "copynumber/CNReferenceMethodAdditionalWaves.h"
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
 * @class CNReferenceMethodAdditionalWavesTest
 * @brief cppunit class for testing CNReferenceMethodAdditionalWaves functions.
 * last change by vliber on 10/28/09
 */

class CNReferenceMethodAdditionalWavesTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNReferenceMethodAdditionalWavesTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
 
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNReferenceMethodAdditionalWavesTest );

void CNReferenceMethodAdditionalWavesTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNReferenceMethodAdditionalWavesTest::functionsTest****");
	// constructor && m_vParams()
	CNReferenceMethodAdditionalWaves cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="additional-waves-reference-method");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="CopyNumber AdditionalWaves");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="wave-correction-reference-method");
	
	vector<SelfDoc::Opt> sv=CNReferenceMethodAdditionalWaves::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==11);
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
	CPPUNIT_ASSERT(sv[2].name=="additional-wave-count");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="1");
	CPPUNIT_ASSERT(sv[2].defaultVal=="1");
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
	CPPUNIT_ASSERT(sv[4].name=="cn-qc-cutoff");
	CPPUNIT_ASSERT(sv[4].type==1);
	CPPUNIT_ASSERT(sv[4].value=="0.27");
	CPPUNIT_ASSERT(sv[4].defaultVal=="0.27");
	CPPUNIT_ASSERT(sv[4].minVal=="NA");
	CPPUNIT_ASSERT(sv[4].maxVal=="NA");
	CPPUNIT_ASSERT(sv[4].descript=="If the CN QC values if over this cutoff then the sample fails QC.");
	CPPUNIT_ASSERT(sv[5].name=="snp-qc-cutoff");
	CPPUNIT_ASSERT(sv[5].type==1);
	CPPUNIT_ASSERT(sv[5].value=="1.1");
	CPPUNIT_ASSERT(sv[5].defaultVal=="1.1");
	CPPUNIT_ASSERT(sv[5].minVal=="NA");
	CPPUNIT_ASSERT(sv[5].maxVal=="NA");
	CPPUNIT_ASSERT(sv[5].descript=="If the SNP QC values if below this cutoff then the sample fails QC.");
	CPPUNIT_ASSERT(sv[6].name=="waviness-seg-count-cutoff");
	CPPUNIT_ASSERT(sv[6].type==3);
	CPPUNIT_ASSERT(sv[6].value=="100");
	CPPUNIT_ASSERT(sv[6].defaultVal=="100");
	CPPUNIT_ASSERT(sv[6].minVal=="0");
	CPPUNIT_ASSERT(sv[6].maxVal=="NA");
	CPPUNIT_ASSERT(sv[6].descript=="The waviness seg count cutoff.");
	CPPUNIT_ASSERT(sv[7].name=="use-high-waviness-seg-count");
	CPPUNIT_ASSERT(sv[7].type==4);
	CPPUNIT_ASSERT(sv[7].value=="true");
	CPPUNIT_ASSERT(sv[7].defaultVal=="true");
	CPPUNIT_ASSERT(sv[7].minVal=="NA");
	CPPUNIT_ASSERT(sv[7].maxVal=="NA");
	CPPUNIT_ASSERT(sv[7].descript=="Use only those cychp files that have a waviness-seg-count > the cutoff if true, else use only those cychp files <= the cutoff.");
	CPPUNIT_ASSERT(sv[8].name=="force");
	CPPUNIT_ASSERT(sv[8].type==4);
	CPPUNIT_ASSERT(sv[8].value=="false");
	CPPUNIT_ASSERT(sv[8].defaultVal=="false");
	CPPUNIT_ASSERT(sv[8].minVal=="NA");
	CPPUNIT_ASSERT(sv[8].maxVal=="NA");
	CPPUNIT_ASSERT(sv[8].descript=="Force the job to run even if there is a mismatch between the cychp files and the input CN reference.");
	CPPUNIT_ASSERT(sv[9].name=="keep-temp-data");
	CPPUNIT_ASSERT(sv[9].type==4);
	CPPUNIT_ASSERT(sv[9].value=="false");
	CPPUNIT_ASSERT(sv[9].defaultVal=="false");
	CPPUNIT_ASSERT(sv[9].minVal=="NA");
	CPPUNIT_ASSERT(sv[9].maxVal=="NA");
	CPPUNIT_ASSERT(sv[9].descript=="If true, then do not delete the temporary data files used by the module.");
	CPPUNIT_ASSERT(sv[10].type==0);
	CPPUNIT_ASSERT(sv[10].value=="snp-qc");
	CPPUNIT_ASSERT(sv[10].defaultVal=="snp-qc");
	CPPUNIT_ASSERT(sv[10].minVal=="");
	CPPUNIT_ASSERT(sv[10].maxVal=="");
	CPPUNIT_ASSERT(sv[10].descript=="Choose the QC option to use.  Available values: snp-qc, contrast-qc, and contrast-qc-nsp");
   
	
	//explainSelf()
	SelfDoc sd=CNReferenceMethodAdditionalWaves::explainSelf();
	CPPUNIT_ASSERT(sd.getState()=="additional-waves-reference-method.trim=2.0.percentile=0.75.additional-wave-count=1.demean=false.cn-qc-cutoff=0.27.snp-qc-cutoff=1.1.waviness-seg-count-cutoff=100.use-high-waviness-seg-count=true.force=false.keep-temp-data=false.selected-qc=snp-qc");
    CPPUNIT_ASSERT(sd.getDocName()=="additional-waves-reference-method");
	CPPUNIT_ASSERT(sd.getDocDescription()=="CopyNumber AdditionalWaves");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==11);
	CPPUNIT_ASSERT(v[0].asString()=="2.0");
	CPPUNIT_ASSERT(v[1].asString()=="0.75");
	CPPUNIT_ASSERT(v[2].asString()=="1");
	CPPUNIT_ASSERT(v[3].asString()=="false");
	CPPUNIT_ASSERT(v[4].asString()=="0.27");
	CPPUNIT_ASSERT(v[5].asString()=="1.1");
	CPPUNIT_ASSERT(v[6].asString()=="100");
	CPPUNIT_ASSERT(v[7].asString()=="true");
	CPPUNIT_ASSERT(v[8].asString()=="false");
	CPPUNIT_ASSERT(v[9].asString()=="false");
	CPPUNIT_ASSERT(v[10].asString()=="snp-qc");
	
	CPPUNIT_ASSERT(sd.getDocOption("trim").asString()=="2.0");
	CPPUNIT_ASSERT(sd.getDocOption("percentile").asString()=="0.75");
	CPPUNIT_ASSERT(sd.getDocOption("additional-wave-count").asString()=="1");
	CPPUNIT_ASSERT(sd.getDocOption("demean").asString()=="false");
	CPPUNIT_ASSERT(sd.getDocOption("cn-qc-cutoff").asString()=="0.27");
	CPPUNIT_ASSERT(sd.getDocOption("snp-qc-cutoff").asString()=="1.1");
	CPPUNIT_ASSERT(sd.getDocOption("waviness-seg-count-cutoff").asString()=="100");
	CPPUNIT_ASSERT(sd.getDocOption("use-high-waviness-seg-count").asString()=="true");
	CPPUNIT_ASSERT(sd.getDocOption("force").asString()=="false");
	CPPUNIT_ASSERT(sd.getDocOption("keep-temp-data").asString()=="false");
	CPPUNIT_ASSERT(sd.getDocOption("selected-qc").asString()=="snp-qc");
	
}

void CNReferenceMethodAdditionalWavesTest::newObjectTest()
{
	
	Verbose::out(1, "****CNReferenceMethodAdditionalWavesTest::newObjectTest****");

     //create default object from factory by string and verify default parameters
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNReferenceMethodAdditionalWaves *cn2=(CNReferenceMethodAdditionalWaves*)obj_CNAnalysisMethodFactory.CNAnalysisMethodForString("additional-waves-reference-method");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj1=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1->size()==11);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5a = obj1->at(4);
	affymetrix_calvin_parameter::ParameterNameValueType param6a = obj1->at(5);	
	affymetrix_calvin_parameter::ParameterNameValueType param7a = obj1->at(6);
	affymetrix_calvin_parameter::ParameterNameValueType param8a = obj1->at(7);
	affymetrix_calvin_parameter::ParameterNameValueType param9a = obj1->at(8);	
	affymetrix_calvin_parameter::ParameterNameValueType param10a = obj1->at(9);
	affymetrix_calvin_parameter::ParameterNameValueType param11a = obj1->at(10);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-trim");
	CPPUNIT_ASSERT(param1a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1a.GetValueFloat()==2.0f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-percentile");
	CPPUNIT_ASSERT(param2a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2a.GetValueFloat()==0.75f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-additional-wave-count");
	CPPUNIT_ASSERT(param3a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-demean");
	CPPUNIT_ASSERT(param4a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param4a.GetValueInt8()==0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5a.GetName())=="affymetrix-algorithm-param-cn-qc-cutoff");
	CPPUNIT_ASSERT(param5a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5a.GetValueFloat()==0.27f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6a.GetName())=="affymetrix-algorithm-param-snp-qc-cutoff");
	CPPUNIT_ASSERT(param6a.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param6a.GetValueFloat()==1.1f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7a.GetName())=="affymetrix-algorithm-param-waviness-seg-count-cutoff");
	CPPUNIT_ASSERT(param7a.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param7a.GetValueInt32()==100);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param8a.GetName())=="affymetrix-algorithm-param-use-high-waviness-seg-count");
	CPPUNIT_ASSERT(param8a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param8a.GetValueInt8()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param9a.GetName())=="affymetrix-algorithm-param-force");
	CPPUNIT_ASSERT(param9a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param9a.GetValueInt8()==0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param10a.GetName())=="affymetrix-algorithm-param-keep-temp-data");
	CPPUNIT_ASSERT(param10a.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param10a.GetValueInt8()==0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param11a.GetName())=="affymetrix-algorithm-param-selected-qc");
    CPPUNIT_ASSERT(param11a.GetParameterType()==ParameterNameValueType::AsciiType);
    CPPUNIT_ASSERT(param11a.GetValueAscii()=="snp-qc");
	delete cn2;

    //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["trim"]="1.0";
	params["percentile"]="1.3";
	params["additional-wave-count"]="2";
	params["demean"]="true";
	params["cn-qc-cutoff"]="1.3";
	params["snp-qc-cutoff"]="2.3";
	params["waviness-seg-count-cutoff"]="20000";
	params["use-high-waviness-seg-count"]="false";
	params["force"]="true";
	params["keep-temp-data"]="true";
    params["selected-qc"]="contrast-qc";
	

	CNAnalysisMethod::getParams()->clear();
	SelfCreate *sc=CNReferenceMethodAdditionalWaves::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> *obj=CNAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj->size()==11);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj->at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj->at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj->at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj->at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5 = obj->at(4);
	affymetrix_calvin_parameter::ParameterNameValueType param6 = obj->at(5);	
	affymetrix_calvin_parameter::ParameterNameValueType param7 = obj->at(6);
	affymetrix_calvin_parameter::ParameterNameValueType param8 = obj->at(7);
	affymetrix_calvin_parameter::ParameterNameValueType param9 = obj->at(8);	
	affymetrix_calvin_parameter::ParameterNameValueType param10 = obj->at(9);
	affymetrix_calvin_parameter::ParameterNameValueType param11 = obj->at(10);
	
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-trim");
	CPPUNIT_ASSERT(param1.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param1.GetValueFloat()==1.0f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-percentile");
	CPPUNIT_ASSERT(param2.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param2.GetValueFloat()==1.3f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-additional-wave-count");
	CPPUNIT_ASSERT(param3.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param3.GetValueInt32()==2);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-demean");
	CPPUNIT_ASSERT(param4.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param4.GetValueInt8()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5.GetName())=="affymetrix-algorithm-param-cn-qc-cutoff");
	CPPUNIT_ASSERT(param5.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param5.GetValueFloat()==1.3f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param6.GetName())=="affymetrix-algorithm-param-snp-qc-cutoff");
	CPPUNIT_ASSERT(param6.GetParameterType()==ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT(param6.GetValueFloat()==2.3f);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param7.GetName())=="affymetrix-algorithm-param-waviness-seg-count-cutoff");
	CPPUNIT_ASSERT(param7.GetParameterType()==ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT(param7.GetValueInt32()==20000);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param8.GetName())=="affymetrix-algorithm-param-use-high-waviness-seg-count");
	CPPUNIT_ASSERT(param8.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param8.GetValueInt8()==0);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param9.GetName())=="affymetrix-algorithm-param-force");
	CPPUNIT_ASSERT(param9.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param9.GetValueInt8()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param10.GetName())=="affymetrix-algorithm-param-keep-temp-data");
	CPPUNIT_ASSERT(param10.GetParameterType()==ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT(param10.GetValueInt8()==1);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param11.GetName())=="affymetrix-algorithm-param-selected-qc");
    CPPUNIT_ASSERT(param11.GetParameterType()==ParameterNameValueType::AsciiType);
    CPPUNIT_ASSERT(param11.GetValueAscii()=="contrast-qc");
		
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNReferenceMethodAdditionalWaves::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="additional-waves-reference-method.trim=2.0.percentile=0.75.additional-wave-count=1.demean=false.cn-qc-cutoff=0.27.snp-qc-cutoff=1.1.waviness-seg-count-cutoff=100.use-high-waviness-seg-count=true.force=false.keep-temp-data=false.selected-qc=snp-qc");
    delete sc;	

	
}


