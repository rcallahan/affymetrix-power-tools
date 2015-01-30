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
#include "copynumber/CNFamilialAnalysisMethodHeteroUPD.h"
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
 * @class CNFamilialAnalysisMethodHeteroUPDTest
 * @brief cppunit class for testing CNFamilialAnalysisMethodHeteroUPD functions.
 * last change by vliber on 10/08/09
  */

class CNFamilialAnalysisMethodHeteroUPDTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNFamilialAnalysisMethodHeteroUPDTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
 
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNFamilialAnalysisMethodHeteroUPDTest );

void CNFamilialAnalysisMethodHeteroUPDTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNFamilialAnalysisMethodHeteroUPDTest::functionsTest****");
	// constructor && m_vParams()
	CNFamilialAnalysisMethodHeteroUPD cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="hetero-upd");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Familial HeteroUPD");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="hetero-upd");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNFamilialAnalysisMethodHeteroUPD::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==5);
	CPPUNIT_ASSERT(sv[0].name=="hetero-upd-marker-count-cutoff");
	CPPUNIT_ASSERT(sv[0].type==3);
	CPPUNIT_ASSERT(sv[0].value=="21");
	CPPUNIT_ASSERT(sv[0].defaultVal=="21");
	CPPUNIT_ASSERT(sv[0].minVal=="1");
	CPPUNIT_ASSERT(sv[0].maxVal=="NA");
	CPPUNIT_ASSERT(sv[0].descript=="HeteroUPD Marker Count Cutoff");
	CPPUNIT_ASSERT(sv[1].name=="hetero-upd-call-cutoff");
	CPPUNIT_ASSERT(sv[1].type==3);
	CPPUNIT_ASSERT(sv[1].value=="4");
	CPPUNIT_ASSERT(sv[1].defaultVal=="4");
	CPPUNIT_ASSERT(sv[1].minVal=="0");
	CPPUNIT_ASSERT(sv[1].maxVal=="NA");
	CPPUNIT_ASSERT(sv[1].descript=="HeteroUPD Call Cutoff");
	CPPUNIT_ASSERT(sv[2].name=="hetero-upd-min-genomic-span");
	CPPUNIT_ASSERT(sv[2].type==3);
	CPPUNIT_ASSERT(sv[2].value=="1000000");
	CPPUNIT_ASSERT(sv[2].defaultVal=="1000000");
	CPPUNIT_ASSERT(sv[2].minVal=="1");
	CPPUNIT_ASSERT(sv[2].maxVal=="NA");
	CPPUNIT_ASSERT(sv[2].descript=="HeteroUPD Minimum Genomic Span");
	CPPUNIT_ASSERT(sv[3].name=="hetero-upd-xChromosome");
	CPPUNIT_ASSERT(sv[3].type==3);
	CPPUNIT_ASSERT(sv[3].value=="24");
	CPPUNIT_ASSERT(sv[3].defaultVal=="24");
	CPPUNIT_ASSERT(sv[3].minVal=="1");
	CPPUNIT_ASSERT(sv[3].maxVal=="NA");
	CPPUNIT_ASSERT(sv[3].descript=="HeteroUPD X Chromosome");
	CPPUNIT_ASSERT(sv[4].name=="hetero-upd-yChromosome");
	CPPUNIT_ASSERT(sv[4].type==3);
	CPPUNIT_ASSERT(sv[4].value=="25");
	CPPUNIT_ASSERT(sv[4].defaultVal=="25");
	CPPUNIT_ASSERT(sv[4].minVal=="1");
	CPPUNIT_ASSERT(sv[4].maxVal=="NA");
	CPPUNIT_ASSERT(sv[4].descript=="HeteroUPD Y Chromosome");
	
	//explainSelf()
	SelfDoc sd=CNFamilialAnalysisMethodHeteroUPD::explainSelf();
	//std::cout<<sd.getState()<<std::endl;
	CPPUNIT_ASSERT(sd.getState()=="hetero-upd.hetero-upd-marker-count-cutoff=21.hetero-upd-call-cutoff=4.hetero-upd-min-genomic-span=1000000.hetero-upd-xChromosome=24.hetero-upd-yChromosome=25");
	CPPUNIT_ASSERT(sd.getDocName()=="hetero-upd");
	CPPUNIT_ASSERT(sd.getDocDescription()=="Familial HeteroUPD");
	vector<SelfDoc::Opt> v=sd.getDocOptions();
	CPPUNIT_ASSERT(v.size()==5);
	CPPUNIT_ASSERT(v[0].asString()=="21");
	CPPUNIT_ASSERT(v[1].asString()=="4");
	CPPUNIT_ASSERT(v[2].asString()=="1000000");
	CPPUNIT_ASSERT(v[3].asString()=="24");
	CPPUNIT_ASSERT(v[4].asString()=="25");
	CPPUNIT_ASSERT(sd.getDocOption("hetero-upd-marker-count-cutoff").asString()=="21");
	CPPUNIT_ASSERT(sd.getDocOption("hetero-upd-call-cutoff").asString()=="4");
	CPPUNIT_ASSERT(sd.getDocOption("hetero-upd-min-genomic-span").asString()=="1000000");
	CPPUNIT_ASSERT(sd.getDocOption("hetero-upd-xChromosome").asString()=="24");
	CPPUNIT_ASSERT(sd.getDocOption("hetero-upd-yChromosome").asString()=="25");
	
}
void CNFamilialAnalysisMethodHeteroUPDTest::newObjectTest()
{
	
	Verbose::out(1, "****CNFamilialAnalysisMethodHeteroUPDTest::newObjectTest****");

    //create default object from factory by string and verify default parameters
	CNFamilialAnalysisMethod::getParams().clear();
	CNFamilialAnalysisMethodFactory obj_CNAnalysisMethodFactory;
	CNFamilialAnalysisMethodHeteroUPD *cn2=(CNFamilialAnalysisMethodHeteroUPD*)obj_CNAnalysisMethodFactory.CNFamilialAnalysisMethodForString("hetero-upd");
    vector<affymetrix_calvin_parameter::ParameterNameValueType> obj1=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj1.size()==5);
	affymetrix_calvin_parameter::ParameterNameValueType param1a = obj1.at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2a = obj1.at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3a = obj1.at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4a = obj1.at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5a = obj1.at(4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1a.GetName())=="affymetrix-algorithm-param-hetero-upd-marker-count-cutoff");
	CPPUNIT_ASSERT(param1a.GetParameterType()==4);
    CPPUNIT_ASSERT(param1a.GetValueInt32()==21);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2a.GetName())=="affymetrix-algorithm-param-hetero-upd-call-cutoff");
	CPPUNIT_ASSERT(param2a.GetParameterType()==4);
	CPPUNIT_ASSERT(param2a.GetValueInt32()==4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3a.GetName())=="affymetrix-algorithm-param-hetero-upd-min-genomic-span");
	CPPUNIT_ASSERT(param3a.GetParameterType()==4);
	CPPUNIT_ASSERT(param3a.GetValueInt32()==1000000);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4a.GetName())=="affymetrix-algorithm-param-hetero-upd-xChromosome");
	CPPUNIT_ASSERT(param4a.GetParameterType()==4);
	CPPUNIT_ASSERT(param4a.GetValueInt32()==24);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5a.GetName())=="affymetrix-algorithm-param-hetero-upd-yChromosome");
	CPPUNIT_ASSERT(param5a.GetParameterType()==4);
	CPPUNIT_ASSERT(param5a.GetValueInt32()==25);
	delete cn2;


   //create newObject with different parameters and test them
	std::map<std::string,std::string> params;
	params["hetero-upd-marker-count-cutoff"]="7";
	params["hetero-upd-call-cutoff"]="34";
	params["hetero-upd-min-genomic-span"]="666";
	params["hetero-upd-xChromosome"]="77";
	params["hetero-upd-yChromosome"]="99";
	CNFamilialAnalysisMethod::getParams().clear();
	SelfCreate *sc=CNFamilialAnalysisMethodHeteroUPD::newObject(params);
    	
	vector<affymetrix_calvin_parameter::ParameterNameValueType> obj=CNFamilialAnalysisMethod::getParams();
	CPPUNIT_ASSERT(obj.size()==5);
	affymetrix_calvin_parameter::ParameterNameValueType param1 = obj.at(0);
	affymetrix_calvin_parameter::ParameterNameValueType param2 = obj.at(1);
	affymetrix_calvin_parameter::ParameterNameValueType param3 = obj.at(2);
	affymetrix_calvin_parameter::ParameterNameValueType param4 = obj.at(3);
	affymetrix_calvin_parameter::ParameterNameValueType param5 = obj.at(4);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param1.GetName())=="affymetrix-algorithm-param-hetero-upd-marker-count-cutoff");
	CPPUNIT_ASSERT(param1.GetParameterType()==4);
    CPPUNIT_ASSERT(param1.GetValueInt32()==7);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param2.GetName())=="affymetrix-algorithm-param-hetero-upd-call-cutoff");
	CPPUNIT_ASSERT(param2.GetParameterType()==4);
	CPPUNIT_ASSERT(param2.GetValueInt32()==34);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param3.GetName())=="affymetrix-algorithm-param-hetero-upd-min-genomic-span");
	CPPUNIT_ASSERT(param3.GetParameterType()==4);
	CPPUNIT_ASSERT(param3.GetValueInt32()==666);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param4.GetName())=="affymetrix-algorithm-param-hetero-upd-xChromosome");
	CPPUNIT_ASSERT(param4.GetParameterType()==4);
	CPPUNIT_ASSERT(param4.GetValueInt32()==77);
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(param5.GetName())=="affymetrix-algorithm-param-hetero-upd-yChromosome");
	CPPUNIT_ASSERT(param5.GetParameterType()==4);
	CPPUNIT_ASSERT(param5.GetValueInt32()==99);
	
	//always return default set of options and value even different object has been created by newObject
	SelfDoc sd1=CNFamilialAnalysisMethodHeteroUPD::explainSelf();
	CPPUNIT_ASSERT(sd1.getState()=="hetero-upd.hetero-upd-marker-count-cutoff=21.hetero-upd-call-cutoff=4.hetero-upd-min-genomic-span=1000000.hetero-upd-xChromosome=24.hetero-upd-yChromosome=25");
    delete sc;	
}

