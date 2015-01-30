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

#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNLog2RatioData.h"
#include "copynumber/CNLog2RatioEngine.h"
#include "copynumber/CNReferenceEngine.h"
#include "copynumber/CNWorkflowEngine.h"
#include "copynumber/CPPTest/Setup.h" 
//
#include "util/AffxArray.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <iostream>
//

/**
 * @class CNLog2RatioDataTest
 * @brief cppunit class for testing CNLog2RatioData functions.
 * last change by vliber on 03/23/09
 */

class CNLog2RatioDataTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(CNLog2RatioDataTest);
	CPPUNIT_TEST(loadGenotypeReportTest);
	CPPUNIT_TEST(loadExperimentsTest);
	CPPUNIT_TEST(loadCyto2ModelFileTest);
	CPPUNIT_TEST_SUITE_END();

public:  
   
   void loadGenotypeReportTest();
   void loadExperimentsTest();
   void loadCyto2ModelFileTest();
   
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNLog2RatioDataTest );


void CNLog2RatioDataTest::loadGenotypeReportTest()
{
	cout<<endl;
	Verbose::out(1, "****CNLog2RatioDataTest::loadGenotypeReportTest****");
	CNLog2RatioData cn1, cn2, cn3;
	CNWorkflowEngine engine;
	cn1.setEngine(&engine);
	//FATAL ERROR: Genotype report file not specfied.
	NEGATIVE_TEST(cn1.loadGenotypeReport(""),Except);
	//FATAL ERROR: Genotype report file must contain the computed_gender column.
	NEGATIVE_TEST(cn1.loadGenotypeReport(INPUT + "/Test.reportBad1.txt"),Except);
	//FATAL ERROR: Cannot open file: ./input/Test.reportNotExist.txt
	NEGATIVE_TEST(cn1.loadGenotypeReport(INPUT + "/Test.reportNotExist.txt"),Except);
    
	engine.setOpt("set-analysis-name", "Test");
	cn2.setEngine(&engine);	
	cn3.setEngine(&engine);	
	//WARNING: Cannot determine sex for experiment:
	// TODO: Fix this test.
//	POSITIVE_TEST(cn2.loadGenotypeReport(INPUT + "/Test.report2.txt"));
}

void CNLog2RatioDataTest::loadExperimentsTest()
{
	Verbose::out(1, "****CNLog2RatioDataTest::loadExperimentsTest****");
	CNLog2RatioData cn1;
	CNWorkflowEngine m_objCNWorkflowEngine;
	m_objCNWorkflowEngine.setOpt("set-analysis-name", "Test");

	///@todo AW: m_objCNWorkflowEngine.setState("create-reference", "true");
	cn1.setEngine(&m_objCNWorkflowEngine);
	
	CPPUNIT_ASSERT(cn1.loadExperiments(INPUT + "/Test.plier.summary.txt")==5);
	CPPUNIT_ASSERT(cn1.getExperiments()->getAt(0)->getExperimentName()=="NA06985_GW6_C.CEL"); 
	CPPUNIT_ASSERT(cn1.getExperiments()->getAt(1)->getExperimentName()=="NA06991_GW6_C.CEL");
	CPPUNIT_ASSERT(cn1.getExperiments()->getAt(2)->getExperimentName()=="NA06993_GW6_C.CEL");
	CPPUNIT_ASSERT(cn1.getExperiments()->getAt(3)->getExperimentName()=="NA06994_GW6_C.CEL");
	CPPUNIT_ASSERT(cn1.getExperiments()->getAt(4)->getExperimentName()=="NA07000_GW6_C.CEL");

	///@todo AW: m_objCNWorkflowEngine.setState("create-reference", "false");
	CPPUNIT_ASSERT(cn1.loadExperiments(INPUT + "/Test.plier.summary.txt")==5); 
	CPPUNIT_ASSERT(cn1.getExperiments()->getAt(0)->getExperimentName()=="NA06985_GW6_C.CEL");
	CPPUNIT_ASSERT(cn1.getExperiments()->getAt(4)->getExperimentName()=="NA07000_GW6_C.CEL");
}
void CNLog2RatioDataTest::loadCyto2ModelFileTest()
{
	
	Verbose::out(1, "****CNLog2RatioDataTest::loadCyto2ModelFileTest****");
	CNLog2RatioData cn;
	CNReferenceEngine en;
	POSITIVE_TEST(cn.loadCyto2ModelFile(INPUT+"/Cyto/ref.a5"));
	CNProbeSetArray * ar=cn.getProbeSets();
	CPPUNIT_ASSERT(ar->getCount()==0);

	en.setOpt("reference-file",INPUT+"/Cyto/ref.a5");
	en.setOpt("annotation-file",INPUT+"/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cn.setEngine(&en);
    ///@todo AW: Fixme. Busted as part of options refactor
    /*
	POSITIVE_TEST(cn.loadAnnotation(false));
	ar=cn.getProbeSets();
	CPPUNIT_ASSERT(ar->getCount()==315666);
	CPPUNIT_ASSERT(ar->at(0)->getProbeSetName()=="C-00IGF");
	CPPUNIT_ASSERT(ar->at(1)->getProbeSetName()=="C-00IGG");
	CPPUNIT_ASSERT(ar->at(ar->getCount()-1)->getProbeSetName()=="S-2XQMG");	
    */
}


