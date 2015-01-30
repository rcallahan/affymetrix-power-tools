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

#include "copynumber/CNReferenceEngine.h"
#include "copynumber/CPPTest/Setup.h"
//

//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <iostream>
//

/**
 * @class CNReferenceEngineTest
 * @brief cppunit class for testing CNReferenceEngine functions.
 * last change by vliber on 02/10/09
 */

class CNReferenceEngineTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(CNReferenceEngineTest);
	CPPUNIT_TEST(defineOptionsTest);
	CPPUNIT_TEST(checkOptionsTest);
	CPPUNIT_TEST(runTest);
    CPPUNIT_TEST_SUITE_END();

public: 
   void defineOptionsTest();
   void checkOptionsTest();
   void runTest();
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNReferenceEngineTest );

void CNReferenceEngineTest::defineOptionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNReferenceEngineTest::defineOptionsTest****");
	CNReferenceEngine m_objCNReferenceEngine;
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOptBool("help")==false); 
	//CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("explain")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOptInt("verbose")==1);
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOptBool("version")==false);
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("command-line")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("exec-guid")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("program-name")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("program-company")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("program-version")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("program-cvs-id")=="");
	
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("expr-summary-file")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("genotype-calls-file")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("genotype-confidences-file")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("genotype-report-file")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("reference-file")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("out-dir")==".");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOptBool("reference-text-output")==false);
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOptBool("log2-input")==false);
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOptBool("adapter-type-normalization")==true);
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("set-analysis-name")=="");
    CPPUNIT_ASSERT(m_objCNReferenceEngine.getOptInt("mem-usage")==0);

	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("probeset-ids")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOpt("annotation-file")=="");
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOptInt("xChromosome")==24);
	CPPUNIT_ASSERT(m_objCNReferenceEngine.getOptInt("yChromosome")==25);
}
void CNReferenceEngineTest::checkOptionsTest()
{
	Verbose::out(1, "****CNReferenceEngineTest::checkOptionsTest****");

	CNReferenceEngine m_objCNReferenceEngine1;
    CPPUNIT_ASSERT(m_objCNReferenceEngine1.isOptDefined("set-analysis-name"));
	CPPUNIT_ASSERT(m_objCNReferenceEngine1.getOpt("set-analysis-name")=="");
	//Must specify a expr-summary-file
	NEGATIVE_TEST(m_objCNReferenceEngine1.checkOptions(),Except);
    
	//Must specify a genotype-calls-file
    CNReferenceEngine m_objCNReferenceEngine2;
    m_objCNReferenceEngine2.setOpt("expr-summary-file",INPUT + "Test.plier.summary.txt");
	NEGATIVE_TEST(m_objCNReferenceEngine2.checkOptions(),Except);

    //Must specify a genotype-confidences-file
    CNReferenceEngine m_objCNReferenceEngine3;
	m_objCNReferenceEngine3.setOpt("expr-summary-file",INPUT + "Test.plier.summary.txt");
    m_objCNReferenceEngine3.setOpt("genotype-calls-file",INPUT + "Test.plier.summary.txt");
	NEGATIVE_TEST(m_objCNReferenceEngine3.checkOptions(), Except);

	//Must specify a genoptype-report-file
    CNReferenceEngine m_objCNReferenceEngine4;
	m_objCNReferenceEngine4.setOpt("expr-summary-file",INPUT + "Test.plier.summary.txt");
    m_objCNReferenceEngine4.setOpt("genotype-calls-file",INPUT + "Test.plier.summary.txt");
	m_objCNReferenceEngine4.setOpt("genotype-confidences-file",INPUT + "Test.plier.summary.txt");
	NEGATIVE_TEST(m_objCNReferenceEngine4.checkOptions(), Except);

	//Must specify a reference-file
    CNReferenceEngine m_objCNReferenceEngine5;
	m_objCNReferenceEngine5.setOpt("expr-summary-file",INPUT + "Test.plier.summary.txt");
    m_objCNReferenceEngine5.setOpt("genotype-calls-file",INPUT + "Test.plier.summary.txt");
	m_objCNReferenceEngine5.setOpt("genotype-confidences-file",INPUT + "Test.plier.summary.txt");
	m_objCNReferenceEngine5.setOpt("genotype-report-file",INPUT + "Test.plier.summary.txt");
	NEGATIVE_TEST(m_objCNReferenceEngine5.checkOptions(),Except);

	//Must specify a netaffx-annotation-file
    CNReferenceEngine m_objCNReferenceEngine6;
	m_objCNReferenceEngine6.setOpt("expr-summary-file",INPUT + "Test.plier.summary.txt");
    m_objCNReferenceEngine6.setOpt("genotype-calls-file",INPUT + "Test.plier.summary.txt");
	m_objCNReferenceEngine6.setOpt("genotype-confidences-file",INPUT + "Test.plier.summary.txt");
	m_objCNReferenceEngine6.setOpt("genotype-report-file",INPUT + "Test.report1.txt");
	m_objCNReferenceEngine6.setOpt("reference-file",INPUT + "/CNReference.a5");
	NEGATIVE_TEST(m_objCNReferenceEngine6.checkOptions(),Except);
    
    
	CNReferenceEngine m_objCNReferenceEngine7;
	m_objCNReferenceEngine7.setOpt("expr-summary-file",INPUT + "Test.plier.summary.txt");
    m_objCNReferenceEngine7.setOpt("genotype-calls-file",INPUT + "Test.plier.summary.txt");
	m_objCNReferenceEngine7.setOpt("genotype-confidences-file",INPUT + "Test.plier.summary.txt");
	m_objCNReferenceEngine7.setOpt("genotype-report-file",INPUT + "Test.report1.txt");
	m_objCNReferenceEngine7.setOpt("reference-file",INPUT + "/CNReference.a5");
	m_objCNReferenceEngine7.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	POSITIVE_TEST(m_objCNReferenceEngine7.checkOptions());
	CPPUNIT_ASSERT(m_objCNReferenceEngine7.getOptBool("create-reference")==true);
	CPPUNIT_ASSERT(m_objCNReferenceEngine7.getOpt("reference-file")==INPUT + "/CNReference.a5");
}

void CNReferenceEngineTest::runTest()
{
	Verbose::out(1, "****CNReferenceEngineTest::runTest****");
	
	CNReferenceEngine m_objCNReferenceEngine1;
	m_objCNReferenceEngine1.setOpt("expr-summary-file",INPUT + "Test.plier.summary.txt");
    m_objCNReferenceEngine1.setOpt("genotype-calls-file",INPUT + "Test.plier.summary.txt");
	m_objCNReferenceEngine1.setOpt("genotype-confidences-file",INPUT + "Test.plier.summary.txt");
	m_objCNReferenceEngine1.setOpt("genotype-report-file",INPUT + "Test.report1.txt");
	m_objCNReferenceEngine1.setOpt("reference-file",INPUT + "/CNReference.a5");
	m_objCNReferenceEngine1.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	//ERROR: SQLiteCode: 26, Message: Failed to prepare SQL statement.
    NEGATIVE_TEST(m_objCNReferenceEngine1.run(),Except);   
}
