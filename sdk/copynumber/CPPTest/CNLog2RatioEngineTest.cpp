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

#include "copynumber/CNLog2RatioEngine.h"
#include "copynumber/CPPTest/Setup.h"
//

#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <iostream>
//
/**
 * @class CNLog2RatioEngineTest
 * @brief cppunit class for testing CNLog2RatioEngine functions.
 * last change by vliber on 01/26/09
 */

class CNLog2RatioEngineTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(CNLog2RatioEngineTest);
	CPPUNIT_TEST(defineOptionsTest);
    CPPUNIT_TEST(checkOptionsTest);
	CPPUNIT_TEST_SUITE_END();

public:  
   void defineOptionsTest();
   void checkOptionsTest();
   void runTest();
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNLog2RatioEngineTest );


void CNLog2RatioEngineTest::defineOptionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNLog2RatioEngineTest::defineOptionsTest****");
	CNLog2RatioEngine m_objCNLog2RatioEngine;
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptBool("help")==false);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptInt("verbose")==1);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptBool("version")==false);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("out-dir")==".");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("command-line")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("exec-guid")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("program-name")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("program-company")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("program-version")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("program-cvs-id")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("version-to-report")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("probeset-ids")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("annotation-file")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("expr-summary-file")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("genotype-calls-file")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("genotype-confidences-file")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("genotype-report-file")=="");
    CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptInt("xChromosome")==24);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptInt("yChromosome")==25);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("reference-file")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptBool("call-copynumber-engine")==true);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptBool("log2ratio-hdf5-output")==false);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptBool("log2ratio-text-output")==false);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptInt("mem-usage")==0);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("analysis")=="");
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptInt("gc-correction-bin-count")==25);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptBool("delete-files")==false);
    CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOptBool("log2-input")==false);
    CPPUNIT_ASSERT(m_objCNLog2RatioEngine.getOpt("gc-content-override-file")=="");
}
void CNLog2RatioEngineTest::checkOptionsTest()
{
	Verbose::out(1, "****CNLog2RatioEngineTest::checkOptionsTest****");
	
	CNLog2RatioEngine m_objCNLog2RatioEngine1;
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine1.isOptDefined("set-analysis-name"));
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine1.getOpt("set-analysis-name")=="");
	//Must specify a reference-file.
	NEGATIVE_TEST(m_objCNLog2RatioEngine1.checkOptions(),Except);
    
	//Must specify a expr-summary-file
    CNLog2RatioEngine m_objCNLog2RatioEngine2a;
    m_objCNLog2RatioEngine2a.setOpt("reference-file",INPUT + "/CNReference.a5");
	NEGATIVE_TEST(m_objCNLog2RatioEngine2a.checkOptions(),Except);

	//Must specify a genotype-calls-file
    CNLog2RatioEngine m_objCNLog2RatioEngine2;
	m_objCNLog2RatioEngine2.setOpt("reference-file",INPUT + "/CNReference.a5");
    m_objCNLog2RatioEngine2.setOpt("expr-summary-file",INPUT + "/Test.plier.summary.txt");
	NEGATIVE_TEST(m_objCNLog2RatioEngine2.checkOptions(), Except);

    //Must specify a genotype-confidences-file
    CNLog2RatioEngine m_objCNLog2RatioEngine3;
	m_objCNLog2RatioEngine3.setOpt("reference-file",INPUT + "/CNReference.a5");
	m_objCNLog2RatioEngine3.setOpt("expr-summary-file",INPUT + "/Test.plier.summary.txt");
    m_objCNLog2RatioEngine3.setOpt("genotype-calls-file",INPUT + "/Test.report1.txt");
	NEGATIVE_TEST(m_objCNLog2RatioEngine3.checkOptions(), Except);

	//Must specify a genoptype-report-file
    CNLog2RatioEngine m_objCNLog2RatioEngine4;
	m_objCNLog2RatioEngine4.setOpt("reference-file",INPUT + "/CNReference.a5");
	m_objCNLog2RatioEngine4.setOpt("expr-summary-file",INPUT + "/Test.plier.summary.txt");
    m_objCNLog2RatioEngine4.setOpt("genotype-calls-file",INPUT + "/Test.report1.txt");
	m_objCNLog2RatioEngine4.setOpt("genotype-confidences-file",INPUT + "/Test.report2.txt");
	NEGATIVE_TEST(m_objCNLog2RatioEngine4.checkOptions(), Except);

	//Must specify a netaffx annotation-file.
    CNLog2RatioEngine m_objCNLog2RatioEngine5;
	m_objCNLog2RatioEngine5.setOpt("reference-file",INPUT + "/CNReference.a5");
	m_objCNLog2RatioEngine5.setOpt("expr-summary-file",INPUT + "/Test.plier.summary.txt");
    m_objCNLog2RatioEngine5.setOpt("genotype-calls-file",INPUT + "/Test.report1.txt");
	m_objCNLog2RatioEngine5.setOpt("genotype-confidences-file",INPUT + "/Test.report2.txt");
	m_objCNLog2RatioEngine5.setOpt("genotype-report-file",INPUT + "/Test.report.txt");
	NEGATIVE_TEST(m_objCNLog2RatioEngine5.checkOptions(), Except);

	
	CNLog2RatioEngine m_objCNLog2RatioEngine6;
	m_objCNLog2RatioEngine6.setOpt("reference-file",INPUT + "/CNReference.a5");
	m_objCNLog2RatioEngine6.setOpt("expr-summary-file",INPUT + "/Test.plier.summary.txt");
    m_objCNLog2RatioEngine6.setOpt("genotype-calls-file",INPUT + "/Test.report1.txt");
	m_objCNLog2RatioEngine6.setOpt("genotype-confidences-file",INPUT + "/Test.report2.txt");
	m_objCNLog2RatioEngine6.setOpt("genotype-report-file",INPUT + "Test.report.txt");
	m_objCNLog2RatioEngine6.setOpt("reference-file",INPUT + "/CNReference.a5");
	m_objCNLog2RatioEngine6.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNLog2RatioEngine6.defineOption("", "cn-calibrate-parameters", PgOpt::STRING_OPT,
                    "SmoothSignal calibration parameters", "");
	POSITIVE_TEST(m_objCNLog2RatioEngine6.setOpt("log2ratio-text-output","true"));
	POSITIVE_TEST(m_objCNLog2RatioEngine6.checkOptions());
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine6.getOptBool("log2ratio-hdf5-output")==true);
	CPPUNIT_ASSERT(m_objCNLog2RatioEngine6.getOptBool("log2ratio-text-output")==false);	
}



