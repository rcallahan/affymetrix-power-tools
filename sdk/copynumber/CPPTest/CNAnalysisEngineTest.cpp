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

#include "copynumber/CNAnalysisEngine.h"
#include "copynumber/CPPTest/Setup.h"
//
#include "util/PgOptions.h"
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
//


/**
 * @class CNAnalysisEngineTest
 * @brief cppunit class for testing CNAnalysisEngine functions.
 * last change by vliber on 03/20/09
 */

class CNAnalysisEngineTest : public CppUnit::TestFixture  
{
  CPPUNIT_TEST_SUITE(CNAnalysisEngineTest);
  CPPUNIT_TEST(defineOptionsTest);  
  CPPUNIT_TEST(checkOptionsTest);
  CPPUNIT_TEST(createAnalysisTest);
  CPPUNIT_TEST_SUITE_END();

public:   
  void defineOptionsTest();
  void checkOptionsTest();
  void createAnalysisTest();
 
  
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisEngineTest );

void CNAnalysisEngineTest::defineOptionsTest()
{
	cout<<endl;
	Util::PrintTextFunctionTitle("CNAnalysisEngineTest","defineOptionsTest");
	CNAnalysisEngine m_objCNAnalysisEngine;
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOptBool("help")==false); 
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOptInt("verbose")==1);
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOptBool("version")==false);
	
	
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("command-line")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("exec-guid")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("program-name")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("program-company")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("program-version")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("program-cvs-id")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("version-to-report")=="");


    CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("log2ratio-file")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("out-dir")==".");
    CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("analysis")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOptInt("gc-correction-bin-count")==25);
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOptInt("xChromosome")==24);
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOptInt("yChromosome")==25);
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOptInt("mem-usage")==0);

	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("geno-qc-file")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOptBool("cyto2")==false);
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("array-name")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOpt("set-analysis-name")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine.getOptBool("text-output")==false);
}
void CNAnalysisEngineTest::checkOptionsTest()
{
	Util::PrintTextFunctionTitle("CNAnalysisEngineTest","checkOptionsTest");

	CNAnalysisEngine m_objCNAnalysisEngine1; 
    CPPUNIT_ASSERT(m_objCNAnalysisEngine1.isOptDefined("set-analysis-name"));
	CPPUNIT_ASSERT(m_objCNAnalysisEngine1.getOpt("set-analysis-name")=="");
	CPPUNIT_ASSERT(m_objCNAnalysisEngine1.isOptDefined("geno-qc-file"));
	//Invalid (missing) log2ratio-file specification
	NEGATIVE_TEST(m_objCNAnalysisEngine1.checkOptions(),Except);

	//Invalid (blank) log2ratio-file specification
	CNAnalysisEngine m_objCNAnalysisEngine2;
	CPPUNIT_ASSERT(m_objCNAnalysisEngine2.isOptDefined("log2ratio-file"));
	std::vector<std::string> vLog2RatioFileNames;
    vLog2RatioFileNames.push_back("");//add nothing
    m_objCNAnalysisEngine2.setOpt("log2ratio-file",vLog2RatioFileNames);
	NEGATIVE_TEST(m_objCNAnalysisEngine2.checkOptions(),Except);
}

void CNAnalysisEngineTest::createAnalysisTest() 
{
	Util::PrintTextFunctionTitle("CNAnalysisEngineTest","createAnalysisTest");

    // cannot call createAnalysis until checkOptions is called
    /*
	CNAnalysisEngine m_objCNAnalysisEngine1;
    POSITIVE_TEST(m_objCNAnalysisEngine1.createAnalysis());
	m_objCNAnalysisEngine1.setOpt("analysis","");
	//Invalid (blank) analysis specification
	NEGATIVE_TEST(m_objCNAnalysisEngine1.createAnalysis(),Except);
    */

	CNAnalysisEngine m_objCNAnalysisEngine2;
	vector<std::string> analysis;
	//FATAL ERROR: Not values specified for option analysis
    ///@todo AW: Broken on merge with options refactor
	//NEGATIVE_TEST(m_objCNAnalysisEngine2.setOpt("analysis",analysis),Except);	

	CNAnalysisEngine m_objCNAnalysisEngine3;
	analysis.push_back("cyto2");//wrong string
	m_objCNAnalysisEngine3.setOpt("analysis",analysis);
	//FATAL ERROR: Don't know how to make object for name: cyto2 of type: CNAnalysisMethod
	NEGATIVE_TEST(m_objCNAnalysisEngine3.createAnalysis(),Except);

	CNAnalysisEngine m_objCNAnalysisEngine4;
	analysis.clear();
	analysis.push_back("cn-cyto2");
	m_objCNAnalysisEngine4.setOpt("analysis",analysis);
	//FATAL ERROR: State yTarget cannot be found in the states for this engine.
	NEGATIVE_TEST(m_objCNAnalysisEngine4.createAnalysis(),Except);

	CNAnalysisEngine m_objCNAnalysisEngine5;
	//Cannot calculate CN calibration parameters.
	NEGATIVE_TEST(m_objCNAnalysisEngine5.calculateCalibrationParameters(),Except);
}

