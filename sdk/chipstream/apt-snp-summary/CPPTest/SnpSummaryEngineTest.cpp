////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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

//
#include "chipstream/apt-snp-summary/SnpSummaryEngine.h"
#include "chipstream/apt-snp-summary/SnpSummaryStats.h"
//
#include "util/Util.h"
#include "util/CPPTest/Setup.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <iostream>
#include <string> 
#include <vector> 

using namespace std;

/**
 * @class SnpSummaryEngineTest
 * @brief cppunit class for testing SnpSummaryEngineTest functions.
 */

class SnpSummaryEngineTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(SnpSummaryEngineTest);   
  CPPUNIT_TEST ( CheckOptionsTest );
  CPPUNIT_TEST( CreateReportFromXMLTest );
  CPPUNIT_TEST(MixedProbeArrayTypeTest);
  CPPUNIT_TEST( SummaryStatsTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void CheckOptionsTest();
  void CreateReportFromXMLTest();
  void MixedProbeArrayTypeTest();
  void SummaryStatsTest();
  void runTest();
  
};

CPPUNIT_TEST_SUITE_REGISTRATION(SnpSummaryEngineTest);

// check the options
void SnpSummaryEngineTest::CheckOptionsTest()
{
	cout<<endl;
	Util::PrintTextFunctionTitle("SnpSummaryEngineTest","checkOptionsTest");
	SnpSummaryEngine engine;
	CPPUNIT_ASSERT(engine.getOptBool("h")==false); 
	string a = engine.getOpt("summary-out-file");
	CPPUNIT_ASSERT(engine.getOpt("summary-out-file") != "");
	string b = engine.getOpt("output-format");
	CPPUNIT_ASSERT(engine.getOpt("output-format") != "");
}

//test calculating metrics with options specified in XML file
void SnpSummaryEngineTest::CreateReportFromXMLTest()
{
	cout<<endl;
	Util::PrintTextFunctionTitle("SnpSummaryEngineTest","addCelTest");	
    Err::setThrowStatus(true);
    bool caughtException = false;
	try
	{
		SnpSummaryEngine engine;
		engine.setOpt("xml-file","SnpSummaryTestAddCel.job");
		engine.run();	
	}
	catch(...)
	{
        caughtException = true;
	}
    CPPUNIT_ASSERT (caughtException == false);
}

// test for failure if two different array types are in list of CEL files
void SnpSummaryEngineTest::MixedProbeArrayTypeTest()
{
	cout<<endl;
	Util::PrintTextFunctionTitle("SnpSummaryEngineTest","mixedProbeArrayTypeTest");	
	cout<<"The fatal error below is expected. ('Unable to extract SNP summaries')" << endl;
	Err::setThrowStatus(true);
        SnpSummaryEngine engine;
        engine.setOpt("xml-file","SnpSummaryTestMixProbeArray.job");
        NEGATIVE_TEST(engine.run(), Except);
}
	
static string ConvertFloatToString(float value)
{
	std::stringstream data;
	data.setf(ios::fixed, ios::floatfield);
	data.precision(8);
	data << value; 
	return data.str();
}

// test calculation of metrics
void SnpSummaryEngineTest::SummaryStatsTest()
{
	cout<<endl;
	Util::PrintTextFunctionTitle("SnpSummaryEngineTest","SummaryStatsTest");	
	
	int nSets = 2;
	vector<vector<u_int8_t> > calls(nSets);
	vector<vector<float> > metrics(nSets);
	for (int iset=0; iset<nSets; iset++)
		calls[iset].resize(2);
	calls[0][0] = 6;
	calls[0][1] = 7;
	calls[1][0] = 8;
	calls[1][1] = 11;
	for (vector<vector<float> >::iterator metricIt=metrics.begin(); metricIt!=metrics.end(); metricIt++)
		metricIt->resize(6);
	SnpSummaryStats stats;
	stats.CalculateMetrics(calls, nSets, metrics);
	CPPUNIT_ASSERT (metrics[0][0] == 100.0);
	CPPUNIT_ASSERT (metrics[0][1] == 50.0);
	CPPUNIT_ASSERT (metrics[0][2] == 0.0);
	CPPUNIT_ASSERT (metrics[0][3] == 50.0);
	CPPUNIT_ASSERT (metrics[0][4] == 0.50);
	string value = ConvertFloatToString(metrics[0][5]);
	CPPUNIT_ASSERT(value == "0.15729921");
	CPPUNIT_ASSERT (metrics[1][0] == 50.0);
	CPPUNIT_ASSERT (metrics[1][1] == 0.0);
	CPPUNIT_ASSERT (metrics[1][2] == 50.0);
	CPPUNIT_ASSERT (metrics[1][3] == 00.0);
	CPPUNIT_ASSERT (metrics[1][4] == 0.50);
	value = ConvertFloatToString(metrics[1][5]);
	CPPUNIT_ASSERT (value == "0.31731051");
}
