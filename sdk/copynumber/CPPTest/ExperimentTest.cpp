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

#include "copynumber/CNExperiment.h"
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
#include <string>
//
 
using namespace std;

/**
 * @class ExperimentTest
 * @brief cppunit class for testing Experiment functions.
 * last change by vliber on 01/17/09
 */

class ExperimentTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(ExperimentTest);
  CPPUNIT_TEST(constructorTest);
  CPPUNIT_TEST(get_setTest);
  CPPUNIT_TEST(get_setGenderTest);
  CPPUNIT_TEST(compareToTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void constructorTest();
  void get_setTest(); 
  void get_setGenderTest(); 
  void compareToTest(); 
};
CPPUNIT_TEST_SUITE_REGISTRATION(ExperimentTest );

void ExperimentTest::constructorTest()
{
	cout<<endl;
	Verbose::out(1, "****ExperimentTest::constructorTest****");
	CNExperiment ex1;
	CPPUNIT_ASSERT(ex1.getExperimentName()=="");
	CPPUNIT_ASSERT(ex1.getCNCallGender()=="unknown");
	CPPUNIT_ASSERT(ex1.hasXX()==false);
	CPPUNIT_ASSERT(ex1.hasY()==false);
	CPPUNIT_ASSERT(ex1.getMedianAutosomeMedian()==0);
	CPPUNIT_ASSERT(ex1.getMadDiffCN()==0);
	CPPUNIT_ASSERT(ex1.getIqr()==0);
	CPPUNIT_ASSERT(ex1.getMeanAbsRle()==0);
	CPPUNIT_ASSERT(CNExperiment::getQCMetricColumnNames()->getCount()==0);
	CPPUNIT_ASSERT(ex1.getQCMetricColumnValues()->getCount()==0);
	CPPUNIT_ASSERT(ex1.getChrXMean()==0);
	CPPUNIT_ASSERT(ex1.getChrYMean()==0);
	CPPUNIT_ASSERT(ex1.getMedianCnState()==0);
	CPPUNIT_ASSERT(ex1.getHomFrequency()==0);
	CPPUNIT_ASSERT(ex1.getHetFrequency()==0);
	CPPUNIT_ASSERT(ex1.getCNCallGenderConfidence()==0);
	CPPUNIT_ASSERT(ex1.getMedianRawIntensity()==0);
	CPPUNIT_ASSERT(ex1.getSNPQC()==0);
	CPPUNIT_ASSERT(ex1.getAntigenomicRatio()==0);
	CPPUNIT_ASSERT(ex1.getGenomeLOH()==0);
	CPPUNIT_ASSERT(ex1.getNumberOfChromosomesToReport()==0);
}
void ExperimentTest::get_setTest()
{
	//set, get
	Verbose::out(1, "****ExperimentTest::get_setTest****");
	CNExperiment ex1;
	ex1.setExperimentName("myTest");
	CPPUNIT_ASSERT(ex1.getExperimentName()=="myTest");
	ex1.setXX(true);
	CPPUNIT_ASSERT(ex1.hasXX()==true);
	ex1.setY(true);
	CPPUNIT_ASSERT(ex1.hasY()==true);
	ex1.setMedianAutosomeMedian(0.3456f);
	CPPUNIT_ASSERT(ex1.getMedianAutosomeMedian()==0.3456f);
	ex1.setMadDiffCN(3.876f);
	CPPUNIT_ASSERT(ex1.getMadDiffCN()==3.876f);
	ex1.setIqr(4.876f);
	CPPUNIT_ASSERT(ex1.getIqr()==4.876f);
	ex1.setMeanAbsRle(3.876f);
	CPPUNIT_ASSERT(ex1.getMeanAbsRle()==3.876f);
}
void ExperimentTest::get_setGenderTest()
{
	//setGender
	CNExperiment ex1;
	Util::PrintTextFunctionTitle("ExperimentTest","get_setGenderTest");
	ex1.setY(false);
	ex1.setXX(false);
	ex1.setRawIntensityRatioGenderFromString("male");
	CPPUNIT_ASSERT(ex1.getRawIntensityRatioGenderAsInt()==1);
	ex1.setRawIntensityRatioGenderFromString("female");
	CPPUNIT_ASSERT(ex1.getRawIntensityRatioGenderAsInt()==0);
	ex1.setRawIntensityRatioGenderFromString("test");
	CPPUNIT_ASSERT(ex1.getRawIntensityRatioGenderAsInt()==2);
}
void ExperimentTest::compareToTest()
{
	//compareTo
	Util::PrintTextFunctionTitle("ExperimentTest","comareToTest");
    CNExperiment ex2,ex3, ex4;
	ex2.setExperimentName("myTest");
	ex3.setExperimentName("MyTest");
	ex4.setExperimentName("MyTest");
	CPPUNIT_ASSERT(ex2.compareTo(ex3,0)==1);
	CPPUNIT_ASSERT(ex3.compareTo(ex2,0)==-1);
	CPPUNIT_ASSERT(ex4.compareTo(ex3,0)==0);
}
