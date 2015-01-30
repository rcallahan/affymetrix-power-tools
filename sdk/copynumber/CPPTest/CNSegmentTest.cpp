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

#include "copynumber/CNSegment.h"
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
#include <string>
//
 
using namespace std;

/**
 * @class CNSegmentTest
 * @brief cppunit class for testing CNSegment functions.
 * last change by vliber on 03/27/09
 */

class CNSegmentTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNSegmentTest);
  CPPUNIT_TEST(constructorTest);
  CPPUNIT_TEST(get_setTest);
  CPPUNIT_TEST(comareToTest);
  CPPUNIT_TEST_SUITE_END();

public:  
  void constructorTest();
  void get_setTest(); 
  void get_setGenderTest(); 
  void comareToTest(); 
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNSegmentTest );

void CNSegmentTest::constructorTest()
{
	cout<<endl;
	Verbose::out(1, "****CNSegmentTest::constructorTest****");
	CNSegment seg1;
	CPPUNIT_ASSERT(seg1.getSegmentType()==0);
	CPPUNIT_ASSERT(seg1.getChromosome()==0);
	CPPUNIT_ASSERT(seg1.getStartPosition()==0);
	CPPUNIT_ASSERT(seg1.getEndPosition()==0);
	CPPUNIT_ASSERT(seg1.getCall()==0);
	CPPUNIT_ASSERT(seg1.getConfidence()==0);
	CPPUNIT_ASSERT(seg1.getMarkerCount()==0);
	CPPUNIT_ASSERT(seg1.getMixture()==MCLASS_NONE);
	CPPUNIT_ASSERT(seg1.getCalibratedCN()==0);
	CPPUNIT_ASSERT(seg1.getFamilialSampleKey()==0);
	CPPUNIT_ASSERT(seg1.getHeterozygosity()==0);
	CPPUNIT_ASSERT(seg1.getHomozygosity()==0);
	CPPUNIT_ASSERT(seg1.getMeanMarkerDistance()==0);
	//Cannot find segment type for analysis name test
	NEGATIVE_TEST(CNSegment::getSegmentType("test"),Except);

}
void CNSegmentTest::get_setTest()
{
	//set, get
	Verbose::out(1, "****CNSegmentTest::get_setTest****");
	CNSegment seg1;
    seg1.setSegmentType(2);
	CPPUNIT_ASSERT(seg1.getSegmentType()==2);
    seg1.setChromosome(3);
	CPPUNIT_ASSERT(seg1.getChromosome()==3);
    seg1.setStartPosition(1000);
	CPPUNIT_ASSERT(seg1.getStartPosition()==1000);
    seg1.setEndPosition(2000);
	CPPUNIT_ASSERT(seg1.getEndPosition()==2000);
    seg1.setCall('a');
	CPPUNIT_ASSERT(seg1.getCall()==97);
    seg1.setConfidence(0.999f);
	CPPUNIT_ASSERT(seg1.getConfidence()==0.999f);
    seg1.setMarkerCount(10);
	CPPUNIT_ASSERT(seg1.getMarkerCount()==10);
//     seg1.setMixture(1.88f);
// 	CPPUNIT_ASSERT(seg1.getMixture()==1.88f);
  seg1.setMixture(MCLASS_n10);
  CPPUNIT_ASSERT(seg1.getMixture()==MCLASS_n10);
    seg1.setCalibratedCN(0.55f);
	CPPUNIT_ASSERT(seg1.getCalibratedCN()==0.55f);
    seg1.setFamilialSampleKey(22);
	CPPUNIT_ASSERT(seg1.getFamilialSampleKey()==22);
    seg1.setHeterozygosity(0.111f);
	CPPUNIT_ASSERT(seg1.getHeterozygosity()==0.111f);
    seg1.setHomozygosity(0.22f);
	CPPUNIT_ASSERT(seg1.getHomozygosity()==0.22f);
    seg1.setMeanMarkerDistance(0.33f);
	CPPUNIT_ASSERT(seg1.getMeanMarkerDistance()==0.33f);
}

void CNSegmentTest::comareToTest()
{
	//compareTo by name (0)
	Verbose::out(1, "****CNSegmentTest::comareToTest****");
    CNSegment seg1, seg2, seg3, seg4;
	seg1.setSegmentName("test");
	seg2.setSegmentName("test");
	CPPUNIT_ASSERT(seg1.compareTo(seg2,0)==0);
	seg2.setSegmentName("Test");
	CPPUNIT_ASSERT(seg1.compareTo(seg2,0)==1);
	seg2.setSegmentName("Test");
	CPPUNIT_ASSERT(seg1.compareTo(seg2,0)==1);
	seg2.setSegmentName("vest");
	CPPUNIT_ASSERT(seg1.compareTo(seg2,0)==-1);

    //compareTo by chromosome, start, end, name (1)
    seg3.setChromosome(2);
	seg4.setChromosome(2);
	CPPUNIT_ASSERT(seg3.compareTo(seg4,1)==0);
	seg3.setChromosome(1);
	seg4.setChromosome(2);
	CPPUNIT_ASSERT(seg3.compareTo(seg4,1)==-1);
	seg3.setChromosome(3);
	seg4.setChromosome(2);
	CPPUNIT_ASSERT(seg3.compareTo(seg4,1)==1);
	seg4.setChromosome(3);
	//StartPosition
	seg3.setStartPosition(1000);
	seg4.setStartPosition(1000);
	CPPUNIT_ASSERT(seg3.compareTo(seg4,1)==0);
	seg3.setStartPosition(100);
	seg4.setStartPosition(1000);
	CPPUNIT_ASSERT(seg3.compareTo(seg4,1)==-1);
	seg3.setStartPosition(2000);
	seg4.setStartPosition(1000);
	CPPUNIT_ASSERT(seg3.compareTo(seg4,1)==1);
	seg4.setStartPosition(2000);
	//EndPosition
	seg3.setEndPosition(1000);
	seg4.setEndPosition(1000);
	CPPUNIT_ASSERT(seg3.compareTo(seg4,1)==0);
	seg3.setEndPosition(100);
	seg4.setEndPosition(1000);
	CPPUNIT_ASSERT(seg3.compareTo(seg4,1)==-1);
	seg3.setEndPosition(2000);
	seg4.setEndPosition(1000);
	CPPUNIT_ASSERT(seg3.compareTo(seg4,1)==1);
	seg4.setEndPosition(2000);
	//SegmentName
	seg3.setSegmentName("test");
	seg4.setSegmentName("test");
	CPPUNIT_ASSERT(seg3.compareTo(seg4,0)==0);
	seg3.setSegmentName("Test");
	CPPUNIT_ASSERT(seg3.compareTo(seg4,0)==-1);
	seg3.setSegmentName("vest");
	CPPUNIT_ASSERT(seg3.compareTo(seg4,0)==1);
}

