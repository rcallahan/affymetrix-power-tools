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

#include "copynumber/CNProbe.h"
#include "copynumber/CPPTest/Setup.h"
//
#include "util/AffxArray.h"
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <iostream>
//
using namespace std;

/**
 * @class CNProbeTest
 * @brief cppunit class for testing CNProbe functions.
 * last change by vliber on 03/24/09
 */

class CNProbeTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(CNProbeTest);
	CPPUNIT_TEST(compareToTest);
	CPPUNIT_TEST_SUITE_END();

public:  
   
   void compareToTest();
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNProbeTest );


void CNProbeTest::compareToTest()
{
	cout<<endl;
	Verbose::out(1, "****CNProbeTest::compareToTest****");
	CNProbe cnP1, cnP2, cnP3, cnP4;
	cnP1.setProbeSetIndex(2);
	cnP1.setAllele('a');
	cnP1.setProbeID(123);
	
	//
	cnP2.setProbeSetIndex(1);
	cnP2.setAllele('c');
	cnP2.setProbeID(129);
	cnP2.setIntensity(0.998f);
	cnP2.setResidual(0.96f);
	cnP2.setProbeEffect(0.94);
	cnP2.setUseForSketch(false);
	cnP2.setMedianIntensity(0.92f);
	cnP2.setPredictedIntensity(0.99f);
	//
	cnP3.setProbeSetIndex(1);
	cnP3.setAllele('g');
	cnP3.setProbeID(126);

	cnP4.setProbeSetIndex(2);
	cnP4.setAllele('a');
	cnP4.setProbeID(126);
	
	

	CNProbeArray m_arProbeSets;
	m_arProbeSets.add(&cnP1);
	m_arProbeSets.add(&cnP2);
	m_arProbeSets.add(&cnP3);
	m_arProbeSets.add(&cnP4);


	m_arProbeSets.quickSort(0);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getProbeSetIndex()==1);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getProbeSetIndex()==1);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getProbeSetIndex()==2);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(3)->getProbeSetIndex()==2);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getAllele()=='c');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getAllele()=='g');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getAllele()=='a');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(3)->getAllele()=='a');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getProbeID()==129);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getProbeID()==126);
	m_arProbeSets.quickSort(1);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getProbeID()==123);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getProbeID()==126);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getProbeID()==126);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(3)->getProbeID()==129);
	m_arProbeSets.quickSort(2);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getProbeSetIndex()==1);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getProbeSetIndex()==1);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getProbeSetIndex()==2);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(3)->getProbeSetIndex()==2);
	m_arProbeSets.quickSort(3);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getProbeSetIndex()==1);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getProbeSetIndex()==1);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getProbeSetIndex()==2);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(3)->getProbeSetIndex()==2);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getAllele()=='c');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getAllele()=='g');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getAllele()=='a');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(3)->getAllele()=='a');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getProbeID()==129);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getProbeID()==126);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getProbeID()==123);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(3)->getProbeID()==126);

    
	//get methods
	CPPUNIT_ASSERT(cnP2.getProbeSetIndex()==1);
	CPPUNIT_ASSERT(cnP2.getAllele()=='c');
	CPPUNIT_ASSERT(cnP2.getProbeID()==129);
	CPPUNIT_ASSERT(cnP2.getIntensity()==0.998f);
	CPPUNIT_ASSERT(cnP2.getResidual()==0.96f);
	CPPUNIT_ASSERT(cnP2.getProbeEffect()==0.94);
	CPPUNIT_ASSERT(cnP2.isUseForSketch()==false);
	CPPUNIT_ASSERT(cnP2.getMedianIntensity()==0.92f);
	CPPUNIT_ASSERT(cnP2.getPredictedIntensity()==0.99f);
	m_arProbeSets.nullAll();
}
