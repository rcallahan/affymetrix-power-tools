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

#include "copynumber/CNProbeSet.h"
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
 * @class CNProbeSetTest
 * @brief cppunit class for testing CNProbeSet functions.
 * last change by vliber on 03/23/09
 */

class CNProbeSetTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(CNProbeSetTest);
	CPPUNIT_TEST(compareToTest);
	CPPUNIT_TEST(getSignalStrengthTest);
	CPPUNIT_TEST(getIntensityContrastTest);
	CPPUNIT_TEST(getSignalContrastTest);
	CPPUNIT_TEST(getIntensityStrengthTest);
	CPPUNIT_TEST(getCNGainTest);
	CPPUNIT_TEST(getCNLossTest);	
	CPPUNIT_TEST_SUITE_END();

public:  
   
   void compareToTest();
   void getSignalStrengthTest();
    void getIntensityContrastTest();
   void getSignalContrastTest();
   void getIntensityStrengthTest();
  
   void getCNGainTest();
   void getCNLossTest();
     
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNProbeSetTest );


void CNProbeSetTest::compareToTest()
{
	cout<<endl;
	Verbose::out(1, "****CNProbeSetTest::compareToTest****");
	CNProbeSet cnPS1, cnPS2, cnPS3;
	cnPS1.setProbeSetName("test2");
	cnPS1.setChromosome('1');
	cnPS1.setPosition(123);
	cnPS1.setMedianSignal(0.999f);
	cnPS1.setXXMedianSignal(0.95f);
	cnPS1.setYMedianSignal(0.96f);
	cnPS1.setAAMedianSignal(0.94f);
	cnPS1.setABMedianSignal(0.93f);
	cnPS1.setBBMedianSignal(0.97f);
	//
	cnPS2.setProbeSetName("test1");
	cnPS2.setChromosome('x');
	cnPS2.setPosition(245);
	cnPS2.setMedianSignal(0.989f);
	cnPS2.setXXMedianSignal(0.85f);
	cnPS2.setYMedianSignal(0.86f);
	cnPS2.setAAMedianSignal(0.84f);
	cnPS2.setABMedianSignal(0.83f);
	cnPS2.setBBMedianSignal(0.87f);
    //
	cnPS3.setProbeSetName("pest1");
	cnPS3.setChromosome('x');
	cnPS3.setPosition(245);

	//
	CPPUNIT_ASSERT(cnPS2.getProbeSetName()=="test1");
	CPPUNIT_ASSERT(cnPS2.getChromosome()=='x');
	CPPUNIT_ASSERT(cnPS2.getPosition()==245);
	CPPUNIT_ASSERT(cnPS2.getMedianSignal()==0.989f);
	CPPUNIT_ASSERT(cnPS2.getXXMedianSignal()==0.85f);
	CPPUNIT_ASSERT(cnPS2.getYMedianSignal()==0.86f);
	CPPUNIT_ASSERT(cnPS2.getAAMedianSignal()==0.84f);
	CPPUNIT_ASSERT(cnPS2.getABMedianSignal()==0.83f);
	CPPUNIT_ASSERT(cnPS2.getBBMedianSignal()==0.87f);

	CNProbeSetArray m_arProbeSets;

	m_arProbeSets.add(&cnPS1);
	m_arProbeSets.add(&cnPS2);
	m_arProbeSets.add(&cnPS3);
	
    m_arProbeSets.quickSort(0);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getProbeSetName()=="pest1");
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getProbeSetName()=="test1");
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getProbeSetName()=="test2");
	m_arProbeSets.quickSort(1);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getChromosome()=='1');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getChromosome()=='x');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getChromosome()=='x');
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getPosition()==123);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getPosition()==245);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getPosition()==245);
	CPPUNIT_ASSERT(m_arProbeSets.getAt(0)->getProbeSetName()=="test2");
	CPPUNIT_ASSERT(m_arProbeSets.getAt(1)->getProbeSetName()=="pest1");
	CPPUNIT_ASSERT(m_arProbeSets.getAt(2)->getProbeSetName()=="test1");

	m_arProbeSets.nullAll();
}

void CNProbeSetTest::getSignalStrengthTest(){
  Verbose::out(1, "****CNProbeSetTest::getSignalStrengthTest****");
  CNProbeSet *cnPS = new CNProbeSet();
  cnPS->setAAlleleSignal(0.0);
  cnPS->setBAlleleSignal(0.0);
  //Zero median signals found. CEL file may be corrupted
  NEGATIVE_TEST(cnPS->getSignalStrength(),Except);
  delete cnPS;
}
void CNProbeSetTest::getIntensityContrastTest(){
  Verbose::out(1, "****CNProbeSetTest::getIntensityContrastTest****");
  CNProbeSet *cnPS = new CNProbeSet();
  cnPS->setAMedianIntensity(0.0);
  cnPS->setBMedianIntensity(0.0);
  //Zero median intensities found. CEL file may be corrupted
  NEGATIVE_TEST(cnPS->getIntensityStrength(),Except);
  delete cnPS;
}
void CNProbeSetTest::getSignalContrastTest(){
  Verbose::out(1, "****CNProbeSetTest::getSignalContrastTest****");
  CNProbeSet *cnPS = new CNProbeSet();
  cnPS->setAAlleleSignal(0.0);
  cnPS->setBAlleleSignal(0.0);
  //Zero median signals found. CEL file may be corrupted
  NEGATIVE_TEST(cnPS->getSignalContrast(2.0),Except);
  delete cnPS;
}
void CNProbeSetTest::getIntensityStrengthTest(){
  Verbose::out(1, "****CNProbeSetTest::getIntensityStrengthTest****");
  CNProbeSet *cnPS = new CNProbeSet();
  cnPS->setAMedianIntensity(0.0);
  cnPS->setBMedianIntensity(0.0);
  //Zero median intensities found. CEL file may be corrupted
  NEGATIVE_TEST(cnPS->getIntensityStrength(),Except);
  delete cnPS;
}

void CNProbeSetTest::getCNGainTest()
{
	Verbose::out(1, "****CNProbeSetTest::getCNGainTest****");
	CNProbeSet *cnPS1 = new CNProbeSet();
	CPPUNIT_ASSERT(cnPS1->getChromosome()==0);
	CPPUNIT_ASSERT(cnPS1->getCNState()==-1);
	CPPUNIT_ASSERT(cnPS1->getCNGain(0,1,2)==0);
    
	CNProbeSet *cnPS2 = new CNProbeSet();
    cnPS2->setChromosome(1);
	cnPS2->setCNState(2);
	CPPUNIT_ASSERT(cnPS2->getChromosome()==1);
	CPPUNIT_ASSERT(cnPS2->getCNState()==2);
	CPPUNIT_ASSERT(cnPS2->getCNGain(1,1,2)==1);
    CPPUNIT_ASSERT(cnPS2->getCNGain(1,2,1)==1);
	cnPS2->setCNState(3);
	CPPUNIT_ASSERT(cnPS2->getCNGain(0,1,2)==1);
	cnPS2->setCNState(1);
    CPPUNIT_ASSERT(cnPS2->getCNGain(0,2,1)==1);

    cnPS2->setCNState(3);
	CPPUNIT_ASSERT(cnPS2->getCNGain(0,3,3)==1);
	CPPUNIT_ASSERT(cnPS2->getCNGain(1,3,3)==1);
	delete cnPS1;
	delete cnPS2;
	
}
void CNProbeSetTest::getCNLossTest()
{
	Verbose::out(1, "****CNProbeSetTest::getCNLossTest****");
	CNProbeSet *cnPS1 = new CNProbeSet();
    CPPUNIT_ASSERT(cnPS1->getChromosome()==0);
	CPPUNIT_ASSERT(cnPS1->getCNState()==-1);
	CPPUNIT_ASSERT(cnPS1->getCNLoss(1,1,2)==1);
	cnPS1->setCNState(3);
	CPPUNIT_ASSERT(cnPS1->getCNState()==3);
	CPPUNIT_ASSERT(cnPS1->getCNLoss(1,1,2)==0);
	
	
	CNProbeSet *cnPS2 = new CNProbeSet();
    cnPS2->setChromosome(1);
	cnPS2->setCNState(0);
	CPPUNIT_ASSERT(cnPS2->getChromosome()==1);
	CPPUNIT_ASSERT(cnPS2->getCNState()==0);
	CPPUNIT_ASSERT(cnPS2->getCNLoss(1,1,2)==1);
    CPPUNIT_ASSERT(cnPS2->getCNLoss(1,2,1)==1);
	cnPS2->setCNState(1);
	CPPUNIT_ASSERT(cnPS2->getCNLoss(0,1,2)==1);
	CPPUNIT_ASSERT(cnPS2->getCNLoss(0,2,1)==0);
	delete cnPS1;
	delete cnPS2;
}
