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

#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNAnalysisMethodReference.h"
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
 * @class CNAnalysisMethodReferenceTest
 * @brief cppunit class for testing CNAnalysisMethodReference functions.
 * last change by vliber on 03/20/09
 */

class CNAnalysisMethodReferenceTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodReferenceTest);
  CPPUNIT_TEST(functionsTest);
  CPPUNIT_TEST(newObjectTest); 
  CPPUNIT_TEST(runPart1Test);
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  void newObjectTest();
  void runPart1Test();
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodReferenceTest );

void CNAnalysisMethodReferenceTest::functionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodReferenceTest::functionsTest****");
	// constructor && m_vParams()
	CNAnalysisMethodReference cnrfCNam;
    CPPUNIT_ASSERT(cnrfCNam.getVersion()=="1.0");
	CPPUNIT_ASSERT(cnrfCNam.getType()=="reference");
	CPPUNIT_ASSERT(cnrfCNam.getDescription()=="Copynumber Reference");
	CPPUNIT_ASSERT(cnrfCNam.getName()=="reference");
	
	//getDefaultDocOptions()
	vector<SelfDoc::Opt> sv=CNAnalysisMethodReference::getDefaultDocOptions();
	CPPUNIT_ASSERT(sv.size()==0);
	
}

void CNAnalysisMethodReferenceTest::newObjectTest()
{
	Verbose::out(1, "****CNAnalysisMethodReferenceTest::newObjectTest****");
     //Don't know how to make object for name: reference of type: CNAnalysisMethod
	CNAnalysisMethod::getParams()->clear();
	CNAnalysisMethodFactory obj_CNAnalysisMethodReference;
	NEGATIVE_TEST((CNAnalysisMethodReference*)obj_CNAnalysisMethodReference.CNAnalysisMethodForString("reference"),Except);
	 
}	

void CNAnalysisMethodReferenceTest::runPart1Test()
{
    Verbose::out(1, "****CNAnalysisMethodReferenceTest::runPart1Test****");
    CNAnalysisMethodReference cnrfCNam;
    //FATAL ERROR: CNAnalysisMethod reference is not setup properly.
	NEGATIVE_TEST(cnrfCNam.runPart1(),Except);	
}
