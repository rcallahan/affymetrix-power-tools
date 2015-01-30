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

#include "copynumber/CNReporterCnchp.h"
#include "copynumber/CNReporterCychp.h"
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
 * @class CNReporterTest
 * @brief cppunit class for testing CNReporter functions.
 * last change by vliber on 01/22/09
 */

class CNReporterTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNReporterTest);
  CPPUNIT_TEST(CNReporterConstructorTest); 
  CPPUNIT_TEST(CNReporterDefineOptionTest); 
  CPPUNIT_TEST(CNReporterRunTest); 
  CPPUNIT_TEST_SUITE_END();

public: 
  void CNReporterConstructorTest();
  void CNReporterDefineOptionTest();
  void CNReporterRunTest(); 
  
 
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNReporterTest );

void CNReporterTest::CNReporterConstructorTest()
{
	cout<<endl;
	Verbose::out(1, "****CNReporterTest::CNReporterConstructorTest****");
	CNReporterCychp cnCychp;
	//FATAL ERROR: CNReport is not setup properly.
   	NEGATIVE_TEST(cnCychp.isSetup(),Except);
		
	CNReporterCnchp cnCnchp;
	//FATAL ERROR: CNReport is not setup properly.
   	NEGATIVE_TEST(cnCnchp.isSetup(),Except);
	
}
void CNReporterTest::CNReporterDefineOptionTest() 
{
    ///@todo assumes old style options handling
    /*
	Verbose::out(1, "****CNReporterTest::CNReporterDefineOptionTest****");
	CNReporterCychp cnCychp;
	PgOptions *opts = new PgOptions;
   	cnCychp.defineOptions(*opts);
	PgOpt opt;
	opt=opts->findOpt("cychp-output");
    CPPUNIT_ASSERT(&opt!=NULL);
    CPPUNIT_ASSERT(opt.m_shortName=="");
    CPPUNIT_ASSERT(opt.m_type==1);
    CPPUNIT_ASSERT(opt.m_help=="Report CYCHP files");
    CPPUNIT_ASSERT(opt.getValue()=="false");
		
	CNReporterCnchp cnCnchp;
	cnCnchp.defineOptions(*opts);
	opt=opts->findOpt("cnchp-output");
    CPPUNIT_ASSERT(&opt!=NULL);
    CPPUNIT_ASSERT(opt.m_shortName=="");
    CPPUNIT_ASSERT(opt.m_type==1);
    CPPUNIT_ASSERT(opt.m_help=="Report CNCHP files");
    CPPUNIT_ASSERT(opt.getValue()=="true");
    */	
}

void CNReporterTest::CNReporterRunTest()
{
	Verbose::out(1, "****CNReporterTest::CNReporterRunTest****");
	CNReporterCychp cnCychp;
	//>FATAL ERROR: CNReport is not setup properly.
   	NEGATIVE_TEST(cnCychp.isSetup(),Except);
	
	
	CNReporterCnchp cnCnchp;
	//>FATAL ERROR: CNReport is not setup properly.
   	NEGATIVE_TEST(cnCnchp.isSetup(),Except);
	
}

