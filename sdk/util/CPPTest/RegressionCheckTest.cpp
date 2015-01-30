////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

#include "util/AffxString.h"
#include "util/RegressionCheck.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
//
#include "util/CPPTest/Setup.h"

using namespace std;

#define NaN numeric_limits<double>::quiet_NaN()
#define INF numeric_limits<double>::infinity()

/**
 * @class RegressionCheckTest
 * @brief cppunit class for testing functions from Util/egressionCheck class.
 *last change by rsatin on 09/15/09
 */

class RegressionCheckTest : public CppUnit::TestFixture
{
  // create nested class with protected members of abstract class exposed
  class RegressionCheckImp : RegressionCheck {
  public:
  virtual bool check(std::string &msg) {return true;}
  void checkFloatTest(float gold, float gen, double eps, bool &success, double &maxDiff, bool log = false, double frac = 0) {
    checkFloat(gold, gen, eps, success, maxDiff, log, frac);
    }
  void setMaxErrorTest(int max) { 
	setMaxError(max); 
    }
  void checkMsgTest(bool condition, const std::string &msg, std::string &summary) {
    checkMsg(condition, msg, summary);
    }
  void reportErrorTest(const std::string &err) {
    reportError(err);
    }
  };

  CPPUNIT_TEST_SUITE( RegressionCheckTest );
  CPPUNIT_TEST( test_checkFloat_finite );
  CPPUNIT_TEST( test_checkFloat_inf_nan );
  CPPUNIT_TEST( test_checkFloat_maxDiff );
  CPPUNIT_TEST( test_setMaxError );
  CPPUNIT_TEST_SUITE_END();

public:
  
  virtual void tearDown();
  void test_checkFloat_finite();
  void test_checkFloat_inf_nan();
  void test_checkFloat_maxDiff();
  void test_setMaxError();

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RegressionCheckTest );

void RegressionCheckTest::tearDown () 
{
  // clean up after failed test cases (pops any dead objects from stack)
  while( Verbose::getParam().m_MsgHandler.size() > 1 )
    Verbose::popMsgHandler();
}

void RegressionCheckTest::test_checkFloat_finite()
{
  cout<<endl;
  Verbose::out(1, "***RegressionCheck testcase***");
  Verbose::out(1, "RegressionCheckTest::test_checkFloat_finite");

  bool success;
  double maxDiff = 0;
  bool log = false;

  // positive test case, fractional tolerance allows all differences
  RegressionCheckImp regressionCheck1;

  success = true;
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0029, success, maxDiff, log );
  CPPUNIT_ASSERT( success );

  success = true;
  regressionCheck1.checkFloatTest( (float)4.17954, (float)4.17671, 0.0029, success, maxDiff, log );
  CPPUNIT_ASSERT( success );

  success = true;
  regressionCheck1.checkFloatTest( (float)-4.17671, (float)-4.17954, 0.0029, success, maxDiff, log );
  CPPUNIT_ASSERT( success );

  success = true;
  regressionCheck1.checkFloatTest( (float)-4.17671, (float)4.17954, 0.0029, success, maxDiff, log );
  CPPUNIT_ASSERT( !success );

  success = true;
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0028, success, maxDiff, log );
  CPPUNIT_ASSERT( !success );

  success = true;
  regressionCheck1.checkFloatTest( (float)-4.17671, (float)-4.17954, 0.0028, success, maxDiff, log );
  CPPUNIT_ASSERT( !success );

  success = true;
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0001, success, maxDiff, log );
  CPPUNIT_ASSERT( !success );

  success = true;
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0001, success, maxDiff, log, 0.00001 );
  CPPUNIT_ASSERT( !success );

  success = true;
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0001, success, maxDiff, log, 0.00065 );
  CPPUNIT_ASSERT( !success );

  success = true;
  regressionCheck1.checkFloatTest( (float)-4.17954, (float)4.17671, 0.0001, success, maxDiff, log, 0.00065 );
  CPPUNIT_ASSERT( !success );

  success = true;
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0001, success, maxDiff, log, 0.0007 );
  CPPUNIT_ASSERT( success );

  success = true;
  regressionCheck1.checkFloatTest( (float)-4.17671, (float)-4.17954, 0.0001, success, maxDiff, log, 0.0007 );
  CPPUNIT_ASSERT( success );

  success = true;
  regressionCheck1.checkFloatTest( (float)4.17671, (float)-4.17954, 0.0001, success, maxDiff, log, 0.0007 );
  CPPUNIT_ASSERT( !success );

  return;
}

void RegressionCheckTest::test_checkFloat_inf_nan()
{
  cout<<endl;
  Verbose::out(1, "RegressionCheckTest::test_checkFloat_inf_nan");

  bool success = true;
  double maxDiff = 0;
  bool log = false;

  RegressionCheckImp regressionCheck1;

  // positive test cases: same non-finite values

  success = true;
  regressionCheck1.checkFloatTest( (float)NaN, (float)NaN, 0.001, success, maxDiff, log );
  CPPUNIT_ASSERT( success );

  success = true;
  regressionCheck1.checkFloatTest( (float)INF, (float)INF, 0.001, success, maxDiff, log );
  CPPUNIT_ASSERT( success );

  success = true;
  regressionCheck1.checkFloatTest( (float)-INF, (float)-INF, 0.001, success, maxDiff, log );
  CPPUNIT_ASSERT( success );

  // negative test cases: different non-finite values

  success = true;
  regressionCheck1.checkFloatTest( (float)NaN, (float)INF, 0.001, success, maxDiff, log );
  CPPUNIT_ASSERT( !success );

  success = true;
  regressionCheck1.checkFloatTest( (float)-INF, (float)INF, 0.001, success, maxDiff, log );
  CPPUNIT_ASSERT( !success );

  // negative test cases: non-finite and finite values cannot be equivalent

  success = true;
  regressionCheck1.checkFloatTest( (float)NaN, (float)0, 0.001, success, maxDiff, log );
  CPPUNIT_ASSERT( !success );

  success = true;
  regressionCheck1.checkFloatTest( (float)100, (float)INF, 0.001, success, maxDiff, log );
  CPPUNIT_ASSERT( !success );

  return;
}

void RegressionCheckTest::test_checkFloat_maxDiff()
{
  cout<<endl;
  Verbose::out(1, "RegressionCheckTest::test_checkFloat_maxDiff");

  bool success = true;
  double maxDiff = 0;
  bool log = false;

  RegressionCheckImp regressionCheck1;

  // test case: small difference between positive values
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0029, success, maxDiff, log=false );
  CPPUNIT_ASSERT( fabs( maxDiff - (4.17954-4.17671) ) < 0.001 );   // allow small rounding errors 

  // test case: larger difference between negative values
  // APT-991
#ifdef _WIN32
  regressionCheck1.checkFloatTest( (float)-100, (float)-101, 0.0029, success, maxDiff, log=false );
#else
  regressionCheck1.checkFloatTest( (float)-100, (float)-101, 0.0029, success, maxDiff, log=true );
#endif
  CPPUNIT_ASSERT( fabs( maxDiff - (101-100) ) < 0.001 );           // allow small rounding errors 

  // test case: larger difference between values
  regressionCheck1.checkFloatTest( (float)-10, (float)0, 0.0029, success, maxDiff, log=false );
  CPPUNIT_ASSERT( fabs( maxDiff - 10 ) < 0.001 );                // allow small rounding errors 

  // test case: smaller difference does not affect maxDiff
  double maxDiff_prev = maxDiff;
  // APT-991
#ifdef _WIN32
  regressionCheck1.checkFloatTest( (float)0, (float)0, 0.0029, success, maxDiff, log=false );
#else
  regressionCheck1.checkFloatTest( (float)0, (float)0, 0.0029, success, maxDiff, log=true );
#endif
  CPPUNIT_ASSERT( maxDiff == maxDiff_prev );

  // test case: non-finite difference does not affect maxDiff (TODO rsatin Is this ok?)
  maxDiff_prev = maxDiff;
  regressionCheck1.checkFloatTest( (float)0, (float)INF, 0.001, success, maxDiff, log=false );
  CPPUNIT_ASSERT( maxDiff == maxDiff_prev );

  return;
}

void RegressionCheckTest::test_setMaxError()
{
  cout<<endl;
  Verbose::out(1, "RegressionCheckTest::test_setMaxError");

  // message output streamed to string object
  ostringstream MessageHandler_oss;
  MsgStream msgHandler1(2, &MessageHandler_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  double maxDiff = 0;
  bool log;

  RegressionCheckImp regressionCheck1;

  int error_count_max = 3;
  regressionCheck1.setMaxErrorTest( error_count_max );

  // test case: insure no messages and no effect on error count if log==false
  bool success = true;
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0001, success, maxDiff, log=false );
  regressionCheck1.checkFloatTest( (float)NaN,     (float)INF,     0.0001, success, maxDiff, log=false );

  CPPUNIT_ASSERT( !success );
  CPPUNIT_ASSERT( MessageHandler_oss.str().find("Error encountered:") == string::npos ); 

  // APT-991
#ifndef _WIN32
  // test case: insure messages and effect on error count if log==true and error_count<error_count_max
  success = true;
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0001, success, maxDiff, log=true );
  regressionCheck1.checkFloatTest( (float)4.17671, (float)4.17954, 0.0001, success, maxDiff, log=true );

  CPPUNIT_ASSERT( !success );
  CPPUNIT_ASSERT( MessageHandler_oss.str().find("Error encountered:") != string::npos ); 
  CPPUNIT_ASSERT( MessageHandler_oss.str().find("Maximum number of errors reported.") == string::npos ); 

  // test case: insure messages stop and single warning after error_count>error_count_max
  success = true;
  regressionCheck1.checkFloatTest( (float)NaN, (float)INF, 0.0001, success, maxDiff, log=true );
  regressionCheck1.checkFloatTest( (float)NaN, (float)INF, 0.0001, success, maxDiff, log=true );
  regressionCheck1.checkFloatTest( (float)NaN, (float)INF, 0.0001, success, maxDiff, log=true );

  // validate number of errors reported matches maximum allowed
  int find_count = 0;
  size_t fpos = 0;
  while( fpos != string::npos ) {
    fpos = MessageHandler_oss.str().find("Error encountered:",fpos);
	if(fpos != string::npos) {
      fpos++;
      find_count++;
      }
    }
  CPPUNIT_ASSERT( find_count == error_count_max );

  // validate warning is reported exactly once
  CPPUNIT_ASSERT( !success );
  CPPUNIT_ASSERT( (fpos = MessageHandler_oss.str().find("Maximum number of errors reported."),0) != string::npos );
  CPPUNIT_ASSERT( MessageHandler_oss.str().find("Maximum number of errors reported.",fpos+1) == string::npos );
#endif

  Verbose::popMsgHandler();

  return;
}
