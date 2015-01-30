////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

///
/// @file   ProbeSetTest.cpp
/// @brief  Testing the ProbeSet functions.

//
#include "chipstream/ProbeSet.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @class ProbeSetTest
 * @brief cppunit class for testing conversion functions.
 */
class ProbeSetTest : public CppUnit::TestFixture {
public:
  ProbeSetTest() {};

  CPPUNIT_TEST_SUITE( ProbeSetTest );
  CPPUNIT_TEST( probeset_test_zerobased );
  CPPUNIT_TEST( probeset_test_2 );
  CPPUNIT_TEST_SUITE_END();

  ///
  void probeset_test_zerobased();
  void probeset_test_2();
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( ProbeSetTest );

#define TZB(_x) printf("TZB(%d)=%d\n",_x,TO_ZEROBASED(_x))
#define FZB(_x) printf("FZB(%d)=%d\n",_x,FROM_ZEROBASED(_x))

// TZB(-3) = -3
// TZB(-2) = -2
// TZB(-1) = -1
// TZB(0)  = -1
// TZB(1)  = 0
// TZB(2)  = 1
// TZB(3)  = 2
// FZB(-3) = -3
// FZB(-2) = -2
// FZB(-1) = -1
// FZB(0)  = 1
// FZB(1)  = 2
// FZB(2)  = 3
// FZB(3)  = 4

void
ProbeSetTest::probeset_test_zerobased()
{
  // for (int i=-3;i<4;i++) { TZB(i); }
  // for (int i=-3;i<4;i++) { FZB(i); }

  // convert to 1-based (external) to 0-based.
  CPPUNIT_ASSERT(TO_ZEROBASED(-2)==-2);
  CPPUNIT_ASSERT(TO_ZEROBASED(-1)==-1);
  CPPUNIT_ASSERT(TO_ZEROBASED(0) ==-1);
  CPPUNIT_ASSERT(TO_ZEROBASED(1) == 0);
  CPPUNIT_ASSERT(TO_ZEROBASED(9) == 8);
  // convert from 0-based (internal) to 1-based (external)
  CPPUNIT_ASSERT(FROM_ZEROBASED(-2)==-2);
  CPPUNIT_ASSERT(FROM_ZEROBASED(-1)==-1);
  CPPUNIT_ASSERT(FROM_ZEROBASED(0) == 1);
  CPPUNIT_ASSERT(FROM_ZEROBASED(1) == 2);
  CPPUNIT_ASSERT(FROM_ZEROBASED(8) == 9);
  //
  printf("probeset_test_zerobased: ok.\n");
}

void
ProbeSetTest::probeset_test_2()
{
  //printf("ok.\n");
}
