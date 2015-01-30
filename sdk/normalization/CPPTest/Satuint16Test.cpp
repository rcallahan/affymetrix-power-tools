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

#include "normalization/CPPTest/Satuint16Test.h"
//
#include "normalization/sat_uint16.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

CPPUNIT_TEST_SUITE_REGISTRATION( Satuint16Test );

void Satuint16Test::setUp()
{
}

void Satuint16Test::tearDown()
{
}

void Satuint16Test::test_Satuint16 () {
  sat_uint16 si0,si1,si2,si3;
  sat_uint16* psi0;
  int i0;
  double d0;

  //printf("Satuint16Test::test_Satuint16()...\n");
  
  // should be zero at start (extra tests to use the vars)
  CPPUNIT_ASSERT(si0==0);
  CPPUNIT_ASSERT(si1==0);
  CPPUNIT_ASSERT(si2==0);
  CPPUNIT_ASSERT(si3==0);

  si0=1;si1=1;
  CPPUNIT_ASSERT(si0==si1);

  si0=-1; si1=si0; si1++;
  CPPUNIT_ASSERT(si0==0);
  CPPUNIT_ASSERT(si0!=si1);
  CPPUNIT_ASSERT(si1==1);

  si0=65535; si1=si0; si1++;
  CPPUNIT_ASSERT(si0==65535);
  CPPUNIT_ASSERT(si0==si1);

  //
  psi0=&si0;
  *psi0=10;
  CPPUNIT_ASSERT(si0==10);

  // does casting work?
  i0=100; si0=i0;
  CPPUNIT_ASSERT(si0==i0);
  // doubles too?
  d0=100.0; si0=d0;
  CPPUNIT_ASSERT(si0==i0);

  // test the test
  //CPPUNIT_ASSERT(0==1);
}

