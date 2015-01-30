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

#include "portability/CPPTest/PortabilityTest.h"
//
#include "portability/affy-base-types.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

CPPUNIT_TEST_SUITE_REGISTRATION( PortabilityTest );

void 
PortabilityTest::setUp()
{
}

void 
PortabilityTest::tearDown()
{
}

void 
PortabilityTest::test_Portability () 
{
  // sizes?
  CPPUNIT_ASSERT(sizeof(  int8_t)==1);
  CPPUNIT_ASSERT(sizeof( uint8_t)==1);
  CPPUNIT_ASSERT(sizeof( int16_t)==2);
  CPPUNIT_ASSERT(sizeof(uint16_t)==2);
  CPPUNIT_ASSERT(sizeof( int32_t)==4);
  CPPUNIT_ASSERT(sizeof(uint32_t)==4);
  CPPUNIT_ASSERT(sizeof( int64_t)==8);
  CPPUNIT_ASSERT(sizeof(uint64_t)==8);
}


