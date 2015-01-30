////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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
#include "calvin_files/portability/test/AffymetrixBaseTypesTest.h"
//
#include "calvin_files/portability/src/AffymetrixBaseTypes.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( AffymetrixBaseTypesTest );

void AffymetrixBaseTypesTest::setUp()
{
}

void AffymetrixBaseTypesTest::tearDown()
{
}

void AffymetrixBaseTypesTest::testdefine_SIZES()
{
	CPPUNIT_ASSERT( sizeof(int8_t) == 1 );
	CPPUNIT_ASSERT( sizeof(int16_t) == 2 );
	CPPUNIT_ASSERT( sizeof(int32_t) == 4 );
	CPPUNIT_ASSERT( sizeof(int64_t) == 8 );
	CPPUNIT_ASSERT( sizeof(u_int8_t) == 1 );
	CPPUNIT_ASSERT( sizeof(u_int16_t) == 2 );
	CPPUNIT_ASSERT( sizeof(u_int32_t) == 4 );
	CPPUNIT_ASSERT( sizeof(u_int64_t) == 8 );
}
