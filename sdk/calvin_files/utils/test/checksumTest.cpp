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

#include "calvin_files/utils/test/checksumTest.h"
//
#include "calvin_files/utils/src/checksum.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( CheckSumTest );

using namespace affymetrix_calvin_utilities;

void CheckSumTest::setUp()
{
}

void CheckSumTest::tearDown()
{
}

void CheckSumTest::testmethod_CheckSumTest_even_number_of_points()
{
	u_int16_t data1[] = {10, 20, 30, 40};
	size_t size = sizeof(data1);
	u_int16_t csum = CheckSum::OnesComplementCheckSum(data1, size);
	CPPUNIT_ASSERT( csum == 65435 );

	u_int16_t data2[] = {10, 20, 30, 40, 65435};
	size = sizeof(data2);
	csum = CheckSum::OnesComplementCheckSum(data2, size);
	CPPUNIT_ASSERT( csum == 0 );
}

void CheckSumTest::testmethod_CheckSumTest_odd_number_of_points()
{
	u_int16_t data1[] = {10, 20, 30, 40, 50};
	size_t size = sizeof(data1);
	u_int16_t csum = CheckSum::OnesComplementCheckSum(data1, size);
	CPPUNIT_ASSERT( csum == 65385 );

	u_int16_t data2[] = {10, 20, 30, 40, 50, 65385};
	size = sizeof(data2);
	csum = CheckSum::OnesComplementCheckSum(data2, size);
	CPPUNIT_ASSERT( csum == 0 );
}

void CheckSumTest::testmethod_CheckSumTest_overflow()
{
	u_int16_t data1[] = {10000, 20000, 30000, 40000, 50000};
	size_t size = sizeof(data1);
	u_int16_t csum = CheckSum::OnesComplementCheckSum(data1, size);
	CPPUNIT_ASSERT( csum == 46605 );

	u_int16_t data2[] = {10000, 20000, 30000, 40000, 50000, 46605};
	size = sizeof(data2);
	csum = CheckSum::OnesComplementCheckSum(data2, size);
	CPPUNIT_ASSERT( csum == 0 );
}

