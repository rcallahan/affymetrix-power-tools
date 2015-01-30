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

//
#include "calvin_files/parameter/test/MIMEValueTest.h"
//
#include <cstring>
#include <string.h>
#include <string>
//

CPPUNIT_TEST_SUITE_REGISTRATION( MIMEValueTest );

using namespace affymetrix_calvin_parameter;

void MIMEValueTest::setUp()
{
	int32_t value;
	value = 76406;
	memset(arr, 0, 10);
	memcpy(arr, &value, 4);
	mimeValue = new MIMEValue(arr, 10);
}

void MIMEValueTest::tearDown()
{
	delete mimeValue;
}

void MIMEValueTest::CreationTest()
{
	CPPUNIT_ASSERT(mimeValue);
}

void MIMEValueTest::SizeTest()
{
	CPPUNIT_ASSERT(mimeValue->Size() == 10);
}

void MIMEValueTest::GetValueTest()
{
	u_int32_t size;
	const void* p = mimeValue->GetValue(size);
	CPPUNIT_ASSERT(memcmp(p, arr, 10) == 0);
	CPPUNIT_ASSERT(size == 10);
}

void MIMEValueTest::CopyConstructorTest()
{
	MIMEValue newValue(*mimeValue);
	u_int32_t size1, size2;
	const void* p1 = mimeValue->GetValue(size1);
	const void* p2 = newValue.GetValue(size2);
	CPPUNIT_ASSERT(memcmp(p1, p2, size1) == 0);
	CPPUNIT_ASSERT(size2 == 10);
	CPPUNIT_ASSERT(memcmp(mimeValue->GetValue(size1), newValue.GetValue(size2), size1) == 0);
}

void MIMEValueTest::AssignmentOperatorTest()
{
	char x[5];
	memset(x, 0, 5);
	x[0] = 4;
	MIMEValue v(x, 5);
	
	CPPUNIT_ASSERT_NO_THROW(v = *mimeValue);
	
	u_int32_t size;
	const void* p = v.GetValue(size);

	CPPUNIT_ASSERT(size == 10);
	CPPUNIT_ASSERT(memcmp(p, arr, 10) == 0);
}

void MIMEValueTest::OperatorEqualsTest()
{
	MIMEValue v(*mimeValue);
	CPPUNIT_ASSERT(v == *mimeValue);
}

void MIMEValueTest::OperatorNotEqualsTest()
{
	char x[5];
	memset(x, 0, 5);
	x[0] = 4;
	MIMEValue v(x, 5);

	CPPUNIT_ASSERT(v != *mimeValue);
}

void MIMEValueTest::SetValueTest()
{
	int32_t value;
	value = 76406;
	memset(arr, 0, 10);
	memcpy(arr, &value, 4);

	MIMEValue v;
	v.SetValue(arr, 10);

	CPPUNIT_ASSERT(v == *mimeValue);
}
