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

#include "calvin_files/parameter/test/ParameterNameValuePairTest.h"
//
#include "calvin_files/parameter/src/Parameter.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( ParameterNameValuePairTest );

using namespace affymetrix_calvin_parameter;
using namespace std;

void ParameterNameValuePairTest::setUp()
{
}

void ParameterNameValuePairTest::tearDown()
{
}


void ParameterNameValuePairTest::test_assignment()
{
	ParameterNameValuePair param1;
	ParameterNameValuePair param2;
	wstring name = L"name";
	wstring value = L"value";

	param1.Name = name;
	param1.Value = value;

	param2 = param1;

	CPPUNIT_ASSERT( param2.Name == param1.Name );
	CPPUNIT_ASSERT( param2.Value == param1.Value );
}

void ParameterNameValuePairTest::test_equality_to_other_object()
{
	ParameterNameValuePair param1;
	ParameterNameValuePair param2;
	wstring name = L"name";
	wstring value = L"value";

	param1.Name = name;
	param1.Value = value;

	param2.Name = name;
	param2.Value = value;

	CPPUNIT_ASSERT( param2 == param1 );
}

void ParameterNameValuePairTest::test_equality_to_string()
{
	ParameterNameValuePair param;
	wstring name = L"name";
	wstring value = L"value";

	param.Name = name;
	param.Value = value;

	CPPUNIT_ASSERT( param == name );
}

void ParameterNameValuePairTest::test_accessibility()
{
	ParameterNameValuePair param;
	wstring name = L"name";
	wstring value = L"value";

	param.Name = name;
	param.Value = value;

	CPPUNIT_ASSERT( param.Name == name );
	CPPUNIT_ASSERT( param.Value == value );
}


// copy ctor
