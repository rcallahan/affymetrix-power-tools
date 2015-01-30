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

#include "calvin_files/parameter/test/ParameterNameValueControlVocabularyTest.h"
//
#include "calvin_files/parameter/src/Parameter.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( ParameterNameValueControlVocabularyTest );

using namespace affymetrix_calvin_parameter;
using namespace std;

void ParameterNameValueControlVocabularyTest::setUp()
{
}

void ParameterNameValueControlVocabularyTest::tearDown()
{
}

void ParameterNameValueControlVocabularyTest::test_assignment()
{
	ParameterNameValueControlVocabulary param1;
	ParameterNameValueControlVocabulary param2;
	wstring name = L"name";
	wstring value = L"value";
	wstring control1 = L"control1";
	wstring control2 = L"control2";

	param1.Name = name;
	param1.Value = value;
	param1.ControlledVocabulary.push_back(control1);
	param1.ControlledVocabulary.push_back(control2);

	param2 = param1;

	CPPUNIT_ASSERT( param2.Name == param1.Name );
	CPPUNIT_ASSERT( param2.Value == param1.Value );
	CPPUNIT_ASSERT( param2.ControlledVocabulary[0] == param1.ControlledVocabulary[0] );
	CPPUNIT_ASSERT( param2.ControlledVocabulary[1] == param1.ControlledVocabulary[1] );
	CPPUNIT_ASSERT( param2.ControlledVocabulary.size() == param1.ControlledVocabulary.size() );
}

void ParameterNameValueControlVocabularyTest::test_equality_to_other_object()
{
	ParameterNameValueControlVocabulary param1;
	ParameterNameValueControlVocabulary param2;
	wstring name = L"name";
	wstring value = L"value";

	param1.Name = name;
	param1.Value = value;

	param2.Name = name;
	param2.Value = value;

	CPPUNIT_ASSERT( param2 == param1 );
}

void ParameterNameValueControlVocabularyTest::test_equality_to_string()
{
	ParameterNameValueControlVocabulary param;
	wstring name = L"name";
	wstring value = L"value";

	param.Name = name;
	param.Value = value;

	CPPUNIT_ASSERT( param == name );
}

void ParameterNameValueControlVocabularyTest::test_accessibility()
{
	ParameterNameValueControlVocabulary param;
	wstring name = L"name";
	wstring value = L"value";
	wstring control1 = L"control1";
	wstring control2 = L"control2";

	param.Name = name;
	param.Value = value;
	param.ControlledVocabulary.push_back(control1);
	param.ControlledVocabulary.push_back(control2);

	CPPUNIT_ASSERT( param.Name == name );
	CPPUNIT_ASSERT( param.Value == value );
	CPPUNIT_ASSERT( param.ControlledVocabulary[0] == control1 );
	CPPUNIT_ASSERT( param.ControlledVocabulary[1] == control2 );

}
