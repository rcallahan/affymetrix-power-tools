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
#include "calvin_files/array/test/ArrayDataTest.h"
//
#include "calvin_files/array/src/ArrayData.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_array;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( ArrayDataTest );

void ArrayDataTest::setUp()
{
}

void ArrayDataTest::tearDown()
{
}

void ArrayDataTest::testCreation()
{
	ArrayData array;
	CPPUNIT_ASSERT(1);
}

void ArrayDataTest::testproperty_Properties()
{
	string id = "file_id";
	ArrayData array;
	array.ArraySetFileIdentifier() = id;
	CPPUNIT_ASSERT(array.ArraySetFileIdentifier() == id);

	array.DataTypeIdentifier() = "affymetrix-calvin-arraysetfile";
	CPPUNIT_ASSERT(array.DataTypeIdentifier() == "affymetrix-calvin-arraysetfile");
	array.CreatedStep() = ArrayRegistrationStep;
	CPPUNIT_ASSERT(array.CreatedStep() == ArrayRegistrationStep);
	array.InitialProject() = L"none";
	CPPUNIT_ASSERT(array.InitialProject() == L"none");
	array.CreationDateTime() = L"now";
	CPPUNIT_ASSERT(array.CreationDateTime() == L"now");
	array.CreatedBy() = L"me";
	CPPUNIT_ASSERT(array.CreatedBy() == L"me");
}

void ArrayDataTest::testmethod_Clear()
{
	string id = "file_id";
	ArrayData array;
	array.ArraySetFileIdentifier() = id;
	array.PhysicalArraysAttributes().resize(2);
	array.UserAttributes().resize(4);
	array.DataTypeIdentifier() = "affymetrix-calvin-arraysetfile";
	array.CreatedStep() = ArrayRegistrationStep;
	array.InitialProject() = L"none";
	array.CreationDateTime() = L"now";
	array.CreatedBy() = L"me";
	array.Clear();

	CPPUNIT_ASSERT(array.ArraySetFileIdentifier() == "");
	CPPUNIT_ASSERT(array.PhysicalArraysAttributes().size() == 0);
	CPPUNIT_ASSERT(array.UserAttributes().size() == 0);
	CPPUNIT_ASSERT(array.DataTypeIdentifier() == "");
	CPPUNIT_ASSERT(array.CreatedStep() == NoStep);
	CPPUNIT_ASSERT(array.InitialProject() == L"");
	CPPUNIT_ASSERT(array.CreationDateTime() == L"");
	CPPUNIT_ASSERT(array.CreatedBy() == L"");
}

void ArrayDataTest::testproperty_PhysicalArraysAttributes()
{
	ArrayData array;
	array.PhysicalArraysAttributes().resize(1);
	array.PhysicalArraysAttributes()[0].Identifier() = "array1-id";
	ParameterNameValuePairVector &params = array.PhysicalArraysAttributes()[0].Attributes();
	ParameterNameValuePair param;

	param.Name = L"array-att-name-1";
	param.Value = L"array-att-value-1";
	params.push_back(param);

	param.Name = L"array-att-name-2";
	param.Value = L"array-att-value-2";
	params.push_back(param);

	CPPUNIT_ASSERT( params.size() == 2 );

	CPPUNIT_ASSERT( array.PhysicalArraysAttributes().size() == 1 );
	CPPUNIT_ASSERT( array.PhysicalArraysAttributes()[0].Identifier() == "array1-id");

	ParameterNameValuePairVector::iterator it = array.PhysicalArraysAttributes()[0].Attributes().begin();
	CPPUNIT_ASSERT( (*it).Name == L"array-att-name-1" );
	CPPUNIT_ASSERT( (*it).Value == L"array-att-value-1" );
	++it;
	CPPUNIT_ASSERT( (*it).Name == L"array-att-name-2" );
	CPPUNIT_ASSERT( (*it).Value == L"array-att-value-2" );
	++it;
	CPPUNIT_ASSERT (it == array.PhysicalArraysAttributes()[0].Attributes().end() );

	CPPUNIT_ASSERT( array.PhysicalArraysAttributes()[0].Attributes()[0].Name == L"array-att-name-1" );
	CPPUNIT_ASSERT( array.PhysicalArraysAttributes()[0].Attributes()[0].Value == L"array-att-value-1" );
	CPPUNIT_ASSERT( array.PhysicalArraysAttributes()[0].Attributes()[1].Name == L"array-att-name-2" );
	CPPUNIT_ASSERT( array.PhysicalArraysAttributes()[0].Attributes()[1].Value == L"array-att-value-2" );
}

void ArrayDataTest::testproperty_UserAttributes()
{

	ArrayData array;
	ParameterNameValueDefaultRequiredTypeList &params = array.UserAttributes();
	ParameterNameValueDefaultRequiredType vparam;

	vparam.SetName(L"user-att-name-1");
	vparam.SetValueAscii("user-att-value-1");
	vparam.SetDefaultValueAscii("user-att-defaultvalue-1");
	vparam.RequiredFlag() = false;
	vparam.ControlledVocabulary().clear() ;
	params.push_back(vparam);

	vparam.SetName(L"user-att-name-2");
	vparam.SetValueText(L"user-att-value-2");
	vparam.SetDefaultValueText(L"user-att-defaultvalue-2");
	vparam.RequiredFlag() = true;
	vparam.ControlledVocabulary().clear()  ;
	params.push_back(vparam);

	vparam.SetName(L"user-att-name-3");
	vparam.SetValueText(L"user-att-value-3");
	vparam.SetDefaultValueText(L"user-att-defaultvalue-3");
	vparam.RequiredFlag() = false;
	vparam.ControlledVocabulary().push_back(L"control1");
	vparam.ControlledVocabulary().push_back(L"control2");
	params.push_back(vparam);

	CPPUNIT_ASSERT( array.UserAttributes().size() == 3 );

	ParameterNameValueDefaultRequiredTypeList::iterator it = array.UserAttributes().begin();

	CPPUNIT_ASSERT( (*it).GetName() == L"user-att-name-1" );
	CPPUNIT_ASSERT( (*it).GetValueAscii() == "user-att-value-1" );
	CPPUNIT_ASSERT( (*it).GetDefaultValueAscii() == "user-att-defaultvalue-1" );
	CPPUNIT_ASSERT( (*it).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*it).RequiredFlag() == false);
	++it;
	CPPUNIT_ASSERT( (*it).GetName() == L"user-att-name-2" );
	CPPUNIT_ASSERT( (*it).GetValueText() == L"user-att-value-2" );
	CPPUNIT_ASSERT( (*it).GetDefaultValueText() == L"user-att-defaultvalue-2" );
	CPPUNIT_ASSERT( (*it).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*it).RequiredFlag() == true);
	++it;
	CPPUNIT_ASSERT( (*it).GetName() == L"user-att-name-3" );
	CPPUNIT_ASSERT( (*it).GetValueText() == L"user-att-value-3" );
	CPPUNIT_ASSERT( (*it).GetDefaultValueText() == L"user-att-defaultvalue-3" );
	CPPUNIT_ASSERT( (*it).ControlledVocabulary().size() == 2 );
	CPPUNIT_ASSERT( (*it).RequiredFlag() == false);

	std::list<std::wstring>::iterator controllIt = (*it).ControlledVocabulary().begin();
	CPPUNIT_ASSERT( (*controllIt) == L"control1" );
	++controllIt;
	CPPUNIT_ASSERT( (*controllIt) == L"control2" );
	++it;
	CPPUNIT_ASSERT (it == array.UserAttributes().end() );

}
