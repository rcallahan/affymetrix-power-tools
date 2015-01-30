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
#include "calvin_files/array/test/ArrayAuditTest.h"
//
#include "calvin_files/array/src/ArrayAudit.h"
#include "calvin_files/array/src/ArrayAuditActionTypes.h"
#include "calvin_files/utils/src/DateTime.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_array;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( ArrayAuditEntryTest );

void ArrayAuditEntryTest::setUp()
{
}

void ArrayAuditEntryTest::tearDown()
{
}

void ArrayAuditEntryTest::testCreation()
{
	ArrayAuditEntryTest entry;
	CPPUNIT_ASSERT(1);
}

void ArrayAuditEntryTest::testproperty_UserName()
{
	ArrayAuditEntry entry;
	wstring user = L"user_name";
	entry.UserName() = user;
	CPPUNIT_ASSERT( entry.UserName() == user);
}

void ArrayAuditEntryTest::testproperty_DateTime()
{
	ArrayAuditEntry entry;
	DateTime dt = DateTime::GetCurrentDateTime();
	entry.DateTime() = dt;
	CPPUNIT_ASSERT( entry.DateTime().Date() == dt.Date());
	CPPUNIT_ASSERT( entry.DateTime().Time() == dt.Time());
}

void ArrayAuditEntryTest::testproperty_ArrayGuid()
{
	ArrayAuditEntry entry;
	entry.ArrayGuid() = "guid";
	CPPUNIT_ASSERT( entry.ArrayGuid() == "guid");
}

void ArrayAuditEntryTest::testmethod_Clear()
{
	ArrayAuditEntry entry;
	AffymetrixGuidType guid; 	 

	entry.ArrayGuid() = "guid";

	guid = "guid1"; 	 
	entry.InputFileGuids().push_back(guid); 	 
	guid = "guid2"; 	 
	entry.InputFileGuids().push_back(guid);

	wstring user = L"user_name";
	entry.UserName() = user;

	DateTime dt = DateTime::GetCurrentDateTime();
	entry.DateTime() = dt;

	guid = "guid1"; 	 
	entry.OutputFileGuids().push_back(guid); 	 
	guid = "guid2"; 	 
	entry.OutputFileGuids().push_back(guid);

	int32_t actionType = 123;
	entry.ActionType() = actionType;

	ParameterNameValuePair param;
	ParameterNameValuePairList::iterator it;

	param.Name = L"name1";
	param.Value = L"value1";
	entry.ActionParameters().push_back(param);
	param.Name = L"name2";
	param.Value = L"value2";
	entry.ActionParameters().push_back(param);

	entry.Clear();

	CPPUNIT_ASSERT( entry.UserName() == L"" );
	CPPUNIT_ASSERT( entry.ArrayGuid() == "" );
	CPPUNIT_ASSERT( entry.DateTime().Date() == L"" );
	CPPUNIT_ASSERT( entry.DateTime().Time() == L"" );
	CPPUNIT_ASSERT( entry.InputFileGuids().size() == 0 ); 	 
	CPPUNIT_ASSERT( entry.OutputFileGuids().size() == 0 );
	CPPUNIT_ASSERT( entry.ActionParameters().size() == 0 );
	CPPUNIT_ASSERT( entry.ActionType() == "" );
}
	 
void ArrayAuditEntryTest::testproperty_InputFileGuids() 	 
{ 	 
	ArrayAuditEntry entry; 	 
	AffymetrixGuidType guid; 	 
	AffymetrixGuidTypeList::iterator it; 	 

	guid = "guid1"; 	 
	entry.InputFileGuids().push_back(guid); 	 
	guid = "guid2"; 	 
	entry.InputFileGuids().push_back(guid); 	 

	CPPUNIT_ASSERT( entry.InputFileGuids().size() == 2 ); 	 
	it = entry.InputFileGuids().begin(); 	 
	guid = *it; 	 
	CPPUNIT_ASSERT( guid == AffymetrixGuidType("guid1") ); 	 
	++it; 	 
	guid = *it; 	 
	CPPUNIT_ASSERT( guid == AffymetrixGuidType("guid2") ); 	 
	++it; 	 
	CPPUNIT_ASSERT( it == entry.InputFileGuids().end() ); 	 
} 	 

void ArrayAuditEntryTest::testproperty_OutputFileGuids() 	 
{ 	 
	ArrayAuditEntry entry; 	 
	AffymetrixGuidType guid; 	 
	AffymetrixGuidTypeList::iterator it; 	 

	guid = "guid1"; 	 
	entry.OutputFileGuids().push_back(guid); 	 
	guid = "guid2"; 	 
	entry.OutputFileGuids().push_back(guid); 	 

	CPPUNIT_ASSERT( entry.OutputFileGuids().size() == 2 ); 	 
	it = entry.OutputFileGuids().begin(); 	 
	guid = *it; 	 
	CPPUNIT_ASSERT( guid == AffymetrixGuidType("guid1") ); 	 
	++it; 	 
	guid = *it; 	 
	CPPUNIT_ASSERT( guid == AffymetrixGuidType("guid2") ); 	 
	++it; 	 
	CPPUNIT_ASSERT( it == entry.OutputFileGuids().end() ); 	 
} 	 

void ArrayAuditEntryTest::testproperty_ActionType()
{
	ArrayAuditEntry entry;
	string actionType = "some_action";
	entry.ActionType() = actionType;
	CPPUNIT_ASSERT( entry.ActionType() == actionType );
}

void ArrayAuditEntryTest::testproperty_ActionParameters()
{
	ArrayAuditEntry entry;
	ParameterNameValuePair param;
	ParameterNameValuePairList::iterator it;

	param.Name = L"name1";
	param.Value = L"value1";
	entry.ActionParameters().push_back(param);
	param.Name = L"name2";
	param.Value = L"value2";
	entry.ActionParameters().push_back(param);

	CPPUNIT_ASSERT( entry.ActionParameters().size() == 2 );
	it = entry.ActionParameters().begin();
	param = *it;
	CPPUNIT_ASSERT( param.Name == L"name1" );
	CPPUNIT_ASSERT( param.Value == L"value1" );
	++it;
	param = *it;
	CPPUNIT_ASSERT( param.Name == L"name2" );
	CPPUNIT_ASSERT( param.Value == L"value2" );
	++it;
	CPPUNIT_ASSERT( it == entry.ActionParameters().end() );
}

