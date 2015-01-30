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
#include "calvin_files/parsers/test/AuditFileReaderTest.h"
//
#include "calvin_files/parsers/src/AuditFileReader.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_array;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

#define NO_FILE "./no_file"
#define AUDIT_FILE "../data/audit.log"

CPPUNIT_TEST_SUITE_REGISTRATION( AuditFileReaderTest );

void AuditFileReaderTest::setUp()
{
}

void AuditFileReaderTest::tearDown()
{
}

void AuditFileReaderTest::testCreation()
{
	AuditFileReader reader;
	CPPUNIT_ASSERT(1);
}

void AuditFileReaderTest::testmethod_Read()
{
	string file;
	AuditFileReader reader;
	ArrayAuditEntryList audit;
	ArrayAuditEntryList::iterator it;
	ArrayAuditEntry entry;
	AffymetrixGuidTypeList::iterator guidIt;
	ParameterNameValuePairList::iterator paramIt;
	ParameterNameValuePair param;

	file = NO_FILE;
	CPPUNIT_ASSERT(reader.Read(file, audit) == false);

	file = AUDIT_FILE;
	CPPUNIT_ASSERT(reader.Read(file, audit) == true);
	CPPUNIT_ASSERT(audit.size() == 2);
	it = audit.begin();
	entry = *it;
	CPPUNIT_ASSERT(entry.UserName() == L"ljevon");
	CPPUNIT_ASSERT(entry.DateTime().Date() == L"10-10-2005");
	CPPUNIT_ASSERT(entry.DateTime().Time() == L"12:00:00Z");
	CPPUNIT_ASSERT(entry.ActionType() == "test.id");
	CPPUNIT_ASSERT(entry.ArrayGuid() == "123-123-123-123");

	guidIt = entry.InputFileGuids().begin();
	CPPUNIT_ASSERT( (*guidIt) == "321-321-321-321" );
	++guidIt;
	CPPUNIT_ASSERT( (*guidIt) == "345-345-345-345" );
	++guidIt;
	CPPUNIT_ASSERT( guidIt ==  entry.InputFileGuids().end() );

	guidIt = entry.OutputFileGuids().begin();
	CPPUNIT_ASSERT( (*guidIt) == "999-999-999-999" );
	++guidIt;
	CPPUNIT_ASSERT( guidIt ==  entry.OutputFileGuids().end() );

	paramIt = entry.ActionParameters().begin();
	param = *paramIt;
	CPPUNIT_ASSERT(param.Name == L"Test");
	CPPUNIT_ASSERT(param.Value == L"true");
	++paramIt;
	param = *paramIt;
	CPPUNIT_ASSERT(param.Name == L"PixelSize");
	CPPUNIT_ASSERT(param.Value == L"3");
	++paramIt;
	CPPUNIT_ASSERT(paramIt == entry.ActionParameters().end());


	++it;
	entry = *it;
	CPPUNIT_ASSERT(entry.UserName() == L"ljevon1");
	CPPUNIT_ASSERT(entry.DateTime().Date() == L"10-11-2005");
	CPPUNIT_ASSERT(entry.DateTime().Time() == L"12:01:00Z");
	CPPUNIT_ASSERT(entry.ActionType() == "test2.id");
	CPPUNIT_ASSERT(entry.ArrayGuid() == "1123-123-123-123");

	guidIt = entry.InputFileGuids().begin();
	CPPUNIT_ASSERT( (*guidIt) == "1321-321-321-321" );
	++guidIt;
	CPPUNIT_ASSERT( (*guidIt) == "1345-345-345-345" );
	++guidIt;
	CPPUNIT_ASSERT( guidIt ==  entry.InputFileGuids().end() );

	guidIt = entry.OutputFileGuids().begin();
	CPPUNIT_ASSERT( (*guidIt) == "1999-999-999-999" );
	++guidIt;
	CPPUNIT_ASSERT( guidIt ==  entry.OutputFileGuids().end() );

	paramIt = entry.ActionParameters().begin();
	param = *paramIt;
	CPPUNIT_ASSERT(param.Name == L"Test");
	CPPUNIT_ASSERT(param.Value == L"false");
	++paramIt;
	param = *paramIt;
	CPPUNIT_ASSERT(param.Name == L"PixelSize");
	CPPUNIT_ASSERT(param.Value == L"1");
	++paramIt;
	param = *paramIt;
	CPPUNIT_ASSERT(param.Name == L"Temp");
	CPPUNIT_ASSERT(param.Value == L"200");
	++paramIt;
	CPPUNIT_ASSERT(paramIt == entry.ActionParameters().end());
}
