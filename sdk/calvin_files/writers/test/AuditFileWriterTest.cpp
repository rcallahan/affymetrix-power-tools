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
#include "calvin_files/writers/test/AuditFileWriterTest.h"
//
#include "calvin_files/parsers/src/AuditFileReader.h"
#include "calvin_files/writers/src/AuditFileWriter.h"
//
#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_array;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( AuditFileWriterTest );

void AuditFileWriterTest::setUp()
{
}

void AuditFileWriterTest::tearDown()
{
}

void AuditFileWriterTest::testCreation()
{
	AuditFileWriter writer;
	CPPUNIT_ASSERT(1);
}

void AuditFileWriterTest::testmethod_WriteFile()
{
	AuditFileWriter writer;
	ArrayAuditEntry entry;
	ParameterNameValuePair param;
	string file = "./audit.log";

	remove(file.c_str());

	entry.UserName() = L"ljevon";
	entry.DateTime().Date(L"10-10-2005");
	entry.DateTime().Time(L"12:00:00Z");
	entry.ActionType() = "test.id";
	entry.ArrayGuid() = "123-123-123-123";
	entry.InputFileGuids().push_back("321-321-321-321");
	entry.InputFileGuids().push_back("345-345-345-345");
	entry.OutputFileGuids().push_back("999-999-999-999");
	param.Name = L"Test";
	param.Value = L"true";
	entry.ActionParameters().push_back(param);
	param.Name = L"PixelSize";
	param.Value = L"3";
	entry.ActionParameters().push_back(param);
	CPPUNIT_ASSERT(writer.Write(file, entry) == true);

	entry.Clear();
	entry.UserName() = L"ljevon1";
	entry.DateTime().Date( L"10-11-2005");
	entry.DateTime().Time( L"12:01:00Z");
	entry.ActionType() = "test2.id";
	entry.ArrayGuid() = "1123-123-123-123";
	entry.InputFileGuids().push_back("1321-321-321-321");
	entry.InputFileGuids().push_back("1345-345-345-345");
	entry.OutputFileGuids().push_back("1999-999-999-999");
	param.Name = L"Test";
	param.Value = L"false";
	entry.ActionParameters().push_back(param);
	param.Name = L"PixelSize";
	param.Value = L"1";
	entry.ActionParameters().push_back(param);
	param.Name = L"Temp";
	param.Value = L"200";
	entry.ActionParameters().push_back(param);
	CPPUNIT_ASSERT(writer.Write(file, entry) == true);





	AuditFileReader reader;
	ArrayAuditEntryList audit;
	ArrayAuditEntryList::iterator it;
	AffymetrixGuidTypeList::iterator guidIt;
	ParameterNameValuePairList::iterator paramIt;

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
