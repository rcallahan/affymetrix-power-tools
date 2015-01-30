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
#include "calvin_files/data/test/FileHeaderTest.h"
//
#include "calvin_files/data/src/FileHeader.h"
#include "calvin_files/utils/src/AffyStlCollectionTypes.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( FileHeaderTest );

void FileHeaderTest::setUp()
{
	header = new FileHeader();
}

void FileHeaderTest::tearDown()
{
	delete header;
}

void FileHeaderTest::testCreation()
{
	FileHeader hdr;
	CPPUNIT_ASSERT(1);
}

void FileHeaderTest::FilenameTest()
{
	std::string p1 = "yo mama";
	header->SetFilename(p1);
	std::string p2 = header->GetFilename();
	CPPUNIT_ASSERT(p1 == p2);
}

void FileHeaderTest::MagicNumberTest()
{
	u_int8_t p1 = header->GetMagicNumber();
	CPPUNIT_ASSERT(p1 == MAGIC_NUM);
}

void FileHeaderTest::VersionTest()
{
	u_int8_t p1 = header->GetVersion();
	CPPUNIT_ASSERT(p1 == CALVINIOVERSION);
}

void FileHeaderTest::GenericDataHdrTest()
{
	GenericDataHeader p1;
	std::string s1 = "claudia schiffer";
	p1.SetFileId(s1);
	header->SetGenericDataHdr(p1);
	GenericDataHeader* p2 = header->GetGenericDataHdr();
	std::string s2 = p2->GetFileId();
	CPPUNIT_ASSERT(s1 == s2);
}

void FileHeaderTest::FindDataGroupHeaderByNameTest()
{
	// Create DataGroupHeaders
	DataGroupHeader d;
	d.SetName(L"Default");
	header->AddDataGroupHdr(d);

	const DataGroupHeader* dch = header->FindDataGroupHeader(L"none");
	CPPUNIT_ASSERT(0 == dch);
	dch = header->FindDataGroupHeader(L"default");
	CPPUNIT_ASSERT(0 == dch);
	dch = header->FindDataGroupHeader(L"Default");
	CPPUNIT_ASSERT(0 != dch);
	CPPUNIT_ASSERT(dch->GetName() == L"Default");
}

