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

#include "calvin_files/data/test/GenericDataHeaderTest.h"
//
#include "calvin_files/data/src/GenericDataHeader.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( GenericDataHeaderTest );

void GenericDataHeaderTest::setUp()
{
	header = new GenericDataHeader();
}

void GenericDataHeaderTest::tearDown()
{
	delete header;
}

void GenericDataHeaderTest::testCreation()
{
	GenericDataHeader hdr;
	CPPUNIT_ASSERT(1);
}

void GenericDataHeaderTest::FileTypeIdTest()
{
	std::string p1 = "file type ID";
	header->SetFileTypeId(p1);
	std::string p2 = header->GetFileTypeId();
	CPPUNIT_ASSERT(p1 == p2);
}

void GenericDataHeaderTest::FileIdTest()
{
	std::string p1 = "unique file ID";
	header->SetFileId(p1);
	std::string p2 = header->GetFileId();
	CPPUNIT_ASSERT(p1 == p2);
}

void GenericDataHeaderTest::FileCreationTimeTest()
{
	std::wstring p1 = L"20050301T12:42:11Z";
	header->SetFileCreationTime(p1);
	std::wstring p2 = header->GetFileCreationTime();
	CPPUNIT_ASSERT(p1 == p2);
}

void GenericDataHeaderTest::LocaleTest()
{
	std::wstring p1 = L"locale";
	header->SetLocale(p1);
	std::wstring p2 = header->GetLocale();
	CPPUNIT_ASSERT(p1 == p2);
}

void GenericDataHeaderTest::NameValTest()
{
	std::wstring p[] = {L"Uma",L"Thurman",L"Angelina",L"Jolie",L"Isabella ",L"Rossellini"};
	ParameterNameValueType t1;
	t1.SetName(p[0]);
	t1.SetValueText(p[1]);
	header->AddNameValParam(t1);

	ParameterNameValueType t2;
	t2.SetName(p[2]);
	t2.SetValueText(p[3]);
	header->AddNameValParam(t2);

	ParameterNameValueType t3;
	t3.SetName(p[4]);
	t3.SetValueText(p[5]);
	header->AddNameValParam(t3);
	ParameterNameValueTypeIt begin;
	ParameterNameValueTypeIt end;
	header->GetNameValIterators(begin, end);
	if(begin != end)
	{
		CPPUNIT_ASSERT(begin->GetName() == p[0]);
		CPPUNIT_ASSERT(begin->GetValueText() == p[1]);
		begin++;
	}
	if(begin != end)
	{
		CPPUNIT_ASSERT(begin->GetName() == p[2]);
		CPPUNIT_ASSERT(begin->GetValueText() == p[3]);
		begin++;
	}
	if(begin != end)
	{
		CPPUNIT_ASSERT(begin->GetName() == p[4]);
		CPPUNIT_ASSERT(begin->GetValueText() == p[5]);
	}
}

void GenericDataHeaderTest::AddParentEntryTest()
{
	GenericDataHeader header1;
	std::string s1 = "Uma Thurman";
	header1.SetFileId(s1);
	GenericDataHeader header2;
	std::string s2 = "Angelina Jolie";
	header2.SetFileId(s2);
	GenericDataHeader header3;
	std::string s3 = "Isabella Rossellini";
	header3.SetFileId(s3);
	header->AddParent(header1);
	header->AddParent(header2);
	header->AddParent(header3);
	GenDataHdrVectorIt begin;
	GenDataHdrVectorIt end;
	header->GetParentIterators(begin, end);
	if(begin != end)
	{
		CPPUNIT_ASSERT(begin->GetFileId() == s1);
		begin++;
	}
	if(begin != end)
	{
		CPPUNIT_ASSERT(begin->GetFileId() == s2);
		begin++;
	}
	if(begin != end)
	{
		CPPUNIT_ASSERT(begin->GetFileId() == s3);
	}
}

void GenericDataHeaderTest::GetNumParentsTest()
{
	GenericDataHeader header1;
	GenericDataHeader header2;
	GenericDataHeader header3;
	header->AddParent(header1);
	header->AddParent(header2);
	header->AddParent(header3);
	int sz = header->GetParentCnt();
	CPPUNIT_ASSERT(sz == 3);
}

void GenericDataHeaderTest::GetNameValPairCntTest()
{
	std::wstring p[] = {L"Uma",L"Thurman",L"Angelina",L"Jolie",L"Isabella ",L"Rossellini"};
	ParameterNameValueType t1;
	t1.SetName(p[0]);
	t1.SetValueText(p[1]);
	header->AddNameValParam(t1);

	ParameterNameValueType t2;
	t2.SetName(p[2]);
	t2.SetValueText(p[3]);
	header->AddNameValParam(t2);

	ParameterNameValueType t3;
	t3.SetName(p[4]);
	t3.SetValueText(p[5]);
	header->AddNameValParam(t3);
	CPPUNIT_ASSERT(header->GetNameValParamCnt() == 3);
}

void GenericDataHeaderTest::GetNameValPairTest()
{
	std::wstring p[] = {L"Uma",L"Thurman",L"Angelina",L"Jolie",L"Isabella ",L"Rossellini"};
	ParameterNameValueType t1;
	t1.SetName(p[0]);
	t1.SetValueText(p[1]);
	header->AddNameValParam(t1);

	ParameterNameValueType t2;
	t2.SetName(p[2]);
	t2.SetValueText(p[3]);
	header->AddNameValParam(t2);

	ParameterNameValueType t3;
	t3.SetName(p[4]);
	t3.SetValueText(p[5]);
	header->AddNameValParam(t3);
	CPPUNIT_ASSERT(header->GetNameValParam(0).GetName() == p[0]);
	CPPUNIT_ASSERT(header->GetNameValParam(0).GetValueText() == p[1]);
	CPPUNIT_ASSERT(header->GetNameValParam(1).GetName() == p[2]);
	CPPUNIT_ASSERT(header->GetNameValParam(1).GetValueText() == p[3]);
	CPPUNIT_ASSERT(header->GetNameValParam(2).GetName() == p[4]);
	CPPUNIT_ASSERT(header->GetNameValParam(2).GetValueText() == p[5]);
}

void GenericDataHeaderTest::UpdateNameValPairTest()
{
	std::wstring p[] = {L"Uma",L"Thurman",L"Angelina",L"Jolie",L"Isabella ",L"Rossellini"};
	ParameterNameValueType t1;
	t1.SetName(p[0]);
	t1.SetValueText(p[1]);
	header->AddNameValParam(t1);

	ParameterNameValueType t2;
	t2.SetName(p[2]);
	t2.SetValueText(p[3]);
	header->AddNameValParam(t2);

	ParameterNameValueType t3;
	t3.SetName(p[4]);
	t3.SetValueText(p[5]);
	header->AddNameValParam(t3);
	CPPUNIT_ASSERT(header->GetNameValParam(0).GetValueText() == p[1]);
}

void GenericDataHeaderTest::FindNameValPairTest()
{
	std::wstring p[] = {L"Uma",L"Thurman",L"Angelina",L"Jolie",L"Isabella ",L"Rossellini"};
	ParameterNameValueType t1;
	t1.SetName(p[0]);
	t1.SetValueText(p[1]);
	header->AddNameValParam(t1);

	ParameterNameValueType t2;
	t2.SetName(p[2]);
	t2.SetValueText(p[3]);
	header->AddNameValParam(t2);

	ParameterNameValueType t3;
	t3.SetName(p[4]);
	t3.SetValueText(p[5]);
	header->AddNameValParam(t3);

	ParameterNameValueType t4;
	t4.SetName(L"");
	t4.SetValueText(L"");
	CPPUNIT_ASSERT(header->FindNameValParam(p[0], t4));
	CPPUNIT_ASSERT(t4.GetName() == p[0]);
	CPPUNIT_ASSERT(t4.GetValueText() == p[1]);
	CPPUNIT_ASSERT(header->FindNameValParam(L"Oingo Boingo", t4) == false);
}

void GenericDataHeaderTest::GetNameValParamsBeginsWithTest()
{
	std::wstring names[] = {L"prefix-one", L"two", L"prefix-three"};
	ParameterNameValueType t;

	for (int32_t i = 0; i < 3; ++i)
	{
		t.SetName(names[i]);
		t.SetValueInt32(i);
		header->AddNameValParam(t);
	}
	
	ParameterNameValueTypeVector v;
	header->GetNameValParamsBeginsWith(L"prefix-", v);
	CPPUNIT_ASSERT(v.size() == 2);
	CPPUNIT_ASSERT(v.at(0).GetName() == names[0]);
	CPPUNIT_ASSERT(v.at(1).GetName() == names[2]);
}

void GenericDataHeaderTest::FindParentByFileTypeIdTest()
{
	GenericDataHeader header1;
	header1.SetFileTypeId("Fred");
	GenericDataHeader header2;
	header2.SetFileTypeId("Wilma");
	GenericDataHeader header3;
	header3.SetFileTypeId("Betty");
	header->AddParent(header1);
	header->AddParent(header2);
	header->AddParent(header3);

	GenericDataHeader* gdh = header->FindParent("Wilma");
	CPPUNIT_ASSERT(gdh != 0);
	CPPUNIT_ASSERT(gdh->GetFileTypeId() == "Wilma");
	gdh = header->FindParent("nonexistant");
	CPPUNIT_ASSERT(gdh == 0);
}
