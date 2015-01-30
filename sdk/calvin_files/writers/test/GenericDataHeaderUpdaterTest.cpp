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

#include "calvin_files/writers/test/GenericDataHeaderUpdaterTest.h"
//
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/writers/src/GenericFileWriter.h"
//
#include "util/Fs.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( GenericDataHeaderUpdaterTest );

void GenericDataHeaderUpdaterTest::setUp()
{
}

void GenericDataHeaderUpdaterTest::tearDown()
{
}

void GenericDataHeaderUpdaterTest::CreationTest()
{
	GenericDataHeaderUpdater updater;
	CPPUNIT_ASSERT(1);
}

void GenericDataHeaderUpdaterTest::UpdateFileIdTest()
{
	std::string filename = "generic_update_fileid_test";
	CreateTestFile(filename, "");

	// Read in the file
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename(filename);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data));

	// Get the file Id
	AffymetrixGuidType fileId = data.Header().GetGenericDataHdr()->GetFileId();

	CPPUNIT_ASSERT(fileId.length() > 0);

	// Update the 
	std::ofstream os;
        Fs::aptOpen(os, filename, std::ios::out|std::ios::binary|std::ios::in);

	AffymetrixGuidType newFileId = AffymetrixGuid::GenerateNewGuid();
	GenericDataHeader updateHdr;
	updateHdr.SetFileId(newFileId);

	GenericDataHeaderUpdater updater;
	CPPUNIT_ASSERT(updater.Update(os, updateHdr, *data.Header().GetGenericDataHdr()));

	os.close();

	// Now read the file again.
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data));

	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->GetFileId() == newFileId);
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->GetFileTypeId() == "funky file");
}

void GenericDataHeaderUpdaterTest::CreateTestFile(const std::string& name, const std::string& fileId)
{
	FileHeader fh;
	fh.SetFilename(name);
	fh.GetGenericDataHdr()->SetFileTypeId("funky file");

	if (fileId.length() != 0)
		fh.GetGenericDataHdr()->SetFileId(fileId);

	GenericFileWriter writer(&fh);
	writer.WriteHeader();
}

void GenericDataHeaderUpdaterTest::UpdateFileIdFailTest()
{
	std::string filename = "generic_update_fileid_test";
	CreateTestFile(filename, "fake file id");

	// Read in the file
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename(filename);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data));

	// Get the file Id
	AffymetrixGuidType fileId = data.Header().GetGenericDataHdr()->GetFileId();

	CPPUNIT_ASSERT(fileId == "fake file id");

	// Update the 
	std::ofstream os;
        Fs::aptOpen(os, filename, std::ios::out|std::ios::binary|std::ios::in);

	AffymetrixGuidType newFileId = AffymetrixGuid::GenerateNewGuid();
	GenericDataHeader updateHdr;
	updateHdr.SetFileId(newFileId);

	// Should not update unless the current and new file ids are the same length.
	GenericDataHeaderUpdater updater;
	CPPUNIT_ASSERT(updater.Update(os, updateHdr, *data.Header().GetGenericDataHdr()) == false);

	os.close();
}

void GenericDataHeaderUpdaterTest::CreateTestFileWithParameters(const std::string& name)
{
	FileHeader fh;
	fh.SetFilename(name);
	fh.GetGenericDataHdr()->SetFileTypeId("funky file");
	fh.GetGenericDataHdr()->SetFileCreationTime(DateTime::GetCurrentDateTime().ToString());

	ParameterNameValueType nvt;
	nvt.SetName(L"Float");
	nvt.SetValueFloat(4.5f);
	fh.GetGenericDataHdr()->AddNameValParam(nvt);
	nvt.SetName(L"String");
	nvt.SetValueText(L"not updatable");
	fh.GetGenericDataHdr()->AddNameValParam(nvt);
	nvt.SetName(L"Integer");
	nvt.SetValueInt32(56);
	fh.GetGenericDataHdr()->AddNameValParam(nvt);

	GenericFileWriter writer(&fh);
	writer.WriteHeader();
}

void GenericDataHeaderUpdaterTest::UpdateMatchingParameterTest()
{
	std::string filename = "generic_update_parameter_test1";
	CreateTestFileWithParameters(filename);

	// Read in the file
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename(filename);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data));

	ParameterNameValueType nvt1;
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"String", nvt1));
	CPPUNIT_ASSERT(nvt1.GetValueText() == L"not updatable");
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"Float", nvt1));
	CPPUNIT_ASSERT(nvt1.GetValueFloat() == 4.5f);
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"Integer", nvt1));
	CPPUNIT_ASSERT(nvt1.GetValueInt32() == 56);

	// Update the parameter list
	std::ofstream os;
        Fs::aptOpen(os, filename, std::ios::out|std::ios::binary|std::ios::in);

	AffymetrixGuidType newFileId = AffymetrixGuid::GenerateNewGuid();
	GenericDataHeader updateHdr;
	updateHdr.SetFileId(newFileId);
	ParameterNameValueType nvt2;
	nvt2.SetName(L"Float");
	nvt2.SetValueFloat(5.7f);
	updateHdr.AddNameValParam(nvt2);
	nvt2.SetName(L"String");
	nvt2.SetValueText(L"updatable");
	updateHdr.AddNameValParam(nvt2);
	nvt2.SetName(L"Integer");
	nvt2.SetValueInt32(112);
	updateHdr.AddNameValParam(nvt2);

	GenericDataHeaderUpdater updater;
	CPPUNIT_ASSERT(updater.Update(os, updateHdr, *data.Header().GetGenericDataHdr()));

	os.close();

	// Now read the file again.
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data));

	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->GetFileId() == newFileId);

	ParameterNameValueType nvt3;
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"String", nvt3));
	CPPUNIT_ASSERT(nvt3.GetValueText() == L"updatable");	// the string is shorter so it can be updated.
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"Float", nvt3));
	CPPUNIT_ASSERT(nvt3.GetValueFloat() == 5.7f);
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"Integer", nvt3));
	CPPUNIT_ASSERT(nvt3.GetValueInt32() == 112);
}

void GenericDataHeaderUpdaterTest::UpdateNonMatchingParamterTest()
{
	std::string filename = "generic_update_parameter_test2";
	CreateTestFileWithParameters(filename);

	// Read in the file
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename(filename);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data));

	ParameterNameValueType nvt1;
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"String", nvt1));
	CPPUNIT_ASSERT(nvt1.GetValueText() == L"not updatable");
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"Float", nvt1));
	CPPUNIT_ASSERT(nvt1.GetValueFloat() == 4.5f);
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"Integer", nvt1));
	CPPUNIT_ASSERT(nvt1.GetValueInt32() == 56);

	// Update the parameter list
	std::ofstream os;
        Fs::aptOpen(os, filename, std::ios::out|std::ios::binary|std::ios::in);

	AffymetrixGuidType newFileId = AffymetrixGuid::GenerateNewGuid();
	GenericDataHeader updateHdr;
	updateHdr.SetFileId(newFileId);
	ParameterNameValueType nvt2;
	nvt2.SetName(L"Float");
	nvt2.SetValueInt32(57);
	updateHdr.AddNameValParam(nvt2);
	nvt2.SetName(L"String");
	nvt2.SetValueText(L"really not updatable");
	updateHdr.AddNameValParam(nvt2);
	nvt2.SetName(L"Integer");
	nvt2.SetValueInt16(112);
	updateHdr.AddNameValParam(nvt2);

	GenericDataHeaderUpdater updater;
	CPPUNIT_ASSERT(updater.Update(os, updateHdr, *data.Header().GetGenericDataHdr()));

	os.close();

	// Now read the file again.
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data));

	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->GetFileId() == newFileId);

	ParameterNameValueType nvt3;
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"String", nvt3));
	CPPUNIT_ASSERT(nvt3.GetValueText() == L"not updatable");
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"Float", nvt3));
	CPPUNIT_ASSERT(nvt3.GetValueFloat() == 4.5f);
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"Integer", nvt3));
	CPPUNIT_ASSERT(nvt3.GetValueInt32() == 56);
}
