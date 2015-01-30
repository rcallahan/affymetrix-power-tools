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

#include "calvin_files/writers/test/GenericFileUpdaterTest.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

#define TEST_FILE_NAME "../data/small_cel_file"

CPPUNIT_TEST_SUITE_REGISTRATION( GenericFileUpdaterTest );

void GenericFileUpdaterTest::setUp() {}

void GenericFileUpdaterTest::tearDown() {}

void GenericFileUpdaterTest::testCreation()
{

}

void GenericFileUpdaterTest::UpdaterTest()
{
	string tempFile = TEST_FILE_NAME;
	tempFile.append(".TMP");
	if(Copy(TEST_FILE_NAME, tempFile))
	{
		GenericFileUpdater* updater = new GenericFileUpdater(tempFile);
		DataSetHeader dataSetHdr;
		dataSetHdr.SetName(L"Test");
		dataSetHdr.SetRowCnt(1);
		dataSetHdr.AddAsciiColumn(L"Test", 32);
		DataSetWriter* writer = updater->GetUpdateDataSetWriter(dataSetHdr);
		writer->WriteHeader();
		writer->Write("Hello World", 32);

		writer->UpdateNextDataSetOffset();
		delete writer;
		delete updater;

		GenericData gData;
		GenericFileReader reader;
		reader.SetFilename(tempFile);
		reader.ReadHeader(gData);
		gData.Close();
		CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSetCnt() == 6);
		std::wstring name = gData.Header().GetDataGroup(0).GetDataSet(5).GetName();
		CPPUNIT_ASSERT(name.compare(L"Test") == 0);
	}
	remove(tempFile.c_str());
}

void GenericFileUpdaterTest::Updater2Test()
{
	// Update by adding a DataSetHeader that has zero rows and cols and has items in the parameter list.
	string tempFile = TEST_FILE_NAME;
	tempFile.append("2.TMP");
	if(Copy(TEST_FILE_NAME, tempFile))
	{
		DataSetHeader dataSetHdr;
		ParameterNameValueType nvt;
		nvt.SetName(L"First");
		nvt.SetValueInt32(11);
		dataSetHdr.AddNameValParam(nvt);
		nvt.SetName(L"Deuxieme");
		nvt.SetValueFloat(5.75f);
		dataSetHdr.AddNameValParam(nvt);
		dataSetHdr.SetName(L"Test2");
		dataSetHdr.SetRowCnt(0);
		GenericFileUpdater* updater = new GenericFileUpdater(tempFile);
		DataSetWriter* writer = updater->GetUpdateDataSetWriter(dataSetHdr);
		writer->WriteHeader();

		writer->UpdateNextDataSetOffset();
		delete writer;
		delete updater;

		GenericData gData;
		GenericFileReader reader;
		reader.SetFilename(tempFile);
		reader.ReadHeader(gData);
		gData.Close();
		CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSetCnt() == 6);
		std::wstring name = gData.Header().GetDataGroup(0).GetDataSet(5).GetName();
		CPPUNIT_ASSERT(name.compare(L"Test2") == 0);
		CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSet(5).GetRowCnt() == 0);
		CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSet(5).GetColumnCnt() == 0);
		CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSet(5).GetNameValParamCnt() == 2);
		
		ParameterNameValueTypeConstIt begin, end;
		gData.Header().GetDataGroup(0).GetDataSet(5).GetNameValIterators(begin, end);
		ParameterNameValueTypeConstIt ii = begin;
		CPPUNIT_ASSERT(ii->GetName() == L"First");
		CPPUNIT_ASSERT(ii->GetValueInt32() == 11);
		++ii;
		CPPUNIT_ASSERT(ii->GetName() == L"Deuxieme");
		CPPUNIT_ASSERT(ii->GetValueFloat() == 5.75f);
	}
	// dont' delete the file because it is used in the next test.
}

void GenericFileUpdaterTest::UpdateExistingDataSet1Test()
{
	// Update an existing data set at the end of the file
	string tempFile = TEST_FILE_NAME;
	tempFile.append("2.TMP");

	DataSetHeader dataSetHdr;
	ParameterNameValueType nvt;
	nvt.SetName(L"Happy Day");
	nvt.SetValueAscii("fun in the sun");
	dataSetHdr.AddNameValParam(nvt);
	dataSetHdr.SetName(L"Test2");
	dataSetHdr.SetRowCnt(0);
	GenericFileUpdater* updater = new GenericFileUpdater(tempFile);
	DataSetWriter* writer = updater->GetUpdateDataSetWriter(dataSetHdr);
	writer->WriteHeader();

	writer->UpdateNextDataSetOffset();
	delete writer;
	delete updater;

	GenericData gData;
	GenericFileReader reader;
	reader.SetFilename(tempFile);
	reader.ReadHeader(gData);
	gData.Close();
	int32_t dataSetCnt = gData.Header().GetDataGroup(0).GetDataSetCnt();
	CPPUNIT_ASSERT(dataSetCnt == 6);
	std::wstring name = gData.Header().GetDataGroup(0).GetDataSet(5).GetName();
	CPPUNIT_ASSERT(name.compare(L"Test2") == 0);
	CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSet(5).GetRowCnt() == 0);
	CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSet(5).GetColumnCnt() == 0);
	CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSet(5).GetNameValParamCnt() == 1);
	
	ParameterNameValueTypeConstIt begin, end;
	gData.Header().GetDataGroup(0).GetDataSet(5).GetNameValIterators(begin, end);
	ParameterNameValueTypeConstIt ii = begin;
	CPPUNIT_ASSERT(ii->GetName() == L"Happy Day");
	CPPUNIT_ASSERT(ii->GetValueAscii() == "fun in the sun");

	// uncomment this line when the updater is working
	remove(tempFile.c_str());
}

void GenericFileUpdaterTest::UpdateExistingDataSet2Test()
{
	// Update an existing data set that is not at the end of the file
	string tempFile = TEST_FILE_NAME;
	tempFile.append("4.TMP");
	if(Copy(TEST_FILE_NAME, tempFile))
	{
		DataSetHeader dataSetHdr;
		ParameterNameValueType nvt;
		nvt.SetName(L"First");
		nvt.SetValueInt32(11);
		dataSetHdr.AddNameValParam(nvt);
		nvt.SetName(L"Deuxieme");
		nvt.SetValueFloat(5.75f);
		dataSetHdr.AddNameValParam(nvt);
		dataSetHdr.SetName(L"Intensity");
		dataSetHdr.SetRowCnt(0);
		GenericFileUpdater* updater = new GenericFileUpdater(tempFile);
		DataSetWriter* writer = updater->GetUpdateDataSetWriter(dataSetHdr);
		writer->WriteHeader();

		writer->UpdateNextDataSetOffset();
		delete writer;
		delete updater;

		GenericData gData;
		GenericFileReader reader;
		reader.SetFilename(tempFile);
		reader.ReadHeader(gData);
		gData.Close();
		int32_t dataSetCnt = gData.Header().GetDataGroup(0).GetDataSetCnt();
		CPPUNIT_ASSERT(dataSetCnt == 6);
		std::wstring name = gData.Header().GetDataGroup(0).GetDataSet(5).GetName();
		CPPUNIT_ASSERT(name.compare(L"Intensity") == 0);
		CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSet(5).GetRowCnt() == 0);
		CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSet(5).GetColumnCnt() == 0);
		CPPUNIT_ASSERT(gData.Header().GetDataGroup(0).GetDataSet(5).GetNameValParamCnt() == 2);
		
		ParameterNameValueTypeConstIt begin, end;
		gData.Header().GetDataGroup(0).GetDataSet(5).GetNameValIterators(begin, end);
		ParameterNameValueTypeConstIt ii = begin;
		CPPUNIT_ASSERT(ii->GetName() == L"First");
		CPPUNIT_ASSERT(ii->GetValueInt32() == 11);
		++ii;
		CPPUNIT_ASSERT(ii->GetName() == L"Deuxieme");
		CPPUNIT_ASSERT(ii->GetValueFloat() == 5.75f);
	}
	// uncomment this line when the updater is working
	remove(tempFile.c_str());
}

bool GenericFileUpdaterTest::Copy(std::string filename, std::string copyFilename)
{
	ifstream in(filename.c_str(), ios::in | ios::binary);
	ofstream out(copyFilename.c_str(), ios::out | ios::binary);
	if(in.is_open() && out.is_open())
	{
		char* buf = new char[1];
		while(!in.eof())
		{
			in.read(buf, 1);
			out.write(buf, 1);
		}
		delete[] buf;
		in.close();
		out.close();
		return true;
	}
	return false;
}
