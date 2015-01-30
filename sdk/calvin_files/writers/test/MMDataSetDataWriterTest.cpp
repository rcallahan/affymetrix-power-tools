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
#include <windows.h>	// MemMapFile is windows only
//
#include "calvin_files/writers/test/MMDataSetDataWriterTest.h"
//
#include "calvin_files/data/src/MemMapFile.h"
#include "calvin_files/writers/src/MMDataSetDataWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( MMDataSetDataWriterTest );

void MMDataSetDataWriterTest::setUp()
{
	// setup the DataSetHeader
	hdr = new DataSetHeader;
	hdr->AddUShortColumn(L"");
	hdr->AddFloatColumn(L"");
	hdr->SetRowCnt(50);
	hdr->SetDataStartFilePos(45);

	// Create the file, MMDataSetDataWriter assumes that the file is already created.
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	mmfile.Create(hdr->GetDataSize() + 45 + 100 /*extra*/);
	mmfile.Close();

	writer = new MMDataSetDataWriter(*hdr, "fred_flintstone");
}

void MMDataSetDataWriterTest::tearDown()
{
	delete writer;
	delete hdr;
	::DeleteFile("fred_flintstone");
}

void MMDataSetDataWriterTest::CreationTest()
{
	CPPUNIT_ASSERT(writer);
}


void MMDataSetDataWriterTest::OpenTest()
{
	CPPUNIT_ASSERT(writer->Open());
	CPPUNIT_ASSERT(writer->Close());
}

void MMDataSetDataWriterTest::CloseTest()
{
	CPPUNIT_ASSERT(writer->Open());
	CPPUNIT_ASSERT(writer->Close());
}

void MMDataSetDataWriterTest::MapDataTest()
{
	CPPUNIT_ASSERT(writer->Open());
	CPPUNIT_ASSERT(writer->MapData(40, 10));
	CPPUNIT_ASSERT(writer->Close());
}

void MMDataSetDataWriterTest::GetMappedDataPtrTest()
{
	CPPUNIT_ASSERT(writer->Open());
	CPPUNIT_ASSERT(writer->GetMappedDataPtr());
	CPPUNIT_ASSERT(writer->Close());
}

void MMDataSetDataWriterTest::GetFirstRowMappedTest()
{
	CPPUNIT_ASSERT(writer->Open());
	CPPUNIT_ASSERT(writer->GetFirstRowMapped() == 0);
	CPPUNIT_ASSERT(writer->MapData(40, 10));
	CPPUNIT_ASSERT(writer->GetFirstRowMapped() == 40);
	CPPUNIT_ASSERT(writer->Close());
}

void MMDataSetDataWriterTest::GetRowsMappedTest()
{
	CPPUNIT_ASSERT(writer->Open());
	CPPUNIT_ASSERT(writer->GetRowsMapped() == 50);
	CPPUNIT_ASSERT(writer->MapData(40, 10));
	CPPUNIT_ASSERT(writer->GetRowsMapped() == 10);
	CPPUNIT_ASSERT(writer->Close());
}

void MMDataSetDataWriterTest::GetMaxRowsToMapTest()
{
	u_int32_t bytes = 200*1024*1024;
	int32_t maxRows = bytes/hdr->GetRowSize();
	CPPUNIT_ASSERT(writer->GetMaxRowsToMap() == maxRows);
}
