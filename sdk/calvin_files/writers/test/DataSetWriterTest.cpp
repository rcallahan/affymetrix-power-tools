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

#include "calvin_files/writers/test/DataSetWriterTest.h"
//
#include "calvin_files/data/src/ColumnInfo.h"
#include "calvin_files/writers/src/DataSetWriter.h"
//
#include "util/Fs.h"
//
using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetWriterTest );

void DataSetWriterTest::setUp()
{
	std::string f = "data_dataGroup";
        Fs::aptOpen(os, f, std::ios::out|std::ios::binary|std::ios::trunc);
	if (!os.is_open() && !os.good())
	{
	CPPUNIT_ASSERT(0);
	}
	else 
	{
		hdr.SetName(L"xyz123");
		hdr.AddAsciiColumn(L"xyz", 64);
		hdr.AddIntColumn(L"123");
		hdr.SetRowCnt(3);
		writer = new DataSetWriter(&os, &hdr);
	}
}

void DataSetWriterTest::tearDown()
{
	delete writer;
	os.close();
	hdr.Clear();
}

void DataSetWriterTest::testCreation()
{
std::ofstream os;
	DataSetHeader hdr;
	DataSetWriter w(&os, &hdr);
	CPPUNIT_ASSERT(1);
}

void DataSetWriterTest::WriteTest()
{
writer->WriteHeader();
	CPPUNIT_ASSERT(1);

	//write out three rows of data to the dataGroup
	writer->Write("testing123", 64);
	writer->Write(6);
	writer->Write("testing456", 64);
	writer->Write(7);
	writer->Write("testing789", 64);
	writer->Write(8);
	CPPUNIT_ASSERT(1);
}
