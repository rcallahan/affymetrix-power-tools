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
#ifdef _MSC_VER
#include <windows.h>	// MemMapFile is windows only
#endif
//
#include "calvin_files/data/test/MemMapFileTest.h"
//
#include "calvin_files/data/src/MemMapFile.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( MemMapFileTest );

void MemMapFileTest::setUp()
{
}

void MemMapFileTest::tearDown()
{
	::DeleteFile("fred_flintstone");
}

void MemMapFileTest::CreationTest()
{
	MemMapFile mmfile;
	CPPUNIT_ASSERT(1);
}

void MemMapFileTest::FilenameTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.GetFilename() == "fred_flintstone");
}

void MemMapFileTest::CreateTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.Create(200));
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}

void MemMapFileTest::OpenTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.Create(200));
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE));
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}

void MemMapFileTest::MapDataTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.Create(200));
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE));
	CPPUNIT_ASSERT(mmfile.MapData(0, 200));
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}

void MemMapFileTest::CloseTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.Create(200));
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE));
	CPPUNIT_ASSERT(mmfile.MapData(0, 200));
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}

void MemMapFileTest::GetDataPtrTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.Create(200));
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE));
	CPPUNIT_ASSERT(mmfile.MapData(0, 200));
	char* ptr = mmfile.GetDataPtr();
	CPPUNIT_ASSERT(ptr != 0);
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}

void MemMapFileTest::GetFirstMappedBytePosTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.Create(200));
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE));
	CPPUNIT_ASSERT(mmfile.MapData(50, 150));
	CPPUNIT_ASSERT(mmfile.GetFirstMappedBytePos() == 50);
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}

void MemMapFileTest::GetBytesMappedTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.Create(200));
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE));
	CPPUNIT_ASSERT(mmfile.MapData(100, 200));
	CPPUNIT_ASSERT(mmfile.GetBytesMapped() == 100);
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}

void MemMapFileTest::GetLastErrorTest()
{
	MemMapFile mmfile;
	std::wstring errorMsg = mmfile.GetLastError();
	CPPUNIT_ASSERT(errorMsg.length() == 0);	
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE) == false);
	errorMsg = mmfile.GetLastError();
	CPPUNIT_ASSERT(errorMsg.length() > 0);
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}

void MemMapFileTest::IsViewMappedTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.IsViewMapped() == false);
	CPPUNIT_ASSERT(mmfile.Create(200));
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE));
	CPPUNIT_ASSERT(mmfile.IsViewMapped() == false);
	CPPUNIT_ASSERT(mmfile.MapData(0, 200));
	CPPUNIT_ASSERT(mmfile.IsViewMapped());
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}

void MemMapFileTest::CreateNoSetFilenameTest()
{
	MemMapFile mmfile;
	CPPUNIT_ASSERT(mmfile.Create(200) == false);
}

void MemMapFileTest::OpenNoCreateTest()
{
	MemMapFile mmfile;
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE) == false);
}

void MemMapFileTest::WriteValuesTest()
{
	MemMapFile mmfile;
	mmfile.SetFilename("fred_flintstone");
	CPPUNIT_ASSERT(mmfile.IsViewMapped() == false);
	CPPUNIT_ASSERT(mmfile.Create(200));
	CPPUNIT_ASSERT(mmfile.Open(READ_AND_WRITE, ALLOWREAD_AND_WRITE));
	CPPUNIT_ASSERT(mmfile.MapData(0, 100));
	char* ptr = mmfile.GetDataPtr();
	CPPUNIT_ASSERT(ptr);
	*ptr = 'c';
	ptr++;
	*ptr = 's';
	// Re-map data
	CPPUNIT_ASSERT(mmfile.MapData(100, 100));
	CPPUNIT_ASSERT(mmfile.MapData(0, 200));
	// Check the value
	ptr = mmfile.GetDataPtr();
	CPPUNIT_ASSERT(*ptr == 'c');
	++ptr;
	CPPUNIT_ASSERT(*ptr == 's');
	CPPUNIT_ASSERT_NO_THROW(mmfile.Close());
}
