////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

#include "file/CPPTest/BARFileDataTest.h"
//
#include "file/BARFileData.h"
//
#include <cmath>
//

#ifdef _MSC_VER
#define TEST_BAR_FILE_NAME ".\\data\\test.bar"
#else
#define TEST_BAR_FILE_NAME "./data/test.bar"
#endif

CPPUNIT_TEST_SUITE_REGISTRATION( CBARFileDataTest );

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

void CBARFileDataTest::setUp()
{
}

void CBARFileDataTest::tearDown()
{
}

void CBARFileDataTest::testCreation()
{
	affxbar::CBARFileData bar;
	CPPUNIT_ASSERT( 1 );
}

void CBARFileDataTest::testproperty_FileName()
{
	affxbar::CBARFileData bar;
	std::string path = "test";
	bar.SetFileName(path.c_str());
}

void CBARFileDataTest::testmethod_Exists()
{
	affxbar::CBARFileData bar;
	std::string path = TEST_BAR_FILE_NAME;
	bar.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bar.Exists() == true );
}

void CBARFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	affxbar::CBARFileData bar;
	std::string path = "test";
	bar.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bar.Exists() == false );
}

void CBARFileDataTest::testmethod_Read()
{
	affxbar::CBARFileData bar;
	affxbar::CGDACSequenceResultItem seq;
	affxbar::BarSequenceResultData data;
	std::string path = TEST_BAR_FILE_NAME;
	bar.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", bar.Exists() == true );
	CPPUNIT_ASSERT( bar.Read() );
	CPPUNIT_ASSERT( CompareFloats(bar.GetVersion(), 2.0f) );
	CPPUNIT_ASSERT( bar.GetNumberSequences() == 1);
	CPPUNIT_ASSERT( bar.GetNumberColumns() == 2);
	CPPUNIT_ASSERT( bar.GetNumberParameters() == 1);
	CPPUNIT_ASSERT( bar.GetParameter(0).Tag == std::string("Param1") );
	CPPUNIT_ASSERT( bar.GetParameter(0).Value == std::string("NoVal") );
	CPPUNIT_ASSERT( bar.GetColumnTypes(0) == affxbar::BAR_DATA_INTEGER );
	CPPUNIT_ASSERT( bar.GetColumnTypes(1) == affxbar::BAR_DATA_FLOAT );

	bar.GetResults(0, seq);
	CPPUNIT_ASSERT(seq.GetColumnType(0) == 	affxbar::BAR_DATA_INTEGER);
	CPPUNIT_ASSERT(seq.GetColumnType(1) == affxbar::BAR_DATA_FLOAT);
	CPPUNIT_ASSERT(seq.GetNumberColumns() == 2);
	CPPUNIT_ASSERT(seq.GetName() == std::string("Seq1") );
	CPPUNIT_ASSERT(seq.GetVersion() == std::string("1.0") );
	CPPUNIT_ASSERT(seq.GetGroupName() == std::string("Group1") );
	CPPUNIT_ASSERT(seq.GetNumberDataPoints() == 5);
	CPPUNIT_ASSERT(seq.GetNumberParameters() == 1);
	CPPUNIT_ASSERT(seq.GetParameter(0).Tag == std::string("SeqParam1") );
	CPPUNIT_ASSERT(seq.GetParameter(0).Value == std::string("NoVal") );
	seq.GetData(0, 0, data);
	CPPUNIT_ASSERT(data.iValue == 0);
	seq.GetData(0, 1, data);
	CPPUNIT_ASSERT(data.fValue == 82.0f);
	seq.GetData(1, 0, data);
	CPPUNIT_ASSERT(data.iValue == 1);
	seq.GetData(1, 1, data);
	CPPUNIT_ASSERT(data.fValue == 38.2f);
	seq.GetData(2, 0, data);
	CPPUNIT_ASSERT(data.iValue == 2);
	seq.GetData(2, 1, data);
	CPPUNIT_ASSERT(data.fValue == 9.4f);
	seq.GetData(3, 0, data);
	CPPUNIT_ASSERT(data.iValue == 3);
	seq.GetData(3, 1, data);
	CPPUNIT_ASSERT(data.fValue == 5.6f);
	seq.GetData(4, 0, data);
	CPPUNIT_ASSERT(data.iValue == 4);
	seq.GetData(4, 1, data);
	CPPUNIT_ASSERT(data.fValue == 57.8f);
	bar.Close();
}

void CBARFileDataTest::testmethod_ReadWhenFileNotExists()
{
	affxbar::CBARFileData bar;
	std::string path = "test";
	bar.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bar.Read() == false );
	bar.Close();
}

void CBARFileDataTest::testmethod_ReadHeader()
{
	affxbar::CBARFileData bar;
	
	std::string path = TEST_BAR_FILE_NAME;
	bar.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", bar.Exists() == true );
	CPPUNIT_ASSERT( bar.ReadHeader() );
	CPPUNIT_ASSERT( CompareFloats(bar.GetVersion(), 2.0f) );
	CPPUNIT_ASSERT( bar.GetNumberSequences() == 1);
	CPPUNIT_ASSERT( bar.GetNumberColumns() == 2);
	CPPUNIT_ASSERT( bar.GetNumberParameters() == 1);
	CPPUNIT_ASSERT( bar.GetParameter(0).Tag == std::string("Param1") );
	CPPUNIT_ASSERT( bar.GetParameter(0).Value == std::string("NoVal") );
	CPPUNIT_ASSERT( bar.GetColumnTypes(0) == affxbar::BAR_DATA_INTEGER );
	CPPUNIT_ASSERT( bar.GetColumnTypes(1) == affxbar::BAR_DATA_FLOAT );
	bar.Close();
}

void CBARFileDataTest::testmethod_ReadHeaderWhenFileNotExists()
{
	affxbar::CBARFileData bar;
	std::string path = "test";
	bar.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bar.ReadHeader() == false );
	bar.Close();
}

void CBARFileDataTest::testmethod_GetErrorFromReadError()
{
	affxbar::CBARFileData bar;
	std::string path = "test";
	bar.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bar.Read() == false );
	CPPUNIT_ASSERT( bar.GetError() == std::string("Unable to open the file.") );
	bar.Close();
}

void CBARFileDataTest::testmethod_GetErrorFromReadHeaderError()
{
	affxbar::CBARFileData bar;
	std::string path = "test";
	bar.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bar.ReadHeader() == false );
	CPPUNIT_ASSERT( bar.GetError() == std::string("Unable to open the file.") );
	bar.Close();
}

void CBARFileDataTest::testmethod_GetFullName()
{
	affxbar::CGDACSequenceResultItem seq;
	seq.SetName("name");
	CPPUNIT_ASSERT(seq.GetFullName() == "name");
	seq.SetVersion("version");
	seq.SetGroupName("group");
	CPPUNIT_ASSERT(seq.GetFullName() == "group:version;name");
}
