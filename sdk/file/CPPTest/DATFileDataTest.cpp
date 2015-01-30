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

#include "file/CPPTest/DATFileDataTest.h"
//
#include "file/DATFileData.h"
//
#include <cmath>
//

using namespace affxdat;
using namespace std;

#define TEST_DAT_FILE_NAME "./data/test.dat"
#define MINI_DAT_FILE_NAME "./data/mini.dat"
#define NOT_A_DAT_FILE_NAME "./data/not_a_dat_file.dat"
#define GRID_DAT_FILE_NAME "./data/grid.dat"

CPPUNIT_TEST_SUITE_REGISTRATION( CDATFileDataTest );

void CDATFileDataTest::setUp()
{
}

void CDATFileDataTest::tearDown()
{
}

void CDATFileDataTest::testproperty_Header_Cols ()
{
	int32_t val = 10;
	CDATFileHeaderData header;
	header.SetCols(val);
	CPPUNIT_ASSERT(header.GetCols() == val);
}

void CDATFileDataTest::testproperty_Header_Rows ()
{
	int32_t val = 10;
	CDATFileHeaderData header;
	header.SetRows(val);
	CPPUNIT_ASSERT(header.GetRows() == val);
}

void CDATFileDataTest::testproperty_Header_MinValue ()
{
	int16_t val = 10;
	CDATFileHeaderData header;
	header.SetMinValue(val);
	CPPUNIT_ASSERT(header.GetMinValue() == val);
}

void CDATFileDataTest::testproperty_Header_MaxValue ()
{
	int16_t val = 10;
	CDATFileHeaderData header;
	header.SetMaxValue(val);
	CPPUNIT_ASSERT(header.GetMaxValue() == val);
}

void CDATFileDataTest::testproperty_Header_Grid ()
{
	GridCoordinatesType grid1;
	GridCoordinatesType grid2;

	grid1.upperleft.x  = 1;
	grid1.upperleft.y  = 2;
	grid1.lowerleft.x  = 3;
	grid1.lowerleft.y  = 4;
	grid1.upperright.x = 5;
	grid1.upperright.y = 7;
	grid1.lowerright.x = 8;
	grid1.lowerright.y = 9;

	CDATFileHeaderData header;
	header.SetGridCorners(grid1);
	grid2 = header.GetGridCorners();

	CPPUNIT_ASSERT( grid1.upperleft.x == grid2.upperleft.x );
	CPPUNIT_ASSERT( grid1.upperleft.y == grid2.upperleft.y );

	CPPUNIT_ASSERT( grid1.lowerleft.x == grid2.lowerleft.x );
	CPPUNIT_ASSERT( grid1.lowerleft.y == grid2.lowerleft.y );

	CPPUNIT_ASSERT( grid1.upperright.x == grid2.upperright.x );
	CPPUNIT_ASSERT( grid1.upperright.y == grid2.upperright.y );

	CPPUNIT_ASSERT( grid1.lowerright.x == grid2.lowerright.x );
	CPPUNIT_ASSERT( grid1.lowerright.y == grid2.lowerright.y );
}

void CDATFileDataTest::testproperty_Header_ArrayType()
{
	string arrayType = "array_type";
	CDATFileHeaderData header;
	header.SetArrayType(arrayType.c_str());
	CPPUNIT_ASSERT( arrayType == header.GetArrayType() );
}


void CDATFileDataTest::testproperty_FileName ()
{
	CDATFileData dat;
	string path = "test";
	dat.SetFileName(path.c_str());
	CPPUNIT_ASSERT( dat.GetFileName() == path );
}

void CDATFileDataTest::testmethod_Exists()
{
	CDATFileData dat;
	string path = TEST_DAT_FILE_NAME;
	dat.SetFileName(path.c_str());
	CPPUNIT_ASSERT( dat.Exists() == true );
}

void CDATFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	CDATFileData dat;
	string path = "test";
	dat.SetFileName(path.c_str());
	CPPUNIT_ASSERT( dat.Exists() == false );
}

void CDATFileDataTest::testmethod_Read ()
{
	CDATFileData dat;
	dat.SetFileName(TEST_DAT_FILE_NAME);
	CPPUNIT_ASSERT( dat.Read() == true );

	CDATFileHeaderData &header = dat.GetHeader();
	CPPUNIT_ASSERT(header.GetCols() == 1280);
	CPPUNIT_ASSERT(header.GetRows() == 1024);
	CPPUNIT_ASSERT(header.GetMinValue() == 208);
	CPPUNIT_ASSERT(header.GetMaxValue() == 4122);
	CPPUNIT_ASSERT(header.GetNumPixels() == header.GetCols()*header.GetRows());
	CPPUNIT_ASSERT( header.GetArrayType() == "u133A_8u");
	CPPUNIT_ASSERT( header.GetScannerID() == "");

	GridCoordinatesType grid = header.GetGridCorners();

	CPPUNIT_ASSERT( grid.upperleft.x == 220 );
	CPPUNIT_ASSERT( grid.upperleft.y == 92 );

	CPPUNIT_ASSERT( grid.lowerleft.x == 223 );
	CPPUNIT_ASSERT( grid.lowerleft.y == 967 );

	CPPUNIT_ASSERT( grid.upperright.x == 1094 );
	CPPUNIT_ASSERT( grid.upperright.y == 89 );

	CPPUNIT_ASSERT( grid.lowerright.x == 1098 );
	CPPUNIT_ASSERT( grid.lowerright.y == 963 );
}

void CDATFileDataTest::testmethod_Read_with_scanner_ID()
{
	CDATFileData dat;
	dat.SetFileName(MINI_DAT_FILE_NAME);
	CPPUNIT_ASSERT( dat.Read() == true );

	CDATFileHeaderData &header = dat.GetHeader();
	CPPUNIT_ASSERT( header.GetScannerID() == "7G_Feas_3");
}

void CDATFileDataTest::testmethod_ReadWhenFileNotExists ()
{
	CDATFileData dat;
	string path = "test";
	dat.SetFileName(path.c_str());
	CPPUNIT_ASSERT( dat.Read() == false );
	dat.Close();
}

void CDATFileDataTest::testmethod_GetPixel ()
{
	CDATFileData dat;
	dat.SetFileName(TEST_DAT_FILE_NAME);
	CPPUNIT_ASSERT( dat.Read() == true );

	CPPUNIT_ASSERT(dat.GetPixel(0, 0) == 1651);
	CPPUNIT_ASSERT(dat.GetPixel(0, 1) == 1651);
	CPPUNIT_ASSERT(dat.GetPixel(1, 0) == 481);
	CPPUNIT_ASSERT(dat.GetPixel(1, 3) == 530);
	CPPUNIT_ASSERT(dat.GetPixel(9, 5) == 991);
	CPPUNIT_ASSERT(dat.GetPixel(1279, 1023) == 1629);
}

void CDATFileDataTest::testmethod_GetPixels_many_rows()
{
	CDATFileData dat;
	dat.SetFileName(TEST_DAT_FILE_NAME);
	CPPUNIT_ASSERT( dat.Read() == true );
	uint16_t *pixels = new uint16_t[dat.GetHeader().GetCols()];
	dat.GetPixels(0, 1, pixels);
	CPPUNIT_ASSERT( pixels[0] == 1651 );
	CPPUNIT_ASSERT( pixels[1] == 481 );
	CPPUNIT_ASSERT( pixels[1279] == 810 );
	dat.GetPixels(1023, pixels);
	CPPUNIT_ASSERT( pixels[0] == 1651 );
	CPPUNIT_ASSERT( pixels[1] == 609 );
	CPPUNIT_ASSERT( pixels[1279] == 1629 );
	delete[] pixels;

	pixels = new uint16_t[dat.GetHeader().GetCols() * dat.GetHeader().GetRows()];
	dat.GetPixels(0, dat.GetHeader().GetRows(), pixels);
	CPPUNIT_ASSERT( pixels[0] == 1651 );
	CPPUNIT_ASSERT( pixels[1] == 481 );
	CPPUNIT_ASSERT( pixels[1279] == 810 );
	int lastRow = dat.GetHeader().GetCols()*(dat.GetHeader().GetRows()-1);
	CPPUNIT_ASSERT( pixels[lastRow] == 1651 );
	CPPUNIT_ASSERT( pixels[lastRow+1] == 609 );
	CPPUNIT_ASSERT( pixels[lastRow+1279] == 1629 );
	delete[] pixels;
}

void CDATFileDataTest::testmethod_GetPixels ()
{
	CDATFileData dat;
	dat.SetFileName(TEST_DAT_FILE_NAME);
	CPPUNIT_ASSERT( dat.Read() == true );
	uint16_t *pixels = new uint16_t[dat.GetHeader().GetCols()];
	dat.GetPixels(0, pixels);
	CPPUNIT_ASSERT( pixels[0] == 1651 );
	CPPUNIT_ASSERT( pixels[1] == 481 );
	CPPUNIT_ASSERT( pixels[1279] == 810 );
	dat.GetPixels(1023, pixels);
	CPPUNIT_ASSERT( pixels[0] == 1651 );
	CPPUNIT_ASSERT( pixels[1] == 609 );
	CPPUNIT_ASSERT( pixels[1279] == 1629 );
	delete[] pixels;
}

void CDATFileDataTest::testproperty_Header_ScannerID()
{
	string id = "123123123";
	CDATFileHeaderData header;
	header.SetScannerID(id.c_str());
	CPPUNIT_ASSERT( header.GetScannerID() == id); 
}


void CDATFileDataTest::testmethod_IsGCOSDATFile()
{
	CPPUNIT_ASSERT(CDATFileData::IsGCOSDATFile(TEST_DAT_FILE_NAME) == true);
	CPPUNIT_ASSERT(CDATFileData::IsGCOSDATFile(NOT_A_DAT_FILE_NAME) == false);
}

void CDATFileDataTest::testmethod_UpdateGridCorners()
{
	GridCoordinatesType grid1;

	grid1.upperleft.x  = 1;
	grid1.upperleft.y  = 2;
	grid1.lowerleft.x  = 3;
	grid1.lowerleft.y  = 4;
	grid1.upperright.x = 5;
	grid1.upperright.y = 7;
	grid1.lowerright.x = 8;
	grid1.lowerright.y = 9;

	CPPUNIT_ASSERT(CDATFileData::UpdateGridCorners(GRID_DAT_FILE_NAME, grid1) == true);

	CDATFileData dat;
	dat.SetFileName(GRID_DAT_FILE_NAME);
	CPPUNIT_ASSERT( dat.Read() == true );
	CPPUNIT_ASSERT( dat.GetHeader().GetScannerID() == "7G_Feas_3");
	CPPUNIT_ASSERT( dat.GetHeader().GetRows() == 26702 );
	CPPUNIT_ASSERT( dat.GetHeader().GetCols() == 26702 );
	CPPUNIT_ASSERT( dat.GetHeader().GetMinValue() == 0 );
	CPPUNIT_ASSERT( dat.GetHeader().GetMaxValue() == 65534 );

	GridCoordinatesType grid2 = dat.GetHeader().GetGridCorners();

	CPPUNIT_ASSERT( grid1.upperleft.x == grid2.upperleft.x );
	CPPUNIT_ASSERT( grid1.upperleft.y == grid2.upperleft.y );

	CPPUNIT_ASSERT( grid1.lowerleft.x == grid2.lowerleft.x );
	CPPUNIT_ASSERT( grid1.lowerleft.y == grid2.lowerleft.y );

	CPPUNIT_ASSERT( grid1.upperright.x == grid2.upperright.x );
	CPPUNIT_ASSERT( grid1.upperright.y == grid2.upperright.y );

	CPPUNIT_ASSERT( grid1.lowerright.x == grid2.lowerright.x );
	CPPUNIT_ASSERT( grid1.lowerright.y == grid2.lowerright.y );

	dat.Close();

	grid1.upperleft.x  = 11;
	grid1.upperleft.y  = 12;
	grid1.lowerleft.x  = 13;
	grid1.lowerleft.y  = 14;
	grid1.upperright.x = 15;
	grid1.upperright.y = 17;
	grid1.lowerright.x = 18;
	grid1.lowerright.y = 19;

	CPPUNIT_ASSERT(CDATFileData::UpdateGridCorners(GRID_DAT_FILE_NAME, grid1) == true);

	dat.SetFileName(GRID_DAT_FILE_NAME);
	CPPUNIT_ASSERT( dat.Read() == true );
	CPPUNIT_ASSERT( dat.GetHeader().GetScannerID() == "7G_Feas_3");
	CPPUNIT_ASSERT( dat.GetHeader().GetRows() == 26702 );
	CPPUNIT_ASSERT( dat.GetHeader().GetCols() == 26702 );
	CPPUNIT_ASSERT( dat.GetHeader().GetMinValue() == 0 );
	CPPUNIT_ASSERT( dat.GetHeader().GetMaxValue() == 65534 );

	grid2 = dat.GetHeader().GetGridCorners();

	CPPUNIT_ASSERT( grid1.upperleft.x == grid2.upperleft.x );
	CPPUNIT_ASSERT( grid1.upperleft.y == grid2.upperleft.y );

	CPPUNIT_ASSERT( grid1.lowerleft.x == grid2.lowerleft.x );
	CPPUNIT_ASSERT( grid1.lowerleft.y == grid2.lowerleft.y );

	CPPUNIT_ASSERT( grid1.upperright.x == grid2.upperright.x );
	CPPUNIT_ASSERT( grid1.upperright.y == grid2.upperright.y );

	CPPUNIT_ASSERT( grid1.lowerright.x == grid2.lowerright.x );
	CPPUNIT_ASSERT( grid1.lowerright.y == grid2.lowerright.y );

	dat.Close();

	CPPUNIT_ASSERT(CDATFileData::UpdateGridCorners(NOT_A_DAT_FILE_NAME, grid1) == false);
}
