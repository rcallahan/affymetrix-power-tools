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
#include "calvin_files/parsers/test/DATFileReaderTest.h"
//
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/parsers/src/DATFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_parameter;

CPPUNIT_TEST_SUITE_REGISTRATION( DATFileReaderTest );

const std::string SMALL_DAT_FILE = "../data/small_DAT_file";
const std::string SMALL_DAT_FILE_WITH_SUBGRIDS = "../data/small_DAT_file_with_subgrids";
const std::string SMALL_DAT_FILE_WITH_SUBGRIDS_AND_PARAMETERS = "../data/small_DAT_file_with_subgrids_and_parameters";

void DATFileReaderTest::setUp()
{
}

void DATFileReaderTest::tearDown()
{
}

void DATFileReaderTest::CreationTest()
{
	DATFileReader reader;
	CPPUNIT_ASSERT(1);
}

void DATFileReaderTest::ReadDATWithGenericReaderTest()
{
	GenericFileReader reader;
	GenericData data;

	reader.SetFilename("../data/DAT_file");
	reader.ReadHeader(data);

	int dataGroups = data.DataGroupCnt();
	CPPUNIT_ASSERT(dataGroups == 1);

	int dataSets = data.DataSetCnt(0);
	CPPUNIT_ASSERT(dataSets == 2);

	affymetrix_calvin_io::DataSet* pixels = data.DataSet(0, 0);
	if(pixels->Open())
	{
		CPPUNIT_ASSERT(1);
		int r = pixels->Rows();
		CPPUNIT_ASSERT(r == 10);
		int c = pixels->Cols();
		CPPUNIT_ASSERT(c == 1);
		int16_t pix1;
		pixels->GetData(0, 0, pix1);
		CPPUNIT_ASSERT(pix1 == 36);
		int16_t pix2;
		pixels->GetData(1, 0, pix2);
		CPPUNIT_ASSERT(pix2 == 3);
	}
	else
	{
		CPPUNIT_ASSERT(0);
	}
	pixels->Delete();

	affymetrix_calvin_io::DataSet* stats = data.DataSet(0, 1);
	if(stats->Open())
	{
		CPPUNIT_ASSERT(1);
		int r = stats->Rows();
		CPPUNIT_ASSERT(r == 1);
		int c = stats->Cols();
		CPPUNIT_ASSERT(c == 2);
		int16_t s1;
		stats->GetData(0, 0, s1);
		CPPUNIT_ASSERT(s1 == 16);
		int16_t s2;
		stats->GetData(0, 1, s2);
		CPPUNIT_ASSERT(s2 == 22);
	}
	else
	{
		CPPUNIT_ASSERT(0);
	}
	stats->Delete();
}

void DATFileReaderTest::ReadSmallDATFileTest()
{
	DATData data;
	DATFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_DAT_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == SMALL_DAT_FILE);
	CPPUNIT_ASSERT(data.GetRows() == 4);
	CPPUNIT_ASSERT(data.GetCols() == 5);
	CPPUNIT_ASSERT(data.GetArrayType() == L"Hg-small");
	CPPUNIT_ASSERT(data.GetPixelSize() == 0.71f);
	CPPUNIT_ASSERT(data.GetScannerType() == L"M10");
	CPPUNIT_ASSERT(data.GetScannerID() == L"main");
	DateTime dt = data.GetScanDate();
	CPPUNIT_ASSERT(dt.ToString() == L"12-25-05T11:12:13Z");

	u_int32_t pixelCnt = data.GetRows()*data.GetCols();
	u_int16_t* expectedPixels = new u_int16_t[pixelCnt];
	for (u_int16_t i=0; i<pixelCnt; ++i)
		expectedPixels[i] = i+10;

	u_int16_t* pixels = new u_int16_t[pixelCnt];
	CPPUNIT_ASSERT(data.GetPixels(pixels, 0, data.GetRows()));
	CPPUNIT_ASSERT(memcmp(expectedPixels, pixels, sizeof(u_int16_t)*data.GetRows()*data.GetCols()) == 0);

	u_int16_t min, max;
	CPPUNIT_ASSERT(data.GetRange(min, max));
	CPPUNIT_ASSERT(min == pixels[0]);
	CPPUNIT_ASSERT(max == pixels[pixelCnt-1]);

	delete [] expectedPixels;
	delete [] pixels;
}

void DATFileReaderTest::ReadSmallDATFileWithSubgridsTest()
{
	DATData data;
	DATFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_DAT_FILE_WITH_SUBGRIDS));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == SMALL_DAT_FILE_WITH_SUBGRIDS);
	CPPUNIT_ASSERT(data.GetRows() == 4);
	CPPUNIT_ASSERT(data.GetCols() == 5);
	CPPUNIT_ASSERT(data.GetArrayType() == L"Hg-small");
	CPPUNIT_ASSERT(data.GetPixelSize() == 0.71f);
	CPPUNIT_ASSERT(data.GetScannerType() == L"M10");
	CPPUNIT_ASSERT(data.GetScannerID() == L"main");
	DateTime dt = data.GetScanDate();
	CPPUNIT_ASSERT(dt.ToString() == L"12-25-05T11:12:13Z");

	u_int32_t pixelCnt = data.GetRows()*data.GetCols();
	u_int16_t* expectedPixels = new u_int16_t[pixelCnt];
	for (u_int16_t i=0; i<pixelCnt; ++i)
		expectedPixels[i] = i+10;

	u_int16_t* pixels = new u_int16_t[pixelCnt];
	CPPUNIT_ASSERT(data.GetPixels(pixels, 0, data.GetRows()));
	CPPUNIT_ASSERT(memcmp(expectedPixels, pixels, sizeof(u_int16_t)*data.GetRows()*data.GetCols()) == 0);

	u_int16_t min, max;
	CPPUNIT_ASSERT(data.GetRange(min, max));
	CPPUNIT_ASSERT(min == pixels[0]);
	CPPUNIT_ASSERT(max == pixels[pixelCnt-1]);

	delete [] expectedPixels;
	delete [] pixels;

	// Compute the expected grid
	float increment = 2.5f;
	FRegion ggRegion;
	float ptVal = 0.0;
	for(int i = 0; i < 4; i++)
	{
		FPoint point;
		point.x = ptVal;
		ptVal += increment;
		point.y = ptVal;
		ggRegion.pts.push_back(point);
	}
	u_int32_t ggStatus = DATData::GridOK | DATData::GridManualAdjust;

	// Check the grid
	CPPUNIT_ASSERT(data.GetGlobalGrid() == ggRegion);
	// Check the status
	CPPUNIT_ASSERT(data.GetGlobalGridStatus() == ggStatus);

	int32_t subgridCnt = 10;

	// Compute the expected subgrids
	for(int n = 0; n < subgridCnt; n++)
	{
		FRegion expectedSGRegion;
		float ptVal = 0.0f + (float)n;
		for(int i = 0; i < 4; i++)
		{
			FPoint point;
			point.x = ptVal;
			ptVal += increment;
			point.y = ptVal;
			expectedSGRegion.pts.push_back(point);
		}

		u_int32_t sgStatus = DATData::GridError | DATData::GridManualAdjust;
		if (n % 2 == 0)
			sgStatus = DATData::GridOK | DATData::GridManualAdjust;

		FRegion sgRegion = data.GetSubgrid(n);

		// Check the subgrid
		CPPUNIT_ASSERT(sgRegion == expectedSGRegion);
		// Check the status
		CPPUNIT_ASSERT(sgStatus == data.GetSubgridStatus(n));
	}

}

void DATFileReaderTest::ReadSmallDATFileWithSubgridsUsingSmallMemoryFootprintTest()
{
	DATData data;
	DATFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_DAT_FILE_WITH_SUBGRIDS));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == SMALL_DAT_FILE_WITH_SUBGRIDS);
	CPPUNIT_ASSERT(data.GetRows() == 4);
	CPPUNIT_ASSERT(data.GetCols() == 5);
	CPPUNIT_ASSERT(data.GetArrayType() == L"Hg-small");
	CPPUNIT_ASSERT(data.GetPixelSize() == 0.71f);
	CPPUNIT_ASSERT(data.GetScannerType() == L"M10");
	CPPUNIT_ASSERT(data.GetScannerID() == L"main");
	DateTime dt = data.GetScanDate();
	CPPUNIT_ASSERT(dt.ToString() == L"12-25-05T11:12:13Z");

	// Compute the expected results for the full image
	u_int32_t fullImagePixelCnt = data.GetRows()*data.GetCols();
	u_int16_t* expectedPixelsForFullImage = new u_int16_t[fullImagePixelCnt];
	for (u_int16_t i=0; i<fullImagePixelCnt; ++i)
		expectedPixelsForFullImage[i] = i+10;

	// Copy a 2 by 3 region of the expected results.
	u_int16_t i=0;
	u_int32_t pixelCnt = 2 * 3; // rows * cols
	u_int16_t* expectedPixels = new u_int16_t[pixelCnt];
	for (u_int16_t irow = 1; irow < 3; ++irow)
		for (u_int16_t icol = 2; icol < 5; ++icol, ++i)
			expectedPixels[i] = expectedPixelsForFullImage[irow*data.GetCols() + icol];

	// Retrieve a 2 by 3 region of the 4 by 5 image
	u_int16_t* pixels = new u_int16_t[pixelCnt];
	data.GetPixels(pixels, 1, 2, 2, 3);
	//CPPUNIT_ASSERT(data.GetPixels(pixels, 1, 2, 2, 3));
	CPPUNIT_ASSERT(memcmp(expectedPixels, pixels, sizeof(u_int16_t)*pixelCnt) == 0);

	u_int16_t min, max;
	CPPUNIT_ASSERT(data.GetRange(min, max));
	CPPUNIT_ASSERT(min == expectedPixelsForFullImage[0]);
	CPPUNIT_ASSERT(max == expectedPixelsForFullImage[fullImagePixelCnt-1]);

	delete [] expectedPixelsForFullImage;
	delete [] expectedPixels;
	delete [] pixels;

	// Compute the expected grid
	float increment = 2.5f;
	FRegion ggRegion;
	float ptVal = 0.0;
	for(int i = 0; i < 4; i++)
	{
		FPoint point;
		point.x = ptVal;
		ptVal += increment;
		point.y = ptVal;
		ggRegion.pts.push_back(point);
	}
	u_int32_t ggStatus = DATData::GridOK | DATData::GridManualAdjust;

	// Check the grid
	CPPUNIT_ASSERT(data.GetGlobalGrid() == ggRegion);
	// Check the status
	CPPUNIT_ASSERT(data.GetGlobalGridStatus() == ggStatus);

	int32_t subgridCnt = 10;

	// Compute the expected subgrids
	for(int n = 0; n < subgridCnt; n++)
	{
		FRegion expectedSGRegion;
		float ptVal = 0.0f + (float)n;
		for(int i = 0; i < 4; i++)
		{
			FPoint point;
			point.x = ptVal;
			ptVal += increment;
			point.y = ptVal;
			expectedSGRegion.pts.push_back(point);
		}

		u_int32_t sgStatus = DATData::GridError | DATData::GridManualAdjust;
		if (n % 2 == 0)
			sgStatus = DATData::GridOK | DATData::GridManualAdjust;

		FRegion sgRegion = data.GetSubgrid(n);

		// Check the subgrid
		CPPUNIT_ASSERT(sgRegion == expectedSGRegion);
		// Check the status
		CPPUNIT_ASSERT(sgStatus == data.GetSubgridStatus(n));
	}
}

void DATFileReaderTest::ReadSmallDATFileWithGridParametersTest()
{
	DATData data;
	DATFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_DAT_FILE_WITH_SUBGRIDS_AND_PARAMETERS));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == SMALL_DAT_FILE_WITH_SUBGRIDS_AND_PARAMETERS);
	CPPUNIT_ASSERT(data.GetRows() == 4);
	CPPUNIT_ASSERT(data.GetCols() == 5);

	ParameterNameValueType nvt;
	ParameterNameValueTypeVector params;
	data.GetGridAlignmentAlgorithmParameters(params);
	CPPUNIT_ASSERT(params.size() == 2);
	CPPUNIT_ASSERT(data.FindGridAlignmentAlgorithmParameter(L"County", nvt));
	CPPUNIT_ASSERT(nvt.GetValueText() == L"Santa Clara");
	CPPUNIT_ASSERT(params[0].GetName() == L"Zip Code");
	CPPUNIT_ASSERT(params[0].GetValueInt32() == 95051);
}
