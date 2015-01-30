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
#include "calvin_files/parsers/test/CelFileReaderTest.h"
//
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/parsers/src/CelFileReader.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_parameter;

const std::string SMALL_CEL_FILE = "../data/small_cel_file";
const std::string SMALL_CEL_FILE_NO_OUTLIER_NO_MASK = "../data/small_cel_file_no_outlier_no_mask";
const std::string SMALL_CEL_FILE_NO_STDEV = "../data/small_cel_file_no_stdev";
const std::string LARGE_CEL_FILE = "../data/large_cel_file";
const std::string SMALL_UINT16_CEL_FILE = "../data/small_uint16_cel_file";

static int32_t intenTestValues[4] = {55683, 4568, 2368, 100};
static float stdevTestValues[4] = {2.345f, 56.23f, 1.53f, 3.875f};

CPPUNIT_TEST_SUITE_REGISTRATION( CelFileReaderTest );

void CelFileReaderTest::setUp()
{
}

void CelFileReaderTest::tearDown()
{
}

void CelFileReaderTest::CreationTest()
{
	CelFileReader reader;
	CPPUNIT_ASSERT(1);
}

void CelFileReaderTest::ReadSmallCelFileNoOutlierNoMaskTest()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename(SMALL_CEL_FILE_NO_OUTLIER_NO_MASK);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	// Read some data
	CPPUNIT_ASSERT(data.GetFileHeader()->GetFilename() == SMALL_CEL_FILE_NO_OUTLIER_NO_MASK);
	CPPUNIT_ASSERT(data.GetVersion() == 1);
	CPPUNIT_ASSERT(data.GetArrayType() == L"Hg-small");
	CPPUNIT_ASSERT(data.GetAlgorithmName() == L"Feature Extraction");
	CPPUNIT_ASSERT(data.GetRows() == 2);
	CPPUNIT_ASSERT(data.GetCols() == 5);
	CPPUNIT_ASSERT(data.GetNumCells() == 10);
	CPPUNIT_ASSERT(data.HasStdev() == true);
	CPPUNIT_ASSERT(data.HasNumPixels() == true);

	ParameterNameValueTypeVector params;
	data.GetAlgorithmParameters(params);
	CPPUNIT_ASSERT(params.size() == 2);
	CPPUNIT_ASSERT(params.at(0).GetName() == L"percentile");
	CPPUNIT_ASSERT(params.at(0).GetValueFloat() == 0.75f);
	CPPUNIT_ASSERT(params.at(1).GetName() == L"outlierlow");
	CPPUNIT_ASSERT(params.at(1).GetValueFloat() == 1.004f);

	float intensity;
	int16_t numPixels;
	bool outlier, masked;
	float stdev;

	int32_t cells = data.GetNumCells();

	for (int32_t cell = 0; cell < cells; ++cell)
	{
		CPPUNIT_ASSERT_NO_THROW(data.GetData(cell, intensity, stdev, numPixels, outlier, masked));

		CPPUNIT_ASSERT(intensity == 100*cell);
		CPPUNIT_ASSERT(stdev == .5*cell);
		CPPUNIT_ASSERT(numPixels == 25);

		// Check the values
		CPPUNIT_ASSERT(outlier == false);
		CPPUNIT_ASSERT(masked == false);
	}

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;
	BoolVector vOutliers;
	BoolVector vMasked;

	// Get individual columns
	CPPUNIT_ASSERT(data.GetIntensities(0, cells, vInten));
	CPPUNIT_ASSERT(data.GetStdev(0, cells, vStdev));
	CPPUNIT_ASSERT(data.GetNumPixels(0, cells, vPixels));
	CPPUNIT_ASSERT(data.GetOutliers(0, cells, vOutliers) == false);
	CPPUNIT_ASSERT(data.GetMasked(0, cells, vMasked) == false);

	CPPUNIT_ASSERT(vInten.size() == 10);
	CPPUNIT_ASSERT(vStdev.size() == 10);
	CPPUNIT_ASSERT(vPixels.size() == 10);
	CPPUNIT_ASSERT(vOutliers.size() == 0);
	CPPUNIT_ASSERT(vMasked.size() == 0);
	
	// Check values
	for (int32_t cell = 0; cell < cells; ++cell)
	{
		CPPUNIT_ASSERT(vInten[cell] == 100*cell);
		CPPUNIT_ASSERT(vStdev[cell] == .5*cell);
		CPPUNIT_ASSERT(vPixels[cell] == 25);
	}

	XYCoordVector xyOutliers;
	XYCoordVector xyMasked;
	data.GetOutlierCoords(xyOutliers);
	data.GetMaskedCoords(xyMasked);
	CPPUNIT_ASSERT(xyOutliers.size() == 0);
	CPPUNIT_ASSERT(xyMasked.size() == 0);
}

void CelFileReaderTest::ReadSmallCelFileTest()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename(SMALL_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFileHeader()->GetFilename() == SMALL_CEL_FILE);
	CheckSmallCelFileHeader(data);

	float intensity;
	int16_t numPixels;
	bool outlier, masked;
	float stdev;

	int32_t cells = data.GetNumCells();

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;
	BoolVector vOutliers;
	BoolVector vMasked;

	// Get individual columns
	CPPUNIT_ASSERT(data.GetIntensities(0, cells, vInten));
	CPPUNIT_ASSERT(data.GetStdev(0, cells, vStdev));
	CPPUNIT_ASSERT(data.GetNumPixels(0, cells, vPixels));
	CPPUNIT_ASSERT(data.GetOutliers(0, cells, vOutliers));
	CPPUNIT_ASSERT(data.GetMasked(0, cells, vMasked));

	CPPUNIT_ASSERT(vInten.size() == 25);
	CPPUNIT_ASSERT(vStdev.size() == 25);
	CPPUNIT_ASSERT(vPixels.size() == 25);
	CPPUNIT_ASSERT(vOutliers.size() == 25);
	CPPUNIT_ASSERT(vMasked.size() == 25);

	
	// Check values
	for (int32_t cell = 0; cell < cells; ++cell)
	{
		CPPUNIT_ASSERT(vInten[cell] == 100*cell);
		CPPUNIT_ASSERT(vStdev[cell] == .5*cell);
		CPPUNIT_ASSERT(vPixels[cell] == 25);
		CheckOutlier(cell, vOutliers[cell]);
		CheckMasked(cell, vMasked[cell]);
	}

	XYCoordVector xyOutliers;
	XYCoordVector xyMasked;
	data.GetOutlierCoords(xyOutliers);
	data.GetMaskedCoords(xyMasked);
	CPPUNIT_ASSERT(xyOutliers.size() == 2);
	CPPUNIT_ASSERT(xyMasked.size() == 3);

	// Change the order of operations in the test
	for (int32_t cell = 0; cell < cells; ++cell)
	{
		CPPUNIT_ASSERT_NO_THROW(data.GetData(cell, intensity, stdev, numPixels, outlier, masked));

		CPPUNIT_ASSERT(intensity == 100*cell);
		CPPUNIT_ASSERT(stdev == .5*cell);
		CPPUNIT_ASSERT(numPixels == 25);

		CheckOutlier(cell, outlier);
		CheckMasked(cell, masked);
	}
}

void CelFileReaderTest::CheckSmallCelFileHeader(CelFileData& data)
{
	// Read some data
	CPPUNIT_ASSERT(data.GetVersion() == 1);
	CPPUNIT_ASSERT(data.GetArrayType() == L"Hg-small");
	CPPUNIT_ASSERT(data.GetAlgorithmName() == L"Feature Extraction");
	CPPUNIT_ASSERT(data.GetRows() == 5);
	CPPUNIT_ASSERT(data.GetCols() == 5);
	CPPUNIT_ASSERT(data.GetNumCells() == 25);
	CPPUNIT_ASSERT(data.HasStdev() == true);
	CPPUNIT_ASSERT(data.HasNumPixels() == true);

	ParameterNameValueTypeVector params;
	data.GetAlgorithmParameters(params);
	CPPUNIT_ASSERT(params.size() == 2);
	CPPUNIT_ASSERT(params.at(0).GetName() == L"percentile");
	CPPUNIT_ASSERT(params.at(0).GetValueFloat() == 0.75f);
	CPPUNIT_ASSERT(params.at(1).GetName() == L"outlierlow");
	CPPUNIT_ASSERT(params.at(1).GetValueFloat() == 1.004f);
}

void CelFileReaderTest::ReadSmallCelFileNoStdevTest()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename(SMALL_CEL_FILE_NO_STDEV);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	// Read some data
	CPPUNIT_ASSERT(data.GetFileHeader()->GetFilename() == SMALL_CEL_FILE_NO_STDEV);
	CPPUNIT_ASSERT(data.GetVersion() == 1);
	CPPUNIT_ASSERT(data.GetArrayType() == L"Hg-small");
	CPPUNIT_ASSERT(data.GetAlgorithmName() == L"Feature Extraction");
	CPPUNIT_ASSERT(data.GetRows() == 5);
	CPPUNIT_ASSERT(data.GetCols() == 5);
	CPPUNIT_ASSERT(data.GetNumCells() == 25);
	CPPUNIT_ASSERT(data.HasStdev() == false);
	CPPUNIT_ASSERT(data.HasNumPixels() == true);

	ParameterNameValueTypeVector params;
	data.GetAlgorithmParameters(params);
	CPPUNIT_ASSERT(params.size() == 2);
	CPPUNIT_ASSERT(params.at(0).GetName() == L"percentile");
	CPPUNIT_ASSERT(params.at(0).GetValueFloat() == 0.75f);
	CPPUNIT_ASSERT(params.at(1).GetName() == L"outlierlow");
	CPPUNIT_ASSERT(params.at(1).GetValueFloat() == 1.004f);

	float intensity;
	int16_t numPixels;
	bool outlier, masked;
	float stdev;

	int32_t cells = data.GetNumCells();

	for (int32_t cell = 0; cell < cells; ++cell)
	{
		CPPUNIT_ASSERT_NO_THROW(data.GetData(cell, intensity, stdev, numPixels, outlier, masked));

		CPPUNIT_ASSERT(intensity == 100*cell);
		CPPUNIT_ASSERT(stdev == 0);
		CPPUNIT_ASSERT(numPixels == 25);

		CheckOutlier(cell, outlier);
		CheckMasked(cell, masked);
	}

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;
	BoolVector vOutliers;
	BoolVector vMasked;

	// Get individual columns
	CPPUNIT_ASSERT(data.GetIntensities(0, cells, vInten));
	CPPUNIT_ASSERT(data.GetStdev(0, cells, vStdev) == false);
	CPPUNIT_ASSERT(data.GetNumPixels(0, cells, vPixels));
	CPPUNIT_ASSERT(data.GetOutliers(0, cells, vOutliers));
	CPPUNIT_ASSERT(data.GetMasked(0, cells, vMasked));
	
	// Check values
	for (int32_t cell = 0; cell < cells; ++cell)
	{
		CPPUNIT_ASSERT(vInten[cell] == 100*cell);
		CPPUNIT_ASSERT(vPixels[cell] == 25);

		CheckOutlier(cell, vOutliers[cell]);
		CheckMasked(cell, vMasked[cell]);
	}

	XYCoordVector xyOutliers;
	XYCoordVector xyMasked;
	data.GetOutlierCoords(xyOutliers);
	data.GetMaskedCoords(xyMasked);
	CPPUNIT_ASSERT(xyOutliers.size() == 2);
	CPPUNIT_ASSERT(xyMasked.size() == 3);
}

void CelFileReaderTest::ReadLargeCelFileCheckHeaderTest()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename(LARGE_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	const int32_t cols = 2560;
	const int32_t rows = 2560;
	const int32_t expectedCells = rows*cols;

	// Read some data
	CPPUNIT_ASSERT(data.GetFileHeader()->GetFilename() == LARGE_CEL_FILE);
	CPPUNIT_ASSERT(data.GetRows() == rows);
	CPPUNIT_ASSERT(data.GetCols() == cols);
	CPPUNIT_ASSERT(data.GetNumCells() == expectedCells);
	CPPUNIT_ASSERT(data.HasStdev() == true);
	CPPUNIT_ASSERT(data.HasNumPixels() == true);
}

void CelFileReaderTest::ReadLargeCelFileSingleDataAccessTest()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename(LARGE_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	const int32_t cols = 2560;
	const int32_t rows = 2560;
	//const int32_t expectedCells = rows*cols;

	// Declare buffers
	float* intensities = new float[cols];
	int32_t* expectedIntensities = new int32_t[cols];
	float* stdev = new float[cols];
	float* expectedStdev = new float[cols];
	int16_t* numPixels = new int16_t[cols];
	int16_t* expectedNumPixels = new int16_t[cols];
	bool outlier, masked;

	// This is slow. Spot check
	for (int32_t row = 0, cell = 0, testIdx = 0;
			 row < rows;
			 row+=100, cell+=(100*cols) )
	{
		for (int32_t col = 0; col < cols; ++col, ++cell, ++testIdx)
		{
			// Get the data
			data.GetData(cell, intensities[col], stdev[col], numPixels[col], outlier, masked);

			// Compute the expected values
			if (testIdx >= 4)
				testIdx = 0;

			expectedIntensities[col] = intenTestValues[testIdx];
			expectedStdev[col] = stdevTestValues[testIdx];
			expectedNumPixels[col] = 25;
		}
		CPPUNIT_ASSERT(memcmp(intensities, expectedIntensities, sizeof(int32_t)*cols) == 0);
		CPPUNIT_ASSERT(memcmp(stdev, expectedStdev, sizeof(float)*cols) == 0);
		CPPUNIT_ASSERT(memcmp(numPixels, expectedNumPixels, sizeof(int16_t)*cols) == 0);
	}


	// clean up
	delete [] intensities;
	delete [] expectedIntensities;
	delete [] stdev;
	delete [] expectedStdev;
	delete [] numPixels;
	delete [] expectedNumPixels;

}

void CelFileReaderTest::ReadLargeCelFileVectorDataAccessTest()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename(LARGE_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	const int32_t cols = 2560;
	const int32_t rows = 2560;
	//const int32_t expectedCells = rows*cols;

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;
	BoolVector vOutliers;
	BoolVector vMasked;

	// This needs to be faster. Spot check
	for (int32_t row = 0/*, testIdx = 0*/; row < rows; row+=100)
	{
		data.GetIntensities(row*cols, cols*100, vInten);
		data.GetStdev(row*cols, cols*100, vStdev);
		data.GetNumPixels(row*cols, cols*100, vPixels);

		for (int32_t col = 0, testIdx = 0; col < cols; ++col, ++testIdx)
		{
			// Compute the expected values
			if (testIdx >= 4)
				testIdx = 0;

			if (vInten[col] != intenTestValues[testIdx] ||
					vStdev[col] != stdevTestValues[testIdx] ||
					vPixels[col] != 25)
				CPPUNIT_ASSERT(0);	// this is slow, only execute it when needed.
		}
	}
	

/*
	XYCoordVector xyOutliers;
	XYCoordVector xyMasked;
	data.GetOutlierCoords(xyOutliers);
	data.GetMaskedCoords(xyMasked);
	CPPUNIT_ASSERT(xyOutliers.size() == 2);
	CPPUNIT_ASSERT(xyMasked.size() == 3);
*/
}

void CelFileReaderTest::CheckOutlier(int32_t cell, bool outlier)
{
	if (cell == 0 || cell == 11 /*(1,2)*/)
		CPPUNIT_ASSERT(outlier == true);
	else
		CPPUNIT_ASSERT(outlier == false);
}

void CelFileReaderTest::CheckMasked(int32_t cell, bool masked)
{
	if (cell == 1 /*(1,0)*/ || cell == 7 /*(2,1)*/ || cell == 13/*(3,2)*/)
		CPPUNIT_ASSERT(masked == true);
	else
		CPPUNIT_ASSERT(masked == false);
}

// Test added to reproduce bug reported by beta testers.
void CelFileReaderTest::ReadSmallCelFileGetOutlierCoordsTest()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename(SMALL_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	XYCoordVector xyOutliers;
	data.GetOutlierCoords(xyOutliers);
	CPPUNIT_ASSERT(xyOutliers.size() == 2);
}

// Test added to reproduce bug reported by beta testers.
void CelFileReaderTest::ReadSmallCelFileGetMaskedCoordsTest()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename(SMALL_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	XYCoordVector xyMasked;
	data.GetMaskedCoords(xyMasked);
	CPPUNIT_ASSERT(xyMasked.size() == 3);
}

void CelFileReaderTest::ReadSmallUInt16CelFile()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename(SMALL_UINT16_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	// Read some data
	CPPUNIT_ASSERT(data.GetFileHeader()->GetFilename() == SMALL_UINT16_CEL_FILE);

	CheckSmallCelFileHeader(data);

	float intensity;
	int16_t numPixels;
	bool outlier, masked;
	float stdev;

	int32_t cells = data.GetNumCells();

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;
	BoolVector vOutliers;
	BoolVector vMasked;

	// Get individual columns
	CPPUNIT_ASSERT(data.GetIntensities(0, cells, vInten));
	CPPUNIT_ASSERT(data.GetStdev(0, cells, vStdev));
	CPPUNIT_ASSERT(data.GetNumPixels(0, cells, vPixels));
	CPPUNIT_ASSERT(data.GetOutliers(0, cells, vOutliers));
	CPPUNIT_ASSERT(data.GetMasked(0, cells, vMasked));

	CPPUNIT_ASSERT(vInten.size() == 25);
	CPPUNIT_ASSERT(vStdev.size() == 25);
	CPPUNIT_ASSERT(vPixels.size() == 25);
	CPPUNIT_ASSERT(vOutliers.size() == 25);
	CPPUNIT_ASSERT(vMasked.size() == 25);

	
	// Check values
	for (int32_t cell = 0; cell < cells; ++cell)
	{
		CPPUNIT_ASSERT(vInten[cell] == 100*cell);
		CPPUNIT_ASSERT(vStdev[cell] == .5*cell);
		CPPUNIT_ASSERT(vPixels[cell] == 25);
		CheckOutlier(cell, vOutliers[cell]);
		CheckMasked(cell, vMasked[cell]);
	}

	XYCoordVector xyOutliers;
	XYCoordVector xyMasked;
	data.GetOutlierCoords(xyOutliers);
	data.GetMaskedCoords(xyMasked);
	CPPUNIT_ASSERT(xyOutliers.size() == 2);
	CPPUNIT_ASSERT(xyMasked.size() == 3);

	// Change the order of operations in the test
	for (int32_t cell = 0; cell < cells; ++cell)
	{
		CPPUNIT_ASSERT_NO_THROW(data.GetData(cell, intensity, stdev, numPixels, outlier, masked));

		CPPUNIT_ASSERT(intensity == 100*cell);
		CPPUNIT_ASSERT(stdev == .5*cell);
		CPPUNIT_ASSERT(numPixels == 25);

		CheckOutlier(cell, outlier);
		CheckMasked(cell, masked);
	}
}

void CelFileReaderTest::ReadMultiChannelCelFileTest()
{
	CelFileData data;
	CelFileReader reader;
	reader.SetFilename("../data/multi_channel.CEL");
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	//read intensities in a wavelength (ie channel or data group)
	data.SetActiveChannel(L"531");
	CPPUNIT_ASSERT(data.GetNumCells() == 2214144);
	FloatVector floats;
	CPPUNIT_ASSERT(data.GetIntensities(0, 3, floats));
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[0], 558.0, 0.1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[1], 3953.3, 0.1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[2], 637.5, 0.1);

	//read intensities in the next wavelength
	data.SetActiveChannel(L"609");
	CPPUNIT_ASSERT(data.GetNumCells() == 2214144);
	floats.clear();
	CPPUNIT_ASSERT(data.GetIntensities(3, 3, floats));
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[0], 3983.0, 0.1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[1], 477.0, 0.1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[2], 522.5, 0.1);

	//switch back to previous wavelength and read standard devs
	data.SetActiveChannel(L"531");
	CPPUNIT_ASSERT(data.GetNumCells() == 2214144);
	floats.clear();
	CPPUNIT_ASSERT(data.GetStdev(0, 3, floats));
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[0], 71.4, 0.1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[1], 240.8, 0.1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[2], 64.4, 0.1);

	//switch wavelength again and read standard devs
	data.SetActiveChannel(L"609");
	CPPUNIT_ASSERT(data.GetNumCells() == 2214144);
	floats.clear();
	CPPUNIT_ASSERT(data.GetStdev(3, 3, floats));
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[0], 317.0, 0.1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[1], 61.09, 0.1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(floats[2], 84.3, 0.1);
}
