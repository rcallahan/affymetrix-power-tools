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
#include "calvin_files/writers/test/MMDATFileWriterTest.h"
//
#include "calvin_files/data/src/ColumnInfo.h"
#include "calvin_files/parsers/src/DATFileReader.h"
#include "calvin_files/writers/src/DataSetWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( MMDATFileWriterTest );

void MMDATFileWriterTest::setUp()
{

}

void MMDATFileWriterTest::tearDown()
{

}

void MMDATFileWriterTest::CreationTest()
{
}

void MMDATFileWriterTest::WriteSmallFileTest()
{
	std::string name = "small_DAT_file";
	std::wstring type = L"Hg-small";
	int32_t rows = 4;
	int32_t cols = 5;
	int32_t pixelCount = rows*cols;
	DATData dataOut(name);
	dataOut.SetPixelCount(pixelCount);
	dataOut.SetStatsCount(1);
	dataOut.SetArrayType(type);
	dataOut.SetPixelSize(0.71f);
	dataOut.SetScannerType(L"M10");
	dataOut.SetScannerID(L"main");
	DateTime dtOut = DateTime::Parse(L"12-25-05T11:12:13Z");
	dataOut.SetScanDate(dtOut);
	dataOut.SetRows(rows);
	dataOut.SetCols(cols);

	u_int16_t* expectedPixels = new u_int16_t[rows*cols];

	{	// Force the destruction of the writer before using the reader.

		MMDATFileWriter writer(dataOut);

		// Use mem-map to add pixel intensities
		CPPUNIT_ASSERT(writer.Open());
		u_int16_t* pixel = writer.GetMappedPixelDataPtr();
		CPPUNIT_ASSERT(pixel != 0);

		for (int32_t row = 0, i = 0; row < rows; ++row)
		{
			// don't re-map data
			for (int32_t col = 0; col < cols; ++col, ++pixel, ++i)
			{
				*pixel = (u_int16_t)(1000*row+10*col);
				expectedPixels[i] = *pixel;
			}
		}

		Uint16Vector stats;
		stats.push_back(0);	// min
		stats.push_back(10000);	// max

		writer.WriteStats(stats);

		CPPUNIT_ASSERT(writer.Close());
	}

	// Check that the correct data was written to the file

	DATData dataIn;
	DATFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(name));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(dataIn));

	CPPUNIT_ASSERT(dataIn.GetFilename() == name);
	CPPUNIT_ASSERT(dataIn.GetRows() == rows);
	CPPUNIT_ASSERT(dataIn.GetCols() == cols);
	CPPUNIT_ASSERT(dataIn.GetArrayType() == type);
	CPPUNIT_ASSERT(dataIn.GetPixelSize() == 0.71f);
	CPPUNIT_ASSERT(dataIn.GetScannerType() == L"M10");
	CPPUNIT_ASSERT(dataIn.GetScannerID() == L"main");
	DateTime dtIn = dataIn.GetScanDate();
	CPPUNIT_ASSERT(dtIn.ToString() == dtOut.ToString());

	u_int16_t* pixels = new u_int16_t[rows*cols];
	CPPUNIT_ASSERT(dataIn.GetPixels(pixels, 0, dataIn.GetRows()));
	CPPUNIT_ASSERT(memcmp(expectedPixels, pixels, sizeof(u_int16_t)*dataIn.GetRows()*dataIn.GetCols()) == 0);

	u_int16_t min, max;
	CPPUNIT_ASSERT(dataIn.GetRange(min, max));
	CPPUNIT_ASSERT(min == 0);
	CPPUNIT_ASSERT(max == 10000);

	delete [] expectedPixels;
	delete [] pixels;
}

// Write a file larger than 200MB (1024*1024*200 bytes) to exercise the remapping logic.

void MMDATFileWriterTest::WriteLargeFileTest()
{
	std::string name = "large_DAT_file";
	std::wstring type = L"Hg-large";
	int32_t rows = 11000;
	int32_t cols = 11000;
	int32_t pixelCount = rows*cols;
	DATData dataOut(name);
	dataOut.SetPixelCount(pixelCount);
	dataOut.SetStatsCount(1);
	dataOut.SetArrayType(type);
	dataOut.SetPixelSize(0.5f);
	dataOut.SetScannerType(L"M10");
	dataOut.SetScannerID(L"lab02");
	DateTime dtOut = DateTime::Parse(L"12-25-05T11:12:14Z");
	dataOut.SetScanDate(dtOut);
	dataOut.SetRows(rows);
	dataOut.SetCols(cols);

	{	// Force the destruction of the writer before using the reader.

		MMDATFileWriter writer(dataOut);

		// Use mem-map to add pixel intensities
		CPPUNIT_ASSERT(writer.Open());

		u_int16_t inten = 10;
		bool stop = false;
		do
		{
			int32_t firstRow = writer.GetFirstPixelRowMapped();
			int32_t numRowsMapped = writer.GetPixelRowsMapped();
			u_int16_t* pixel = (u_int16_t*)writer.GetMappedPixelDataPtr();
			CPPUNIT_ASSERT(pixel != 0);
			for (int32_t row = 0; row < numRowsMapped; ++row)
			{
				for (int32_t col = 0; col < cols; ++col, ++pixel)
				{
					if (++inten > 46000)
					{
						inten = 0;
					}
					*pixel = inten;
				}
			}

			if (firstRow+numRowsMapped >= rows)
				stop = true;
			else
				writer.MapPixelData(firstRow+numRowsMapped, writer.GetMaxPixelRowsToMap());
		} while(stop == false);

		Uint16Vector stats;
		stats.push_back(0);	// min
		stats.push_back(46000);	// max

		writer.WriteStats(stats);

		CPPUNIT_ASSERT(writer.Close());
	}

	// Check that the correct data was written to the file

	DATData dataIn;
	DATFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(name));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(dataIn));

	CPPUNIT_ASSERT(dataIn.GetFilename() == name);
	CPPUNIT_ASSERT(dataIn.GetRows() == rows);
	CPPUNIT_ASSERT(dataIn.GetCols() == cols);
	CPPUNIT_ASSERT(dataIn.GetArrayType() == type);
	CPPUNIT_ASSERT(dataIn.GetPixelSize() == 0.5f);
	CPPUNIT_ASSERT(dataIn.GetScannerType() == L"M10");
	CPPUNIT_ASSERT(dataIn.GetScannerID() == L"lab02");
	DateTime dtIn = dataIn.GetScanDate();
	CPPUNIT_ASSERT(dtIn.ToString() == dtOut.ToString());

	int32_t MaxRowsToFetch = 1000;
	u_int16_t* pixels = new u_int16_t[MaxRowsToFetch*cols];
	u_int16_t* expectedPixels = new u_int16_t[MaxRowsToFetch*cols];

	u_int16_t inten = 10;
	int32_t startRow = 0;
	while (startRow < rows)
	{
		int32_t rowsToFetch = MaxRowsToFetch;
		if (startRow + rowsToFetch > rows)
			rowsToFetch = rows - startRow;

		// Compute the expected values
		for (int32_t row = 0, i = 0; row < rowsToFetch; ++row)
		{
			for (int32_t col = 0; col < cols; ++col)
			{
				if (++inten > 46000)
				{
					inten = 0;
				}
				expectedPixels[i++] = inten;
			}
		}

		// Compare results
		CPPUNIT_ASSERT(dataIn.GetPixels(pixels, startRow, rowsToFetch));
		CPPUNIT_ASSERT(memcmp(expectedPixels, pixels, sizeof(u_int16_t)*rowsToFetch*dataIn.GetCols()) == 0);

		startRow += rowsToFetch;
	}

	u_int16_t min, max;
	CPPUNIT_ASSERT(dataIn.GetRange(min, max));
	CPPUNIT_ASSERT(min == 0);
	CPPUNIT_ASSERT(max == 46000);

	delete [] expectedPixels;
	delete [] pixels;
} 
