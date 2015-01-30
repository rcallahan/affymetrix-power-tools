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
#include "calvin_files/writers/test/TestFileGenerator.h"
//
#include "calvin_files/array/src/ArrayId.h"
#include "calvin_files/data/src/ColumnInfo.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/utils/src/DateTime.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CDFCntrlFileWriter.h"
#include "calvin_files/writers/src/CDFFileWriter.h"
#include "calvin_files/writers/src/CalvinCelFileWriter.h"
#include "calvin_files/writers/src/DATFileUpdater.h"
#include "calvin_files/writers/src/DATFileWriter.h"
#include "calvin_files/writers/src/FileOutput.h"
#include "calvin_files/writers/src/GenericDataHeaderWriter.h"
#include "calvin_files/writers/src/GenericFileWriter.h"
//
#include "util/Fs.h"
//
#include <fstream>
#include <iostream>
#include <stdio.h>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_parameter;


CPPUNIT_TEST_SUITE_REGISTRATION( TestFileGenerator );

void TestFileGenerator::GenerateTestFiles()
{
	//
	// UNCOMMENT METHODS HERE TO GENERATE TEST FILES.
	//

//	//WriteOutGenericDATDataFile1UsingRawWrites();
//	//WriteOutGenericDATDataFile1UsingGenericWriter();	// mostly uses generic writer
	
	// uses the Calvin writers
//	WriteOutGenericDATDataFileNoGrid();
//	WriteOutGenericDATDataFileWithGrid();
//	WriteOutGenericDataFileWithAllColumnTypes();
//	WriteSmallCelFileNoOutlierNoMask();
//	WriteSmallCelFileNoStdev();
//	WriteSmallCelFile();
//	WriteLargeCelFile();
//	WriteCelFileWithADataSetWithZeroRows();
//	WriteSmallDatFile();
//	WriteLargeDatFile();
//	WriteActualSizeExpressionCDFFile();
//	WriteSmallExpressionCDFFile();
//	WriteSmallQCCDFFile();
//	WriteSmallDatFileWithGridAndSubgrids();
//	WriteSmallCelFileWithAPartialDatHeaderTest();
//	WriteSmallCelFileWithAFullDatHeaderTest();
//	//WriteSmallUint16CelFile();
//	WriteSmallDatFileWithGridAndSubgridsAndParameters();
//	WriteSmallDatFileWithReservedStringParameters();
}


void TestFileGenerator::WriteOutGenericDATDataFile1UsingRawWrites()
{
	// Open the file.
	ofstream fileStream;
        Fs::aptOpen(fileStream,  "test.file.data_dat.raw", ios_base::out | ios_base::binary | ios_base::trunc);

	//
	// Header
	//

	// magic number
	FileOutput::WriteUInt8(fileStream, 59);

	// version
	FileOutput::WriteInt8(fileStream, 1);

	// Number of data cubes - confirm this
	FileOutput::WriteUInt32(fileStream, 1);

	// offset to the data cube byte offset array
	int offsetLocation = fileStream.tellp();
	FileOutput::WriteUInt32(fileStream, 1);	//????

	//
	// Generic Header
	//

	// dataTypeIdentifier
	FileOutput::WriteString8(fileStream, SCAN_ACQUISITION_DATA_TYPE);

	// unique file identifier
	FileOutput::WriteString8(fileStream, "test-dat-guid");

	// Date and time creation
//	FileOutput::WriteString16(fileStream, DateTime::GetCurrentDateTime().ToString());
	FileOutput::WriteString16(fileStream, L"2004-07-04T11:12:13Z");

	// Locale
	FileOutput::WriteString16(fileStream, L"en-US");

	// Number of parameters
	FileOutput::WriteUInt32(fileStream, 3);

	// Name-value pairs
	// probe array type
	FileOutput::WriteString16(fileStream, ARRAY_TYPE_PARAM_NAME);
	FileOutput::WriteString16(fileStream, L"Hg-U133A");
	
	// barcode
	FileOutput::WriteString16(fileStream, ARRAY_BARCODE_PARAM_NAME);
	FileOutput::WriteString16(fileStream, L"Barcode");

	// Other
	FileOutput::WriteString16(fileStream, L"Parameter1");
	FileOutput::WriteString16(fileStream, L"Value1");

	// Number of parent files
	FileOutput::WriteUInt32(fileStream, 1);

	// Parent GenericDataHeader - array file

	// dataTypeIdentifier
	FileOutput::WriteString8(fileStream, ARRAY_FILE_TYPE_IDENTIFIER);

	// unique file identifier
	FileOutput::WriteString8(fileStream, "test-array-guid");

	// Date and time creation
//	FileOutput::WriteString16(fileStream, DateTime::GetCurrentDateTime().ToString());
	FileOutput::WriteString16(fileStream, L"2004-07-01T13:14:15Z");

	// Locale
	FileOutput::WriteString16(fileStream, L"en-US");

	// Number of parameters
	FileOutput::WriteUInt32(fileStream, 2);

	// Name-value pair
	// probe array type
	FileOutput::WriteString16(fileStream, ARRAY_TYPE_PARAM_NAME);
	FileOutput::WriteString16(fileStream, L"Hg-U133A");
	FileOutput::WriteString16(fileStream, ARRAY_LOT_PARAM_NAME);
	FileOutput::WriteString16(fileStream, L"Thanks alot");

	// Number of parent files
	FileOutput::WriteUInt32(fileStream, 0);

	//
	// Data Cube
	//

	int offset = fileStream.tellp();

	// Data cube header

	// Data cube name
	FileOutput::WriteString16(fileStream, L"acquired data");

	// Number of parameters
	FileOutput::WriteUInt32(fileStream, 2);

	// Name-value pairs
	FileOutput::WriteString16(fileStream, L"Scanner");
	FileOutput::WriteString16(fileStream, L"M10");
	FileOutput::WriteString16(fileStream, L"Pixel Size");
	FileOutput::WriteString16(fileStream, L"0.051");

	// Number of columns - there is one column of unsigned shorts (pixel data)
	FileOutput::WriteUInt32(fileStream, 1);

	// Array of value types
	FileOutput::WriteInt8(fileStream, 3);	// UShort

	// Array of value sizes
	FileOutput::WriteInt32(fileStream, 2);

	// Number of rows in the data cube
	int32_t rows = 100;
	FileOutput::WriteInt32(fileStream, rows);

	// Write out the data
	for( int32_t i=0; i < rows; ++i )
	{
		u_int16_t value = (u_int16_t)(i*10+i);
		FileOutput::WriteUInt16(fileStream, value);
	}

	// write the offset
	fileStream.seekp(offsetLocation);
	FileOutput::WriteUInt32(fileStream, offset);

	fileStream.close();

}

void TestFileGenerator::AddStandardGenericDataHeader(GenericDataHeader& gdh)
{
	// Fill the GenericDataHeader.
	gdh.SetFileTypeId(SCAN_ACQUISITION_DATA_TYPE);
	gdh.SetFileId("test-dat-guid");
	gdh.SetFileCreationTime(L"2004-07-04T11:12:13Z");
	gdh.SetLocale(L"en-US");
	ParameterNameValueType nvt;
	nvt.SetName(ARRAY_TYPE_PARAM_NAME);
	nvt.SetValueText(L"Hg-U133A");
	gdh.AddNameValParam(nvt);
	nvt.SetName(ARRAY_BARCODE_PARAM_NAME);
	nvt.SetValueText(L"Barcode");
	gdh.AddNameValParam(nvt);
	nvt.SetName(L"Parameter1");
	nvt.SetValueText(L"Value1");
	gdh.AddNameValParam(nvt);

	GenericDataHeader gdhParent;
	gdhParent.SetFileTypeId(ARRAY_TYPE_IDENTIFIER);
	gdhParent.SetFileId("test-array-guid");
	gdhParent.SetFileCreationTime(L"2004-07-01T13:14:15Z");
	gdhParent.SetLocale(L"en-US");
	nvt.SetName(ARRAY_TYPE_PARAM_NAME);
	nvt.SetValueText(L"Hg-U133A");
	gdhParent.AddNameValParam(nvt);
	nvt.SetName(ARRAY_LOT_PARAM_NAME);
	nvt.SetValueText(L"Thanks alot");
	gdhParent.AddNameValParam(nvt);

	gdh.AddParent(gdhParent);
}

// Uses a mix of the GenericFileWriter and raw commands
void TestFileGenerator::WriteOutGenericDATDataFile1UsingGenericWriter()
{
	GenericDataHeader gdh;
	AddStandardGenericDataHeader(gdh);

	// Fill the DataSetHeader
	DataSetHeader dph;
	dph.SetName(L"acquired data");
	ParameterNameValueType nvt;
	nvt.SetName(L"Scanner");
	nvt.SetValueText(L"M10");
	dph.AddNameValParam(nvt);
	nvt.SetName(L"Pixel Size");
	nvt.SetValueFloat(0.051f);
	dph.AddNameValParam(nvt);
	dph.AddColumn(UShortColumn(L"Pixel"));

	int32_t rows = 100;
	dph.SetRowCnt(rows);

	// Open the file.
	ofstream fileStream;
        Fs::aptOpen(fileStream, "test.file.data_dat", ios_base::out | ios_base::binary | ios_base::trunc);

	// Write the file header using raw methods until the writer is available.
	// magic number
	FileOutput::WriteUInt8(fileStream, 59);

	// version
	FileOutput::WriteInt8(fileStream, 1);

	// Number of data cubes - confirm this
	FileOutput::WriteUInt32(fileStream, 1);

	// offset to the data cube byte offset array
	int offsetLocation = fileStream.tellp();
	FileOutput::WriteUInt32(fileStream, 1);	//????

	// Write the GenericDataHeader to the file
	GenericDataHeaderWriter gdhWriter;
	gdhWriter.Write(fileStream, gdh);


	//
	// Data Cube
	//

	int offset = fileStream.tellp();

	// Write the DataSetHeader
	DataSetWriter dphWriter(&fileStream, &dph);
	dphWriter.WriteHeader();

	// Write out the data
	for( int32_t i=0; i < rows; ++i )
	{
		u_int16_t value = (u_int16_t)(i*10+i);
		FileOutput::WriteUInt16(fileStream, value);
	}

	// write the offset
	fileStream.seekp(offsetLocation);
	FileOutput::WriteUInt32(fileStream, offset);

	fileStream.close();

}

/*
void TestFileGenerator::WriteOutGenericDATDataFile1UsingGenericWriter_ForReal()
{
	GenericDataHeader gdh;

	// Fill the GenericDataHeader.
	gdh.SetFileTypeId(SCAN_ACQUISITION_DATA_TYPE);
	gdh.SetFileId("test-dat-guid");
	gdh.SetFileCreationTime(L"2004-07-04T11:12:13Z");
	gdh.SetLocale(L"en-US");
	gdh.AddNameValPair(ARRAY_TYPE_PARAM_NAME, L"Hg-U133A");
	gdh.AddNameValPair(ARRAY_BARCODE_PARAM_NAME, L"Barcode");
	gdh.AddNameValPair(L"Parameter1", L"Value1");

	GenericDataHeader gdhParent;
	gdhParent.SetFileTypeId(ARRAY_FILE_TYPE_IDENTIFIER);
	gdhParent.SetFileId("test-array-guid");
	gdhParent.SetFileCreationTime(L"2004-07-01T13:14:15Z");
	gdhParent.SetLocale(L"en-US");
	gdhParent.AddNameValPair(ARRAY_TYPE_PARAM_NAME, L"Hg-U133A");
	gdhParent.AddNameValPair(ARRAY_LOT_PARAM_NAME, L"Thanks alot");

	gdh.AddParent(gdhParent);

	// Fill the DataSetHeader
	DataSetHeader dph;
	dph.SetName(L"acquired data");
	dph.AddNameValPair(L"Scanner", L"M10");
	dph.AddNameValPair(L"Pixel Size", L"0.051");

	dph.AddColumnType(UShortColumnType());

	int32_t rows = 100;
	dph.SetRowCnt(rows);

	// Set the FileHeader
	FileHeader fh;
	fh.SetFilename("test.file.data_dat");
	fh.SetMagicNumber(59);	// remove these
	fh.SetVersion(1);	// remove these
	fh.SetGenericDataHdr(gdh);
	fh.AddDataGroupHdr(dph);

	// Create the generic file writer
	GenericFileWriter gfWriter(fh);
	gfWriter.WriteHeader();

	DataSetWriterIt dcwBegin, dcwEnd;
	gfWriter.GetDataSetWriters(dcwBegin, dcwEnd);

	dcwBegin->WriteHeader();

	// Write out the data
	for( int32_t i=0; i < rows; ++i )
	{
		u_int16_t value = (u_int16_t)(i*10+i);
		dcwBegin->Write(value);
	}
}*/

void TestFileGenerator::WriteOutGenericDATDataFileNoGrid()
{
	GenericDataHeader gdh;
	AddStandardGenericDataHeader(gdh);

	// Fill the DataGroupHeader
	DataGroupHeader dch;
	dch.SetName(L"First Data Cube");

	// Fill the DataSetHeader
	DataSetHeader dph;
	dph.SetName(L"acquired data");
	ParameterNameValueType nvt;
	nvt.SetName(L"Scanner");
	nvt.SetValueText(L"M10");
	dph.AddNameValParam(nvt);
	nvt.SetName(L"Pixel Size");
	nvt.SetValueFloat(0.051f);
	dph.AddNameValParam(nvt);

	dph.AddColumn(UShortColumn(L"Pixel"));

	int32_t rows = 100;
	dph.SetRowCnt(rows);

	dch.AddDataSetHdr(dph);

	// Set the FileHeader
	FileHeader fh;
	fh.SetFilename("test.file.data_dat");
	fh.SetGenericDataHdr(gdh);
	fh.AddDataGroupHdr(dch);

	// Create the generic file writer
	GenericFileWriter gfWriter(&fh);
	gfWriter.WriteHeader();

	DataGroupWriterIt dcwBegin, dcwEnd;
	gfWriter.GetDataGroupWriters(dcwBegin, dcwEnd);

	DataGroupWriter d = *dcwBegin;
	dcwBegin->WriteHeader();

	DataSetWriterIt dpwBegin, dpwEnd;
	dcwBegin->GetDataSetWriters(dpwBegin, dpwEnd);

	dpwBegin->WriteHeader();

	// Write out the data
	for( int32_t i=0; i < rows; ++i )
	{
		u_int16_t value = (u_int16_t)(i*10+i);
		dpwBegin->Write(value);
	}

	dpwBegin->UpdateNextDataSetOffset();
	dcwBegin->Close();
}


void TestFileGenerator::WriteOutGenericDATDataFileWithGrid()
{
	GenericDataHeader gdh;
	AddStandardGenericDataHeader(gdh);

	// Fill the DataGroupHeader
	DataGroupHeader dch;
	dch.SetName(L"");	// unnamed DataGroup

	// Fill the pixel intensity DataSetHeader
	DataSetHeader dphPixel;
	dphPixel.SetName(L"acquired data");
	ParameterNameValueType nvt;
	nvt.SetName(L"Scanner");
	nvt.SetValueText(L"M10");
	dphPixel.AddNameValParam(nvt);
	nvt.SetName(L"Pixel Size");
	nvt.SetValueFloat(0.051f);
	dphPixel.AddNameValParam(nvt);

	dphPixel.AddColumn(UShortColumn(L"Pixel"));

	int32_t rows = 1000;
	dphPixel.SetRowCnt(rows);

	dch.AddDataSetHdr(dphPixel);

	// Fill the grid DataSetHeader
	DataSetHeader dphGrid;
	dphGrid.SetName(L"grid position");
	nvt.SetName(L"GhostGrids");
	nvt.SetValueText(L"True");
	dphGrid.AddNameValParam(nvt);
	nvt.SetName(L"Pixel Size");
	nvt.SetValueFloat(0.051f);
	dphGrid.AddNameValParam(nvt);

	dphGrid.AddColumn(FloatColumn(L"Upper left x"));
	dphGrid.AddColumn(FloatColumn(L"Upper left y"));
	dphGrid.AddColumn(FloatColumn(L"Upper right x"));
	dphGrid.AddColumn(FloatColumn(L"Upper right y"));
	dphGrid.AddColumn(FloatColumn(L"Lower right x"));
	dphGrid.AddColumn(FloatColumn(L"Lower right y"));
	dphGrid.AddColumn(FloatColumn(L"Lower left x"));
	dphGrid.AddColumn(FloatColumn(L"Lower left y"));

	int32_t grids = 5;	// first is the global grid with 4 subgrids
	dphGrid.SetRowCnt(grids);

	dch.AddDataSetHdr(dphGrid);

	// Set the FileHeader
	FileHeader fh;
	fh.SetFilename("test.file.data_dat_with_grid");
	fh.SetGenericDataHdr(gdh);
	fh.AddDataGroupHdr(dch);

	// Create the generic file writer
	GenericFileWriter gfWriter(&fh);
	gfWriter.WriteHeader();

	DataGroupWriterIt dcwBegin, dcwEnd;
	gfWriter.GetDataGroupWriters(dcwBegin, dcwEnd);

	DataGroupWriter d = *dcwBegin;
	dcwBegin->WriteHeader();

	DataSetWriterIt dpwBegin, dpwEnd;
	dcwBegin->GetDataSetWriters(dpwBegin, dpwEnd);

	// Write out the pixel DataSet

	dpwBegin->WriteHeader();

	for( int32_t i=0; i < rows; ++i )
	{
		u_int16_t value = (u_int16_t)(i*10+i);
		dpwBegin->Write(value);
	}

	dpwBegin->UpdateNextDataSetOffset();

	++dpwBegin;

	// Write out the grid DataSet

	dpwBegin->WriteHeader();

	for( int32_t i=0; i < grids; ++i )
	{
		for (int32_t corner = 0; corner < 4; ++corner)
		{
			float value = (float)(i*100 + corner);
			dpwBegin->Write(value);
			dpwBegin->Write(value);
		}
	}

	dpwBegin->UpdateNextDataSetOffset();

	dcwBegin->Close();
}

void TestFileGenerator::WriteOutGenericDataFileWithAllColumnTypes()
{
	GenericDataHeader gdh;
	AddStandardGenericDataHeader(gdh);

	// Fill the DataGroupHeader
	DataGroupHeader dch;
	dch.SetName(L"default");	// default DataGroup

	// Fill the all types DataSetHeader
	DataSetHeader dphAT;
	dphAT.SetName(L"all types");
	ParameterNameValueType nvt;
	nvt.SetName(L"How many types");
	nvt.SetValueText(L"All types");
	dphAT.AddNameValParam(nvt);
	nvt.SetName(L"Powered by");
	nvt.SetValueText(L"Affymetrix");
	dphAT.AddNameValParam(nvt);

	dphAT.AddColumn(ByteColumn(L"Byte type"));
	dphAT.AddColumn(UByteColumn(L"UByte type"));
	dphAT.AddColumn(ASCIIColumn(L"ASCII type", 10));
	dphAT.AddColumn(ShortColumn(L"Short type"));
	dphAT.AddColumn(UShortColumn(L"UShort type"));
	dphAT.AddColumn(IntColumn(L"Int type"));
	dphAT.AddColumn(UIntColumn(L"UInt type"));
	dphAT.AddColumn(UnicodeColumn(L"Unicode type", 15));
	dphAT.AddColumn(FloatColumn(L"Float type"));

	int32_t rows = 2;
	dphAT.SetRowCnt(rows);

	dch.AddDataSetHdr(dphAT);

	// Set the FileHeader
	FileHeader fh;
	fh.SetFilename("test.file.data_all_column_types");
	fh.SetGenericDataHdr(gdh);
	fh.AddDataGroupHdr(dch);

	// Create the generic file writer
	GenericFileWriter gfWriter(&fh);
	gfWriter.WriteHeader();

	DataGroupWriterIt dcwBegin, dcwEnd;
	gfWriter.GetDataGroupWriters(dcwBegin, dcwEnd);

	DataGroupWriter d = *dcwBegin;
	dcwBegin->WriteHeader();

	DataSetWriterIt dpwBegin, dpwEnd;
	dcwBegin->GetDataSetWriters(dpwBegin, dpwEnd);

	// Write out the all types DataSet

	dpwBegin->WriteHeader();

	for( int32_t row = 0; row < rows; ++row )
	{
		char str[10];
		wchar_t wstr[15];
		int8_t b = 1+10*row;
		u_int8_t ub = 2+10*row;
		sprintf(str, "%d", 3+10*row);
		int16_t s = 4+10*row;
		u_int16_t us = 5+10*row;
		int32_t i = 6+10*row;
		u_int32_t ui = 7+10*row;
		FormatString1(wstr, 15, L"%d", 8+10*row);
		float f = 9+10*row;

		dpwBegin->Write(b);	// btye
		dpwBegin->Write(ub);	// unsigned byte
		dpwBegin->Write(str, 10);	// ACSII string
		dpwBegin->Write(s);	// short
		dpwBegin->Write(us);	// unsigned short
		dpwBegin->Write(i);	// int
		dpwBegin->Write(ui);	// unsigned int
		dpwBegin->Write(wstr, 15);	// Unicode string
		dpwBegin->Write(f);	// float
	}

	dpwBegin->UpdateNextDataSetOffset();
	dcwBegin->Close();
}

void TestFileGenerator::WriteSmallCelFile()
{
	CelFileData data("small_cel_file");
	data.SetIntensityCount(25);
	data.SetStdDevCount(25);
	data.SetPixelCount(25);
	data.SetOutlierCount(2);
	data.SetMaskCount(3);
	data.SetRows(5);
	data.SetCols(5);
	data.SetArrayType(L"Hg-small");
	data.SetAlgorithmName(L"Feature Extraction");
	ParameterNameValueType nvt;
	nvt.SetName(L"percentile");
	nvt.SetValueFloat(0.75f);
	data.AddAlgorithmParameter(nvt);
	nvt.SetName(L"outlierlow");
	nvt.SetValueFloat(1.004f);
	data.AddAlgorithmParameter(nvt);

	CelFileWriter* writer = new CelFileWriter(data);

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;

	for (int i=0; i<25; ++i)
	{
		vInten.push_back(100.0f*i);
		vStdev.push_back(.5*i);
		vPixels.push_back(25);
	}

	// Do some writing
	writer->WriteIntensities(vInten);
	writer->WriteStdDevs(vStdev);
	writer->WritePixels(vPixels);

//	XYCoordVector
	XYCoordVector outlier;
	XYCoord xy(0,0);
	outlier.push_back(xy);
	xy.xCoord = 1;
	xy.yCoord = 2;
	outlier.push_back(xy);
	writer->WriteOutlierCoords(outlier);

	XYCoordVector masked;
	xy.xCoord = 1;
	xy.yCoord = 0;
	masked.push_back(xy);
	xy.xCoord = 2;
	xy.yCoord = 1;
	masked.push_back(xy);
	xy.xCoord = 3;
	xy.yCoord = 2;
	masked.push_back(xy);
	writer->WriteMaskCoords(masked);

	delete writer;
}

// Requires changes in CelFileWriter and CELData.
//
//void TestFileGenerator::WriteSmallUint16CelFile()
//{
//	CelFileData data("small_uint16_cel_file");
//	data.SetIntensityCount(25);
//	data.SetStdDevCount(25);
//	data.SetPixelCount(25);
//	data.SetOutlierCount(2);
//	data.SetMaskCount(3);
//	data.SetRows(5);
//	data.SetCols(5);
//	data.SetArrayType(L"Hg-small");
//	data.SetAlgorithmName(L"Feature Extraction");
//	ParameterNameValueType nvt;
//	nvt.SetName(L"percentile");
//	nvt.SetValueFloat(0.75f);
//	data.AddAlgorithmParameter(nvt);
//	nvt.SetName(L"outlierlow");
//	nvt.SetValueFloat(1.004f);
//	data.AddAlgorithmParameter(nvt);
//
//	CelFileWriter* writer = new CelFileWriter(data);
//
//	Uint16Vector vInten;
//	FloatVector vStdev;
//	Int16Vector vPixels;
//
//	for (int i=0; i<25; ++i)
//	{
//		vInten.push_back(100*i);
//		vStdev.push_back(.5*i);
//		vPixels.push_back(25);
//	}
//
//	// Do some writing
//	writer->WriteIntensities(vInten);
//	writer->WriteStdDevs(vStdev);
//	writer->WritePixels(vPixels);
//
//	//XYCoordVector
//	XYCoordVector outlier;
//	XYCoord xy(0,0);
//	outlier.push_back(xy);
//	xy.xCoord = 1;
//	xy.yCoord = 2;
//	outlier.push_back(xy);
//	writer->WriteOutlierCoords(outlier);
//
//	XYCoordVector masked;
//	xy.xCoord = 1;
//	xy.yCoord = 0;
//	masked.push_back(xy);
//	xy.xCoord = 2;
//	xy.yCoord = 1;
//	masked.push_back(xy);
//	xy.xCoord = 3;
//	xy.yCoord = 2;
//	masked.push_back(xy);
//	writer->WriteMaskCoords(masked);
//
//	delete writer;
//}

void TestFileGenerator::WriteRemaingSmallCelFileWithGridParameters(CelFileData& data)
{
	data.SetIntensityCount(25);
	data.SetStdDevCount(25);
	data.SetPixelCount(25);
	data.SetOutlierCount(2);
	data.SetMaskCount(3);
	data.SetRows(5);
	data.SetCols(5);
	data.SetArrayType(L"Hg-small");
    data.SetLibraryPackageName(L"Hg-small-lib-package");
    data.SetMasterFileName(L"Hg-small-master-file");
	data.SetAlgorithmName(L"Feature Extraction");
	ParameterNameValueType nvt;
	nvt.SetName(L"percentile");
	nvt.SetValueFloat(0.75f);
	data.AddAlgorithmParameter(nvt);
	nvt.SetName(L"outlierlow");
	nvt.SetValueFloat(1.004f);
	data.AddAlgorithmParameter(nvt);
	nvt.SetName(L"CellMargin");
	nvt.SetValueInt32(2);
	data.AddAlgorithmParameter(nvt);

	// Add grid
	const wchar_t* gridParams[] = {
    L"GridULX", L"GridULY", L"GridURX", L"GridURY",
    L"GridLRX", L"GridLRY", L"GridLLX", L"GridLLY"};
	for (int32_t i = 0; i < 8; ++i)
	{
		nvt.SetName(gridParams[i]);
		nvt.SetValueFloat(2.0f + (float)i);
		data.AddAlgorithmParameter(nvt);
	}

	CelFileWriter* writer = new CelFileWriter(data);

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;

	for (int i=0; i<25; ++i)
	{
		vInten.push_back(100.0f*i);
		vStdev.push_back(.5*i);
		vPixels.push_back(25);
	}

	// Do some writing
	writer->WriteIntensities(vInten);
	writer->WriteStdDevs(vStdev);
	writer->WritePixels(vPixels);

//	XYCoordVector
	XYCoordVector outlier;
	XYCoord xy(0,0);
	outlier.push_back(xy);
	xy.xCoord = 1;
	xy.yCoord = 2;
	outlier.push_back(xy);
	writer->WriteOutlierCoords(outlier);

	XYCoordVector masked;
	xy.xCoord = 1;
	xy.yCoord = 0;
	masked.push_back(xy);
	xy.xCoord = 2;
	xy.yCoord = 1;
	masked.push_back(xy);
	xy.xCoord = 3;
	xy.yCoord = 2;
	masked.push_back(xy);
	writer->WriteMaskCoords(masked);

	delete writer;
}

void TestFileGenerator::WriteSmallCelFileWithAPartialDatHeaderTest()
{
	CelFileData data("small_cel_file_partial_datheader");

	// Write 
	ParameterNameValueType nvt;
	GenericDataHeader datHdr;
	datHdr.SetFileId(AffymetrixGuid::GenerateNewGuid());
	datHdr.SetFileTypeId("affymetrix-calvin-scan-acquisition");
	datHdr.SetFileCreationTime(L"2004-07-01T13:14:15Z");
	nvt.SetName(L"affymetrix-partial-dat-header");
	std::wstring datHeaderString = L"  small_cel_file_partial_datheader:CLS=25   RWS=25   XIN=1  YIN=1  VE=0         0   05/19/05 02:45:59 ScannerID:  ScannerTyp   \x14  \x14 Hg-Small.1sq \x14  \x14  \x14  \x14  \x14 570 \x14 45.200001 \x14 0.340000 \x14 1.0900 \x14 3";
	nvt.SetValueText(datHeaderString);
	datHdr.AddNameValParam(nvt);
	nvt.SetName(L"affymetrix-max-pixel-intensity");
	nvt.SetValueUInt16(46001);
	datHdr.AddNameValParam(nvt);
	nvt.SetName(L"affymetrix-min-pixel-intensity");
	nvt.SetValueUInt16(1);
	datHdr.AddNameValParam(nvt);

	// Add DAT GenericDataHeader as parent.
	data.GetFileHeader()->GetGenericDataHdr()->AddParent(datHdr);

	WriteRemaingSmallCelFileWithGridParameters(data);
}

void TestFileGenerator::WriteSmallCelFileWithAFullDatHeaderTest()	// Files converted from GCOS will have a full DatHeader
{
	CelFileData data("small_cel_file_full_datheader");

	// Write 
	ParameterNameValueType nvt;
	GenericDataHeader datHdr;
	datHdr.SetFileId(AffymetrixGuid::GenerateNewGuid());
	datHdr.SetFileTypeId("affymetrix-calvin-scan-acquisition");
	datHdr.SetFileCreationTime(L"2004-07-01T13:14:15Z");
	nvt.SetName(L"affymetrix-dat-header");
	std::wstring datHeaderString = L"[45..56789]  small_cel_file_full_datheader:CLS=25   RWS=25   XIN=1  YIN=1  VE=0         0   05/19/05 02:45:59 ScannerID:  ScannerTyp   \x14  \x14 Hg-Small.1sq \x14  \x14  \x14  \x14  \x14 570 \x14 45.200001 \x14 0.340000 \x14 1.0900 \x14 3";
	nvt.SetValueText(datHeaderString);
	datHdr.AddNameValParam(nvt);

	// Add DAT GenericDataHeader as parent.
	data.GetFileHeader()->GetGenericDataHdr()->AddParent(datHdr);

	WriteRemaingSmallCelFileWithGridParameters(data);
}

void TestFileGenerator::WriteSmallCelFileNoOutlierNoMask()
{
	CelFileData data("small_cel_file_no_outlier_no_mask");
	data.SetIntensityCount(10);
	data.SetStdDevCount(10);
	data.SetPixelCount(10);
	data.SetOutlierCount(0);
	data.SetMaskCount(0);
	data.SetRows(2);
	data.SetCols(5);
	data.SetArrayType(L"Hg-small");
	data.SetAlgorithmName(L"Feature Extraction");
	ParameterNameValueType nvt;
	nvt.SetName(L"percentile");
	nvt.SetValueFloat(0.75f);
	data.AddAlgorithmParameter(nvt);
	nvt.SetName(L"outlierlow");
	nvt.SetValueFloat(1.004f);
	data.AddAlgorithmParameter(nvt);

	CelFileWriter* writer = new CelFileWriter(data);

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;

	for (int i=0; i<10; ++i)
	{
		vInten.push_back(100.0f*i);
		vStdev.push_back(.5*i);
		vPixels.push_back(25);
	}

	// Do some writing
	writer->WriteIntensities(vInten);
	writer->WriteStdDevs(vStdev);
	writer->WritePixels(vPixels);

	delete writer;
}

void TestFileGenerator::WriteSmallCelFileNoStdev()
{
	CelFileData data("small_cel_file_no_stdev");
	data.SetIntensityCount(25);
	data.SetStdDevCount(0);
	data.SetPixelCount(25);
	data.SetOutlierCount(2);
	data.SetMaskCount(3);
	data.SetRows(5);
	data.SetCols(5);
	data.SetArrayType(L"Hg-small");
	data.SetAlgorithmName(L"Feature Extraction");
	ParameterNameValueType nvt;
	nvt.SetName(L"percentile");
	nvt.SetValueFloat(0.75f);
	data.AddAlgorithmParameter(nvt);
	nvt.SetName(L"outlierlow");
	nvt.SetValueFloat(1.004f);
	data.AddAlgorithmParameter(nvt);

	CelFileWriter* writer = new CelFileWriter(data);

	FloatVector vInten;
	Int16Vector vPixels;

	for (int i=0; i<25; ++i)
	{
		vInten.push_back(100.0f*i);
		vPixels.push_back(25);
	}

	// Do some writing
	writer->WriteIntensities(vInten);
	writer->WritePixels(vPixels);

//	XYCoordVector
	XYCoordVector outlier;
	XYCoord xy(0,0);
	outlier.push_back(xy);
	xy.xCoord = 1;
	xy.yCoord = 2;
	outlier.push_back(xy);
	writer->WriteOutlierCoords(outlier);

	XYCoordVector masked;
	xy.xCoord = 1;
	xy.yCoord = 0;
	masked.push_back(xy);
	xy.xCoord = 2;
	xy.yCoord = 1;
	masked.push_back(xy);
	xy.xCoord = 3;
	xy.yCoord = 2;
	masked.push_back(xy);
	writer->WriteMaskCoords(masked);


	delete writer;
}

void TestFileGenerator::WriteLargeCelFile()
{
	CelFileData data("large_cel_file");
	data.SetRows(2560);
	data.SetCols(2560);
	int32_t cells = 2560*2560;
	data.SetIntensityCount(cells);
	data.SetStdDevCount(cells);
	data.SetPixelCount(cells);
	data.SetOutlierCount(2000);
	data.SetMaskCount(3000);
	ParameterNameValueType nvt;
	nvt.SetName(L"percentile");
	nvt.SetValueFloat(0.75f);
	data.AddAlgorithmParameter(nvt);
	nvt.SetName(L"outlierlow");
	nvt.SetValueFloat(1.004f);
	data.AddAlgorithmParameter(nvt);

	CelFileWriter* writer = new CelFileWriter(data);

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;

	vInten.resize(cells);
	vStdev.resize(cells);
	vPixels.resize(cells);

	int32_t intenTestValues[4] = {55683.0f, 4568.0f, 2368.0f, 100.0f};
	float stdevTestValues[4] = {2.345f, 56.23f, 1.53f, 3.875f};

	//int32_t testIdx = 0;
	for (int32_t i=0, testIdx=0; i<cells; ++i, ++testIdx)
	{
		if (testIdx >= 4)
			testIdx = 0;

		vInten[i] = intenTestValues[testIdx];
		vStdev[i] = stdevTestValues[testIdx];
		vPixels[i] = 25;
	}

	// Do some writing
	writer->WriteIntensities(vInten);
	writer->WriteStdDevs(vStdev);
	writer->WritePixels(vPixels);

	XYCoord xy(0,0);
	XYCoordVector outlier;
	for (int32_t i=0, col=2000; i<2000; ++i, --col)
	{
		xy.xCoord = i;
		xy.yCoord = col;
		outlier.push_back(xy);
	}
	writer->WriteOutlierCoords(outlier);

	XYCoordVector masked;
	for (int32_t i=0; i<1500; ++i)
	{
		xy.xCoord = 1200;
		xy.yCoord = i;
		masked.push_back(xy);
		xy.xCoord = 2400;
		xy.yCoord = i;
		masked.push_back(xy);
	}
	writer->WriteMaskCoords(masked);

	delete writer;
}

void TestFileGenerator::WriteCelFileWithADataSetWithZeroRows()
{
	CelFileData data("small_cel_file_with_dataset_of_zero_rows");
	data.SetIntensityCount(25);
	data.SetStdDevCount(25);
	data.SetPixelCount(25);
	data.SetOutlierCount(0);
	data.SetMaskCount(0);
	data.SetRows(5);
	data.SetCols(5);
	data.SetArrayType(L"Hg-small");
	data.SetAlgorithmName(L"Feature Extraction");
	ParameterNameValueType nvt;
	nvt.SetName(L"percentile");
	nvt.SetValueFloat(0.75f);
	data.AddAlgorithmParameter(nvt);
	nvt.SetName(L"outlierlow");
	nvt.SetValueFloat(1.004f);
	data.AddAlgorithmParameter(nvt);

	CelFileWriter* writer = new CelFileWriter(data);

	FloatVector vInten;
	FloatVector vStdev;
	Int16Vector vPixels;

	for (int i=0; i<25; ++i)
	{
		vInten.push_back(100.0f*i);
		vStdev.push_back(.123f);
		vPixels.push_back(25);
	}

	// Do some writing
	writer->WriteIntensities(vInten);
	writer->WriteStdDevs(vStdev);
	writer->WritePixels(vPixels);

	XYCoordVector outlier;
	writer->WriteOutlierCoords(outlier);

	XYCoordVector masked;
	writer->WriteMaskCoords(masked);

	delete writer;
}

void TestFileGenerator::WriteDatFile(const std::string& name, const std::wstring& type, int32_t rows, int32_t cols, bool showProgress)
{
	cout << "Started writing " << name.c_str() << endl;

	int32_t pixelCount = rows*cols;
	DATData data(name);
	data.SetPixelCount(pixelCount);
	data.SetStatsCount(1);
	data.SetArrayType(type);
	data.SetPixelSize(0.71f);
	data.SetScannerType(L"M10");
	data.SetScannerID(L"main");
	DateTime dt = DateTime::Parse(L"2005-12-25T11:12:13Z");
	data.SetScanDate(dt);
	data.SetRows(rows);
	data.SetCols(cols);
	DATFileWriter writer(data);

	Uint16Vector pixels;
	pixels.reserve(pixelCount);

	u_int16_t inten = 10;
	for (int32_t i=0; i<pixelCount; ++i, ++inten)
	{
		if (inten > 46000)
		{
			inten = 0;
			if (showProgress)
				cout << "Written " << i << " pixels so far..." << endl;
		}
		pixels.push_back(inten);
	}

	Uint16Vector stats;
	stats.push_back(pixels[0]);	// min
	stats.push_back(pixels[19]);	// max

	writer.WriteStats(stats);
	writer.WritePixels(pixels);

	cout << "Completed writing " << name.c_str() << endl;
}

void TestFileGenerator::AddGridAndSubgrids(DATData& data, float increment, int32_t subgridCnt)
{
	// Write the grid
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
	data.SetGlobalGrid(ggStatus, ggRegion);

	for(int n = 0; n < subgridCnt; n++)
	{
		FRegion sgRegion;
		float ptVal = 0.0 + (float)n;
		for(int i = 0; i < 4; i++)
		{
			FPoint point;
			point.x = ptVal;
			ptVal += increment;
			point.y = ptVal;
			sgRegion.pts.push_back(point);
		}

		u_int32_t sgStatus = DATData::GridError | DATData::GridManualAdjust;
		if (n % 2 == 0)
			sgStatus = DATData::GridOK | DATData::GridManualAdjust;

		data.AddSubgrid(sgStatus, sgRegion);
	}
}

void TestFileGenerator::WriteSmallDatFile()
{
	WriteDatFile("small_DAT_file", L"Hg-small", 4, 5, true);
}

void TestFileGenerator::WriteLargeDatFile()
{
	WriteDatFile("large_DAT_file", L"Hg-large", 25100, 25100, true);
}

void TestFileGenerator::WriteSmallDatFileWithGridAndSubgrids()
{
	std::string filename = "small_DAT_file_with_subgrids";
	WriteDatFile(filename, L"Hg-small", 4, 5, true);

	DATData data(filename);
	AddGridAndSubgrids(data, 2.5, 10);

	DATFileUpdater update(data);
	update.Update();
}

void TestFileGenerator::WriteExpressionCDFFile(const std::string& filename, u_int32_t probeSetCnt)
{
	wchar_t name[100];

	u_int8_t unitType = 3;
	u_int8_t direction = 1;
	u_int32_t atoms = 11;
	u_int8_t cellsPerAtom = 2;
	u_int32_t cells = atoms*cellsPerAtom;

	CDFData data(filename);
	data.SetProbeSetCnt(probeSetCnt, Expression);
	data.SetArrayCols(atoms);
	data.SetArrayRows(probeSetCnt*cellsPerAtom);
	CDFFileWriter writer(data);

	for (u_int32_t i = 0; i < probeSetCnt; ++i)
	{
		// make the unit name
		FormatString1(name, 100, L"biob_%d", i);

		writer.OpenDataGroup(name, 1);
		CDFProbeSetWriter* probeWriter = writer.CreateProbeSetWriter(name,unitType,direction,atoms,cells,i,cellsPerAtom);
		probeWriter->WriteHeader();

		for (u_int32_t atom = 0; atom < atoms; ++atom)
		{
			for (u_int32_t cellInAtom = 0; cellInAtom < cellsPerAtom; ++cellInAtom)
			{
				probeWriter->Write(atom,i*cellsPerAtom+cellInAtom,atom,atom,'C','G');
			}
		}
		probeWriter->Close();
		delete probeWriter;
		writer.CloseDataGroup();	
	}
}

void TestFileGenerator::WriteActualSizeExpressionCDFFile()
{
	WriteExpressionCDFFile("actualsize_CDF_file", 50000);
}

void TestFileGenerator::WriteSmallExpressionCDFFile()
{
	WriteExpressionCDFFile("small_CDF_file", 10);
}

void TestFileGenerator::WriteQCCDFFile(const std::string& filename, u_int32_t probeSetCnt)
{
	wchar_t name[100];

	u_int32_t atoms = 11;
	u_int8_t cellsPerAtom = 2;
	//u_int32_t cells = atoms*cellsPerAtom;

	u_int8_t probeLength = 25;
	u_int8_t perfectMatchFlag = 0;
	u_int8_t backgroundProbeFlag = 0;

	CDFData data(filename);
	data.SetProbeSetCnt(probeSetCnt, Control);
	data.SetArrayCols(atoms);
	data.SetArrayRows(probeSetCnt*cellsPerAtom);
	CDFCntrlFileWriter writer(data);

	for (u_int32_t i = 0; i < probeSetCnt; ++i)
	{
		// make the unit name
		FormatString1(name, 100, L"control_%d", i);

		writer.OpenDataGroup(name, 1);
		CDFCntrlProbeSetWriter* probeWriter = writer.GetCntrlProbeSetWriter(name);
		probeWriter->WriteHeader();

		for (u_int32_t atom = 0; atom < atoms; ++atom)
		{
			for (u_int32_t cellInAtom = 0; cellInAtom < cellsPerAtom; ++cellInAtom)
			{
				perfectMatchFlag = (cellInAtom % 2)? 1 : 0;
				backgroundProbeFlag = (atom % 2) ? 1 : 0;
				probeWriter->Write(atom,i*cellsPerAtom+cellInAtom, probeLength, perfectMatchFlag, backgroundProbeFlag);
			}
		}
		probeWriter->Close();
		delete probeWriter;
		writer.CloseDataGroup();	
	}
}

void TestFileGenerator::WriteSmallQCCDFFile()
{
	WriteQCCDFFile("small_QCCDF_file", 10);
}

void TestFileGenerator::WriteSmallDatFileWithGridAndSubgridsAndParameters()
{
	std::string filename = "small_DAT_file_with_subgrids_and_parameters";
	WriteDatFile(filename, L"Hg-small", 4, 5, true);

	DATData data(filename);
	AddGridAndSubgrids(data, 2.5, 10);

	ParameterNameValueType nvt;
	nvt.SetName(L"Zip Code");
	nvt.SetValueInt32(95051);
	data.AddGridAlignmentAlgorithmParameter(nvt);
	nvt.SetName(L"County");
	nvt.SetValueText(L"Santa Clara");
	data.AddGridAlignmentAlgorithmParameter(nvt);

	DATFileUpdater update(data);
	update.Update();
}

void TestFileGenerator::WriteSmallDatFileWithReservedStringParameters()
{
	int32_t rows = 4;
	int32_t cols = 5;
	int32_t pixelCount = rows*cols;
	DATData data("small_dat_file_with_reserved_string_parameters");
	data.SetPixelCount(pixelCount);
	data.SetStatsCount(1);
	data.SetArrayType(L"Hg_small");
	data.SetPixelSize(0.71f);
	data.SetScannerType(L"M10");
	data.SetScannerID(L"main");
	DateTime dt = DateTime::Parse(L"2005-12-25T11:12:13Z");
	data.SetScanDate(dt);
	data.SetRows(rows);
	data.SetCols(cols);
	std::string arrayId = "smellsliketeenspirit";
	data.SetArrayId(arrayId);
	
	ParameterNameValueType nvt;
	nvt.SetName(L"fixedlen");
	nvt.SetValueText(L"twenty-five", 25);
	data.GetFileHeader()->GetGenericDataHdr()->AddNameValParam(nvt);

	DATFileWriter writer(data);

	Uint16Vector pixels;
	pixels.reserve(pixelCount);

	u_int16_t inten = 10;
	for (int32_t i=0; i<pixelCount; ++i, ++inten)
	{
		if (inten > 46000)
		{
			inten = 0;
		}
		pixels.push_back(inten);
	}

	Uint16Vector stats;
	stats.push_back(pixels[0]);	// min
	stats.push_back(pixels[19]);	// max

	writer.WriteStats(stats);
	writer.WritePixels(pixels);
}
