////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
#include "calvin_files/writers/test/CalvinCHPMultiDataFileBufferWriterTest.h"
//
#include "calvin_files/parsers/src/CHPMultiDataFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileBufferWriter.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
//
#include <stdio.h>

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;

CPPUNIT_TEST_SUITE_REGISTRATION( CalvinCHPMultiDataFileBufferWriterTest );

#ifdef _MSC_VER
#pragma warning(disable: 4996) // ignore deprecated functions warning
#endif

#define TEST_FILE1 "multi_data_file1"
#define TEST_FILE2a "multi_data_file2a"
#define TEST_FILE2b "multi_data_file2b"
#define TEST_FILE3a "multi_data_file3a"
#define TEST_FILE3b "multi_data_file3b"
#define DMET_TEST_FILE "dmet_multi_data"

static string IntToString(int i)
{
   char buf[64];
   sprintf(buf, "%d", i);
   return buf;
}

void CalvinCHPMultiDataFileBufferWriterTest::setUp()
{
}

void CalvinCHPMultiDataFileBufferWriterTest::tearDown()
{
}

void CalvinCHPMultiDataFileBufferWriterTest::CreateReferenceFile1()
{
	int ng = 10000;
	int nx = 5000;
	CHPMultiDataData data(TEST_FILE1);
	data.SetEntryCount(GenotypeMultiDataType, ng, 10);
	data.SetEntryCount(ExpressionMultiDataType, nx, 10);
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	ProbeSetMultiDataGenotypeData gn;
	ProbeSetMultiDataExpressionData ex;

	writer->SeekToDataSet(GenotypeMultiDataType);
	gn.call = 0;
	gn.confidence = 0.0f;
	for (int i=0; i<ng; i++)
	{
		gn.name = IntToString(i);
		writer->WriteEntry(gn);
	}
	writer->SeekToDataSet(ExpressionMultiDataType);
	ex.quantification = 0;
	for (int i=0; i<nx; i++)
	{
		ex.name = IntToString(i);
		writer->WriteEntry(ex);
	}

	delete writer;
}

void CalvinCHPMultiDataFileBufferWriterTest::testBuffer1()
{
	CreateReferenceFile1();
	UpdateFile1();

	ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;
	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	reader.SetFilename(TEST_FILE1);
	reader.Read(data);

	int ng = 10000;
	int nx = 5000;

	CPPUNIT_ASSERT(data.GetEntryCount(ExpressionMultiDataType) == nx);
	CPPUNIT_ASSERT(data.GetEntryCount(GenotypeMultiDataType) == ng);

	for (int i=0; i<nx; i++)
	{
		data.GetExpressionEntry(ExpressionMultiDataType, i, ex);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, i, 0.0001f);
		CPPUNIT_ASSERT(ex.name == IntToString(i));
	}
	for (int i=0; i<ng; i++)
	{
		data.GetGenotypeEntry(GenotypeMultiDataType, i, gn);
		CPPUNIT_ASSERT(gn.call == i%4);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, i, 0.0001f);
		CPPUNIT_ASSERT(gn.name == IntToString(i));
	}
}

void CalvinCHPMultiDataFileBufferWriterTest::UpdateFile1()
{
	CHPMultiDataFileBufferWriter writer;
	vector<string> chps;
	vector<MultiDataType> dataTypes;
	map<MultiDataType, int> maxNameLengths;
	dataTypes.push_back(GenotypeMultiDataType);
	dataTypes.push_back(ExpressionMultiDataType);
	chps.push_back(TEST_FILE1);
	maxNameLengths[GenotypeMultiDataType] = 10;
	maxNameLengths[ExpressionMultiDataType] = 10;
	writer.Initialize(&chps, dataTypes, maxNameLengths);
	writer.SetMaxBufferSize(102400);
	ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;

	int ng = 10000;
	int nx = 5000;

	for (int i=0; i<nx; i++)
	{
		ex.quantification = (float)i;
		ex.name = IntToString(i);
		writer.WriteMultiDataExpressionEntry(ExpressionMultiDataType, 0, ex);
	}
	for (int i=0; i<ng; i++)
	{
		gn.name = IntToString(i);
		gn.call = i%4;
		gn.confidence = (float) i;
		writer.WriteMultiDataGenotypeEntry(GenotypeMultiDataType, 0, gn);
	}
	writer.FlushBuffer();
}

void CalvinCHPMultiDataFileBufferWriterTest::CreateReferenceFile2(const char *fileName, int offset)
{
	vector<ColumnInfo> cols;
	ParameterNameValueType nv;

	IntColumn icol(L"int");
	cols.push_back(icol);

	FloatColumn fcol(L"float");
	cols.push_back(fcol);

	int ng = 10000;
	int nx = 5000;
	CHPMultiDataData data(fileName);
	data.SetEntryCount(ExpressionMultiDataType, nx, 10, cols);
	data.SetEntryCount(GenotypeMultiDataType, ng, 10, cols);
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	ProbeSetMultiDataGenotypeData gn;
	ProbeSetMultiDataExpressionData ex;

	nv.SetName(L"int");
	nv.SetValueInt32(0);
	ex.metrics.push_back(nv);
	gn.metrics.push_back(nv);
	nv.SetName(L"float");
	nv.SetValueFloat(0);
	ex.metrics.push_back(nv);
	gn.metrics.push_back(nv);

	writer->SeekToDataSet(GenotypeMultiDataType);
	gn.call = 0;
	gn.confidence = 0.0f;
	for (int i=0; i<ng; i++)
	{
		gn.name = IntToString(i+offset);
		writer->WriteEntry(gn);
	}

	writer->SeekToDataSet(ExpressionMultiDataType);
	ex.quantification = 0;
	for (int i=0; i<nx; i++)
	{
		ex.name = IntToString(i+offset);
		writer->WriteEntry(ex);
	}
	delete writer;
}

void CalvinCHPMultiDataFileBufferWriterTest::testBuffer2()
{
	CreateReferenceFile2(TEST_FILE2a, 0);
	CreateReferenceFile2(TEST_FILE2b, 1001);
	vector<string> chps;
	vector<int> offset;
	chps.push_back(TEST_FILE2a);
	offset.push_back(0);
	chps.push_back(TEST_FILE2b);
	offset.push_back(1001);
	UpdateFile2(chps, offset);
	CheckFile2(TEST_FILE2a, 0);
	CheckFile2(TEST_FILE2b, 1001);
}

void CalvinCHPMultiDataFileBufferWriterTest::CheckFile2(const char *fileName, int offset)
{
	ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;
	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	reader.SetFilename(fileName);
	reader.Read(data);

	int ng = 10000;
	int nx = 5000;

	CPPUNIT_ASSERT(data.GetEntryCount(ExpressionMultiDataType) == nx);
	CPPUNIT_ASSERT(data.GetEntryCount(GenotypeMultiDataType) == ng);

	for (int i=0; i<nx; i++)
	{
		data.GetExpressionEntry(ExpressionMultiDataType, i, ex);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, i+offset, 0.0001f);
		CPPUNIT_ASSERT(ex.name == IntToString(i+offset));

		CPPUNIT_ASSERT(ex.metrics[0].GetValueInt32() == i+offset);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
	}
	for (int i=0; i<ng; i++)
	{
		data.GetGenotypeEntry(GenotypeMultiDataType, i, gn);
		CPPUNIT_ASSERT(gn.call == (i+offset)%4);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, i+offset, 0.0001f);
		CPPUNIT_ASSERT(gn.name == IntToString(i+offset));

		CPPUNIT_ASSERT(gn.metrics[0].GetValueInt32() == i+offset);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
	}
}

void CalvinCHPMultiDataFileBufferWriterTest::UpdateFile2(std::vector<std::string> &fileNames, vector<int> offset)
{
	ParameterNameValueType nv;
	CHPMultiDataFileBufferWriter writer;
	vector<MultiDataType> dataTypes;
	dataTypes.push_back(ExpressionMultiDataType);
	dataTypes.push_back(GenotypeMultiDataType);
	map<MultiDataType, int> maxNameLengths;
	maxNameLengths[ExpressionMultiDataType]=10;
	maxNameLengths[GenotypeMultiDataType]=10;
	writer.Initialize(&fileNames, dataTypes, maxNameLengths);
	writer.SetMaxBufferSize(102400);

	int ng = 10000;
	int nx = 5000;

	ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;

	nv.SetName(L"int");
	nv.SetValueInt32(0);
	ex.metrics.push_back(nv);
	gn.metrics.push_back(nv);
	nv.SetName(L"float");
	nv.SetValueFloat(0);
	ex.metrics.push_back(nv);
	gn.metrics.push_back(nv);
	for (int i=0; i<nx; i++)
	{
		for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
		{
			ex.name = IntToString(i+offset[ifile]);
			ex.quantification = (float)(i+offset[ifile]);
			ex.metrics[0].SetValueInt32(i+offset[ifile]);
			ex.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
			writer.WriteMultiDataExpressionEntry(ExpressionMultiDataType, ifile, ex);
		}
	}

	for (int i=0; i<ng; i++)
	{
		for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
		{
			gn.name = IntToString(i+offset[ifile]);
			gn.call = (i+offset[ifile])%4;
			gn.confidence = (float)(i+offset[ifile]);
			gn.metrics[0].SetValueInt32(i+offset[ifile]);
			gn.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
			writer.WriteMultiDataGenotypeEntry(GenotypeMultiDataType, ifile, gn);
		}
	}
	writer.FlushBuffer();
}

void CalvinCHPMultiDataFileBufferWriterTest::CreateReferenceFile3(const char *fileName, int offset)
{
	vector<ColumnInfo> cols;
	ParameterNameValueType nv;

	IntColumn icol(L"int");
	cols.push_back(icol);

	FloatColumn fcol(L"float");
	cols.push_back(fcol);

	int ng = 10000;
	int nc = 10000;
	int nx = 5000;
	CHPMultiDataData data(fileName);
	data.SetEntryCount(ExpressionMultiDataType, nx, 10, cols);
	data.SetEntryCount(GenotypeMultiDataType, ng, 10, cols);
	data.SetEntryCount(CopyNumberMultiDataType, nc, 10, cols);
	data.SetEntryCount(CytoMultiDataType, nc, 10);
	data.SetEntryCount(CopyNumberVariationMultiDataType, nx, 10, cols);

	DataSetHeader *dsh = data.GetDataSetHeader(CopyNumberMultiDataType);
	for (int i=0; i<5; ++i)
	{
		ParameterNameValueType param;
		param.SetName(StringUtils::ConvertMBSToWCS(IntToString(i)));
		param.SetValueAscii(IntToString(i));
		dsh->AddNameValParam(param);
	}

	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	ProbeSetMultiDataGenotypeData gn;
	ProbeSetMultiDataCopyNumberData cn;
	ProbeSetMultiDataCytoRegionData cy;
	ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataCopyNumberVariationRegionData cv;

	nv.SetName(L"int");
	nv.SetValueInt32(0);
	ex.metrics.push_back(nv);
	gn.metrics.push_back(nv);
	cv.metrics.push_back(nv);
	nv.SetName(L"float");
	nv.SetValueFloat(0);
	ex.metrics.push_back(nv);
	gn.metrics.push_back(nv);
	cv.metrics.push_back(nv);

	writer->SeekToDataSet(GenotypeMultiDataType);
	gn.call = 0;
	gn.confidence = 0.0f;
	writer->SeekToDataSet(CopyNumberMultiDataType);
	cn.chr = 0;
	cn.position = 0;
	writer->SeekToDataSet(CytoMultiDataType);
	cy.call = 0;
	cy.confidenceScore = 0;
	cy.chr = 0;
	cy.startPosition = 0;
	cy.stopPosition = 0;
	writer->SeekToDataSet(ExpressionMultiDataType);
	ex.quantification = 0;
	for (int i=0; i<nx; i++)
	{
		ex.name = IntToString(i+offset);
		writer->WriteEntry(ex);
	}

	writer->SeekToDataSet(CopyNumberVariationMultiDataType);
	cv.call = 0;
	cv.confidenceScore = 0.0f;
	cv.signal = 0.0f;    
	delete writer;
}

void CalvinCHPMultiDataFileBufferWriterTest::testBuffer3()
{
	CreateReferenceFile3(TEST_FILE3a, 0);
	CreateReferenceFile3(TEST_FILE3b, 1001);
	vector<string> chps;
	vector<int> offset;
	chps.push_back(TEST_FILE3a);
	offset.push_back(0);
	chps.push_back(TEST_FILE3b);
	offset.push_back(1001);
	UpdateFile3(chps, offset);
	CheckFile3(TEST_FILE3a, 0);
	CheckFile3(TEST_FILE3b, 1001);
}

void CalvinCHPMultiDataFileBufferWriterTest::CheckFile3(const char *fileName, int offset)
{
	ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;
	ProbeSetMultiDataCopyNumberData cn;
	ProbeSetMultiDataCytoRegionData cy;
	ProbeSetMultiDataCopyNumberVariationRegionData cv;
	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	reader.SetFilename(fileName);
	reader.Read(data);

	int ng = 10000;
	int nc = 10000;
	int nx = 5000;

	CPPUNIT_ASSERT(data.GetEntryCount(ExpressionMultiDataType) == nx);
	CPPUNIT_ASSERT(data.GetEntryCount(GenotypeMultiDataType) == ng);
	CPPUNIT_ASSERT(data.GetEntryCount(CopyNumberMultiDataType) == nc);
	CPPUNIT_ASSERT(data.GetEntryCount(CytoMultiDataType) == nc);
	CPPUNIT_ASSERT(data.GetEntryCount(CopyNumberVariationMultiDataType) == nx);

	for (int i=0; i<nx; i++)
	{
		data.GetExpressionEntry(ExpressionMultiDataType, i, ex);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, i+offset, 0.0001f);
		CPPUNIT_ASSERT(ex.name == IntToString(i+offset));

		CPPUNIT_ASSERT(ex.metrics[0].GetValueInt32() == i+offset);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
	}
	for (int i=0; i<ng; i++)
	{
		data.GetGenotypeEntry(GenotypeMultiDataType, i, gn);
		CPPUNIT_ASSERT(gn.call == (i+offset)%4);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, i+offset, 0.0001f);
		CPPUNIT_ASSERT(gn.name == IntToString(i+offset));

		CPPUNIT_ASSERT(gn.metrics[0].GetValueInt32() == i+offset);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
	}
	for (int i=0; i<nc; i++)
	{
		data.GetCopyNumberEntry(CopyNumberMultiDataType, i, cn);
		CPPUNIT_ASSERT(cn.chr == (i+offset)%4);
		CPPUNIT_ASSERT(cn.position == i+offset);
		CPPUNIT_ASSERT(cn.name == IntToString(i+offset));

		CPPUNIT_ASSERT(cn.metrics[0].GetValueInt32() == i+offset);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(cn.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
	}
	for (int i=0; i<nc; i++)
	{
		data.GetCytoEntry(CytoMultiDataType, i, cy);
		CPPUNIT_ASSERT(cy.call == (i+offset)%4);
		CPPUNIT_ASSERT(cy.chr == (i+offset)%10);
		CPPUNIT_ASSERT(cy.startPosition == i);
		CPPUNIT_ASSERT(cy.stopPosition == i+1);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(cy.confidenceScore, i+offset, 0.00001f);
		CPPUNIT_ASSERT(cy.name == IntToString(i+offset));
	}
	for (int i=0; i<nx; i++)
	{
		data.GetCopyNumberVariationEntry(CopyNumberVariationMultiDataType, i, cv);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(cv.confidenceScore, i+offset, 0.00001f);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(cv.signal, i+offset, 0.00001f);
		CPPUNIT_ASSERT(cv.call == (i+offset)%4);
		CPPUNIT_ASSERT(cv.name == IntToString(i+offset));
	}
}

void CalvinCHPMultiDataFileBufferWriterTest::UpdateFile3(std::vector<std::string> &fileNames, vector<int> offset)
{
	ParameterNameValueType nv;
	CHPMultiDataFileBufferWriter writer;
	vector<MultiDataType> dataTypes;
	dataTypes.push_back(ExpressionMultiDataType);
	dataTypes.push_back(GenotypeMultiDataType);
	dataTypes.push_back(CopyNumberMultiDataType);
	dataTypes.push_back(CytoMultiDataType);
	dataTypes.push_back(CopyNumberVariationMultiDataType);
	map<MultiDataType, int> maxNameLengths;
	maxNameLengths[ExpressionMultiDataType]=10;
	maxNameLengths[GenotypeMultiDataType]=10;
	maxNameLengths[CopyNumberMultiDataType]=10;
	maxNameLengths[CytoMultiDataType]=10;
	maxNameLengths[CopyNumberVariationMultiDataType]=10;
	writer.Initialize(&fileNames, dataTypes, maxNameLengths);
	writer.SetMaxBufferSize(102400);

	int ng = 10000;
	int nc = 10000;
	int nx = 5000;

	ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;
	ProbeSetMultiDataCopyNumberData cn;
	ProbeSetMultiDataCytoRegionData cy;
	ProbeSetMultiDataCopyNumberVariationRegionData cv;

	nv.SetName(L"int");
	nv.SetValueInt32(0);
	ex.metrics.push_back(nv);
	gn.metrics.push_back(nv);
	cn.metrics.push_back(nv);
	cv.metrics.push_back(nv);
	nv.SetName(L"float");
	nv.SetValueFloat(0);
	ex.metrics.push_back(nv);
	gn.metrics.push_back(nv);
	cn.metrics.push_back(nv);
	cv.metrics.push_back(nv);

	for (int i=0; i<nx; i++)
	{
		for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
		{
			ex.name = IntToString(i+offset[ifile]);
			ex.quantification = (float)(i+offset[ifile]);
			ex.metrics[0].SetValueInt32(i+offset[ifile]);
			ex.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
			writer.WriteMultiDataExpressionEntry(ExpressionMultiDataType, ifile, ex);
		}
	}

	for (int i=0; i<ng; i++)
	{
		for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
		{
			gn.name = IntToString(i+offset[ifile]);
			gn.call = (i+offset[ifile])%4;
			gn.confidence = (float)(i+offset[ifile]);
			gn.metrics[0].SetValueInt32(i+offset[ifile]);
			gn.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
			writer.WriteMultiDataGenotypeEntry(GenotypeMultiDataType, ifile, gn);
		}
	}
	for (int i=0; i<nc; i++)
	{
		for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
		{
			cn.name = IntToString(i+offset[ifile]);
			cn.chr = ((i+offset[ifile])%4);
			cn.position = (i+offset[ifile]);
			cn.metrics[0].SetValueInt32(i+offset[ifile]);
			cn.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
			writer.WriteMultiDataCopyNumberEntry(CopyNumberMultiDataType, ifile, cn);
		}
	}
	for (int i=0; i<nc; i++)
	{
		for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
		{
			cy.name = IntToString(i+offset[ifile]);
			cy.call = ((i+offset[ifile])%4);
			cy.chr = ((i+offset[ifile])%10);
			cy.startPosition = i;
			cy.stopPosition = i+1;
			cy.confidenceScore = (float)(i+offset[ifile]);
			writer.WriteMultiDataCytoRegionEntry(CytoMultiDataType, ifile, cy);
		}
	}

	for (int i=0; i<nx; i++)
	{
		for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
		{
			cv.name = IntToString(i+offset[ifile]);
			cv.signal = (float)(i+offset[ifile]);
			cv.call = ((i+offset[ifile])%4);
			cv.confidenceScore = (float)(i+offset[ifile]);
			cv.metrics[0].SetValueInt32(i+offset[ifile]);
			cv.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
			writer.WriteMultiDataCopyNumberVariationRegionEntry(CopyNumberVariationMultiDataType,
				ifile, cv);
		}
	}
	writer.FlushBuffer(); 
}

void CalvinCHPMultiDataFileBufferWriterTest::testDMET()
{
	CreateDmetReferenceFile();
	UpdateDmetFile();

	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	reader.SetFilename(DMET_TEST_FILE);
	reader.Read(data);

	int gRows = 100;
	int cRows = 50;
	int mRows = 20;

	CPPUNIT_ASSERT(data.GetEntryCount(DmetBiAllelicMultiDataType) == gRows);
	CPPUNIT_ASSERT(data.GetEntryCount(DmetCopyNumberMultiDataType) == cRows);
	CPPUNIT_ASSERT(data.GetEntryCount(DmetMultiAllelicMultiDataType) == mRows);

	DmetBiAllelicData gData;
	for (int i = 0; i < gRows; i++)
	{
		data.GetEntry(DmetBiAllelicMultiDataType, i, gData);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(gData.confidence, i, 0.0001f);
		CPPUNIT_ASSERT(gData.name == IntToString(i));
	}
	DmetCopyNumberData cData;
	for (int i = 0; i < cRows; i++)
	{
		data.GetEntry(DmetCopyNumberMultiDataType, i, cData);
		CPPUNIT_ASSERT(cData.call == i%4);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(cData.confidence, i, 0.0001f);
		CPPUNIT_ASSERT(cData.name == IntToString(i));
	}
	DmetMultiAllelicData mData;
	for (int i = 0; i < mRows; i++)
	{
		data.GetEntry(DmetMultiAllelicMultiDataType, i, mData);
		CPPUNIT_ASSERT(mData.call == i%4);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(mData.confidence, i, 0.0001f);
		CPPUNIT_ASSERT(mData.name == IntToString(i));
	}
    
}

void CalvinCHPMultiDataFileBufferWriterTest::CreateDmetReferenceFile()
{
	int gRows = 100;
	int cRows = 50;
	int mRows = 20;
	CHPMultiDataData data(DMET_TEST_FILE);
	data.SetEntryCount(DmetBiAllelicMultiDataType, gRows, 10);
	data.SetEntryCount(DmetCopyNumberMultiDataType, cRows, 10);
	data.SetEntryCount(DmetMultiAllelicMultiDataType, mRows, 10);
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	DmetBiAllelicData gData;
	writer->SeekToDataSet(DmetBiAllelicMultiDataType);
	gData.call = 22;
	gData.confidence = 1.1f;
	gData.force = 22;
	gData.signalA = 22;
	gData.signalB = 22;
	gData.contextA = 22;
	gData.contextB = 22;
	for (int i = 0; i < gRows; i++)
	{
		gData.name = IntToString(i);
		writer->WriteEntry(gData);
	}
	DmetCopyNumberData cData;
	writer->SeekToDataSet(DmetCopyNumberMultiDataType);
	cData.call = 33;
	cData.confidence = 33.3f;
	cData.force = 99;
	cData.estimate = 66.6f;
	cData.lower = 44.4f;
	cData.upper = 55.5f;
	for (int i = 0; i < cRows; i++)
	{
		cData.name = IntToString(i);
		writer->WriteEntry(cData);
	}
	DmetMultiAllelicData mData;
	writer->SeekToDataSet(DmetMultiAllelicMultiDataType);
	mData.call = 22;
	mData.confidence = 77.7f;
	mData.force = 55;
	mData.alleleCount = 0;
	mData.signalA = 0;
	mData.signalB = 0;
	mData.signalC = 0;
	mData.signalD = 0;
	mData.signalE = 0;
	mData.signalF = 0;
	mData.contextA = 0;
	mData.contextB = 0;
	mData.contextC = 0;
	mData.contextD = 0;
	mData.contextE = 0;
	mData.contextF = 0;
	for (int i = 0; i < mRows; i++)
	{
		mData.name = IntToString(i);
		writer->WriteEntry(mData);
	}
	delete writer;
}

void CalvinCHPMultiDataFileBufferWriterTest::UpdateDmetFile()
{
	CHPMultiDataFileBufferWriter writer;
	vector<string> chps;
	vector<MultiDataType> dataTypes;
	map<MultiDataType, int> maxNameLengths;
	dataTypes.push_back(DmetBiAllelicMultiDataType);
	dataTypes.push_back(DmetCopyNumberMultiDataType);
	dataTypes.push_back(DmetMultiAllelicMultiDataType);
	chps.push_back(DMET_TEST_FILE);
	maxNameLengths[DmetBiAllelicMultiDataType] = 10;
	maxNameLengths[DmetCopyNumberMultiDataType] = 10;
	maxNameLengths[DmetMultiAllelicMultiDataType] = 10;
	writer.Initialize(&chps, dataTypes, maxNameLengths);
	writer.SetMaxBufferSize(102400);

	DmetBiAllelicData gData;
	int gRows = 100;
	for (int i = 0; i < gRows; i++)
	{
		gData.call = 22;
		gData.confidence = (float)i;
		gData.force = 22;
		gData.signalA = 22;
		gData.signalB = 22;
		gData.contextA = 22;
		gData.contextB = 22;
		gData.name = IntToString(i);
		writer.WriteEntry(DmetBiAllelicMultiDataType, 0, gData);
	}
	DmetCopyNumberData cData;
	int cRows = 50;
	for (int i = 0; i < cRows; i++)
	{
		cData.name = IntToString(i);
		cData.call = i%4;
		cData.confidence = (float)i;
		cData.force = 99;
		cData.estimate = 66.6f;
		cData.lower = 44.4f;
		cData.upper = 55.5f;
		writer.WriteEntry(DmetCopyNumberMultiDataType, 0, cData);
	}

	DmetMultiAllelicData mData;
	int mRows = 20;
	for (int i = 0; i < mRows; i++)
	{
		mData.name = IntToString(i);
		mData.call = i%4;
		mData.confidence = (float)i;
		mData.force = 55;
		mData.alleleCount = 0;
		mData.signalA = 0;
		mData.signalB = 0;
		mData.signalC = 0;
		mData.signalD = 0;
		mData.signalE = 0;
		mData.signalF = 0;
		mData.contextA = 0;
		mData.contextB = 0;
		mData.contextC = 0;
		mData.contextD = 0;
		mData.contextE = 0;
		mData.contextF = 0;
		writer.WriteEntry(DmetMultiAllelicMultiDataType, 0, mData);
	}
	writer.FlushBuffer();
}

static void CreateChrSumFile(const char *file)
{
    std::list<std::wstring> groupNames;
    groupNames.push_back(L"g1");
	CHPMultiDataData data(file, &groupNames);
	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(ChromosomeSummaryMultiDataType, 3, 10, L"g1");
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
    delete writer;
}

void CalvinCHPMultiDataFileBufferWriterTest::WriteTestChrSummary()
{
    CreateChrSumFile("chrsum1");
    CreateChrSumFile("chrsum2");
    affymetrix_calvin_data::ChromosomeMultiDataSummaryData p;

    CHPMultiDataFileBufferWriter writer;
	vector<string> chps;
	vector<MultiDataType> dataTypes;

    chps.push_back("chrsum1");
	chps.push_back("chrsum2");
	dataTypes.push_back(ChromosomeSummaryMultiDataType);

    std::map<MultiDataType, int> maxNameLn;
    maxNameLn[ChromosomeSummaryMultiDataType] = 10;

    writer.Initialize(&chps, dataTypes, maxNameLn);
    float fval = 1.0f;
	//char buf[16];
    for (int i=0; i<3; i++)
    {
        p.chr=i;
		p.display=IntToString(i);
		p.startIndex=(u_int32_t)fval++;
		p.markerCount=(u_int32_t)fval++;
        p.minSignal=fval++;
        p.maxSignal=fval++;
        p.medianCnState=fval++;
        p.homFrequency=fval++;
        p.hetFrequency=fval++;
        //p.mosaicism=fval++;
        writer.WriteEntry(ChromosomeSummaryMultiDataType, 0, p);
        p.chr=i+1;
		p.display=IntToString(i+1);
		p.startIndex=(u_int32_t)fval++;
		p.markerCount=(u_int32_t)fval++;
        p.minSignal=fval++;
        p.maxSignal=fval++;
        p.medianCnState=fval++;
        p.homFrequency=fval++;
        p.hetFrequency=fval++;
        //p.mosaicism=fval++;
        writer.WriteEntry(ChromosomeSummaryMultiDataType, 1, p);
    }
    writer.FlushBuffer();



    fval = 1.0f;
	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("chrsum2");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(ChromosomeSummaryMultiDataType) == 3);

    for (int i=0; i<3; i++)
    {
        fval += 7;
        data2.GetChromosomeSummaryEntry(ChromosomeSummaryMultiDataType, i, p);
   	    CPPUNIT_ASSERT(p.chr == i+1);
   	    CPPUNIT_ASSERT(strcmp(p.display.c_str(), IntToString(i+1).c_str()) == 0);
		CPPUNIT_ASSERT(p.startIndex == (u_int32_t)fval++);
		CPPUNIT_ASSERT(p.markerCount == (u_int32_t)fval++);
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.minSignal, fval++, 0.0001f);
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.maxSignal, fval++, 0.0001f);
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.medianCnState, fval++, 0.0001f);
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.homFrequency, fval++, 0.0001f);
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.hetFrequency, fval++, 0.0001f);
	    //CPPUNIT_ASSERT_DOUBLES_EQUAL(p.mosaicism, fval++, 0.0001f);
    }
}

static void CreateChrSegFileEx(const char *file)
{
    std::list<std::wstring> groupNames;
    groupNames.push_back(L"g1");
    groupNames.push_back(L"g2");

	CHPMultiDataData data(file, &groupNames);

	ParameterNameValueType param;

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");

    const MultiDataType dataTypes[] = {
        SegmentGenotypeConcordanceMultiDataType,
        SegmentGenotypeDiscordanceMultiDataType,
        SegmentCNLossLOHConcordanceMultiDataType,
        SegmentCNNeutralLOHConcordanceMultiDataType,
        SegmentHeteroUPDMultiDataType,
        SegmentIsoUPDMultiDataType,
        SegmentDenovoCopyNumberMultiDataType,
        SegmentHemizygousParentOfOriginMultiDataType
    };
    int nTypes = sizeof(dataTypes) / sizeof(MultiDataType);

    for (int itype=0; itype<nTypes; itype++)
    {
        if (dataTypes[itype] == SegmentGenotypeDiscordanceMultiDataType)
            data.SetEntryCount(dataTypes[itype], 3+itype, 10, cols, L"g1");
        else
            data.SetEntryCount(dataTypes[itype], 3+itype, 10, L"g2");
    }
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
    delete writer;
}

void CalvinCHPMultiDataFileBufferWriterTest::WriteTestChrSegmentEx()
{
    CreateChrSegFileEx("chrseg1ex");

    affymetrix_calvin_data::ChromosomeSegmentDataEx p;
	ParameterNameValueType param;
    CHPMultiDataFileBufferWriter writer;
	vector<string> chps;
	vector<MultiDataType> dataTypesV;
    std::map<MultiDataType, int> maxNameLn;

    const MultiDataType dataTypes[] = {
        SegmentGenotypeConcordanceMultiDataType,
        SegmentGenotypeDiscordanceMultiDataType,
        SegmentCNLossLOHConcordanceMultiDataType,
        SegmentCNNeutralLOHConcordanceMultiDataType,
        SegmentHeteroUPDMultiDataType,
        SegmentIsoUPDMultiDataType,
        SegmentDenovoCopyNumberMultiDataType,
        SegmentHemizygousParentOfOriginMultiDataType
    };
    int nTypes = sizeof(dataTypes) / sizeof(MultiDataType);

    chps.push_back("chrseg1ex");
    for (int i=0; i<nTypes; i++)
    {
        dataTypesV.push_back(dataTypes[i]);
        maxNameLn[dataTypes[i]] = 10;
    }
    writer.Initialize(&chps, dataTypesV, maxNameLn);
    int ival=0;
    //char buf[10];
    for (int itype=0; itype<nTypes; itype++)
    {
        for (int i=0; i<3+itype; i++)
        {
			p.referenceSampleKey = ival++;
			p.familialSampleKey = ival++;
            p.chr = ival % 200;
            ++ival;
            p.call = ival % 200;
            ++ival;
            p.confidence = ival++;
            p.startPosition = ival++;
            p.stopPosition = ival++;
            p.markerCount = ival++;
			p.homozygosity = ival++;
			p.heterozygosity = ival++;
            p.segmentId = ival++;
            p.metrics.clear();
	        if (dataTypes[itype] == SegmentGenotypeDiscordanceMultiDataType)
            {
                param.SetValueInt32(ival++);
                p.metrics.push_back(param);
                param.SetValueFloat(ival++);
                p.metrics.push_back(param);
            }
            writer.WriteEntry(dataTypes[itype], 0, p);
            p.metrics.clear();
        }
    }
    writer.FlushBuffer();



	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("chrseg1ex");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");

    ival = 0;
    for (int itype=0; itype<nTypes; itype++)
    {
    	CPPUNIT_ASSERT(data2.GetEntryCount(dataTypes[itype]) == 3+itype);
    
        if (dataTypes[itype] == SegmentGenotypeDiscordanceMultiDataType)
        {
            CPPUNIT_ASSERT(data2.GetMetricColumnName(SegmentGenotypeDiscordanceMultiDataType, 0) == L"int");
            CPPUNIT_ASSERT(data2.GetMetricColumnName(SegmentGenotypeDiscordanceMultiDataType, 1) == L"float");
        }
        for (int i=0; i<3+itype; i++)
        {
            data2.GetChromosomeSegmentEntry(dataTypes[itype], i, p);
			CPPUNIT_ASSERT(p.referenceSampleKey == ival++);
			CPPUNIT_ASSERT(p.familialSampleKey == ival++);
   	        CPPUNIT_ASSERT(p.chr == ival%200);
            ++ival;
	        CPPUNIT_ASSERT(p.call == ival%200);
            ++ival;
            CPPUNIT_ASSERT_DOUBLES_EQUAL(p.confidence, ival++, 0.0001f);
	        CPPUNIT_ASSERT(p.startPosition == ival++);
            CPPUNIT_ASSERT(p.stopPosition == ival++);
            CPPUNIT_ASSERT(p.markerCount == ival++);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(p.homozygosity, ival++, 0.0001f);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(p.heterozygosity, ival++, 0.0001f);
            CPPUNIT_ASSERT(p.segmentId == ival++);
            if (dataTypes[itype] == SegmentGenotypeDiscordanceMultiDataType)
            {
                param = p.metrics[0];
                CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
                CPPUNIT_ASSERT(param.GetValueInt32() == ival++);
                param = p.metrics[1];
                CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), ival++, 0.00001f);
            }
        }
    }

}

static void CreateChrSegFile(const char *file)
{
    std::list<std::wstring> groupNames;
    groupNames.push_back(L"g1");
    groupNames.push_back(L"g2");

	CHPMultiDataData data(file, &groupNames);

	ParameterNameValueType param;

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");

    const MultiDataType dataTypes[] = {
        SegmentCNMultiDataType,
        SegmentLOHMultiDataType,
        SegmentCNNeutralLOHMultiDataType,
        SegmentNormalDiploidMultiDataType,
        SegmentMosaicismMultiDataType,
        SegmentNoCallMultiDataType
    };
    int nTypes = sizeof(dataTypes) / sizeof(MultiDataType);

    for (int itype=0; itype<nTypes; itype++)
    {
        if (dataTypes[itype] == SegmentNormalDiploidMultiDataType)
            data.SetEntryCount(dataTypes[itype], 3+itype, 10, cols, L"g1");
        else
            data.SetEntryCount(dataTypes[itype], 3+itype, 10, L"g2");
    }
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
    delete writer;
}

void CalvinCHPMultiDataFileBufferWriterTest::WriteTestChrSegment()
{
    CreateChrSegFile("chrseg1");

    affymetrix_calvin_data::ChromosomeSegmentData p;
	ParameterNameValueType param;
    CHPMultiDataFileBufferWriter writer;
	vector<string> chps;
	vector<MultiDataType> dataTypesV;
    std::map<MultiDataType, int> maxNameLn;

    const MultiDataType dataTypes[] = {
        SegmentCNMultiDataType,
        SegmentLOHMultiDataType,
        SegmentCNNeutralLOHMultiDataType,
        SegmentNormalDiploidMultiDataType,
        SegmentMosaicismMultiDataType,
        SegmentNoCallMultiDataType
    };
    int nTypes = sizeof(dataTypes) / sizeof(MultiDataType);

    chps.push_back("chrseg1");
    for (int i=0; i<nTypes; i++)
    {
        dataTypesV.push_back(dataTypes[i]);
        maxNameLn[dataTypes[i]] = 10;
    }
    writer.Initialize(&chps, dataTypesV, maxNameLn);
    int ival=0;
    //char buf[10];
    for (int itype=0; itype<nTypes; itype++)
    {
        for (int i=0; i<3+itype; i++)
        {
            p.chr = ival % 200;
            ++ival;
            p.startPosition = ival++;
            p.stopPosition = ival++;
            p.markerCount = ival++;
			p.meanMarkerDistance = ival++;
            p.segmentId = ival++;
            p.metrics.clear();
	        if (dataTypes[itype] == SegmentNormalDiploidMultiDataType)
            {
                param.SetValueInt32(ival++);
                p.metrics.push_back(param);
                param.SetValueFloat(ival++);
                p.metrics.push_back(param);
            }
            writer.WriteEntry(dataTypes[itype], 0, p);
            p.metrics.clear();
        }
    }
    writer.FlushBuffer();



	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("chrseg1");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");

    ival = 0;
    for (int itype=0; itype<nTypes; itype++)
    {
    	CPPUNIT_ASSERT(data2.GetEntryCount(dataTypes[itype]) == 3+itype);
    
        if (dataTypes[itype] == SegmentNormalDiploidMultiDataType)
        {
            CPPUNIT_ASSERT(data2.GetMetricColumnName(SegmentNormalDiploidMultiDataType, 0) == L"int");
            CPPUNIT_ASSERT(data2.GetMetricColumnName(SegmentNormalDiploidMultiDataType, 1) == L"float");
        }
        for (int i=0; i<3+itype; i++)
        {
            data2.GetChromosomeSegmentEntry(dataTypes[itype], i, p);
   	        CPPUNIT_ASSERT(p.chr == ival%200);
            ++ival;
	        CPPUNIT_ASSERT(p.startPosition == ival++);
            CPPUNIT_ASSERT(p.stopPosition == ival++);
            CPPUNIT_ASSERT(p.markerCount == ival++);
			CPPUNIT_ASSERT(p.meanMarkerDistance == ival++);
            CPPUNIT_ASSERT(p.segmentId == ival++);
            if (dataTypes[itype] == SegmentNormalDiploidMultiDataType)
            {
                param = p.metrics[0];
                CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
                CPPUNIT_ASSERT(param.GetValueInt32() == ival++);
                param = p.metrics[1];
                CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), ival++, 0.00001f);
            }
        }
    }

}

static void CreateFamilialFile(const char *file)
{
    affymetrix_calvin_data::FamilialSample s;
    affymetrix_calvin_data::FamilialSegmentOverlap o;

    std::list<std::wstring> groupNames;
    groupNames.push_back(L"g1");
    groupNames.push_back(L"g2");
	CHPMultiDataData data(file, &groupNames);

	ParameterNameValueType param;

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");

    data.SetEntryCount(FamilialSegmentOverlapsMultiDataType, 3, 10, 10, 10, L"g2");
    data.SetEntryCount(FamilialSamplesMultiDataType, 3, 10, 10, 10, 10, L"g1");

    CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
/*
	writer->SeekToDataSet(FamilialSamplesMultiDataType);
    int ival = 0;
    char buf[10];
    for (int i=0; i<3; i++)
    {
        s.sampleKey = ival++;
        sprintf(buf,"%i",ival++);
	    s.arrID = buf;
        sprintf(buf,"%i",ival++);
	    s.chpID = buf;
        //wsprintf(wbuf,L"%i",ival++);
	    s.chpFilename = L"wbuf";
        sprintf(buf,"%i",ival++);
	    s.role = buf;
        s.roleValidity = (ival%2);
        ++ival;
	    s.roleConfidence = ival++;
        writer->WriteEntry(s);
    }
	writer->SeekToDataSet(FamilialSegmentOverlapsMultiDataType);
    for (int i=0; i<3; i++)
    {
        sprintf(buf,"%i",ival++);
        o.segmentType = buf;
        o.referenceSampleKey = ival++;
        sprintf(buf,"%i",ival++);
        o.referenceSegmentID = buf;
        o.familialSampleKey = ival++;
        sprintf(buf,"%i",ival++);
        o.familialSegmentID = buf;
        writer->WriteEntry(o);
    }
*/
    delete writer;
}

void CalvinCHPMultiDataFileBufferWriterTest::WriteTestFamilial()
{
    CreateFamilialFile("familial1");
    CreateFamilialFile("familial2");

    CHPMultiDataFileBufferWriter writer;
	vector<string> chps;
	vector<MultiDataType> dataTypes;
    std::map<MultiDataType, int> maxSegmentTypeLn;
    std::map<MultiDataType, int> maxReferenceSampleKeyLn;
    std::map<MultiDataType, int> maxFamilialSegmentIDLn;	  
    std::map<MultiDataType, int> maxFamilialARRIDLn;
    std::map<MultiDataType, int> maxFamilialCHPIDLn;
    std::map<MultiDataType, int> maxFamilialCHPFilenameLn;
    std::map<MultiDataType, int> maxFamilialRoleLn;

    chps.push_back("familial1");
	chps.push_back("familial2");
	dataTypes.push_back(FamilialSegmentOverlapsMultiDataType);
	dataTypes.push_back(FamilialSamplesMultiDataType);

    maxSegmentTypeLn[FamilialSegmentOverlapsMultiDataType] = 10;
    maxReferenceSampleKeyLn[FamilialSegmentOverlapsMultiDataType] = 10;
    maxFamilialSegmentIDLn[FamilialSegmentOverlapsMultiDataType] = 10;	  
    maxFamilialARRIDLn[FamilialSegmentOverlapsMultiDataType] = 10;
    maxFamilialCHPIDLn[FamilialSegmentOverlapsMultiDataType] = 10;
    maxFamilialCHPFilenameLn[FamilialSegmentOverlapsMultiDataType] = 10;
    maxFamilialRoleLn[FamilialSegmentOverlapsMultiDataType] = 10;

    maxSegmentTypeLn[FamilialSamplesMultiDataType] = 10;
    maxReferenceSampleKeyLn[FamilialSamplesMultiDataType] = 10;
    maxFamilialSegmentIDLn[FamilialSamplesMultiDataType] = 10;	  
    maxFamilialARRIDLn[FamilialSamplesMultiDataType] = 10;
    maxFamilialCHPIDLn[FamilialSamplesMultiDataType] = 10;
    maxFamilialCHPFilenameLn[FamilialSamplesMultiDataType] = 10;
    maxFamilialRoleLn[FamilialSamplesMultiDataType] = 10;

    writer.Initialize(&chps, dataTypes, 
                    maxSegmentTypeLn,
                    maxReferenceSampleKeyLn,
                    maxFamilialSegmentIDLn,	  
                    maxFamilialARRIDLn,
                    maxFamilialCHPIDLn,
                    maxFamilialCHPFilenameLn,
                    maxFamilialRoleLn);

    affymetrix_calvin_data::FamilialSample s;
    affymetrix_calvin_data::FamilialSegmentOverlap o;

    char buf[10];
    int ival = 0;
    for (int i=0; i<3; i++)
    {
        s.sampleKey = ival++;
        sprintf(buf,"%i",ival++);
	    s.arrID = buf;
        sprintf(buf,"%i",ival++);
	    s.chpID = buf;
        //wsprintf(wbuf,L"%i",ival++);
	    s.chpFilename = L"wbuf";
        sprintf(buf,"%i",ival++);
	    s.role = buf;
        s.roleValidity = (ival%2);
        ++ival;
	    s.roleConfidence = ival++;
        writer.WriteEntry(FamilialSamplesMultiDataType, 0, s);
                
        s.sampleKey = ival++;
        sprintf(buf,"%i",ival++);
	    s.arrID = buf;
        sprintf(buf,"%i",ival++);
	    s.chpID = buf;
        //wsprintf(wbuf,L"%i",ival++);
	    s.chpFilename = L"wbuf";
        sprintf(buf,"%i",ival++);
	    s.role = buf;
        s.roleValidity = (ival%2);
        ++ival;
	    s.roleConfidence = ival++;
        writer.WriteEntry(FamilialSamplesMultiDataType, 1, s);
    }
    for (int i=0; i<3; i++)
    {
        sprintf(buf,"%i",ival++);
        o.segmentType = buf;
        o.referenceSampleKey = ival++;
        sprintf(buf,"%i",ival++);
        o.referenceSegmentID = buf;
        o.familialSampleKey = ival++;
        sprintf(buf,"%i",ival++);
        o.familialSegmentID = buf;
        writer.WriteEntry(FamilialSegmentOverlapsMultiDataType, 0, o);
        
        sprintf(buf,"%i",ival++);
        o.segmentType = buf;
        o.referenceSampleKey = ival++;
        sprintf(buf,"%i",ival++);
        o.referenceSegmentID = buf;
        o.familialSampleKey = ival++;
        sprintf(buf,"%i",ival++);
        o.familialSegmentID = buf;
        writer.WriteEntry(FamilialSegmentOverlapsMultiDataType, 1, o);
    }
    writer.FlushBuffer();

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("familial1");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");

	CPPUNIT_ASSERT(data2.GetEntryCount(FamilialSamplesMultiDataType) == 3);
	CPPUNIT_ASSERT(data2.GetEntryCount(FamilialSegmentOverlapsMultiDataType) == 3);
    
    ival = 0;
    for (int i=0; i<3; i++)
    {
        data2.GetFamilialSampleEntry(FamilialSamplesMultiDataType, i, s);
        CPPUNIT_ASSERT(s.sampleKey == ival++);
        sprintf(buf,"%i",ival++);
	    CPPUNIT_ASSERT(s.arrID == buf);
        sprintf(buf,"%i",ival++);
	    CPPUNIT_ASSERT(s.chpID == buf);
        //wsprintf(wbuf,L"%i",ival++);
	    CPPUNIT_ASSERT(s.chpFilename == L"wbuf");
        sprintf(buf,"%i",ival++);
	    CPPUNIT_ASSERT(s.role == buf);
        CPPUNIT_ASSERT(s.roleValidity == (ival%2 == 0 ? false : true));
        ++ival;
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(s.roleConfidence, ival++, 0.0001f);

        ival += 6;
    }
    for (int i=0; i<3; i++)
    {
        data2.GetFamilialSegmentOverlapEntry(FamilialSegmentOverlapsMultiDataType, i, o);
        sprintf(buf,"%i",ival++);
        CPPUNIT_ASSERT(o.segmentType == buf);
        CPPUNIT_ASSERT(o.referenceSampleKey == ival++);
        sprintf(buf,"%i",ival++);
        CPPUNIT_ASSERT(o.referenceSegmentID == buf);
        CPPUNIT_ASSERT(o.familialSampleKey == ival++);
        sprintf(buf,"%i",ival++);
        CPPUNIT_ASSERT(o.familialSegmentID == buf);

        ival += 5;
    }
}

static void CreateAllelePeakFile()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;
    affymetrix_calvin_data::AllelePeaks s;

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	CHPMultiDataData data("CHP_MultiData_allele_peaks_buffer");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(AllelePeaksMultiDataType, 2, 10, cols);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::AllelePeaks e;

	writer->SeekToDataSet(AllelePeaksMultiDataType);
	e.name = "abc";
    e.chr = 0;
    e.position = 1;
	e.peaks.clear();
    param.SetValueInt32(0);
    e.peaks.push_back(param);
    param.SetValueFloat(0.0f);
    e.peaks.push_back(param);
	writer->WriteEntry(e);

	e.name = "xyz";
    e.chr = 0;
    e.position = 0;
	e.peaks.clear();
    param.SetValueInt32(0);
    e.peaks.push_back(param);
    param.SetValueFloat(0.0f);
    e.peaks.push_back(param);
	writer->WriteEntry(e);

	delete writer;
}

static void UpdateAllelePeakFile()
{
	ParameterNameValueType nv;
	CHPMultiDataFileBufferWriter writer;
	vector<MultiDataType> dataTypes;
	dataTypes.push_back(AllelePeaksMultiDataType);
	map<MultiDataType, int> maxNameLengths;
	maxNameLengths[AllelePeaksMultiDataType]=10;
	std::vector<std::string> fileNames(1);
	fileNames[0] = "CHP_MultiData_allele_peaks_buffer";
	writer.Initialize(&fileNames, dataTypes, maxNameLengths);

	affymetrix_calvin_data::AllelePeaks e;

	nv.SetName(L"int");
	nv.SetValueInt32(0);
	e.peaks.push_back(nv);
	nv.SetName(L"float");
	nv.SetValueFloat(0);
	e.peaks.push_back(nv);
	for (int i=0; i<2; i++)
	{
		e.name = IntToString(i);
		e.chr = i+1;
		e.position = i+2;
		e.peaks[0].SetValueInt32(i+3);
		e.peaks[1].SetValueFloat(i+4);
		writer.WriteEntry(AllelePeaksMultiDataType, 0, e);
	}
	writer.FlushBuffer();
}

void CalvinCHPMultiDataFileBufferWriterTest::WriteTestAllelePeaks()
{
	CreateAllelePeakFile();
	UpdateAllelePeakFile();

    /////////////////////////////

	ParameterNameValueTypeList params;
	ParameterNameValueType param;
    affymetrix_calvin_data::AllelePeaks e;

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_allele_peaks_buffer");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetNumMetricColumns(AllelePeaksMultiDataType) == 2);
    CPPUNIT_ASSERT(data2.GetMetricColumnName(AllelePeaksMultiDataType, 0) == L"int");
    CPPUNIT_ASSERT(data2.GetMetricColumnName(AllelePeaksMultiDataType, 1) == L"float");

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(AllelePeaksMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

	ParameterNameValueTypeList p = data2.GetAlgParams();
	ParameterNameValueTypeList::iterator it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	p = data2.GetSummaryParams();
	it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");

	data2.GetEntry(AllelePeaksMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.name == "0");
	CPPUNIT_ASSERT(e.chr == 1);
	CPPUNIT_ASSERT(e.position == 2);
	CPPUNIT_ASSERT(e.peaks.size() == 2);
    param = e.peaks[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 3);
    param = e.peaks[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 4.0f, 0.00001f);


	data2.GetEntry(AllelePeaksMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.name == "1");
	CPPUNIT_ASSERT(e.chr == 2);
	CPPUNIT_ASSERT(e.position == 3);
	CPPUNIT_ASSERT(e.peaks.size() == 2);
    param = e.peaks[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 4);
    param = e.peaks[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 5.0f, 0.00001f);

}

// Marker AB Signals


static void CreateMarkerABSignalFile()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	CHPMultiDataData data("CHP_MultiData_marker_signals_buffer");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(MarkerABSignalsMultiDataType, 2, 10, cols);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::MarkerABSignals e;

	writer->SeekToDataSet(MarkerABSignalsMultiDataType);
	e.index = 0;
    e.metrics.clear();
    param.SetValueInt32(0);
    e.metrics.push_back(param);
    param.SetValueFloat(0.0f);
    e.metrics.push_back(param);
	writer->WriteEntry(e);
	e.index = 0;
    e.metrics.clear();
    param.SetValueInt32(0);
    e.metrics.push_back(param);
    param.SetValueFloat(0.0f);
    e.metrics.push_back(param);
	writer->WriteEntry(e);

	delete writer;
}

static void UpdateMarkerABSignalFile()
{
	ParameterNameValueType nv;
	CHPMultiDataFileBufferWriter writer;
	vector<MultiDataType> dataTypes;
	dataTypes.push_back(MarkerABSignalsMultiDataType);
	map<MultiDataType, int> maxNameLengths;
	maxNameLengths[MarkerABSignalsMultiDataType]=10;
	std::vector<std::string> fileNames(1);
	fileNames[0] = "CHP_MultiData_marker_signals_buffer";
	writer.Initialize(&fileNames, dataTypes, maxNameLengths);

	affymetrix_calvin_data::MarkerABSignals e;

	nv.SetName(L"int");
	nv.SetValueInt32(0);
	e.metrics.push_back(nv);
	nv.SetName(L"float");
	nv.SetValueFloat(0);
	e.metrics.push_back(nv);

	for (int i=0; i<2; i++)
	{
		e.index = i;
		e.metrics[0].SetValueInt32(i+3);
		e.metrics[1].SetValueFloat(i+4);
		writer.WriteEntry(MarkerABSignalsMultiDataType, 0, e);
	}
	writer.FlushBuffer();
}

void CalvinCHPMultiDataFileBufferWriterTest::WriteTestMarkerABSignals()
{
	CreateMarkerABSignalFile();
	UpdateMarkerABSignalFile();

    /////////////////////////////

	ParameterNameValueTypeList params;
	ParameterNameValueType param;
    affymetrix_calvin_data::MarkerABSignals e;

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_marker_signals_buffer");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetNumMetricColumns(MarkerABSignalsMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(MarkerABSignalsMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

	ParameterNameValueTypeList p = data2.GetAlgParams();
	ParameterNameValueTypeList::iterator it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	p = data2.GetSummaryParams();
	it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");

	data2.GetEntry(MarkerABSignalsMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.index == 0);
	CPPUNIT_ASSERT(e.metrics.size() == 2);
    param = e.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 3);
    param = e.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 4.0f, 0.00001f);


	data2.GetEntry(MarkerABSignalsMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.index == 1);
	CPPUNIT_ASSERT(e.metrics.size() == 2);
    param = e.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 4);
    param = e.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 5.0f, 0.00001f);

}



static void CreateCytoGenotypeCallDataFile()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	CHPMultiDataData data("CHP_MultiData_cytogenotype_data_buffer");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(CytoGenotypeCallMultiDataType, 2, 10, cols);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::CytoGenotypeCallData e;

	writer->SeekToDataSet(CytoGenotypeCallMultiDataType);
	e.index = 0;
	e.call = 1;
	e.aSignal = 1.0f;
	e.bSignal = 1.0f;
	e.confidence = 1.0f;
	e.signalStrength = 1.0f;
	e.forcedCall = 1;
	e.contrast = 1.0f;

    e.metrics.clear();
    param.SetValueInt32(0);
    e.metrics.push_back(param);
    param.SetValueFloat(0.0f);
    e.metrics.push_back(param);
	writer->WriteEntry(e);
	e.index = 0;
    e.metrics.clear();
    param.SetValueInt32(0);
    e.metrics.push_back(param);
    param.SetValueFloat(0.0f);
    e.metrics.push_back(param);
	writer->WriteEntry(e);

	delete writer;
}

static void UpdateCytoGenotypeCallDataFile()
{
	ParameterNameValueType nv;
	CHPMultiDataFileBufferWriter writer;
	vector<MultiDataType> dataTypes;
	dataTypes.push_back(CytoGenotypeCallMultiDataType);
	map<MultiDataType, int> maxNameLengths;
	maxNameLengths[CytoGenotypeCallMultiDataType]=10;
	std::vector<std::string> fileNames(1);
	fileNames[0] = "CHP_MultiData_cytogenotype_data_buffer";
	writer.Initialize(&fileNames, dataTypes, maxNameLengths);

	affymetrix_calvin_data::CytoGenotypeCallData e;

	nv.SetName(L"int");
	nv.SetValueInt32(4);
	e.metrics.push_back(nv);
	nv.SetName(L"float");
	nv.SetValueFloat(5.0f);
	e.metrics.push_back(nv);

	float testFloat = 0.0f;
	for (int i=0; i<2; i++)
	{
		testFloat++;
		e.index = i;
		e.call = i+1;
		e.forcedCall = i+1;
		e.aSignal = testFloat;
		e.bSignal = testFloat;
		e.confidence = testFloat;
		e.signalStrength = testFloat;
		e.contrast = testFloat;
		e.metrics[0].SetValueInt32(i+3);
		e.metrics[1].SetValueFloat(i+4);
		writer.WriteEntry(CytoGenotypeCallMultiDataType, 0, e);
	}
	writer.FlushBuffer();

}

void CalvinCHPMultiDataFileBufferWriterTest::WriteTestCytoGenotypeCallData()
{
	CreateCytoGenotypeCallDataFile();
	UpdateCytoGenotypeCallDataFile();

    /////////////////////////////

	ParameterNameValueTypeList params;
	ParameterNameValueType param;
    affymetrix_calvin_data::CytoGenotypeCallData e;

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_cytogenotype_data_buffer");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(CytoGenotypeCallMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

	ParameterNameValueTypeList p = data2.GetAlgParams();
	ParameterNameValueTypeList::iterator it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	p = data2.GetSummaryParams();
	it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");

	data2.GetEntry(CytoGenotypeCallMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.index == 0);
	CPPUNIT_ASSERT(e.call == 1);
	CPPUNIT_ASSERT(e.confidence == 1);
	CPPUNIT_ASSERT(e.forcedCall == 1);
	CPPUNIT_ASSERT(e.aSignal == 1.0f);
	CPPUNIT_ASSERT(e.bSignal == 1.0f);
	CPPUNIT_ASSERT(e.signalStrength == 1.0f);
	CPPUNIT_ASSERT(e.contrast == 1.0f);
	CPPUNIT_ASSERT(e.metrics.size() == 2);
    param = e.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 3);
    param = e.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 4.0f, 0.00001f);


	data2.GetEntry(CytoGenotypeCallMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.index == 1);
	CPPUNIT_ASSERT(e.call == 2);
	CPPUNIT_ASSERT(e.confidence == 2);
	CPPUNIT_ASSERT(e.forcedCall == 2);
	CPPUNIT_ASSERT(e.aSignal == 2.0f);
	CPPUNIT_ASSERT(e.bSignal == 2.0f);
	CPPUNIT_ASSERT(e.signalStrength == 2.0f);
	CPPUNIT_ASSERT(e.contrast == 2.0f);
	CPPUNIT_ASSERT(e.metrics.size() == 2);
    param = e.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 4);
    param = e.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 5.0f, 0.00001f);

}





