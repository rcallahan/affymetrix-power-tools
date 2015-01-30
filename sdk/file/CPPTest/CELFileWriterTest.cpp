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

#include "file/CPPTest/CELFileWriterTest.h"
//
#include "file/CELFileWriter.h"
//
#include <cmath>
//

#ifdef _MSC_VER
#define snprintf _snprintf
#endif


// Delimiter character in DAT header 
#define DELIMCHAR 0x14

CPPUNIT_TEST_SUITE_REGISTRATION( CCELFileWriterTest );

#ifdef _MSC_VER
#define V3_OUT_FILE ".\\data\\test.out.v3.CEL"
#define V4_OUT_FILE ".\\data\\test.out.v4.CEL"
#define BCEL_OUT_FILE ".\\data\\test.out.bcel.CEL"
#define CCEL_OUT_FILE ".\\data\\test.out.ccel.CEL"
#else
#define V3_OUT_FILE "./data/test.out.v3.CEL"
#define V4_OUT_FILE "./data/test.out.v4.CEL"
#define BCEL_OUT_FILE "./data/test.out.bcel.CEL"
#define CCEL_OUT_FILE "./data/test.out.ccel.CEL"
#endif

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

void CCELFileWriterTest::setUp()
{
}

void CCELFileWriterTest::tearDown()
{
}

void CCELFileWriterTest::testWrite_v3()
{
	CCELFileWriter cel;
	cel.SetFileName(V3_OUT_FILE);
	cel.SetFileFormat(CCELFileData::TEXT_CEL);
	CPPUNIT_ASSERT(cel.GetFileFormat() == CCELFileData::TEXT_CEL);
	setData(cel);
	CPPUNIT_ASSERT(cel.WriteTextCel() == true);

	cel.Clear();
	cel.SetFileName(V3_OUT_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);

	CPPUNIT_ASSERT(cel.GetVersion() == 3);
	char scanInfo[512];
	snprintf(scanInfo,sizeof(scanInfo), " %c %c chip type.1sq %c %c %c %c %c %c %c %c %c ",
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR,
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR);
	std::string strHeader = "Cols=5\nRows=5\nTotalX=5\nTotalY=5\nOffsetX=0\nOffsetY=0\n";
	strHeader += "GridCornerUL=1 2\nGridCornerUR=3 4\nGridCornerLR=5 6\nGridCornerLL=7 8\n";
	strHeader += "Axis-invertX=0\nAxisInvertY=0\nswapXY=0\nDatHeader=";
	strHeader += scanInfo;
	strHeader += "\nAlgorithm=alg\nAlgorithmParameters=tag1:value1;tag2:value2;CellMargin:2\n";
	CPPUNIT_ASSERT( cel.GetHeaderString() == strHeader);

	for (int ir=0; ir<5; ir++)
	{
		for (int ic=0; ic<5; ic++)
		{
			CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(ic, ir), (float) (ir+ic+0.5)));
			CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(ic, ir), 1.3f));
			CPPUNIT_ASSERT(cel.GetPixels(ic, ir) == 10);
			if ((ic % 2) == 0) CPPUNIT_ASSERT(cel.IsMasked(ic,ir) == true);
			else CPPUNIT_ASSERT(cel.IsMasked(ic,ir) == false);
			if ((ic % 2) == 1) CPPUNIT_ASSERT(cel.IsOutlier(ic,ir) == true);
			else CPPUNIT_ASSERT(cel.IsOutlier(ic,ir) == false);
		}
	}
}

void CCELFileWriterTest::testWrite_v4()
{
	CCELFileWriter cel;
	cel.SetFileName(V4_OUT_FILE);
	cel.SetFileFormat(CCELFileData::XDA_BCEL);
	CPPUNIT_ASSERT(cel.GetFileFormat() == CCELFileData::XDA_BCEL);
	setData(cel);
	CPPUNIT_ASSERT(cel.WriteXDABCel() == true);

	cel.Clear();
	cel.SetFileName(V4_OUT_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);

	CPPUNIT_ASSERT(cel.GetVersion() == 4);
	char scanInfo[512];
	snprintf(scanInfo,sizeof(scanInfo), " %c %c chip type.1sq %c %c %c %c %c %c %c %c %c ",
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR,
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR);
	std::string strHeader = "Cols=5\nRows=5\nTotalX=5\nTotalY=5\nOffsetX=0\nOffsetY=0\n";
	strHeader += "GridCornerUL=1 2\nGridCornerUR=3 4\nGridCornerLR=5 6\nGridCornerLL=7 8\n";
	strHeader += "Axis-invertX=0\nAxisInvertY=0\nswapXY=0\nDatHeader=";
	strHeader += scanInfo;
	strHeader += "\nAlgorithm=alg\nAlgorithmParameters=tag1:value1;tag2:value2;CellMargin:2\n";
	CPPUNIT_ASSERT( cel.GetHeaderString() == strHeader);

	for (int ir=0; ir<5; ir++)
	{
		for (int ic=0; ic<5; ic++)
		{
			CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(ic, ir), (float) (ir+ic+0.5)));
			CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(ic, ir), 1.3f));
			CPPUNIT_ASSERT(cel.GetPixels(ic, ir) == 10);
			if ((ic % 2) == 0) CPPUNIT_ASSERT(cel.IsMasked(ic,ir) == true);
			else CPPUNIT_ASSERT(cel.IsMasked(ic,ir) == false);
			if ((ic % 2) == 1) CPPUNIT_ASSERT(cel.IsOutlier(ic,ir) == true);
			else CPPUNIT_ASSERT(cel.IsOutlier(ic,ir) == false);
		}
	}
}

void CCELFileWriterTest::testWrite_bcel()
{
	CCELFileWriter cel;
	cel.SetFileName(BCEL_OUT_FILE);
	cel.SetFileFormat(CCELFileData::TRANSCRIPTOME_BCEL);
	CPPUNIT_ASSERT(cel.GetFileFormat() == CCELFileData::TRANSCRIPTOME_BCEL);
	setData(cel);
	CPPUNIT_ASSERT(cel.WriteTranscriptomeBCel() == true);

	cel.Clear();
	cel.SetFileName(BCEL_OUT_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);

	CPPUNIT_ASSERT(cel.GetVersion() == 4);
	char scanInfo[512];
	snprintf(scanInfo,sizeof(scanInfo), " %c %c chip type.1sq %c %c %c %c %c %c %c %c %c ",
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR,
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR);
	std::string strHeader = "Cols=5\nRows=5\nTotalX=5\nTotalY=5\nOffsetX=0\nOffsetY=0\n";
	strHeader += "GridCornerUL=1 2\nGridCornerUR=3 4\nGridCornerLR=5 6\nGridCornerLL=7 8\n";
	strHeader += "Axis-invertX=0\nAxisInvertY=0\nswapXY=0\nDatHeader=";
	strHeader += scanInfo;
	strHeader += "\nAlgorithm=alg\nAlgorithmParameters=tag1:value1;tag2:value2;CellMargin:2\n";
	CPPUNIT_ASSERT( cel.GetHeaderString() == strHeader);

	for (int ir=0; ir<5; ir++)
	{
		for (int ic=0; ic<5; ic++)
		{
			CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(ic, ir), (float) (ir+ic+1)));
			CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(ic, ir), 1.0f));
			CPPUNIT_ASSERT(cel.GetPixels(ic, ir) == 10);
			if ((ic % 2) == 0) CPPUNIT_ASSERT(cel.IsMasked(ic,ir) == true);
			else CPPUNIT_ASSERT(cel.IsMasked(ic,ir) == false);
			if ((ic % 2) == 1) CPPUNIT_ASSERT(cel.IsOutlier(ic,ir) == true);
			else CPPUNIT_ASSERT(cel.IsOutlier(ic,ir) == false);
		}
	}
}

void CCELFileWriterTest::testWrite_ccel()
{
	CCELFileWriter cel;
	cel.SetFileName(CCEL_OUT_FILE);
	cel.SetFileFormat(CCELFileData::COMPACT_BCEL);
	CPPUNIT_ASSERT(cel.GetFileFormat() == CCELFileData::COMPACT_BCEL);
	setData(cel);
	CPPUNIT_ASSERT(cel.WriteCompactBCel() == true);

	cel.Clear();
	cel.SetFileName(CCEL_OUT_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);

	CPPUNIT_ASSERT(cel.GetVersion() == 4);
	char scanInfo[512];
	snprintf(scanInfo,sizeof(scanInfo), " %c %c chip type.1sq %c %c %c %c %c %c %c %c %c ",
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR,
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR);
	std::string strHeader = "Cols=5\nRows=5\nTotalX=5\nTotalY=5\nOffsetX=0\nOffsetY=0\n";
	strHeader += "GridCornerUL=1 2\nGridCornerUR=3 4\nGridCornerLR=5 6\nGridCornerLL=7 8\n";
	strHeader += "Axis-invertX=0\nAxisInvertY=0\nswapXY=0\nDatHeader=";
	strHeader += scanInfo;
	strHeader += "\nAlgorithm=alg\nAlgorithmParameters=tag1:value1;tag2:value2;CellMargin:2\n";
	CPPUNIT_ASSERT(cel.GetHeaderString() == strHeader);

	for (int ir=0; ir<5; ir++)
	{
		for (int ic=0; ic<5; ic++)
		{
			CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(ic, ir), (float) (ir+ic+1)));
			CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(ic, ir), 0.0f));
			CPPUNIT_ASSERT(cel.GetPixels(ic, ir) == 0);
			if ((ic % 2) == 0) CPPUNIT_ASSERT(cel.IsMasked(ic,ir) == true);
			else CPPUNIT_ASSERT(cel.IsMasked(ic,ir) == false);
			if ((ic % 2) == 1) CPPUNIT_ASSERT(cel.IsOutlier(ic,ir) == false);
			else CPPUNIT_ASSERT(cel.IsOutlier(ic,ir) == false);
		}
	}
}

void CCELFileWriterTest::setData(CCELFileWriter& cel)
{
	cel.SetDimensions(5, 5);

	GridCoordinatesType corners;
	corners.upperleft.x = 1;
	corners.upperleft.y = 2;
	corners.upperright.x = 3;
	corners.upperright.y = 4;
    corners.lowerright.x = 5;
    corners.lowerright.y = 6;
    corners.lowerleft.x = 7;
    corners.lowerleft.y = 8;
	cel.SetGridCorners(corners);

	cel.SetChipType("chip type");
	cel.SetAlgorithmName("alg");
	cel.AddAlgorithmParameter("tag1", "value1");
	cel.AddAlgorithmParameter("tag2", "value2");
	cel.SetMargin(2);

	for (int ir=0; ir<5; ir++)
	{
		for (int ic=0; ic<5; ic++)
		{
			cel.SetIntensity(ic, ir, (float) (ir+ic+0.5));
			cel.SetStdv(ic, ir, 1.3f);
			cel.SetPixels(ic, ir, 10);
			if ((ic % 2) == 0) cel.SetMask(ic, ir, 1);
			if ((ic % 2) == 1) cel.SetOutlier(ic, ir, 1);
		}
	}
}
