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

#include "file/CPPTest/CELFileDataTest.h"
//
#include "file/CELFileData.h"
//
#include <cmath>
#include <iostream>
//

#ifdef _MSC_VER
#pragma warning(disable: 4996) // don't show deprecated warnings.
#endif
#ifdef _MSC_VER
#define snprintf _snprintf
#endif

CPPUNIT_TEST_SUITE_REGISTRATION( CCELFileDataTest );

using namespace affxcel;

// Delimiter character in DAT header 
#define DELIMCHAR 0x14

#ifdef _MSC_VER
#define V3_FILE ".\\data\\test.v3.CEL"
#define V4_FILE ".\\data\\test.v4.CEL"
#define NO_FILE ".\\data\\no_file.CEL"
#define BCEL_FILE ".\\data\\test.bcel.CEL"
#define CCEL_FILE ".\\data\\test.ccel.CEL"
#else
#define V3_FILE "./data/test.v3.CEL"
#define V4_FILE "./data/test.v4.CEL"
#define NO_FILE "./data/no_file.CEL"
#define BCEL_FILE "./data/test.bcel.CEL"
#define CCEL_FILE "./data/test.ccel.CEL"
#endif

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

void CCELFileDataTest::setUp()
{
}

void CCELFileDataTest::tearDown()
{
}

void CCELFileDataTest::testmethod_IsVersion3CompatibleFile()
{
	affxcel::CCELFileData cel;
	cel.SetFileName(NO_FILE);
	CPPUNIT_ASSERT(cel.IsVersion3CompatibleFile() == false);
	cel.SetFileName(V3_FILE);
	CPPUNIT_ASSERT(cel.IsVersion3CompatibleFile() == true);
	cel.SetFileName(V4_FILE);
	CPPUNIT_ASSERT(cel.IsVersion3CompatibleFile() == false);
	cel.SetFileName(BCEL_FILE);
	CPPUNIT_ASSERT(cel.IsVersion3CompatibleFile() == false);
	cel.SetFileName(CCEL_FILE);
	CPPUNIT_ASSERT(cel.IsVersion3CompatibleFile() == false);
}

void CCELFileDataTest::testHeader_Version()
{
	affxcel::CCELFileHeaderData header;
	header.SetVersion(3);
	CPPUNIT_ASSERT( header.GetVersion() == 3);
}

void CCELFileDataTest::testHeader_AlgParams()
{
	affxcel::CCELFileHeaderData header;

	header.AddAlgorithmParameter("tag1", "value1");
	header.AddAlgorithmParameter("tag2", "value2");
	CPPUNIT_ASSERT(header.GetAlgorithmParameter("tag1") == "value1");
	CPPUNIT_ASSERT(header.GetAlgorithmParameter("tag2") == "value2");
	CPPUNIT_ASSERT(header.GetAlgorithmParameters() == "tag1:value1;tag2:value2");
	header.SetAlgorithmParameter("tag2", "value2_new");
	header.SetAlgorithmParameter("tag3", "value3");
	CPPUNIT_ASSERT(header.GetAlgorithmParameter("tag2") == "value2_new");
	CPPUNIT_ASSERT(header.GetAlgorithmParameter("tag3") == "");
	CPPUNIT_ASSERT(header.GetAlgorithmParameters() == "tag1:value1;tag2:value2_new");
	CPPUNIT_ASSERT(header.GetParams() == "");
	header.SetParams("tag4:value4;tag5:value5");
	CPPUNIT_ASSERT(header.GetParams() == "tag4:value4;tag5:value5");
	header.ParseAlgorithmParameters();
	CPPUNIT_ASSERT(header.GetAlgorithmParameter("tag4") == "value4");
	CPPUNIT_ASSERT(header.GetAlgorithmParameter("tag5") == "value5");
	CPPUNIT_ASSERT(header.GetAlgorithmParameters() == "tag1:value1;tag2:value2_new;tag4:value4;tag5:value5");
}

void CCELFileDataTest::testHeader_SetGridCorners()
{
	affxcel::CCELFileHeaderData header;
	GridCoordinatesType corners;
	corners.upperleft.x = 1;
	corners.upperleft.y = 2;
	corners.upperright.x = 3;
	corners.upperright.y = 4;
    corners.lowerright.x = 5;
    corners.lowerright.y = 6;
    corners.lowerleft.x = 7;
    corners.lowerleft.y = 8;
	header.SetGridCorners(corners);
	corners = header.GetGridCorners();
	CPPUNIT_ASSERT(corners.upperleft.x == 1);
	CPPUNIT_ASSERT(corners.upperleft.y == 2);
	CPPUNIT_ASSERT(corners.upperright.x == 3);
	CPPUNIT_ASSERT(corners.upperright.y == 4);
	CPPUNIT_ASSERT(corners.lowerright.x == 5);
	CPPUNIT_ASSERT(corners.lowerright.y == 6);
	CPPUNIT_ASSERT(corners.lowerleft.x == 7);
	CPPUNIT_ASSERT(corners.lowerleft.y == 8);
}

void CCELFileDataTest::testHeader_Dimensions()
{
	affxcel::CCELFileHeaderData header;
	header.SetRows(1);
	header.SetCols(2);
	CPPUNIT_ASSERT( header.GetRows() == 1);
	CPPUNIT_ASSERT( header.GetCols() == 2);
}

void CCELFileDataTest::testHeader_AlgName()
{
	affxcel::CCELFileHeaderData header;
	header.SetAlg("alg");
	CPPUNIT_ASSERT( header.GetAlg() == "alg");
}

void CCELFileDataTest::testHeader_ChipType()
{
	affxcel::CCELFileHeaderData header;
	header.SetChipType("chip type");
	CPPUNIT_ASSERT( header.GetChipType() == "chip type");
}

void CCELFileDataTest::testHeader_Margin()
{
	affxcel::CCELFileHeaderData header;
	header.SetMargin(23);
	CPPUNIT_ASSERT( header.GetMargin() == 23);
}

void CCELFileDataTest::testHeader_DatHeader()
{
	affxcel::CCELFileHeaderData header;
	header.SetDatHeader("datheader");
	CPPUNIT_ASSERT( header.GetDatHeader() == "datheader");
}

void CCELFileDataTest::testHeader_Header()
{
	affxcel::CCELFileHeaderData header;

	GridCoordinatesType corners;
	corners.upperleft.x = 1;
	corners.upperleft.y = 2;
	corners.upperright.x = 3;
	corners.upperright.y = 4;
    corners.lowerright.x = 5;
    corners.lowerright.y = 6;
    corners.lowerleft.x = 7;
    corners.lowerleft.y = 8;
	header.SetGridCorners(corners);

	header.SetRows(1);
	header.SetCols(2);
	header.SetChipType("chip type");
	header.SetAlg("alg");
	header.AddAlgorithmParameter("tag1", "value1");
	header.AddAlgorithmParameter("tag2", "value2");
	header.SetDatHeader();
	char scanInfo[1024];
	snprintf(scanInfo,sizeof(scanInfo), " %c %c chip type.1sq %c %c %c %c %c %c %c %c %c ",
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR,
		DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR, DELIMCHAR);
	std::string strHeader = "Cols=2\nRows=1\nTotalX=2\nTotalY=1\nOffsetX=0\nOffsetY=0\n";
	strHeader += "GridCornerUL=1 2\nGridCornerUR=3 4\nGridCornerLR=5 6\nGridCornerLL=7 8\n";
	strHeader += "Axis-invertX=0\nAxisInvertY=0\nswapXY=0\nDatHeader=";
	strHeader += scanInfo;
	strHeader += "\nAlgorithm=alg\nAlgorithmParameters=tag1:value1;tag2:value2\n";
	CPPUNIT_ASSERT( header.GetHeader() == strHeader);
}

void CCELFileDataTest::testHeader_Cells()
{
	affxcel::CCELFileHeaderData header;
	header.SetCells(300);
	CPPUNIT_ASSERT( header.GetCells() == 300);
}

void CCELFileDataTest::testHeader_MaskedCells()
{
	affxcel::CCELFileHeaderData header;
	header.SetMasked(300);
	CPPUNIT_ASSERT( header.GetMasked() == 300);
	header.IncrementMasked();
	CPPUNIT_ASSERT( header.GetMasked() == 301);
	header.DecrementMasked();
	CPPUNIT_ASSERT( header.GetMasked() == 300);
}

void CCELFileDataTest::testHeader_Outliers()
{
	affxcel::CCELFileHeaderData header;
	header.SetOutliers(300);
	CPPUNIT_ASSERT( header.GetOutliers() == 300);
	header.IncrementOutliers();
	CPPUNIT_ASSERT( header.GetOutliers() == 301);
	header.DecrementOutliers();
	CPPUNIT_ASSERT( header.GetOutliers() == 300);
}

void CCELFileDataTest::testCELFileEntryType_size()
{
	CPPUNIT_ASSERT( sizeof(CELFileEntryType) == 10 );
}

void CCELFileDataTest::testCELFileTranscriptomeEntryType_size()
{
	CPPUNIT_ASSERT( sizeof(CELFileTranscriptomeEntryType) == 5 );
}

void CCELFileDataTest::testPaths()
{
	CCELFileData cel;
	cel.SetFileName(V3_FILE);
	CPPUNIT_ASSERT(cel.GetFileName() == V3_FILE);
}

void CCELFileDataTest::testReadHeader_fail()
{
	CCELFileData cel;
	cel.SetFileName(NO_FILE);
	CPPUNIT_ASSERT(cel.ReadHeader() == false);
}

void CCELFileDataTest::testRead_fail()
{
	CCELFileData cel;
	cel.SetFileName(NO_FILE);
	CPPUNIT_ASSERT(cel.Read() == false);
}

void CCELFileDataTest::testRead_v3()
{
	CCELFileData cel;
	cel.SetFileName(V3_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);
	CPPUNIT_ASSERT(cel.GetFileFormat() == CCELFileData::TEXT_CEL);

	std::string header = "Cols=126\n";
	header += "Rows=126\n";
	header += "TotalX=126\n";
	header += "TotalY=126\n";
	header += "OffsetX=0\n";
	header += "OffsetY=0\n";
	header += "GridCornerUL=239 234\n";
	header += "GridCornerUR=4504 232\n";
	header += "GridCornerLR=4505 4497\n";
	header += "GridCornerLL=240 4500\n";
	header += "Axis-invertX=0\n";
	header += "AxisInvertY=0\n";
	header += "swapXY=0\n";
	header += "DatHeader=[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6\n";
	header += "Algorithm=Percentile\n";
	header += "AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000\n";

	// Accessors for header information.
	CPPUNIT_ASSERT( cel.GetVersion() == 3 );
	CPPUNIT_ASSERT( cel.GetCols() == 126 );
	CPPUNIT_ASSERT( cel.GetRows() == 126 );
	CPPUNIT_ASSERT( cel.GetNumCells() == 126*126 );
	CPPUNIT_ASSERT( cel.GetHeaderString() == header );
	CPPUNIT_ASSERT( cel.GetAlg() == "Percentile");
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetChipType() == "Test3");
	CPPUNIT_ASSERT( cel.GetCellMargin() == 2);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	
	CPPUNIT_ASSERT( cel.GetHeaderKey("HEADER") == header);
	CPPUNIT_ASSERT( cel.GetHeaderKey("VERSION") == "3");
	CPPUNIT_ASSERT( cel.GetHeaderKey("COLS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ROWS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALX") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALY") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUL") == "(239, 234)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUR") == "(4504, 232)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLR") == "(4505, 4497)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLL") == "(240, 4500)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXIS-INVERTX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXISINVERTY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("SWAPXY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHM") == "Percentile");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHMPARAMETERS") == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetHeaderKey("DATHEADER") == "[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERCELLS") == "15876");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBEROUTLIERCELLS") == "3917");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERMASKEDCELLS") == "1");

	GridCoordinatesType corners = cel.GetGridCorners();
	CPPUNIT_ASSERT(corners.upperleft.x == 239);
	CPPUNIT_ASSERT(corners.upperleft.y == 234);
	CPPUNIT_ASSERT(corners.upperright.x == 4504);
	CPPUNIT_ASSERT(corners.upperright.y == 232);
	CPPUNIT_ASSERT(corners.lowerright.x == 4505);
	CPPUNIT_ASSERT(corners.lowerright.y == 4497);
	CPPUNIT_ASSERT(corners.lowerleft.x == 240);
	CPPUNIT_ASSERT(corners.lowerleft.y == 4500);

	CPPUNIT_ASSERT(cel.IsMasked(0,0) == true);
	for (int i=1; i<cel.GetNumCells(); i++)
	{
		CPPUNIT_ASSERT(cel.IsMasked(i % cel.GetCols(), i / cel.GetCols()) == false);
	}

	CPPUNIT_ASSERT(cel.IsOutlier(0,0) == true);
	CPPUNIT_ASSERT(cel.IsOutlier(5,0) == true);
	CPPUNIT_ASSERT(cel.IsOutlier(125,125) == true);
	CPPUNIT_ASSERT(cel.IsOutlier(0,1) == false);

	CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(0,0), 17047.0f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(122,125), 577.0f));

	CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(0,0), 10412.4f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(122,125), 6930.0f));

	CPPUNIT_ASSERT(cel.GetPixels(0,0) == 1024);
	CPPUNIT_ASSERT(cel.GetPixels(122,125) == 992);

	CELFileEntryType entry;
	cel.GetEntry(0, 0, entry);
	CPPUNIT_ASSERT(CompareFloats(entry.Intensity, 17047.0f));
	CPPUNIT_ASSERT(CompareFloats(entry.Stdv, 10412.4f));
	CPPUNIT_ASSERT(entry.Pixels == 1024);

	cel.GetEntry(122, 125, entry);
	CPPUNIT_ASSERT(CompareFloats(entry.Intensity, 577.0f));
	CPPUNIT_ASSERT(CompareFloats(entry.Stdv, 6930.0f));
	CPPUNIT_ASSERT(entry.Pixels == 992);
}

void CCELFileDataTest::testRead_v4()
{
	CCELFileData cel;
	cel.SetFileName(V4_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);
	CPPUNIT_ASSERT(cel.GetFileFormat() == CCELFileData::XDA_BCEL);

	std::string header = "Cols=126\n";
	header += "Rows=126\n";
	header += "TotalX=126\n";
	header += "TotalY=126\n";
	header += "OffsetX=0\n";
	header += "OffsetY=0\n";
	header += "GridCornerUL=239 234\n";
	header += "GridCornerUR=4504 232\n";
	header += "GridCornerLR=4505 4497\n";
	header += "GridCornerLL=240 4500\n";
	header += "Axis-invertX=0\n";
	header += "AxisInvertY=0\n";
	header += "swapXY=0\n";
	header += "DatHeader=[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6\n";
	header += "Algorithm=Percentile\n";
	header += "AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000\n";

	// Accessors for header information.
	CPPUNIT_ASSERT( cel.GetVersion() == 4 );
	CPPUNIT_ASSERT( cel.GetCols() == 126 );
	CPPUNIT_ASSERT( cel.GetRows() == 126 );
	CPPUNIT_ASSERT( cel.GetNumCells() == 126*126 );
	CPPUNIT_ASSERT( cel.GetHeaderString() == header );
	CPPUNIT_ASSERT( cel.GetAlg() == "Percentile");
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetChipType() == "Test3");
	CPPUNIT_ASSERT( cel.GetCellMargin() == 2);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	
	CPPUNIT_ASSERT( cel.GetHeaderKey("HEADER") == header);
	CPPUNIT_ASSERT( cel.GetHeaderKey("VERSION") == "4");
	CPPUNIT_ASSERT( cel.GetHeaderKey("COLS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ROWS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALX") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALY") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUL") == "(239, 234)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUR") == "(4504, 232)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLR") == "(4505, 4497)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLL") == "(240, 4500)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXIS-INVERTX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXISINVERTY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("SWAPXY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHM") == "Percentile");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHMPARAMETERS") == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetHeaderKey("DATHEADER") == "[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERCELLS") == "15876");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBEROUTLIERCELLS") == "3917");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERMASKEDCELLS") == "1");

	GridCoordinatesType corners = cel.GetGridCorners();
	CPPUNIT_ASSERT(corners.upperleft.x == 239);
	CPPUNIT_ASSERT(corners.upperleft.y == 234);
	CPPUNIT_ASSERT(corners.upperright.x == 4504);
	CPPUNIT_ASSERT(corners.upperright.y == 232);
	CPPUNIT_ASSERT(corners.lowerright.x == 4505);
	CPPUNIT_ASSERT(corners.lowerright.y == 4497);
	CPPUNIT_ASSERT(corners.lowerleft.x == 240);
	CPPUNIT_ASSERT(corners.lowerleft.y == 4500);

	CPPUNIT_ASSERT(cel.IsMasked(0,0) == true);
	for (int i=1; i<cel.GetNumCells(); i++)
	{
		CPPUNIT_ASSERT(cel.IsMasked(i % cel.GetCols(), i / cel.GetCols()) == false);
	}

	CPPUNIT_ASSERT(cel.IsOutlier(0,0) == true);
	CPPUNIT_ASSERT(cel.IsOutlier(5,0) == true);
	CPPUNIT_ASSERT(cel.IsOutlier(125,125) == true);
	CPPUNIT_ASSERT(cel.IsOutlier(0,1) == false);

	CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(0,0), 17047.0f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(122,125), 577.0f));

	CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(0,0), 10412.4f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(122,125), 6930.0f));

	CPPUNIT_ASSERT(cel.GetPixels(0,0) == 1024);
	CPPUNIT_ASSERT(cel.GetPixels(122,125) == 992);

	CELFileEntryType entry;
	cel.GetEntry(0, 0, entry);
	CPPUNIT_ASSERT(CompareFloats(entry.Intensity, 17047.0f));
	CPPUNIT_ASSERT(CompareFloats(entry.Stdv, 10412.4f));
	CPPUNIT_ASSERT(entry.Pixels == 1024);

	cel.GetEntry(122, 125, entry);
	CPPUNIT_ASSERT(CompareFloats(entry.Intensity, 577.0f));
	CPPUNIT_ASSERT(CompareFloats(entry.Stdv, 6930.0f));
	CPPUNIT_ASSERT(entry.Pixels == 992);
}

void CCELFileDataTest::testRead_bcel()
{
	CCELFileData cel;
	cel.SetFileName(BCEL_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);
	CPPUNIT_ASSERT(cel.GetFileFormat() == CCELFileData::TRANSCRIPTOME_BCEL);

	std::string header = "Cols=126\n";
	header += "Rows=126\n";
	header += "TotalX=126\n";
	header += "TotalY=126\n";
	header += "OffsetX=0\n";
	header += "OffsetY=0\n";
	header += "GridCornerUL=239 234\n";
	header += "GridCornerUR=4504 232\n";
	header += "GridCornerLR=4505 4497\n";
	header += "GridCornerLL=240 4500\n";
	header += "Axis-invertX=0\n";
	header += "AxisInvertY=0\n";
	header += "swapXY=0\n";
	header += "DatHeader=[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6\n";
	header += "Algorithm=Percentile\n";
	header += "AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000\n";

	// Accessors for header information.
	CPPUNIT_ASSERT( cel.GetVersion() == 4 );
	CPPUNIT_ASSERT( cel.GetCols() == 126 );
	CPPUNIT_ASSERT( cel.GetRows() == 126 );
	CPPUNIT_ASSERT( cel.GetNumCells() == 126*126 );
	CPPUNIT_ASSERT( cel.GetHeaderString() == header );
	CPPUNIT_ASSERT( cel.GetAlg() == "Percentile");
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetChipType() == "Test3");
	CPPUNIT_ASSERT( cel.GetCellMargin() == 2);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	
	CPPUNIT_ASSERT( cel.GetHeaderKey("HEADER") == header);
	CPPUNIT_ASSERT( cel.GetHeaderKey("VERSION") == "4");
	CPPUNIT_ASSERT( cel.GetHeaderKey("COLS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ROWS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALX") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALY") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUL") == "(239, 234)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUR") == "(4504, 232)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLR") == "(4505, 4497)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLL") == "(240, 4500)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXIS-INVERTX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXISINVERTY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("SWAPXY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHM") == "Percentile");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHMPARAMETERS") == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetHeaderKey("DATHEADER") == "[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERCELLS") == "15876");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBEROUTLIERCELLS") == "3917");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERMASKEDCELLS") == "1");

	GridCoordinatesType corners = cel.GetGridCorners();
	CPPUNIT_ASSERT(corners.upperleft.x == 239);
	CPPUNIT_ASSERT(corners.upperleft.y == 234);
	CPPUNIT_ASSERT(corners.upperright.x == 4504);
	CPPUNIT_ASSERT(corners.upperright.y == 232);
	CPPUNIT_ASSERT(corners.lowerright.x == 4505);
	CPPUNIT_ASSERT(corners.lowerright.y == 4497);
	CPPUNIT_ASSERT(corners.lowerleft.x == 240);
	CPPUNIT_ASSERT(corners.lowerleft.y == 4500);

	CPPUNIT_ASSERT(cel.IsMasked(0,0) == true);
	for (int i=1; i<cel.GetNumCells(); i++)
	{
		CPPUNIT_ASSERT(cel.IsMasked(i % cel.GetCols(), i / cel.GetCols()) == false);
	}

	CPPUNIT_ASSERT(cel.IsOutlier(0,0) == true);
	CPPUNIT_ASSERT(cel.IsOutlier(5,0) == true);
	CPPUNIT_ASSERT(cel.IsOutlier(125,125) == true);
	CPPUNIT_ASSERT(cel.IsOutlier(0,1) == false);

	CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(0,0), 17047.0f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(6,0), 1418.0f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(122,125), 577.0f));

	CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(0,0), 10412.0f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(6,0), 3392.0f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(122,125), 6930.0f));

	CPPUNIT_ASSERT(cel.GetPixels(0,0) == 0);
	CPPUNIT_ASSERT(cel.GetPixels(122,125) == 224);

	CELFileTranscriptomeEntryType entry;
	cel.GetTranscriptomeEntry(0, 0, entry);
	CPPUNIT_ASSERT(CompareFloats(entry.Intensity, 17047.0f));
	CPPUNIT_ASSERT(CompareFloats(entry.Stdv, 10412.0f));
	CPPUNIT_ASSERT(entry.Pixels == 0);

	cel.GetTranscriptomeEntry(6, 0, entry);
	CPPUNIT_ASSERT(CompareFloats(entry.Intensity, 1418.0f));
	CPPUNIT_ASSERT(CompareFloats(entry.Stdv, 3392.0f));

	cel.GetTranscriptomeEntry(122, 125, entry);
	CPPUNIT_ASSERT(CompareFloats(entry.Intensity, 577.0f));
	CPPUNIT_ASSERT(CompareFloats(entry.Stdv, 6930.0f));
	CPPUNIT_ASSERT(entry.Pixels == 224);
}

void CCELFileDataTest::testRead_ccel()
{
	CCELFileData cel;
	cel.SetFileName(CCEL_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);
	CPPUNIT_ASSERT(cel.GetFileFormat() == CCELFileData::COMPACT_BCEL);

	std::string header = "Cols=126\n";
	header += "Rows=126\n";
	header += "TotalX=126\n";
	header += "TotalY=126\n";
	header += "OffsetX=0\n";
	header += "OffsetY=0\n";
	header += "GridCornerUL=239 234\n";
	header += "GridCornerUR=4504 232\n";
	header += "GridCornerLR=4505 4497\n";
	header += "GridCornerLL=240 4500\n";
	header += "Axis-invertX=0\n";
	header += "AxisInvertY=0\n";
	header += "swapXY=0\n";
	header += "DatHeader=[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6\n";
	header += "Algorithm=Percentile\n";
	header += "AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000\n";

	// Accessors for header information.
	CPPUNIT_ASSERT( cel.GetVersion() == 4 );
	CPPUNIT_ASSERT( cel.GetCols() == 126 );
	CPPUNIT_ASSERT( cel.GetRows() == 126 );
	CPPUNIT_ASSERT( cel.GetNumCells() == 126*126 );
	CPPUNIT_ASSERT( cel.GetHeaderString() == header );
	CPPUNIT_ASSERT( cel.GetAlg() == "Percentile");
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetChipType() == "Test3");
	CPPUNIT_ASSERT( cel.GetCellMargin() == 2);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	
	CPPUNIT_ASSERT( cel.GetHeaderKey("HEADER") == header);
	CPPUNIT_ASSERT( cel.GetHeaderKey("VERSION") == "4");
	CPPUNIT_ASSERT( cel.GetHeaderKey("COLS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ROWS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALX") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALY") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUL") == "(239, 234)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUR") == "(4504, 232)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLR") == "(4505, 4497)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLL") == "(240, 4500)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXIS-INVERTX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXISINVERTY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("SWAPXY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHM") == "Percentile");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHMPARAMETERS") == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetHeaderKey("DATHEADER") == "[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERCELLS") == "15876");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBEROUTLIERCELLS") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERMASKEDCELLS") == "1");

	GridCoordinatesType corners = cel.GetGridCorners();
	CPPUNIT_ASSERT(corners.upperleft.x == 239);
	CPPUNIT_ASSERT(corners.upperleft.y == 234);
	CPPUNIT_ASSERT(corners.upperright.x == 4504);
	CPPUNIT_ASSERT(corners.upperright.y == 232);
	CPPUNIT_ASSERT(corners.lowerright.x == 4505);
	CPPUNIT_ASSERT(corners.lowerright.y == 4497);
	CPPUNIT_ASSERT(corners.lowerleft.x == 240);
	CPPUNIT_ASSERT(corners.lowerleft.y == 4500);

	CPPUNIT_ASSERT(cel.IsMasked(0,0) == true);
	for (int i=1; i<cel.GetNumCells(); i++)
	{
		CPPUNIT_ASSERT(cel.IsMasked(i % cel.GetCols(), i / cel.GetCols()) == false);
	}

	CPPUNIT_ASSERT(cel.IsOutlier(0,0) == false);
	CPPUNIT_ASSERT(cel.IsOutlier(5,0) == false);
	CPPUNIT_ASSERT(cel.IsOutlier(125,125) == false);
	CPPUNIT_ASSERT(cel.IsOutlier(0,1) == false);

	CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(0,0), 17047.0f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetIntensity(122,125), 577.0f));

	CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(0,0), 0.0f));
	CPPUNIT_ASSERT(CompareFloats(cel.GetStdv(122,125), 0.0f));

	CPPUNIT_ASSERT(cel.GetPixels(0,0) == 0);
	CPPUNIT_ASSERT(cel.GetPixels(122,125) == 0);
}

void CCELFileDataTest::testRead_v3_no_mask_outlier()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(V3_FILE, CCELFileData::CEL_DATA) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_DATA);
}

void CCELFileDataTest::testRead_v4_no_mask_outlier()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(V4_FILE, CCELFileData::CEL_DATA) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_DATA);
}

void CCELFileDataTest::testRead_bcel_no_mask_outlier()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(BCEL_FILE, CCELFileData::CEL_DATA) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_DATA);
}

void CCELFileDataTest::testRead_ccel_no_mask_outlier()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(CCEL_FILE, CCELFileData::CEL_DATA) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_DATA);
}

void CCELFileDataTest::testRead_v3_no_mask()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(V3_FILE, CCELFileData::CEL_OUTLIER) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_OUTLIER);
}

void CCELFileDataTest::testRead_v4_no_mask()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(V4_FILE, CCELFileData::CEL_OUTLIER) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_OUTLIER);
}

void CCELFileDataTest::testRead_bcel_no_mask()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(BCEL_FILE, CCELFileData::CEL_OUTLIER) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_OUTLIER);
}

void CCELFileDataTest::testRead_ccel_no_mask()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(CCEL_FILE, CCELFileData::CEL_OUTLIER) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_OUTLIER);
}

void CCELFileDataTest::testRead_v3_no_outlier()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(V3_FILE, CCELFileData::CEL_MASK) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_MASK);
}

void CCELFileDataTest::testRead_v4_no_outlier()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(V4_FILE, CCELFileData::CEL_MASK) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_MASK);
}

void CCELFileDataTest::testRead_bcel_no_outlier()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(BCEL_FILE, CCELFileData::CEL_MASK) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_MASK);
}

void CCELFileDataTest::testRead_ccel_no_outlier()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(CCEL_FILE, CCELFileData::CEL_MASK) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_MASK);
}

void CCELFileDataTest::testRead_v3_all()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(V3_FILE, CCELFileData::CEL_ALL) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_ALL);
}

void CCELFileDataTest::testRead_v4_all()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(V4_FILE, CCELFileData::CEL_ALL) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_ALL);
}

void CCELFileDataTest::testRead_bcel_all()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(BCEL_FILE, CCELFileData::CEL_ALL) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_ALL);
}

void CCELFileDataTest::testRead_ccel_all()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.ReadEx(CCEL_FILE, CCELFileData::CEL_ALL) == true);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetReadState() == CCELFileData::CEL_ALL);
}

void CCELFileDataTest::testReadHeader_v3()
{
	CCELFileData cel;
	cel.SetFileName(V3_FILE);
	CPPUNIT_ASSERT(cel.ReadHeader() == true);

	std::string header = "Cols=126\n";
	header += "Rows=126\n";
	header += "TotalX=126\n";
	header += "TotalY=126\n";
	header += "OffsetX=0\n";
	header += "OffsetY=0\n";
	header += "GridCornerUL=239 234\n";
	header += "GridCornerUR=4504 232\n";
	header += "GridCornerLR=4505 4497\n";
	header += "GridCornerLL=240 4500\n";
	header += "Axis-invertX=0\n";
	header += "AxisInvertY=0\n";
	header += "swapXY=0\n";
	header += "DatHeader=[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6\n";
	header += "Algorithm=Percentile\n";
	header += "AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000\n";

	// Accessors for header information.
	CPPUNIT_ASSERT( cel.GetVersion() == 3 );
	CPPUNIT_ASSERT( cel.GetCols() == 126 );
	CPPUNIT_ASSERT( cel.GetRows() == 126 );
	CPPUNIT_ASSERT( cel.GetNumCells() == cel.GetCols() * cel.GetRows() );
	CPPUNIT_ASSERT( cel.GetHeaderString() == header );
	CPPUNIT_ASSERT( cel.GetAlg() == "Percentile");
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetChipType() == "Test3");
	CPPUNIT_ASSERT( cel.GetCellMargin() == 2);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	
	CPPUNIT_ASSERT( cel.GetHeaderKey("HEADER") == header);
	CPPUNIT_ASSERT( cel.GetHeaderKey("VERSION") == "3");
	CPPUNIT_ASSERT( cel.GetHeaderKey("COLS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ROWS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALX") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALY") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUL") == "(239, 234)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUR") == "(4504, 232)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLR") == "(4505, 4497)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLL") == "(240, 4500)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXIS-INVERTX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXISINVERTY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("SWAPXY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHM") == "Percentile");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHMPARAMETERS") == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetHeaderKey("DATHEADER") == "[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERCELLS") == "15876");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBEROUTLIERCELLS") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERMASKEDCELLS") == "0");

	GridCoordinatesType corners = cel.GetGridCorners();
	CPPUNIT_ASSERT(corners.upperleft.x == 239);
	CPPUNIT_ASSERT(corners.upperleft.y == 234);
	CPPUNIT_ASSERT(corners.upperright.x == 4504);
	CPPUNIT_ASSERT(corners.upperright.y == 232);
	CPPUNIT_ASSERT(corners.lowerright.x == 4505);
	CPPUNIT_ASSERT(corners.lowerright.y == 4497);
	CPPUNIT_ASSERT(corners.lowerleft.x == 240);
	CPPUNIT_ASSERT(corners.lowerleft.y == 4500);
}

void CCELFileDataTest::testReadHeader_v4()
{
	CCELFileData cel;
	cel.SetFileName(V4_FILE);
	CPPUNIT_ASSERT(cel.ReadHeader() == true);

	std::string header = "Cols=126\n";
	header += "Rows=126\n";
	header += "TotalX=126\n";
	header += "TotalY=126\n";
	header += "OffsetX=0\n";
	header += "OffsetY=0\n";
	header += "GridCornerUL=239 234\n";
	header += "GridCornerUR=4504 232\n";
	header += "GridCornerLR=4505 4497\n";
	header += "GridCornerLL=240 4500\n";
	header += "Axis-invertX=0\n";
	header += "AxisInvertY=0\n";
	header += "swapXY=0\n";
	header += "DatHeader=[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6\n";
	header += "Algorithm=Percentile\n";
	header += "AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000\n";

	// Accessors for header information.
	CPPUNIT_ASSERT( cel.GetVersion() == 4 );
	CPPUNIT_ASSERT( cel.GetCols() == 126 );
	CPPUNIT_ASSERT( cel.GetRows() == 126 );
	CPPUNIT_ASSERT( cel.GetNumCells() == 126*126 );
	CPPUNIT_ASSERT( cel.GetHeaderString() == header );
	CPPUNIT_ASSERT( cel.GetAlg() == "Percentile");
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetChipType() == "Test3");
	CPPUNIT_ASSERT( cel.GetCellMargin() == 2);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 3917);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	
	CPPUNIT_ASSERT( cel.GetHeaderKey("HEADER") == header);
	CPPUNIT_ASSERT( cel.GetHeaderKey("VERSION") == "4");
	CPPUNIT_ASSERT( cel.GetHeaderKey("COLS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ROWS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALX") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALY") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUL") == "(239, 234)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUR") == "(4504, 232)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLR") == "(4505, 4497)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLL") == "(240, 4500)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXIS-INVERTX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXISINVERTY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("SWAPXY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHM") == "Percentile");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHMPARAMETERS") == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetHeaderKey("DATHEADER") == "[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERCELLS") == "15876");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBEROUTLIERCELLS") == "3917");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERMASKEDCELLS") == "1");

	GridCoordinatesType corners = cel.GetGridCorners();
	CPPUNIT_ASSERT(corners.upperleft.x == 239);
	CPPUNIT_ASSERT(corners.upperleft.y == 234);
	CPPUNIT_ASSERT(corners.upperright.x == 4504);
	CPPUNIT_ASSERT(corners.upperright.y == 232);
	CPPUNIT_ASSERT(corners.lowerright.x == 4505);
	CPPUNIT_ASSERT(corners.lowerright.y == 4497);
	CPPUNIT_ASSERT(corners.lowerleft.x == 240);
	CPPUNIT_ASSERT(corners.lowerleft.y == 4500);
}

void CCELFileDataTest::testReadHeader_bcel()
{
	CCELFileData cel;
	cel.SetFileName(BCEL_FILE);
	CPPUNIT_ASSERT(cel.ReadHeader() == true);

	std::string header = "Cols=126\n";
	header += "Rows=126\n";
	header += "TotalX=126\n";
	header += "TotalY=126\n";
	header += "OffsetX=0\n";
	header += "OffsetY=0\n";
	header += "GridCornerUL=239 234\n";
	header += "GridCornerUR=4504 232\n";
	header += "GridCornerLR=4505 4497\n";
	header += "GridCornerLL=240 4500\n";
	header += "Axis-invertX=0\n";
	header += "AxisInvertY=0\n";
	header += "swapXY=0\n";
	header += "DatHeader=[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6\n";
	header += "Algorithm=Percentile\n";
	header += "AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000\n";

	// Accessors for header information.
	CPPUNIT_ASSERT( cel.GetVersion() == 4 );
	CPPUNIT_ASSERT( cel.GetCols() == 126 );
	CPPUNIT_ASSERT( cel.GetRows() == 126 );
	CPPUNIT_ASSERT( cel.GetNumCells() == 126*126 );
	CPPUNIT_ASSERT( cel.GetHeaderString() == header );
	CPPUNIT_ASSERT( cel.GetAlg() == "Percentile");
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetChipType() == "Test3");
	CPPUNIT_ASSERT( cel.GetCellMargin() == 2);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 0);
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	
	CPPUNIT_ASSERT( cel.GetHeaderKey("HEADER") == header);
	CPPUNIT_ASSERT( cel.GetHeaderKey("VERSION") == "4");
	CPPUNIT_ASSERT( cel.GetHeaderKey("COLS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ROWS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALX") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALY") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUL") == "(239, 234)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUR") == "(4504, 232)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLR") == "(4505, 4497)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLL") == "(240, 4500)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXIS-INVERTX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXISINVERTY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("SWAPXY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHM") == "Percentile");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHMPARAMETERS") == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetHeaderKey("DATHEADER") == "[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERCELLS") == "15876");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBEROUTLIERCELLS") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERMASKEDCELLS") == "0");

	GridCoordinatesType corners = cel.GetGridCorners();
	CPPUNIT_ASSERT(corners.upperleft.x == 239);
	CPPUNIT_ASSERT(corners.upperleft.y == 234);
	CPPUNIT_ASSERT(corners.upperright.x == 4504);
	CPPUNIT_ASSERT(corners.upperright.y == 232);
	CPPUNIT_ASSERT(corners.lowerright.x == 4505);
	CPPUNIT_ASSERT(corners.lowerright.y == 4497);
	CPPUNIT_ASSERT(corners.lowerleft.x == 240);
	CPPUNIT_ASSERT(corners.lowerleft.y == 4500);
}

void CCELFileDataTest::testReadHeader_ccel()
{
	CCELFileData cel;
	cel.SetFileName(CCEL_FILE);
	CPPUNIT_ASSERT(cel.ReadHeader() == true);

	std::string header = "Cols=126\n";
	header += "Rows=126\n";
	header += "TotalX=126\n";
	header += "TotalY=126\n";
	header += "OffsetX=0\n";
	header += "OffsetY=0\n";
	header += "GridCornerUL=239 234\n";
	header += "GridCornerUR=4504 232\n";
	header += "GridCornerLR=4505 4497\n";
	header += "GridCornerLL=240 4500\n";
	header += "Axis-invertX=0\n";
	header += "AxisInvertY=0\n";
	header += "swapXY=0\n";
	header += "DatHeader=[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6\n";
	header += "Algorithm=Percentile\n";
	header += "AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000\n";

	// Accessors for header information.
	CPPUNIT_ASSERT( cel.GetVersion() == 4 );
	CPPUNIT_ASSERT( cel.GetCols() == 126 );
	CPPUNIT_ASSERT( cel.GetRows() == 126 );
	CPPUNIT_ASSERT( cel.GetNumCells() == 126*126 );
	CPPUNIT_ASSERT( cel.GetHeaderString() == header );
	CPPUNIT_ASSERT( cel.GetAlg() == "Percentile");
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetChipType() == "Test3");
	CPPUNIT_ASSERT( cel.GetCellMargin() == 2);
	CPPUNIT_ASSERT( cel.GetNumOutliers() == 0);
	CPPUNIT_ASSERT( cel.GetNumMasked() == 1);
	CPPUNIT_ASSERT( cel.GetParams() == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	
	CPPUNIT_ASSERT( cel.GetHeaderKey("HEADER") == header);
	CPPUNIT_ASSERT( cel.GetHeaderKey("VERSION") == "4");
	CPPUNIT_ASSERT( cel.GetHeaderKey("COLS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ROWS") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALX") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("TOTALY") == "126");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUL") == "(239, 234)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERUR") == "(4504, 232)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLR") == "(4505, 4497)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("GRIDCORNERLL") == "(240, 4500)");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("OFFSETY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXIS-INVERTX") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("AXISINVERTY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("SWAPXY") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHM") == "Percentile");
	CPPUNIT_ASSERT( cel.GetHeaderKey("ALGORITHMPARAMETERS") == "Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:FALSE;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:FALSE;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000");
	CPPUNIT_ASSERT( cel.GetHeaderKey("DATHEADER") == "[81..46222]  Test3:CLS=4733 RWS=4733 XIN=3  YIN=3  VE=17        2.0 02/25/ 0 10:35:25       GridVerify=None  Test3.1sq                  6");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERCELLS") == "15876");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBEROUTLIERCELLS") == "0");
	CPPUNIT_ASSERT( cel.GetHeaderKey("NUMBERMASKEDCELLS") == "1");

	GridCoordinatesType corners = cel.GetGridCorners();
	CPPUNIT_ASSERT(corners.upperleft.x == 239);
	CPPUNIT_ASSERT(corners.upperleft.y == 234);
	CPPUNIT_ASSERT(corners.upperright.x == 4504);
	CPPUNIT_ASSERT(corners.upperright.y == 232);
	CPPUNIT_ASSERT(corners.lowerright.x == 4505);
	CPPUNIT_ASSERT(corners.lowerright.y == 4497);
	CPPUNIT_ASSERT(corners.lowerleft.x == 240);
	CPPUNIT_ASSERT(corners.lowerleft.y == 4500);
}

void CCELFileDataTest::testExists()
{
	CCELFileData cel;
	cel.SetFileName(NO_FILE);
	CPPUNIT_ASSERT(cel.Exists() == false);
	cel.SetFileName(V3_FILE);
	CPPUNIT_ASSERT(cel.Exists() == true);
}

void CCELFileDataTest::testIsXDACompatibleFile()
{
	CCELFileData cel;
	cel.SetFileName(V4_FILE);
	CPPUNIT_ASSERT(cel.IsXDACompatibleFile() == true);
	cel.SetFileName(V3_FILE);
	CPPUNIT_ASSERT(cel.IsXDACompatibleFile() == false);
	cel.SetFileName(BCEL_FILE);
	CPPUNIT_ASSERT(cel.IsXDACompatibleFile() == false);
	cel.SetFileName(CCEL_FILE);
	CPPUNIT_ASSERT(cel.IsXDACompatibleFile() == false);
}

void CCELFileDataTest::testIsTranscriptomeBcelFile()
{
	CCELFileData cel;
	cel.SetFileName(V4_FILE);
	CPPUNIT_ASSERT(cel.IsTranscriptomeBcelFile() == false);
	cel.SetFileName(V3_FILE);
	CPPUNIT_ASSERT(cel.IsTranscriptomeBcelFile() == false);
	cel.SetFileName(BCEL_FILE);
	CPPUNIT_ASSERT(cel.IsTranscriptomeBcelFile() == true);
	cel.SetFileName(CCEL_FILE);
	CPPUNIT_ASSERT(cel.IsTranscriptomeBcelFile() == false);
}

void CCELFileDataTest::testIsCompactCelFile()
{
	CCELFileData cel;
	cel.SetFileName(V4_FILE);
	CPPUNIT_ASSERT(cel.IsCompactCelFile() == false);
	cel.SetFileName(V3_FILE);
	CPPUNIT_ASSERT(cel.IsCompactCelFile() == false);
	cel.SetFileName(BCEL_FILE);
	CPPUNIT_ASSERT(cel.IsCompactCelFile() == false);
	cel.SetFileName(CCEL_FILE);
	CPPUNIT_ASSERT(cel.IsCompactCelFile() == true);
}

void CCELFileDataTest::testXYToIndex()
{
	CCELFileData cel;
	cel.SetFileName(V4_FILE);
	cel.Read();
	CPPUNIT_ASSERT( cel.XYToIndex(0,0) == 0);
	CPPUNIT_ASSERT( cel.XYToIndex(1,0) == 1);
	CPPUNIT_ASSERT( cel.XYToIndex(0,1) == cel.GetCols());
	CPPUNIT_ASSERT( cel.XYToIndex(3,2) == 2*cel.GetCols() + 3);

	CPPUNIT_ASSERT( CCELFileData::XYToIndex(0, 0, 100, 100) == 0);
	CPPUNIT_ASSERT( CCELFileData::XYToIndex(1, 0, 100, 100) == 1);
	CPPUNIT_ASSERT( CCELFileData::XYToIndex(0, 1, 100, 100) == 100);
	CPPUNIT_ASSERT( CCELFileData::XYToIndex(3, 2, 100, 100) == 203);
}

void CCELFileDataTest::testIndexToXY()
{
	CCELFileData cel;
	cel.SetFileName(V4_FILE);
	cel.Read();
	CPPUNIT_ASSERT( cel.IndexToX(0) == 0);
	CPPUNIT_ASSERT( cel.IndexToX(1) == 1);
	CPPUNIT_ASSERT( cel.IndexToX(cel.GetCols()-1) == cel.GetCols()-1);
	CPPUNIT_ASSERT( cel.IndexToX(cel.GetCols()) == 0);

	CPPUNIT_ASSERT( cel.IndexToY(0) == 0);
	CPPUNIT_ASSERT( cel.IndexToY(1) == 0);
	CPPUNIT_ASSERT( cel.IndexToY(cel.GetCols()-1) == 0);
	CPPUNIT_ASSERT( cel.IndexToY(cel.GetCols()) == 1);
}

void CCELFileDataTest::testParameterNames()
{
	affxcel::CCELFileHeaderData header;

	header.AddAlgorithmParameter("tag1", "value1");
	header.AddAlgorithmParameter("tag2", "value2");
	CPPUNIT_ASSERT(header.GetNumberAlgorithmParameters() == 2);
	CPPUNIT_ASSERT(header.GetAlgorithmParameterTag(0) == "tag1");
	CPPUNIT_ASSERT(header.GetAlgorithmParameterTag(1) == "tag2");

	CCELFileData cel;
	cel.SetFileName(V4_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);
	CPPUNIT_ASSERT(cel.GetNumberAlgorithmParameters() == 12);
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(0) == "Percentile");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(1) == "CellMargin");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(2) == "OutlierHigh");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(3) == "OutlierLow");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(4) == "AlgVersion");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(5) == "FixedCellSize");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(6) == "IgnoreOutliersInShiftRows");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(7) == "FeatureExtraction");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(8) == "UseSubgrids");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(9) == "RandomizePixels");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(10) == "ErrorBasis");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(11) == "StdMult");

	cel.Clear();
	CPPUNIT_ASSERT(cel.GetNumberAlgorithmParameters() == 0);
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(0) == "");

	cel.SetFileName(V3_FILE);
	CPPUNIT_ASSERT(cel.Read() == true);
	CPPUNIT_ASSERT(cel.GetNumberAlgorithmParameters() == 12);
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(0) == "Percentile");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(1) == "CellMargin");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(2) == "OutlierHigh");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(3) == "OutlierLow");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(4) == "AlgVersion");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(5) == "FixedCellSize");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(6) == "IgnoreOutliersInShiftRows");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(7) == "FeatureExtraction");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(8) == "UseSubgrids");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(9) == "RandomizePixels");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(10) == "ErrorBasis");
	CPPUNIT_ASSERT(cel.GetAlgorithmParameterTag(11) == "StdMult");
}

void CCELFileDataTest::testmethod_GetHeaderKey_with_invalid_input()
{
	CCELFileData cel;
	CPPUNIT_ASSERT(cel.GetHeaderKey("invalid_key") == "");
}
