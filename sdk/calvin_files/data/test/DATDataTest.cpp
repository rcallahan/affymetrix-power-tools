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
#include "calvin_files/data/test/DATDataTest.h"
//
#include "calvin_files/data/src/DATData.h"
#include "calvin_files/utils/src/AffyStlCollectionTypes.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DATDataTest );

void DATDataTest::setUp()
{
	data = new DATData("DATMetaData");
}

void DATDataTest::tearDown()
{
	delete data;
}

void DATDataTest::testCreation()
{
	DATData hdr("DATMetaData");
	CPPUNIT_ASSERT(1);
}

void DATDataTest::FilenameTest()
{
	std::string p1 = "DATMetaData";
	data->SetFilename(p1);
	std::string p2 = data->GetFilename();
	CPPUNIT_ASSERT(p1 == p2);
}

void DATDataTest::SetPixelCountTest()
{
	u_int32_t p1 = 100;
	data->SetPixelCount(p1);
	CPPUNIT_ASSERT(1);
}

void DATDataTest::SetStatsCountTest()
{
	u_int32_t p1 = 100;
	data->SetStatsCount(p1);
	CPPUNIT_ASSERT(1);
}

void DATDataTest::ArrayTypeTest()
{
	CPPUNIT_ASSERT_NO_THROW(data->SetArrayType(L"coolio"));
	CPPUNIT_ASSERT(data->GetArrayType() == L"coolio");
}

void DATDataTest::PixelSizeTest()
{
	CPPUNIT_ASSERT_NO_THROW(data->SetPixelSize(0.71f));
	CPPUNIT_ASSERT(data->GetPixelSize() == 0.71f);
}

void DATDataTest::ScannerTypeTest()
{
	CPPUNIT_ASSERT_NO_THROW(data->SetScannerType(L"M10"));
	CPPUNIT_ASSERT(data->GetScannerType() == L"M10");
}

void DATDataTest::ScannerIDTest()
{
	CPPUNIT_ASSERT_NO_THROW(data->SetScannerID(L"Lab1"));
	CPPUNIT_ASSERT(data->GetScannerID() == L"Lab1");
}

void DATDataTest::ScanDateTest()
{
	DateTime dt = DateTime::GetCurrentDateTime();
	CPPUNIT_ASSERT_NO_THROW(data->SetScanDate(dt));
	DateTime dt2 = data->GetScanDate();
	CPPUNIT_ASSERT(dt.Time() == dt2.Time());
	CPPUNIT_ASSERT(dt.Date() == dt2.Date());
	CPPUNIT_ASSERT(dt.IsUTC() == dt2.IsUTC());
}

void DATDataTest::RowsTest()
{
	CPPUNIT_ASSERT_NO_THROW(data->SetRows(234));
	CPPUNIT_ASSERT(data->GetRows() == 234);
}

void DATDataTest::ColsTest()
{
	CPPUNIT_ASSERT_NO_THROW(data->SetCols(8760));
	CPPUNIT_ASSERT(data->GetCols() == 8760);
}

void DATDataTest::GlobalGridTest()
{
	FGridCoords grid;
	grid.lowerleft.x = 1.0f;
	grid.lowerleft.y = 1.2f;
	grid.lowerright.x = 11.0f;
	grid.lowerright.y = 14.2f;
	grid.upperright.x = 111.0f;
	grid.upperright.y = 114.2f;
	grid.upperleft.x = 1001.0f;
	grid.upperleft.y = 1200.2f;

	FRegion reg = grid;
	data->SetGlobalGrid(DATData::GridError, reg);
	FGridCoords grid2 = (FGridCoords)data->GetGlobalGrid();

	CPPUNIT_ASSERT(grid2.lowerleft.x == 1.0f);
	CPPUNIT_ASSERT(grid2.lowerleft.y == 1.2f);
	CPPUNIT_ASSERT(grid2.lowerright.x == 11.0f);
	CPPUNIT_ASSERT(grid2.lowerright.y == 14.2f);
	CPPUNIT_ASSERT(grid2.upperright.x == 111.0f);
	CPPUNIT_ASSERT(grid2.upperright.y == 114.2f);
	CPPUNIT_ASSERT(grid2.upperleft.x == 1001.0f);
	CPPUNIT_ASSERT(grid2.upperleft.y == 1200.2f);
	CPPUNIT_ASSERT(data->GetGlobalGridStatus() == DATData::GridError);
}

void DATDataTest::SubgridTest()
{
	FGridCoords grid;
	grid.lowerleft.x = 1.0f;
	grid.lowerleft.y = 1.2f;
	grid.lowerright.x = 11.0f;
	grid.lowerright.y = 14.2f;
	grid.upperright.x = 111.0f;
	grid.upperright.y = 114.2f;
	grid.upperleft.x = 1001.0f;
	grid.upperleft.y = 1200.2f;

	CPPUNIT_ASSERT(data->GetSubgridCnt() == 0);

	data->AddSubgrid(DATData::GridOK, (FRegion)grid);
	FGridCoords grid2 = grid;
	grid2.lowerleft.x = 2.0f;
	data->AddSubgrid(DATData::GridManualAdjust, (FRegion)grid2);

	CPPUNIT_ASSERT(data->GetSubgridCnt() == 2);

	FGridCoords grid3 = (FGridCoords)data->GetSubgrid(0);
	CPPUNIT_ASSERT(grid3.lowerleft.x == 1.0f);
	CPPUNIT_ASSERT(grid3.lowerleft.y == 1.2f);
	CPPUNIT_ASSERT(grid3.lowerright.x == 11.0f);
	CPPUNIT_ASSERT(grid3.lowerright.y == 14.2f);
	CPPUNIT_ASSERT(grid3.upperright.x == 111.0f);
	CPPUNIT_ASSERT(grid3.upperright.y == 114.2f);
	CPPUNIT_ASSERT(grid3.upperleft.x == 1001.0f);
	CPPUNIT_ASSERT(grid3.upperleft.y == 1200.2f);
	CPPUNIT_ASSERT(data->GetSubgridStatus(0) == DATData::GridOK);

	FGridCoords grid4 = (FGridCoords)data->GetSubgrid(1);
	CPPUNIT_ASSERT(grid4.lowerleft.x == 2.0f);
	CPPUNIT_ASSERT(data->GetSubgridStatus(1) == DATData::GridManualAdjust);

	data->ClearSubgrids();

	CPPUNIT_ASSERT(data->GetSubgridCnt() == 0);
}

void DATDataTest::ArrayIdTest()
{
	AffymetrixGuidType typeId = AffymetrixGuid::GenerateNewGuid();
	data->SetArrayId(typeId);
	CPPUNIT_ASSERT(data->GetArrayId() == typeId);
}

void DATDataTest::ArrayBarcodeTest()
{
	std::wstring barcode = L"soncanyouplaymeamemoryimnotreallysurehowitgoes";
	data->SetArrayBarcode(barcode);
	CPPUNIT_ASSERT(data->GetArrayBarcode() == barcode);
}

void DATDataTest::AddGridAlignmentAlgorithmParameterTest()
{
	ParameterNameValueType nvt;
	nvt.SetName(L"Nirvana");
	nvt.SetValueInt32(1234);
	CPPUNIT_ASSERT_NO_THROW(data->AddGridAlignmentAlgorithmParameter(nvt));
}

void DATDataTest::FindGridAlignmentAlgorithmParameterTest()
{
	ParameterNameValueType nvt;
	nvt.SetName(L"Nirvana");
	nvt.SetValueInt32(1234);
	data->AddGridAlignmentAlgorithmParameter(nvt);
	nvt.SetName(L"Hole");
	nvt.SetValueText(L"OU812");
	data->AddGridAlignmentAlgorithmParameter(nvt);

	ParameterNameValueType nvtResult;
	CPPUNIT_ASSERT(data->FindGridAlignmentAlgorithmParameter(L"Hole", nvtResult));
	CPPUNIT_ASSERT(nvtResult.GetValueText() == L"OU812");
	CPPUNIT_ASSERT(data->FindGridAlignmentAlgorithmParameter(L"Nirvana", nvtResult));
	CPPUNIT_ASSERT(nvtResult.GetValueInt32() == 1234);
	CPPUNIT_ASSERT(data->FindGridAlignmentAlgorithmParameter(L"Missing", nvtResult) == false);
}

void DATDataTest::GetGridAlignmentAlgorithmParametersTest()
{
	ParameterNameValueType nvt;
	nvt.SetName(L"Nirvana");
	nvt.SetValueInt32(94061);
	data->AddGridAlignmentAlgorithmParameter(nvt);
	nvt.SetName(L"Pi");
	nvt.SetValueFloat(3.141592654f);
	data->AddGridAlignmentAlgorithmParameter(nvt);

	ParameterNameValueTypeVector params;
	data->GetGridAlignmentAlgorithmParameters(params);
	CPPUNIT_ASSERT(params.size() == 2);
	CPPUNIT_ASSERT(params[0].GetName() == L"Nirvana");
	CPPUNIT_ASSERT(params[0].GetValueInt32() == 94061);
}

void DATDataTest::ClearGridAlignmentAlgorithmParametersTest()
{
	ParameterNameValueType nvt;
	nvt.SetName(L"Nirvana");
	nvt.SetValueInt32(94061);
	data->AddGridAlignmentAlgorithmParameter(nvt);
	nvt.SetName(L"Pi");
	nvt.SetValueFloat(3.141592654f);
	data->AddGridAlignmentAlgorithmParameter(nvt);

	ParameterNameValueTypeVector params;
	data->GetGridAlignmentAlgorithmParameters(params);
	CPPUNIT_ASSERT(params.size() == 2);
	data->ClearGridAlignmentAlgorithmParameters();
	data->GetGridAlignmentAlgorithmParameters(params);
	CPPUNIT_ASSERT(params.size() == 0);
}

