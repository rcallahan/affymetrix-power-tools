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

#include "calvin_files/writers/test/DATFileUpdaterTest.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DATFileUpdaterTest );

void DATFileUpdaterTest::setUp() {}

void DATFileUpdaterTest::tearDown() {}

void DATFileUpdaterTest::testCreation()
{
	std::string filename = "DAT_file";
	WriteDatFile(filename);
	DATData data(filename);
	DATFileUpdater updater(data);
	CPPUNIT_ASSERT(1);
}

void DATFileUpdaterTest::GlobalGridTest()
{
	std::string filename = "DAT_GlobalGrid";
	WriteGlobalGridDatFile(filename, 1.0);
	u_int32_t status1;
	FRegion region1;
	ParameterNameValueTypeVector params1;
	ReadGlobalGridData(filename, region1, status1, params1);
	DATData data(filename);
	AddGlobalGridData(data, 1.0);
	AddParameters(data);
	u_int32_t status2 = data.GetGlobalGridStatus();
	FRegion region2 = data.GetGlobalGrid();
	ParameterNameValueTypeVector params2;
	data.GetGridAlignmentAlgorithmParameters(params2);
	CPPUNIT_ASSERT(region1 == region2);
	CPPUNIT_ASSERT(status1 == status2);
	CompareParameterNameValueTypeVectors(params1, params2);
}

void DATFileUpdaterTest::GlobalGridIITest()
{
	std::string filename = "DAT_GlobalGrid_II";
	WriteGlobalGridDatFile(filename, 1.0);
	WriteGlobalGridDatFile(filename, 3.0);
	u_int32_t status1;
	FRegion region1;
	ParameterNameValueTypeVector params1;
	ReadGlobalGridData(filename, region1, status1, params1);
	DATData data(filename);
	AddGlobalGridData(data, 3.0);
	AddParameters(data);
	FRegion region2 = data.GetGlobalGrid();
	u_int32_t status2 = data.GetGlobalGridStatus();
	ParameterNameValueTypeVector params2;
	data.GetGridAlignmentAlgorithmParameters(params2);
	CPPUNIT_ASSERT(region1 == region2);
	CPPUNIT_ASSERT(status1 == status2);
	CompareParameterNameValueTypeVectors(params1, params2);
}

void DATFileUpdaterTest::GlobalGridIIITest()
{
	std::string filename = "DAT_GlobalGrid_III";
	WriteSubGridDatFile(filename, 1.0);
	WriteGlobalGridDatFile(filename, 3.0);
	u_int32_t status1;
	FRegion region1;
	ParameterNameValueTypeVector params1;
	ReadGlobalGridData(filename, region1, status1, params1);
	DATData data(filename);
	AddGlobalGridData(data, 3.0);
	AddParameters(data);
	FRegion region2 = data.GetGlobalGrid();
	u_int32_t status2 = data.GetGlobalGridStatus();
	ParameterNameValueTypeVector params2;
	data.GetGridAlignmentAlgorithmParameters(params2);
	CPPUNIT_ASSERT(region1 == region2);
	CPPUNIT_ASSERT(status1 == status2);
	CompareParameterNameValueTypeVectors(params1, params2);
}

void DATFileUpdaterTest::SubGridTest()
{
	std::string filename = "DAT_SubGrid";
	WriteSubGridDatFile(filename, 1.0);
	Uint32Vector status1;
	FRegionVector region1;
	ParameterNameValueTypeVector params1;
	ReadSubGridData(filename, region1, status1, params1);
	DATData data(filename);
	AddSubGridData(data, 1.0);
	AddParameters(data);
	Uint32Vector status2;
	FRegionVector region2;
	GetSubGridData(data, region2, status2);
	ParameterNameValueTypeVector params2;
	data.GetGridAlignmentAlgorithmParameters(params2);
	if (region1 == region2)
		CPPUNIT_ASSERT(region1 == region2);
	CPPUNIT_ASSERT(status1 == status2);
	CompareParameterNameValueTypeVectors(params1, params2);
}

void DATFileUpdaterTest::SubGridIITest()
{
	std::string filename = "DAT_SubGrid_II";
	WriteSubGridDatFile(filename, 1.0);
	WriteSubGridDatFile(filename, 3.0);
	Uint32Vector status1;
	FRegionVector region1;
	ParameterNameValueTypeVector params1;
	ReadSubGridData(filename, region1, status1, params1);
	DATData data(filename);
	AddSubGridData(data, 3.0);
	AddParameters(data);
	Uint32Vector status2;
	FRegionVector region2;
	GetSubGridData(data, region2, status2);
	ParameterNameValueTypeVector params2;
	data.GetGridAlignmentAlgorithmParameters(params2);
	CPPUNIT_ASSERT(region1 == region2);
	CPPUNIT_ASSERT(status1 == status2);
	CompareParameterNameValueTypeVectors(params1, params2);
}

void DATFileUpdaterTest::SubGridIIITest()
{
	std::string filename = "DAT_SubGrid_III";
	WriteSubGridDatFile(filename, 1.0);
	DATData data(filename);
	AddSubGridData(data, 7.0);
	DATFileUpdater updater(data);
	updater.Update();
	Uint32Vector status1;
	FRegionVector region1;
	ParameterNameValueTypeVector params1;
	ReadSubGridData(filename, region1, status1, params1);
	DATData data2(filename);
	AddSubGridData(data2, 7.0);
	AddParameters(data);
	Uint32Vector status2;
	FRegionVector region2;
	GetSubGridData(data2, region2, status2);
	ParameterNameValueTypeVector params2;
	data.GetGridAlignmentAlgorithmParameters(params2);
	CPPUNIT_ASSERT(region1 == region2);
	CPPUNIT_ASSERT(status1 == status2);
	CompareParameterNameValueTypeVectors(params1, params2);
}

void DATFileUpdaterTest::WriteDatFile(const std::string& filename)
{
	DATData data(filename);
	data.SetPixelCount(10);
	data.SetStatsCount(10);
	data.SetPixelSize(.7F);
	data.SetScannerID(L"main");
	data.SetScannerType(L"M10");
	data.SetArrayType(L"Hg-small");
	data.SetRows(4);
	data.SetCols(5);
	DATFileWriter* writer = new DATFileWriter(data);
	u_int16_t s[] = { 16, 22, 14, 39, 26, 36 };
	Uint16Vector stats;
	for(int i = 0; i < 6; i++)
	{
		stats.push_back(s[i]);
	}
	writer->WriteStats(stats);

	u_int16_t p[] = { 36, 3, 7, 8, 11, 2 };
	Uint16Vector pixels;
	for(int i = 0; i < 6; i++)
	{
		pixels.push_back(p[i]);
	}
	writer->WritePixels(pixels);
	delete writer;
}

bool DATFileUpdaterTest::IsEqual(const FRegionVector& r1, const FRegionVector& r2)
{
	if(r1.size() == r2.size())
	{
		size_t sz = r1.size();
		for(size_t i = 0; i < sz; i++)
		{
			if(r1[i] == r2[i]) continue;
			return false;
		}
		return true;
	}
	return false;
}

void DATFileUpdaterTest::GetSubGridData(const DATData& data, FRegionVector& regions, Uint32Vector& status)
{
	int32_t cnt = data.GetSubgridCnt();
	for(int i = 0; i < cnt; i++)
	{
		FRegion region = data.GetSubgrid(i);
		regions.push_back(region);
		status.push_back(data.GetSubgridStatus(i));
	}
}

void DATFileUpdaterTest::ReadGlobalGridData(const std::string& filename, FRegion& region, u_int32_t& status, ParameterNameValueTypeVector& params)
{
	DATFileReader reader;
	DATData data(filename);
	reader.Read(data);
	region = data.GetGlobalGrid();
	status = data.GetGlobalGridStatus();
	data.GetGridAlignmentAlgorithmParameters(params);
}

void DATFileUpdaterTest::ReadSubGridData(const std::string& filename, FRegionVector& regions, Uint32Vector& status, ParameterNameValueTypeVector& params)
{
	DATFileReader reader;
	DATData data(filename);
	reader.Read(data);
	GetSubGridData(data, regions, status);
	data.GetGridAlignmentAlgorithmParameters(params);
}

void DATFileUpdaterTest::WriteGlobalGridDatFile(const std::string& filename, float increment)
{
	WriteDatFile(filename);
	DATData data(filename);
	AddGlobalGridData(data, increment);
	AddParameters(data);
	DATFileUpdater updater(data);
	updater.Update();
}

void DATFileUpdaterTest::WriteSubGridDatFile(const std::string& filename, float increment)
{
	WriteGlobalGridDatFile(filename, increment);
	DATData data(filename);
	AddSubGridData(data, increment);
	AddParameters(data);
	DATFileUpdater updater(data);
	updater.Update();
}

void DATFileUpdaterTest::AddSubGridData(DATData& data, float increment)
{
	for(int n = 0; n < 4; n++)
	{
		FRegion region;
		float ptVal = 0.0;
		for(int i = 0; i < 4; i++)
		{
			FPoint point;
			point.x = ptVal;
			ptVal += increment;
			point.y = ptVal;
			region.pts.push_back(point);
		}

		u_int32_t status = DATData::GridError | DATData::GridManualAdjust;
		if (n % 2 == 0)
			status = DATData::GridOK | DATData::GridManualAdjust;

		data.AddSubgrid(status, region);
	}
}

void DATFileUpdaterTest::AddGlobalGridData(DATData& data, float increment)
{
	FRegion region;
	float ptVal = 0.0;
	for(int i = 0; i < 4; i++)
	{
		FPoint point;
		point.x = ptVal;
		ptVal += increment;
		point.y = ptVal;
		region.pts.push_back(point);
	}
	u_int32_t status = DATData::GridOK | DATData::GridManualAdjust;
	data.SetGlobalGrid(status, region);
}

void DATFileUpdaterTest::UpdateFileIdTest()
{
	std::string filename = "DAT_FileIdTest";
	WriteDatFile(filename);
	AffymetrixGuidType fileId = ReadFileId(filename);

	DATData data(filename);
	AddGlobalGridData(data, 1.0);
	
	{ // Need to force DATFileUpdater destructor between being able to read the file properly.
		DATFileUpdater updater(data);
		updater.Update();
	}
	AffymetrixGuidType newFileId = ReadFileId(filename);

	CPPUNIT_ASSERT(fileId != newFileId);
}

AffymetrixGuidType DATFileUpdaterTest::ReadFileId(const std::string& filename)
{
	DATFileReader reader;
	DATData data(filename);
	reader.Read(data);
	
	return data.GetFileHeader()->GetGenericDataHdr()->GetFileId();
}

void DATFileUpdaterTest::AddParameters(DATData& data)
{
	ParameterNameValueType nvt;
	nvt.SetName(L"Santa Clara");
	nvt.SetValueUInt8(56);
	data.AddGridAlignmentAlgorithmParameter(nvt);

	nvt.SetName(L"San Mateo");
	nvt.SetValueFloat(4.56f);
	data.AddGridAlignmentAlgorithmParameter(nvt);
}

void DATFileUpdaterTest::CompareParameterNameValueTypeVectors(ParameterNameValueTypeVector& left, ParameterNameValueTypeVector& right)
{
	CPPUNIT_ASSERT(left.size() == right.size());

	for (ParameterNameValueTypeVector::size_type i = 0; i < left.size(); ++i)
	{
		CPPUNIT_ASSERT(left[i].GetName() == right[i].GetName());
		CPPUNIT_ASSERT(left[i].GetParameterType() == right[i].GetParameterType());

		switch(left[i].GetParameterType())
		{
		case ParameterNameValueType::FloatType:
			CPPUNIT_ASSERT(left[i].GetValueFloat() == right[i].GetValueFloat());
			break;
		case ParameterNameValueType::UInt8Type:
			CPPUNIT_ASSERT(left[i].GetValueUInt8() == right[i].GetValueUInt8());
			break;
    default:
			CPPUNIT_ASSERT(0);
      break;
		}
	}
}
