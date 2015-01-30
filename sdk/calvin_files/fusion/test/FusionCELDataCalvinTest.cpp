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
#include "calvin_files/fusion/test/FusionCELDataCalvinTest.h"
//
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/exception/src/DevelopmentException.h"
#include "calvin_files/parameter/src/CELAlgorithmParameterNames.h"
#include "calvin_files/parsers/src/CelFileReader.h"
#include "calvin_files/parsers/src/FileException.h"
//
#include <cstring>
#include <string>
//
//#include "StringUtils.h"

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_exceptions;

const std::string EMPTY_CEL_FILE = "../data/empty.cel";
const std::string CALVIN_CEL_FILE = "../data/small_cel_file_partial_datheader";
const std::string CALVIN_CEL_FILE_FROM_GCOS = "../data/small_cel_file_full_datheader";
const std::string CALVIN_CEL_FILE_MISSING_GRID_AND_MARGIN_AND_PARENT_DAT = "../data/small_cel_file";
const std::string NON_EXISTANT_FILE = "../data/doesnt_exist";
const std::string QC_FILE = "../data/exp_qc.CEL";

CPPUNIT_TEST_SUITE_REGISTRATION(FusionCELDataCalvinTest);

void FusionCELDataCalvinTest::setUp()
{
	fusionCel = new FusionCELData;
	calvinCel = new CelFileData;
}

void FusionCELDataCalvinTest::tearDown()
{
	delete fusionCel;
	delete calvinCel;
}

void FusionCELDataCalvinTest::GetSetParameters()
{
	fusionCel->SetFileName(QC_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	ParameterNameValueTypeList params;
	
	params = fusionCel->GetDataSetParameters(L"QC Analysis Parameters");
	CPPUNIT_ASSERT(params.size() == 5);

	ParameterNameValueTypeList::iterator it = params.begin();
	it++;
	it++;
	it++;
	it++;
	ParameterNameValueType &param = *it;

	CPPUNIT_ASSERT(param.GetName() == L"Array Metric");
	CPPUNIT_ASSERT(param.GetValueAscii() == "pm_mean");

	params = fusionCel->GetDataSetParameters(L"QC Analysis Results");
	CPPUNIT_ASSERT(params.size() > 0);
}

void FusionCELDataCalvinTest::CreationTest()
{
	CPPUNIT_ASSERT(fusionCel);
}

void FusionCELDataCalvinTest::EmptyTest()
{
	fusionCel->SetFileName(EMPTY_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == false);
}

void FusionCELDataCalvinTest::GetFileIdTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetFileId() == "0000065535-1147280233-0000005844-0000011008-0000006224");
	CPPUNIT_ASSERT(fusionCel->GetGenericData()->FileIdentifier() == "0000065535-1147280233-0000005844-0000011008-0000006224");
}

void FusionCELDataCalvinTest::FileNameTest()
{
	fusionCel->SetFileName("Merced");
	CPPUNIT_ASSERT(fusionCel->GetFileName() == "Merced");
}

void FusionCELDataCalvinTest::ExistsTest()
{
	fusionCel->SetFileName(NON_EXISTANT_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Exists() == false);
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Exists() == true);
}

void FusionCELDataCalvinTest::ReadHeaderTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
}

//void FusionCELDataCalvinTest::DimensionsTest()
//{
//	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
//	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
//	fusionCel->SetDimensions(3,4);
//	CPPUNIT_ASSERT(fusionCel->GetRows() == 3);
//	CPPUNIT_ASSERT(fusionCel->GetCols() == 4);
//}

void FusionCELDataCalvinTest::AlgorithmNameTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetAlg() == L"Feature Extraction");
}

void FusionCELDataCalvinTest::ChipTypeTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetChipType() == L"Hg-small");
}

void FusionCELDataCalvinTest::LibPackageTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetLibraryPackageName() == L"Hg-small-lib-package");
}

void FusionCELDataCalvinTest::MasterFileTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetMasterFileName() == L"Hg-small-master-file");
}

void FusionCELDataCalvinTest::MarginTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetCellMargin() == 2);
}

void FusionCELDataCalvinTest::MarginMissingTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE_MISSING_GRID_AND_MARGIN_AND_PARENT_DAT.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetCellMargin() == 0);
}

void FusionCELDataCalvinTest::HeaderTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetHeader() == L"");	// header is not stored in calvin files
}

void FusionCELDataCalvinTest::ParamsTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	std::wstring s = fusionCel->GetParams();
	CPPUNIT_ASSERT(s.find_last_of(L"GridLLY:9.000000") != std::wstring::npos);
}

void FusionCELDataCalvinTest::GetAlgorithmParametersTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	std::wstring s = fusionCel->GetAlgorithmParameters();
	CPPUNIT_ASSERT(s.find_last_of(L"GridLLY:9.000000") != std::wstring::npos);
}

void FusionCELDataCalvinTest::GetNumberAlgorithmParametersTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	CPPUNIT_ASSERT(fusionCel->GetNumberAlgorithmParameters() == 11);
}

void FusionCELDataCalvinTest::GetAlgorithmParameterTagTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	int num = fusionCel->GetNumberAlgorithmParameters();
	CPPUNIT_ASSERT(fusionCel->GetAlgorithmParameterTag(num-1) == L"GridLLY");
	CPPUNIT_ASSERT(fusionCel->GetAlgorithmParameterTag(num) == L"");	// out-of-bounds test
}

void FusionCELDataCalvinTest::GetAlgorithmParameterTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	CPPUNIT_ASSERT(fusionCel->GetAlgorithmParameter(L"GridLLY") == L"9.000000");
	CPPUNIT_ASSERT(fusionCel->GetAlgorithmParameter(L"bogus") == L"");	// out-of-bounds test
}

void FusionCELDataCalvinTest::ParametersTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	FusionTagValuePairTypeList& list = fusionCel->GetParameters();
	std::list<FusionTagValuePairType>::reverse_iterator ri = list.rbegin();
	CPPUNIT_ASSERT(ri->Tag == L"GridLLY");
	CPPUNIT_ASSERT(ri->Value == L"9.000000");
}

void FusionCELDataCalvinTest::ReadHeaderAndVerifyTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	CompareHeaderData();
	ComparePartialDatHeaders();
}

void FusionCELDataCalvinTest::CompareHeaderData()
{
	// Check all header values
	CPPUNIT_ASSERT(fusionCel->GetVersion() == calvinCel->GetVersion());
	CPPUNIT_ASSERT(fusionCel->GetCols() == calvinCel->GetCols());
	CPPUNIT_ASSERT(fusionCel->GetRows() == calvinCel->GetRows());
	CPPUNIT_ASSERT(fusionCel->GetNumCells() == calvinCel->GetNumCells());

	CPPUNIT_ASSERT(fusionCel->GetHeader() == L"");
	CPPUNIT_ASSERT(fusionCel->GetAlg() == calvinCel->GetAlgorithmName());
	CPPUNIT_ASSERT(fusionCel->GetChipType() == calvinCel->GetArrayType());

	ParameterNameValueTypeVector algParams;
	calvinCel->GetAlgorithmParameters(algParams);

	// Test GetParams and GetAlgorithmParameters
	std::wstring calvinParamString;
	for (ParameterNameValueTypeIt ci = algParams.begin(); ci != algParams.end(); ++ci)
	{
		// This is a repeat of the calvin adapter implementation, oh well, a file change will not break it
		if (algParams.begin() != ci)
			calvinParamString += L";";
		calvinParamString += ci->GetName();
		calvinParamString += L":";
		calvinParamString += ci->ToString();
	}
	CPPUNIT_ASSERT(calvinParamString == fusionCel->GetParams());
	CPPUNIT_ASSERT(calvinParamString == fusionCel->GetAlgorithmParameters());

	CPPUNIT_ASSERT(fusionCel->GetNumberAlgorithmParameters() == (int)algParams.size());
	int paramCount = fusionCel->GetNumberAlgorithmParameters();
	ParameterNameValueTypeIt ci = algParams.begin();
	for (int ip = 0; ip < paramCount && ci != algParams.end(); ++ip, ++ci)
	{
		wstring fusionParamTag = fusionCel->GetAlgorithmParameterTag(ip);
		CPPUNIT_ASSERT(ci->GetName() == fusionParamTag);
		CPPUNIT_ASSERT(ci->ToString() == fusionCel->GetAlgorithmParameter(fusionParamTag.c_str()));
	}

	// Test GetParameters
	FusionTagValuePairTypeList& fusionList = fusionCel->GetParameters();
	std::list<FusionTagValuePairType>::iterator fi = fusionList.begin();

	int32_t pi = 0;
	int32_t algParamCnt = (int32_t)algParams.size();

	for (; fi != fusionList.end() && pi < algParamCnt; ++fi, ++pi)
	{
		std::wstring t1 = fi->Tag;
		std::wstring v1 = fi->Value;
		std::wstring t2 = algParams.at(pi).GetName();
		std::wstring v2 = algParams.at(pi).ToString();
		CPPUNIT_ASSERT(fi->Tag == algParams.at(pi).GetName());
		CPPUNIT_ASSERT(fi->Value == algParams.at(pi).ToString());
		if (fi->DetailedType().GetParameterType() == ParameterNameValueType::FloatType)
			CPPUNIT_ASSERT(fi->DetailedType().GetValueFloat() == algParams.at(pi).GetValueFloat());
		else if (fi->DetailedType().GetParameterType() == ParameterNameValueType::Int32Type)
			CPPUNIT_ASSERT(fi->DetailedType().GetValueInt32() == algParams.at(pi).GetValueInt32());
	}

	CPPUNIT_ASSERT(fusionCel->GetCellMargin() == 2);

	XYCoordVector outlierCoords;
	calvinCel->GetOutlierCoords(outlierCoords);
	CPPUNIT_ASSERT(fusionCel->GetNumOutliers() == (int)outlierCoords.size());

	XYCoordVector maskCoords;
	calvinCel->GetMaskedCoords(maskCoords);
	CPPUNIT_ASSERT(fusionCel->GetNumMasked() == (int)maskCoords.size());

	CompareGrids();
}

void FusionCELDataCalvinTest::ComparePartialDatHeaders()
{
	// Find the parent DAT
	GenDataHdrVectorIt begin, end; 
	calvinCel->GetFileHeader()->GetGenericDataHdr()->GetParentIterators(begin, end);
	std::wstring datHeader;
	u_int16_t maxPixel = 0, minPixel = 0;

	for (GenDataHdrVectorIt ii = begin; ii != end; ++ii)
	{
		if (ii->GetFileTypeId() == SCAN_ACQUISITION_DATA_TYPE)
		{
			// found the right header, now look for the parameter
			ParameterNameValueType nvt;
			if (ii->FindNameValParam(DAT_HEADER_PARAM_NAME, nvt))
			{
				datHeader = nvt.GetValueText();
			}
			if (ii->FindNameValParam(MAX_PIXEL_INTENSITY_PARAM_NAME, nvt))
			{
				if (nvt.GetParameterType() == ParameterNameValueType::UInt16Type)
					maxPixel = nvt.GetValueUInt16();
			}

			if (ii->FindNameValParam(MIN_PIXEL_INTENSITY_PARAM_NAME, nvt))
			{
				if (nvt.GetParameterType() == ParameterNameValueType::UInt16Type)
					minPixel = nvt.GetValueUInt16();
			}
		}
	}

	u_int32_t expectedMax = 0, expectedMin = 0;	// Needs to be a 4 byte int to keep swscanf from messing up the stack frame.
	std::wstring fusionDatHeader = fusionCel->GetDatHeader();
	CPPUNIT_ASSERT(fusionDatHeader.find(datHeader, 0) != wstring::npos);
	CPPUNIT_ASSERT(swscanf(fusionDatHeader.c_str(), L"[%d..%d]", &expectedMin, &expectedMax) == 2);
	CPPUNIT_ASSERT(minPixel == expectedMin);
	CPPUNIT_ASSERT(maxPixel == expectedMax);

}

void FusionCELDataCalvinTest::CompareGrids()
{
	affymetrix_fusion_io::FGridCoords fusionGrid = fusionCel->GetGridCorners();

	ParameterNameValueType nvt;
	CPPUNIT_ASSERT(calvinCel->FindAlgorithmParameter(GRIDULX_PARAM_NAME, nvt));
	CPPUNIT_ASSERT(nvt.GetValueFloat() == fusionGrid.upperleft.x);
	CPPUNIT_ASSERT(calvinCel->FindAlgorithmParameter(GRIDULY_PARAM_NAME, nvt));
	CPPUNIT_ASSERT(nvt.GetValueFloat() == fusionGrid.upperleft.y);

	CPPUNIT_ASSERT(calvinCel->FindAlgorithmParameter(GRIDURX_PARAM_NAME, nvt));
	CPPUNIT_ASSERT(nvt.GetValueFloat() == fusionGrid.upperright.x);
	CPPUNIT_ASSERT(calvinCel->FindAlgorithmParameter(GRIDURY_PARAM_NAME, nvt));
	CPPUNIT_ASSERT(nvt.GetValueFloat() == fusionGrid.upperright.y);

	CPPUNIT_ASSERT(calvinCel->FindAlgorithmParameter(GRIDLRX_PARAM_NAME, nvt));
	CPPUNIT_ASSERT(nvt.GetValueFloat() == fusionGrid.lowerright.x);
	CPPUNIT_ASSERT(calvinCel->FindAlgorithmParameter(GRIDLRY_PARAM_NAME, nvt));
	CPPUNIT_ASSERT(nvt.GetValueFloat() == fusionGrid.lowerright.y);

	CPPUNIT_ASSERT(calvinCel->FindAlgorithmParameter(GRIDLLX_PARAM_NAME, nvt));
	CPPUNIT_ASSERT(nvt.GetValueFloat() == fusionGrid.lowerleft.x);
	CPPUNIT_ASSERT(calvinCel->FindAlgorithmParameter(GRIDLLY_PARAM_NAME, nvt));
	CPPUNIT_ASSERT(nvt.GetValueFloat() == fusionGrid.lowerleft.y);
}

void FusionCELDataCalvinTest::ReadTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	CompareHeaderData();
	ComparePartialDatHeaders();
}

void FusionCELDataCalvinTest::CloseTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);
	CPPUNIT_ASSERT_NO_THROW(fusionCel->Close());
}

void FusionCELDataCalvinTest::GetIntensityIndexTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	int32_t count = calvinCel->GetNumCells();
	FloatVector inten;
	calvinCel->GetIntensities(0, count, inten);

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = 1;

	for (int i = 0; i < cells; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->GetIntensity(i) == inten[i]);
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::GetIntensityXYTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	int32_t count = calvinCel->GetNumCells();
	FloatVector inten;
	calvinCel->GetIntensities(0, count, inten);

	// Spot check the values
	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			int index = fusionCel->XYToIndex(col, row);
			CPPUNIT_ASSERT(fusionCel->GetIntensity(col, row) == inten[index]);
		}
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::GetStdvIndexTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	int32_t count = calvinCel->GetNumCells();
	FloatVector stdev;
	calvinCel->GetStdev(0, count, stdev);

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = 1;

	for (int i = 0; i < cells; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->GetStdv(i) == stdev[i]);
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::GetStdvXYTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	int32_t count = calvinCel->GetNumCells();
	FloatVector stdev;
	calvinCel->GetStdev(0, count, stdev);

	// Spot check the values
	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			int index = fusionCel->XYToIndex(col, row);
			CPPUNIT_ASSERT(fusionCel->GetStdv(col, row) == stdev[index]);
		}
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::GetPixelsIndexTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	int32_t count = calvinCel->GetNumCells();
	Int16Vector pixels;
	calvinCel->GetNumPixels(0, count, pixels);

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = 1;

	for (int i = 0; i < cells; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->GetPixels(i) == pixels[i]);
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::GetPixelsXYTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	int32_t count = calvinCel->GetNumCells();
	Int16Vector pixels;
	calvinCel->GetNumPixels(0, count, pixels);

	// Spot check the values
	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			int index = fusionCel->XYToIndex(col, row);
			CPPUNIT_ASSERT(fusionCel->GetPixels(col, row) == pixels[index]);
		}
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::IsOutlierIndexTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = 1;

	for (int i = 0; i < cells; i+= incr)
	{
		BoolVector v;
		calvinCel->GetOutliers(i, 1, v);
		CPPUNIT_ASSERT(fusionCel->IsOutlier(i) == v.at(0));
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::IsOutlierXYTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	// Spot check the values
	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			BoolVector v;
			int index = fusionCel->XYToIndex(col, row);
			calvinCel->GetOutliers(index, 1, v);
			CPPUNIT_ASSERT(fusionCel->IsOutlier(col, row) == v.at(0));
		}
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::IsMaskedIndexTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = 1;

	for (int i = 0; i < cells; i+= incr)
	{
		BoolVector v;
		calvinCel->GetMasked(i, 1, v);
		CPPUNIT_ASSERT(fusionCel->IsMasked(i) == v.at(0));
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::IsMaskedXYTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	// Spot check the values
	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			BoolVector v;
			int index = fusionCel->XYToIndex(col, row);
			calvinCel->GetMasked(index, 1, v);
			CPPUNIT_ASSERT(fusionCel->IsMasked(col, row) == v.at(0));
		}
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::IndexToXTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	int cells = fusionCel->GetNumCells();
	int incr = 1;

	for (int i = 0; i < cells; i+= incr)
	{
		int x = i % calvinCel->GetCols();
		CPPUNIT_ASSERT(fusionCel->IndexToX(i) == x);
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::IndexToYTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	int cells = fusionCel->GetNumCells();
	int incr = 1;

	for (int i = 0; i < cells; i+= incr)
	{
		int y = i / calvinCel->GetCols();
		CPPUNIT_ASSERT(fusionCel->IndexToY(i) == y);
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::XYToIndexTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			int index = row*calvinCel->GetCols() + col;
			CPPUNIT_ASSERT(fusionCel->XYToIndex(col, row) == index);
		}
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::ClearTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);
	CPPUNIT_ASSERT_NO_THROW(fusionCel->Clear());

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));
	calvinCel->Clear();

	// Check some values
	CPPUNIT_ASSERT(fusionCel->GetRows() == calvinCel->GetRows());
	CPPUNIT_ASSERT(fusionCel->GetNumCells() == calvinCel->GetNumCells());
	CPPUNIT_ASSERT(fusionCel->GetFileName() == calvinCel->GetFilename());	// semantics of GetFileName changed in the fusion layer
}

void FusionCELDataCalvinTest::GetEntryIndexTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = 1;

	FusionCELFileEntryType calvinEntry;

	float intensity;
	float stdev;
	int16_t numPixels; 
	bool outlier, masked;
	for (int i = 0; i < cells; i+= incr)
	{
		fusionCel->GetEntry(i, calvinEntry);
		calvinCel->GetData(i, intensity, stdev, numPixels, outlier, masked);

		CPPUNIT_ASSERT(calvinEntry.Intensity == intensity);
		CPPUNIT_ASSERT(calvinEntry.Stdv == stdev);
		CPPUNIT_ASSERT(calvinEntry.Pixels == numPixels);
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::GetEntryXYTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	// Spot check the values
	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	FusionCELFileEntryType calvinEntry;

	float intensity;
	float stdev;
	int16_t numPixels; 
	bool outlier, masked;
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			int index = row*calvinCel->GetCols() + col;
			fusionCel->GetEntry(col, row, calvinEntry);
			calvinCel->GetData(index, intensity, stdev, numPixels, outlier, masked);

			CPPUNIT_ASSERT(calvinEntry.Intensity == intensity);
			CPPUNIT_ASSERT(calvinEntry.Stdv == stdev);
			CPPUNIT_ASSERT(calvinEntry.Pixels == numPixels);
		}
	}
	
	fusionCel->Close();
}

void FusionCELDataCalvinTest::CallMethodsWhenObjectIsNotReadyTest()
{
//	CPPUNIT_ASSERT_THROW(fusionCel->ReadHeader(), FileNotOpenException);

	CPPUNIT_ASSERT_THROW(fusionCel->GetVersion(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetRows(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetCols(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetNumCells(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetHeader(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetAlg(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetParams(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetChipType(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetCellMargin(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetNumOutliers(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetNumMasked(), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetParameters(), FileNotOpenException);
//	CPPUNIT_ASSERT_THROW(fusionCel->AddAlgorithmParameter(L"", L""), FileNotOpenException);
//	CPPUNIT_ASSERT_THROW(fusionCel->SetDimensions(0,0), FileNotOpenException);
//	CPPUNIT_ASSERT_THROW(fusionCel->SetAlgorithmName(L""), FileNotOpenException);
//	CPPUNIT_ASSERT_THROW(fusionCel->SetChipType(L""), FileNotOpenException);
//	CPPUNIT_ASSERT_THROW(fusionCel->SetMargin(0), FileNotOpenException);

	FusionCELFileEntryType entry;
	CPPUNIT_ASSERT_THROW(fusionCel->GetEntry(0, entry), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetEntry(0,0, entry), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetIntensity(0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetIntensity(0,0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetStdv(0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetStdv(0,0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetPixels(0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->GetPixels(0,0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->IsOutlier(0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->IsOutlier(0,0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->IsMasked(0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->IsMasked(0,0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->IndexToX(0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->IndexToY(0), FileNotOpenException);
	CPPUNIT_ASSERT_THROW(fusionCel->XYToIndex(0,0), FileNotOpenException);
	//CPPUNIT_ASSERT_THROW(fusionCel->Close(), FileNotOpenException);
	CPPUNIT_ASSERT_NO_THROW(fusionCel->Clear());
}

void FusionCELDataCalvinTest::ErrorTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	
	fusionCel->SetError(L"test error");
	CPPUNIT_ASSERT(fusionCel->GetError() == L"");	// calvin will return a blank error.
}

void FusionCELDataCalvinTest::GetFileSizeTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->GetFileSize() == 4190);
}

void FusionCELDataCalvinTest::ReadExTest()
{
	CPPUNIT_ASSERT(fusionCel->ReadEx(CALVIN_CEL_FILE.c_str(), FusionCELData::CEL_DATA));

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	// Check some header values
	CPPUNIT_ASSERT(fusionCel->GetVersion() == calvinCel->GetVersion());
	CPPUNIT_ASSERT(fusionCel->GetCols() == calvinCel->GetCols());
	CPPUNIT_ASSERT(fusionCel->GetRows() == calvinCel->GetRows());
	CPPUNIT_ASSERT(fusionCel->GetNumCells() == calvinCel->GetNumCells());

	CPPUNIT_ASSERT(fusionCel->GetHeader() == L"");
	CPPUNIT_ASSERT(fusionCel->GetAlg() == calvinCel->GetAlgorithmName());
	CPPUNIT_ASSERT(fusionCel->GetChipType() == calvinCel->GetArrayType());
}

void FusionCELDataCalvinTest::GetReadStateTest()
{
	CPPUNIT_ASSERT(fusionCel->ReadEx(CALVIN_CEL_FILE.c_str(), FusionCELData::CEL_DATA));
	CPPUNIT_ASSERT(fusionCel->GetReadState() == FusionCELData::CEL_ALL);	// always CEL_ALL
}

void FusionCELDataCalvinTest::GetHeaderKeyTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	CPPUNIT_ASSERT(fusionCel->GetHeaderKey(L"anything") == L"");
}

void FusionCELDataCalvinTest::GetDatHeaderTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	std::wstring datHeaderString = L"[1..46001]  small_cel_file_partial_datheader:CLS=25   RWS=25   XIN=1  YIN=1  VE=0         0   05/19/05 02:45:59 ScannerID:  ScannerTyp   \x14  \x14 Hg-Small.1sq \x14  \x14  \x14  \x14  \x14 570 \x14 45.200001 \x14 0.340000 \x14 1.0900 \x14 3";
	std::wstring s= fusionCel->GetDatHeader();
	CPPUNIT_ASSERT(fusionCel->GetDatHeader() == datHeaderString);
}

void FusionCELDataCalvinTest::GetDatHeaderFromGCOSMigratedCELTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE_FROM_GCOS.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	std::wstring datHeaderString = L"[45..56789]  small_cel_file_full_datheader:CLS=25   RWS=25   XIN=1  YIN=1  VE=0         0   05/19/05 02:45:59 ScannerID:  ScannerTyp   \x14  \x14 Hg-Small.1sq \x14  \x14  \x14  \x14  \x14 570 \x14 45.200001 \x14 0.340000 \x14 1.0900 \x14 3";
	std::wstring s= fusionCel->GetDatHeader();
	CPPUNIT_ASSERT(fusionCel->GetDatHeader() == datHeaderString);
}

void FusionCELDataCalvinTest::GetDatHeaderFailedTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE_MISSING_GRID_AND_MARGIN_AND_PARENT_DAT.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	CPPUNIT_ASSERT(fusionCel->GetDatHeader() == L"");
}

void FusionCELDataCalvinTest::GetGridCornersTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	CelFileReader reader;
	reader.SetFilename(CALVIN_CEL_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(*calvinCel));

	CompareGrids();
}

void FusionCELDataCalvinTest::GetGridCornersFailedTest()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE_MISSING_GRID_AND_MARGIN_AND_PARENT_DAT.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	affymetrix_fusion_io::FGridCoords emptyRect;
	affymetrix_fusion_io::FGridCoords rect = fusionCel->GetGridCorners();

	CPPUNIT_ASSERT(rect.IsEmpty());
}

void FusionCELDataCalvinTest::TestMultiChannel()
{
	fusionCel->SetFileName("../data/multi-data.cel");
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->IsMultiColor() == true);
	WStringVector c = fusionCel->GetChannels();
	CPPUNIT_ASSERT(c.size() == 2);
	CPPUNIT_ASSERT(c[0] == L"531");
	CPPUNIT_ASSERT(c[1] == L"609");
}

void FusionCELDataCalvinTest::TestSingleChannel()
{
	fusionCel->SetFileName(CALVIN_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->IsMultiColor() == false);
	WStringVector c = fusionCel->GetChannels();
	CPPUNIT_ASSERT(c.size() == 1);
	CPPUNIT_ASSERT(c[0] == L"Default Group");
}
