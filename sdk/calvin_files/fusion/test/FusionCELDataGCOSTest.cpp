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
#include "calvin_files/fusion/test/FusionCELDataGCOSTest.h"
//
#include "calvin_files/exception/src/DevelopmentException.h"
#include "calvin_files/fusion/src/FusionCoords.h"
#include "calvin_files/parsers/src/FileException.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affxcel;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_utilities;

// Underlying CCelFileData class expects backslashes rather than forward slashes on Windows hmmmm
const std::string GCOS_CEL_FILE = "../data/ethan1-1.CEL";
const std::string NON_EXISTANT_FILE = "../data/doesnt_exist";

CPPUNIT_TEST_SUITE_REGISTRATION(FusionCELDataGCOSTest);

void FusionCELDataGCOSTest::setUp()
{
	fusionCel = new FusionCELData;
	gcosCel = new CCELFileData;
}

void FusionCELDataGCOSTest::tearDown()
{
	delete fusionCel;
	delete gcosCel;
}

void FusionCELDataGCOSTest::TestSingleChannel()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->IsMultiColor() == false);
	WStringVector c = fusionCel->GetChannels();
	CPPUNIT_ASSERT(c.size() == 0);
}

void FusionCELDataGCOSTest::CreationTest()
{
	CPPUNIT_ASSERT(fusionCel);
}

void FusionCELDataGCOSTest::FileNameTest()
{
	fusionCel->SetFileName("Merced");
	CPPUNIT_ASSERT(fusionCel->GetFileName() == "Merced");
}

void FusionCELDataGCOSTest::GetFileIdTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetFileId() == "");
	CPPUNIT_ASSERT(fusionCel->GetGenericData() == NULL);
}

void FusionCELDataGCOSTest::ExistsTest()
{
	fusionCel->SetFileName(NON_EXISTANT_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Exists() == false);
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Exists() == true);
}

void FusionCELDataGCOSTest::ReadHeaderTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
}

//void FusionCELDataGCOSTest::DimensionsTest()
//{
//	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
//	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
//	fusionCel->SetDimensions(3,4);
//	CPPUNIT_ASSERT(fusionCel->GetRows() == 3);
//	CPPUNIT_ASSERT(fusionCel->GetCols() == 4);
//}

void FusionCELDataGCOSTest::AlgorithmNameTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetAlg() == L"Percentile");
}

void FusionCELDataGCOSTest::ChipTypeTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetChipType() == L"HG-U133A");
}

void FusionCELDataGCOSTest::LibPackageTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetLibraryPackageName() == L"");
}

void FusionCELDataGCOSTest::MasterFileTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetMasterFileName() == L"");
}

void FusionCELDataGCOSTest::MarginTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	CPPUNIT_ASSERT(fusionCel->GetCellMargin() == 2);
}

void FusionCELDataGCOSTest::ParamsTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	std::wstring s = fusionCel->GetParams();
	CPPUNIT_ASSERT(s.find_last_of(L"StdMult:1.000000") != std::wstring::npos);
}

void FusionCELDataGCOSTest::GetAlgorithmParametersTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	std::wstring s = fusionCel->GetAlgorithmParameters();
	CPPUNIT_ASSERT(s.find_last_of(L"StdMult:1.000000") != std::wstring::npos);
}

void FusionCELDataGCOSTest::GetNumberAlgorithmParametersTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	CPPUNIT_ASSERT(fusionCel->GetNumberAlgorithmParameters() == 12);
}

void FusionCELDataGCOSTest::GetAlgorithmParameterTagTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	int num = fusionCel->GetNumberAlgorithmParameters();
	CPPUNIT_ASSERT(fusionCel->GetAlgorithmParameterTag(num-1) == L"StdMult");
	CPPUNIT_ASSERT(fusionCel->GetAlgorithmParameterTag(num) == L"");	// out-of-bounds test
}

void FusionCELDataGCOSTest::GetAlgorithmParameterTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	CPPUNIT_ASSERT(fusionCel->GetAlgorithmParameter(L"OutlierLow") == L"1.004");
	CPPUNIT_ASSERT(fusionCel->GetAlgorithmParameter(L"bogus") == L"");	// out-of-bounds test
}

void FusionCELDataGCOSTest::ParametersTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	FusionTagValuePairTypeList& list = fusionCel->GetParameters();
	std::list<FusionTagValuePairType>::iterator ii = list.begin();
	CPPUNIT_ASSERT(ii->Tag == L"Percentile");
	CPPUNIT_ASSERT(ii->Value == L"75");
}

void FusionCELDataGCOSTest::ReadHeaderAndVerifyTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->ReadHeader());

	CompareHeaderData();
}

void FusionCELDataGCOSTest::CompareHeaderData()
{
	// Check all header values
	CPPUNIT_ASSERT(fusionCel->GetVersion() == gcosCel->GetVersion());
	CPPUNIT_ASSERT(fusionCel->GetCols() == gcosCel->GetCols());
	CPPUNIT_ASSERT(fusionCel->GetRows() == gcosCel->GetRows());
	CPPUNIT_ASSERT(fusionCel->GetNumCells() == gcosCel->GetNumCells());

	std::string headerString = StringUtils::ConvertWCSToMBS(fusionCel->GetHeader());
	CPPUNIT_ASSERT(headerString == gcosCel->GetHeaderString());
	std::string algString = StringUtils::ConvertWCSToMBS(fusionCel->GetAlg());
	CPPUNIT_ASSERT(algString == gcosCel->GetAlg());
	std::string paramString = StringUtils::ConvertWCSToMBS(fusionCel->GetParams());
	CPPUNIT_ASSERT(paramString == gcosCel->GetParams());
	std::string chipTypetring = StringUtils::ConvertWCSToMBS(fusionCel->GetChipType());
	CPPUNIT_ASSERT(chipTypetring == gcosCel->GetChipType());

	std::string params = StringUtils::ConvertWCSToMBS(fusionCel->GetParams());
	CPPUNIT_ASSERT(gcosCel->GetParams() == params);
	std::string algorithmParameters = StringUtils::ConvertWCSToMBS(fusionCel->GetAlgorithmParameters());
	CPPUNIT_ASSERT(gcosCel->GetAlgorithmParameters() == algorithmParameters);
	CPPUNIT_ASSERT(fusionCel->GetNumberAlgorithmParameters() == gcosCel->GetNumberAlgorithmParameters());
	int paramCount = fusionCel->GetNumberAlgorithmParameters();
	for (int ip = 0; ip < paramCount; ++ip)
	{
		std::wstring fusionParamTag = fusionCel->GetAlgorithmParameterTag(ip);
		std::string gcosParamTag = gcosCel->GetAlgorithmParameterTag(ip);
		CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(fusionParamTag) == gcosParamTag);
		std::string fusionParamValue = StringUtils::ConvertWCSToMBS(fusionCel->GetAlgorithmParameter(fusionParamTag.c_str()));
		CPPUNIT_ASSERT(fusionParamValue == gcosCel->GetAlgorithmParameter(gcosParamTag.c_str()));
	}

	// Test GetParameters
	FusionTagValuePairTypeList& calvinList = fusionCel->GetParameters();
	int gcosParamCount = gcosCel->GetNumberAlgorithmParameters();

	std::list<FusionTagValuePairType>::iterator ci = calvinList.begin();

	for (int gi = 0; ci != calvinList.end() && gi < gcosParamCount; ++ci, ++gi)
	{
		std::string calvinTag = StringUtils::ConvertWCSToMBS(ci->Tag);
		CPPUNIT_ASSERT(calvinTag == gcosCel->GetAlgorithmParameterTag(gi));
		std::string calvinValue = StringUtils::ConvertWCSToMBS(ci->Value);
		CPPUNIT_ASSERT(calvinValue == gcosCel->GetAlgorithmParameter(gcosCel->GetAlgorithmParameterTag(gi).c_str()));
		CPPUNIT_ASSERT(ci->Tag == ci->DetailedType().GetName());
		CPPUNIT_ASSERT(ci->Value == ci->DetailedType().GetValueText());
	}

	CPPUNIT_ASSERT(fusionCel->GetCellMargin() == gcosCel->GetCellMargin());
	CPPUNIT_ASSERT(fusionCel->GetNumOutliers() == gcosCel->GetNumOutliers());
	CPPUNIT_ASSERT(fusionCel->GetNumMasked() == gcosCel->GetNumMasked());
	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(fusionCel->GetDatHeader()) == gcosCel->GetDatHeader());
	
	CompareGridCoordinates();
}

void FusionCELDataGCOSTest::CompareGridCoordinates()
{
	// Check the grid
	affymetrix_fusion_io::FGridCoords fusionGrid = fusionCel->GetGridCorners();
	GridCoordinatesType gcosRect = gcosCel->GetGridCorners();
	CPPUNIT_ASSERT(fusionGrid.upperleft.x = (float)gcosRect.upperleft.x);
	CPPUNIT_ASSERT(fusionGrid.upperleft.y = (float)gcosRect.upperleft.y);
	CPPUNIT_ASSERT(fusionGrid.upperright.x = (float)gcosRect.upperright.x);
	CPPUNIT_ASSERT(fusionGrid.upperright.y = (float)gcosRect.upperright.y);
	CPPUNIT_ASSERT(fusionGrid.lowerright.x = (float)gcosRect.lowerright.x);
	CPPUNIT_ASSERT(fusionGrid.lowerright.y = (float)gcosRect.lowerright.y);
	CPPUNIT_ASSERT(fusionGrid.lowerleft.x = (float)gcosRect.lowerleft.x);
	CPPUNIT_ASSERT(fusionGrid.lowerleft.y = (float)gcosRect.lowerleft.y);
}

void FusionCELDataGCOSTest::ReadTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->ReadHeader());

	CompareHeaderData();
}

void FusionCELDataGCOSTest::CloseTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);
	CPPUNIT_ASSERT_NO_THROW(fusionCel->Close());
}

void FusionCELDataGCOSTest::GetIntensityIndexTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = cells / 10;

	for (int i = 0; i < 10; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->GetIntensity(i) == gcosCel->GetIntensity(i));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::GetIntensityXYTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int rows = fusionCel->GetRows();

	for (int row = 0; row < rows; ++row)
	{
		CPPUNIT_ASSERT(fusionCel->GetIntensity(0, row) == gcosCel->GetIntensity(0, row));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::GetStdvIndexTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = cells / 10;

	for (int i = 0; i < cells; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->GetStdv(i) == gcosCel->GetStdv(i));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::GetStdvXYTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	for (int row = 0; row < rows; ++row)
	{
		CPPUNIT_ASSERT(fusionCel->GetStdv(cols-1, row) == gcosCel->GetStdv(cols-1, row));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::GetPixelsIndexTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = cells / 10;

	for (int i = 0; i < cells; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->GetPixels(i) == gcosCel->GetPixels(i));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::GetPixelsXYTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int rows = fusionCel->GetRows();

	for (int row = 0; row < rows; ++row)
	{
		CPPUNIT_ASSERT(fusionCel->GetPixels(0, row) == gcosCel->GetPixels(0, row));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::IsOutlierIndexTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = cells / 10;

	for (int i = 0; i < cells; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->IsOutlier(i) == gcosCel->IsOutlier(i));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::IsOutlierXYTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	for (int row = 0; row < rows; ++row)
	{
		CPPUNIT_ASSERT(fusionCel->IsOutlier(cols-1, row) == gcosCel->IsOutlier(cols-1, row));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::IsMaskedIndexTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = cells / 10;

	for (int i = 0; i < cells; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->IsMasked(i) == gcosCel->IsMasked(i));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::IsMaskedXYTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int rows = fusionCel->GetRows();

	for (int row = 0; row < rows; ++row)
	{
		CPPUNIT_ASSERT(fusionCel->IsMasked(0, row) == gcosCel->IsMasked(0, row));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::IndexToXTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = cells / 10;

	for (int i = 0; i < cells; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->IndexToX(i) == gcosCel->IndexToX(i));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::IndexToYTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = cells / 10;

	for (int i = 0; i < cells; i+= incr)
	{
		CPPUNIT_ASSERT(fusionCel->IndexToY(i) == gcosCel->IndexToY(i));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::XYToIndexTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int rows = fusionCel->GetRows();

	for (int row = 0; row < rows; ++row)
	{
		CPPUNIT_ASSERT(fusionCel->XYToIndex(0, row) == gcosCel->XYToIndex(0, row));
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::ClearTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);
	CPPUNIT_ASSERT_NO_THROW(fusionCel->Clear());

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());
	gcosCel->Clear();

	// Check some values
	CPPUNIT_ASSERT(fusionCel->GetRows() == gcosCel->GetRows());
	CPPUNIT_ASSERT(fusionCel->GetNumCells() == gcosCel->GetNumCells());
	CPPUNIT_ASSERT(fusionCel->GetFileName() == gcosCel->GetFileName());
}

void FusionCELDataGCOSTest::GetEntryIndexTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int cells = fusionCel->GetNumCells();
	int incr = cells / 10;

	affxcel::CELFileEntryType gcosEntry;
	FusionCELFileEntryType calvinEntry;

	for (int i = 0; i < cells; i+= incr)
	{
		fusionCel->GetEntry(i, calvinEntry);
		gcosCel->GetEntry(i, gcosEntry);

		CPPUNIT_ASSERT(calvinEntry.Intensity == gcosEntry.Intensity);
		CPPUNIT_ASSERT(calvinEntry.Stdv == gcosEntry.Stdv);
		CPPUNIT_ASSERT(calvinEntry.Pixels == gcosEntry.Pixels);
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::GetEntryXYTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->Read() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->Read());

	// Spot check the values
	int rows = fusionCel->GetRows();
	int cols = fusionCel->GetCols();

	affxcel::CELFileEntryType gcosEntry;
	FusionCELFileEntryType calvinEntry;

	for (int row = 0; row < rows; ++row)
	{
		fusionCel->GetEntry(cols-1, row, calvinEntry);
		gcosCel->GetEntry(cols-1, row, gcosEntry);

		CPPUNIT_ASSERT(calvinEntry.Intensity == gcosEntry.Intensity);
		CPPUNIT_ASSERT(calvinEntry.Stdv == gcosEntry.Stdv);
		CPPUNIT_ASSERT(calvinEntry.Pixels == gcosEntry.Pixels);
	}
	
	fusionCel->Close();
	gcosCel->Close();
}

void FusionCELDataGCOSTest::CallMethodsWhenObjectIsNotReadyTest()
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

void FusionCELDataGCOSTest::ErrorTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);
	
	fusionCel->SetError(L"test error");
	CPPUNIT_ASSERT(fusionCel->GetError() == L"test error");
}

void FusionCELDataGCOSTest::GetFileSizeTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->GetFileSize() == 5078705);
}

void FusionCELDataGCOSTest::ReadExTest()
{
	CPPUNIT_ASSERT(fusionCel->ReadEx(GCOS_CEL_FILE.c_str(), FusionCELData::CEL_DATA));

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->ReadHeader());

	// Check some header values
	CPPUNIT_ASSERT(fusionCel->GetVersion() == gcosCel->GetVersion());
	CPPUNIT_ASSERT(fusionCel->GetCols() == gcosCel->GetCols());
	CPPUNIT_ASSERT(fusionCel->GetRows() == gcosCel->GetRows());
	CPPUNIT_ASSERT(fusionCel->GetNumCells() == gcosCel->GetNumCells());

	std::string headerString = StringUtils::ConvertWCSToMBS(fusionCel->GetHeader());
	CPPUNIT_ASSERT(headerString == gcosCel->GetHeaderString());
	std::string algString = StringUtils::ConvertWCSToMBS(fusionCel->GetAlg());
	CPPUNIT_ASSERT(algString == gcosCel->GetAlg());
	std::string paramString = StringUtils::ConvertWCSToMBS(fusionCel->GetParams());
	CPPUNIT_ASSERT(paramString == gcosCel->GetParams());
	std::string chipTypetring = StringUtils::ConvertWCSToMBS(fusionCel->GetChipType());
	CPPUNIT_ASSERT(chipTypetring == gcosCel->GetChipType());
}

void FusionCELDataGCOSTest::GetReadStateTest()
{
	CPPUNIT_ASSERT(fusionCel->ReadEx(GCOS_CEL_FILE.c_str(), FusionCELData::CEL_DATA));
	CPPUNIT_ASSERT(fusionCel->GetReadState() == FusionCELData::CEL_DATA);
}

void FusionCELDataGCOSTest::GetHeaderKeyTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	//std::wstring s = fusionCel->GetHeader();
	CPPUNIT_ASSERT(fusionCel->GetHeaderKey(L"TotalX") == L"712");

	// Don't do this, it causes a crash in GCOS CCELFileData
	//CPPUNIT_ASSERT(fusionCel->GetHeaderKey(L"bogus") == L"");
}

void FusionCELDataGCOSTest::GetDatHeaderTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->ReadHeader());

	CPPUNIT_ASSERT(StringUtils::ConvertWCSToMBS(fusionCel->GetDatHeader()) == gcosCel->GetDatHeader());
}

void FusionCELDataGCOSTest::GetGridCornersTest()
{
	fusionCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(fusionCel->ReadHeader() == true);

	gcosCel->SetFileName(GCOS_CEL_FILE.c_str());
	CPPUNIT_ASSERT(gcosCel->ReadHeader());

	CompareGridCoordinates();
}
