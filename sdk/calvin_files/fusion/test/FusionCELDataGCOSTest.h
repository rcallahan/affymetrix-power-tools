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

#ifndef __FUSIONCELDATAGCOSTEST_H_
#define __FUSIONCELDATAGCOSTEST_H_

//
#include "calvin_files/fusion/src/FusionCELData.h"
//
#include "file/CELFileData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

// This test fixture tests the FusionCELData class
// with GCOS CEL data.

class FusionCELDataGCOSTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE (FusionCELDataGCOSTest);

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (FileNameTest);
	CPPUNIT_TEST (ExistsTest);
	CPPUNIT_TEST (ReadHeaderTest);
//	CPPUNIT_TEST (DimensionsTest);
	CPPUNIT_TEST (AlgorithmNameTest);
	CPPUNIT_TEST (ChipTypeTest);
	CPPUNIT_TEST (LibPackageTest);
	CPPUNIT_TEST (MasterFileTest);
	CPPUNIT_TEST (MarginTest);
	CPPUNIT_TEST (ParamsTest);
	CPPUNIT_TEST (GetAlgorithmParametersTest);
	CPPUNIT_TEST (GetNumberAlgorithmParametersTest);
	CPPUNIT_TEST (GetAlgorithmParameterTagTest);
	CPPUNIT_TEST (GetAlgorithmParameterTest);
	CPPUNIT_TEST (ParametersTest);
	CPPUNIT_TEST (ReadHeaderAndVerifyTest);
	CPPUNIT_TEST (ReadTest);
	CPPUNIT_TEST (CloseTest);
	CPPUNIT_TEST (ClearTest);
	CPPUNIT_TEST (GetIntensityIndexTest);
	CPPUNIT_TEST (GetIntensityXYTest);
	CPPUNIT_TEST (GetStdvIndexTest);
	CPPUNIT_TEST (GetStdvXYTest);
	CPPUNIT_TEST (GetPixelsIndexTest);
	CPPUNIT_TEST (GetPixelsXYTest);
	CPPUNIT_TEST (IsOutlierIndexTest);
	CPPUNIT_TEST (IsOutlierXYTest);
	CPPUNIT_TEST (IsMaskedIndexTest);
	CPPUNIT_TEST (IsMaskedXYTest);
	CPPUNIT_TEST (IndexToXTest);
	CPPUNIT_TEST (IndexToYTest);
	CPPUNIT_TEST (XYToIndexTest);
	CPPUNIT_TEST (GetEntryIndexTest);
	CPPUNIT_TEST (GetEntryXYTest);
	CPPUNIT_TEST (CallMethodsWhenObjectIsNotReadyTest);
	CPPUNIT_TEST (ErrorTest);
	CPPUNIT_TEST (GetFileSizeTest);
	CPPUNIT_TEST (ReadExTest);
	CPPUNIT_TEST (GetReadStateTest);
	CPPUNIT_TEST (GetHeaderKeyTest);
	CPPUNIT_TEST (GetDatHeaderTest);
	CPPUNIT_TEST (GetGridCornersTest);
	CPPUNIT_TEST (GetFileIdTest);
	CPPUNIT_TEST (TestSingleChannel);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void CreationTest();
	void FileNameTest();
	void ExistsTest();
	void ReadHeaderTest();
//	void DimensionsTest();
	void AlgorithmNameTest();
	void ChipTypeTest();
    void LibPackageTest();
	void MasterFileTest();
	void MarginTest();
	void ParamsTest();
	void GetAlgorithmParametersTest();
	void GetNumberAlgorithmParametersTest();
	void GetAlgorithmParameterTagTest();
	void GetAlgorithmParameterTest();
	void ParametersTest();
	void ReadHeaderAndVerifyTest();
	void ReadTest();
	void CloseTest();
	void ClearTest();
	void GetIntensityIndexTest();
	void GetIntensityXYTest();
	void GetStdvIndexTest();
	void GetStdvXYTest();
	void GetPixelsIndexTest();
	void GetPixelsXYTest();
	void IsOutlierIndexTest();
	void IsOutlierXYTest();
	void IsMaskedIndexTest();
	void IsMaskedXYTest();
	void IndexToXTest();
	void IndexToYTest();
	void XYToIndexTest();
	void GetEntryIndexTest();
	void GetEntryXYTest();
	void CallMethodsWhenObjectIsNotReadyTest();
	void ErrorTest();
	void GetFileSizeTest();
	void ReadExTest();
	void GetReadStateTest();
	void GetHeaderKeyTest();
	void GetDatHeaderTest();
	void GetGridCornersTest();
	void GetFileIdTest();
	void TestSingleChannel();

// Helper methods
private:
	void CompareHeaderData();
	void CompareGridCoordinates();

private:
	affymetrix_fusion_io::FusionCELData* fusionCel;
	affxcel::CCELFileData* gcosCel;
};

#endif // __FUSIONCELDATAGCOSTEST_H_
