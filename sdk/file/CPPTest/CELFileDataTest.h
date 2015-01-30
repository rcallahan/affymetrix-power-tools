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

#ifndef __CELFILEDATATEST_H_
#define __CELFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CCELFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CCELFileDataTest );
	
	CPPUNIT_TEST(testHeader_Version);
	CPPUNIT_TEST(testHeader_AlgParams);
	CPPUNIT_TEST(testHeader_Dimensions);
	CPPUNIT_TEST(testHeader_AlgName);
	CPPUNIT_TEST(testHeader_ChipType);
	CPPUNIT_TEST(testHeader_Margin);
	CPPUNIT_TEST(testHeader_DatHeader);
	CPPUNIT_TEST(testHeader_Header);
	CPPUNIT_TEST(testHeader_SetGridCorners);
	CPPUNIT_TEST(testHeader_Cells);
	CPPUNIT_TEST(testHeader_MaskedCells);
	CPPUNIT_TEST(testHeader_Outliers);
	CPPUNIT_TEST(testCELFileEntryType_size);
	CPPUNIT_TEST(testCELFileTranscriptomeEntryType_size);
	CPPUNIT_TEST(testRead_v3);
	CPPUNIT_TEST(testRead_v4);
	CPPUNIT_TEST(testRead_bcel);
	CPPUNIT_TEST(testRead_ccel);
	CPPUNIT_TEST(testRead_v3_no_mask_outlier);
	CPPUNIT_TEST(testRead_v4_no_mask_outlier);
	CPPUNIT_TEST(testRead_bcel_no_mask_outlier);
	CPPUNIT_TEST(testRead_ccel_no_mask_outlier);
	CPPUNIT_TEST(testRead_v3_no_mask);
	CPPUNIT_TEST(testRead_v4_no_mask);
	CPPUNIT_TEST(testRead_bcel_no_mask);
	CPPUNIT_TEST(testRead_ccel_no_mask);
	CPPUNIT_TEST(testRead_v3_no_outlier);
	CPPUNIT_TEST(testRead_v4_no_outlier);
	CPPUNIT_TEST(testRead_bcel_no_outlier);
	CPPUNIT_TEST(testRead_ccel_no_outlier);
	CPPUNIT_TEST(testRead_v3_all);
	CPPUNIT_TEST(testRead_v4_all);
	CPPUNIT_TEST(testRead_bcel_all);
	CPPUNIT_TEST(testRead_ccel_all);
	CPPUNIT_TEST(testReadHeader_v3);
	CPPUNIT_TEST(testReadHeader_v4);
	CPPUNIT_TEST(testReadHeader_bcel);
	CPPUNIT_TEST(testReadHeader_ccel);
	CPPUNIT_TEST(testExists);
	CPPUNIT_TEST(testIsXDACompatibleFile);
	CPPUNIT_TEST(testIsTranscriptomeBcelFile);
	CPPUNIT_TEST(testIsCompactCelFile);
	CPPUNIT_TEST(testReadHeader_fail);
	CPPUNIT_TEST(testRead_fail);
	CPPUNIT_TEST(testPaths);
	CPPUNIT_TEST(testXYToIndex);
	CPPUNIT_TEST(testIndexToXY);
	CPPUNIT_TEST(testParameterNames);
	CPPUNIT_TEST(testmethod_GetHeaderKey_with_invalid_input);
	CPPUNIT_TEST(testmethod_IsVersion3CompatibleFile);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testHeader_Version();
	void testHeader_AlgParams();
	void testHeader_Dimensions();
	void testHeader_AlgName();
	void testHeader_ChipType();
	void testHeader_Margin();
	void testHeader_DatHeader();
	void testHeader_Header();
	void testHeader_SetGridCorners();
	void testHeader_Cells();
	void testHeader_MaskedCells();
	void testHeader_Outliers();
	void testCELFileEntryType_size();
	void testCELFileTranscriptomeEntryType_size();
	void testRead_v3();
	void testRead_v4();
	void testRead_bcel();
	void testRead_ccel();
	void testRead_v3_no_mask_outlier();
	void testRead_v4_no_mask_outlier();
	void testRead_bcel_no_mask_outlier();
	void testRead_ccel_no_mask_outlier();
	void testRead_v3_no_mask();
	void testRead_v4_no_mask();
	void testRead_bcel_no_mask();
	void testRead_ccel_no_mask();
	void testRead_v3_no_outlier();
	void testRead_v4_no_outlier();
	void testRead_bcel_no_outlier();
	void testRead_ccel_no_outlier();
	void testRead_v3_all();
	void testRead_v4_all();
	void testRead_bcel_all();
	void testRead_ccel_all();
	void testReadHeader_v3();
	void testReadHeader_v4();
	void testReadHeader_bcel();
	void testReadHeader_ccel();
	void testExists();
	void testIsXDACompatibleFile();
	void testIsTranscriptomeBcelFile();
	void testIsCompactCelFile();
	void testReadHeader_fail();
	void testRead_fail();
	void testPaths();
	void testXYToIndex();
	void testIndexToXY();
	void testParameterNames();
	void testmethod_GetHeaderKey_with_invalid_input();
	void testmethod_IsVersion3CompatibleFile();
};


#endif // __CELFILEDATATEST_H_
