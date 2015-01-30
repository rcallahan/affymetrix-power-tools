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

#ifndef __DATDATATEST_H_
#define __DATDATATEST_H_

#include "calvin_files/data/src/DATData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_calvin_io;

class DATDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( DATDataTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( FilenameTest );
	CPPUNIT_TEST ( SetPixelCountTest );
	CPPUNIT_TEST ( SetStatsCountTest );
	CPPUNIT_TEST ( ArrayTypeTest );
	CPPUNIT_TEST ( PixelSizeTest );
	CPPUNIT_TEST ( ScannerTypeTest );
	CPPUNIT_TEST ( ScannerIDTest );
	CPPUNIT_TEST ( ScanDateTest );
	CPPUNIT_TEST ( RowsTest );
	CPPUNIT_TEST ( ColsTest );
	CPPUNIT_TEST ( GlobalGridTest );
	CPPUNIT_TEST ( SubgridTest );
	CPPUNIT_TEST ( ArrayIdTest );
	CPPUNIT_TEST ( ArrayBarcodeTest );
	CPPUNIT_TEST ( AddGridAlignmentAlgorithmParameterTest );
	CPPUNIT_TEST ( FindGridAlignmentAlgorithmParameterTest );
	CPPUNIT_TEST ( GetGridAlignmentAlgorithmParametersTest );
	CPPUNIT_TEST ( ClearGridAlignmentAlgorithmParametersTest );

	CPPUNIT_TEST_SUITE_END();

private:

	DATData *data;

public:
	void setUp();
	void tearDown();

	void testCreation();
	void FilenameTest();
	void SetPixelCountTest();
	void SetStatsCountTest();
	void ArrayTypeTest();
	void PixelSizeTest();
	void ScannerTypeTest();
	void ScannerIDTest();
	void ScanDateTest();
	void RowsTest();
	void ColsTest();
	void GlobalGridTest();
	void SubgridTest();
	void ArrayIdTest();
	void ArrayBarcodeTest();
	void AddGridAlignmentAlgorithmParameterTest();
	void FindGridAlignmentAlgorithmParameterTest();
	void GetGridAlignmentAlgorithmParametersTest();
	void ClearGridAlignmentAlgorithmParametersTest();
};

#endif // __DATDATATEST_H_
