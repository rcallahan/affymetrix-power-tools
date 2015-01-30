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

#ifndef __DATFILEDATATEST_H_
#define __DATFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CDATFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CDATFileDataTest );

	CPPUNIT_TEST( testproperty_Header_Cols );
	CPPUNIT_TEST( testproperty_Header_Rows );
	CPPUNIT_TEST( testproperty_Header_MinValue );
	CPPUNIT_TEST( testproperty_Header_MaxValue );
	CPPUNIT_TEST( testproperty_Header_Grid );
	CPPUNIT_TEST( testproperty_Header_ArrayType );
	CPPUNIT_TEST( testproperty_Header_ScannerID );

	CPPUNIT_TEST( testproperty_FileName );
	CPPUNIT_TEST( testmethod_Exists );
	CPPUNIT_TEST( testmethod_ExistsWhenFileNotExists );
	CPPUNIT_TEST( testmethod_Read );
	CPPUNIT_TEST( testmethod_Read_with_scanner_ID );
	CPPUNIT_TEST( testmethod_ReadWhenFileNotExists );
	CPPUNIT_TEST( testmethod_GetPixel );
	CPPUNIT_TEST( testmethod_GetPixels );
	CPPUNIT_TEST( testmethod_GetPixels_many_rows );
	CPPUNIT_TEST( testmethod_IsGCOSDATFile );
	CPPUNIT_TEST( testmethod_UpdateGridCorners );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testproperty_Header_Cols ();
	void testproperty_Header_Rows ();
	void testproperty_Header_MinValue ();
	void testproperty_Header_MaxValue ();
	void testproperty_Header_Grid ();
	void testproperty_Header_ArrayType();
	void testproperty_Header_ScannerID();

	void testproperty_FileName ();
	void testmethod_Exists ();
	void testmethod_ExistsWhenFileNotExists ();
	void testmethod_Read ();
	void testmethod_Read_with_scanner_ID();
	void testmethod_ReadWhenFileNotExists ();
	void testmethod_GetPixel ();
	void testmethod_GetPixels ();
	void testmethod_GetPixels_many_rows ();
	void testmethod_IsGCOSDATFile ();
	void testmethod_UpdateGridCorners ();
};


#endif // __DATFILEDATATEST_H_
