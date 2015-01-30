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

#ifndef __FILEIOTEST_H_
#define __FILEIOTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CFileIOTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CFileIOTest );

	CPPUNIT_TEST( testdefine_SIZES );
	CPPUNIT_TEST( testfunction_affy_swap );
	CPPUNIT_TEST(testfunction_ReadInt32_I);
	CPPUNIT_TEST(testfunction_ReadUInt32_I);
	CPPUNIT_TEST(testfunction_ReadInt16_I);
	CPPUNIT_TEST(testfunction_ReadUInt16_I);
	CPPUNIT_TEST(testfunction_ReadInt8);
	CPPUNIT_TEST(testfunction_ReadUInt8);
	CPPUNIT_TEST(testfunction_ReadFloat_I);
	CPPUNIT_TEST(testfunction_ReadCString_I);
	CPPUNIT_TEST(testfunction_ReadString_I);
	CPPUNIT_TEST(testfunction_ReadUIntLenString_I);
	CPPUNIT_TEST(testfunction_ReadInt32_N);
	CPPUNIT_TEST(testfunction_ReadUInt32_N);
	CPPUNIT_TEST(testfunction_ReadInt16_N);
	CPPUNIT_TEST(testfunction_ReadUInt16_N);
	CPPUNIT_TEST(testfunction_ReadFloat_N);
	CPPUNIT_TEST(testfunction_ReadCString_N);
	CPPUNIT_TEST(testfunction_ReadString_N);
	CPPUNIT_TEST(testfunction_ReadUIntLenString_N);
	CPPUNIT_TEST(testfunction_ReadNextLine_dosfile);
	CPPUNIT_TEST(testfunction_ReadNextLine_unixfile);
	CPPUNIT_TEST(testfunction_ReadFloatFromOldBPMAP);
	CPPUNIT_TEST(testfunction_MmGetFloatFromOldBPMAP);
	CPPUNIT_TEST(testfunction_ReadFixedString);
	CPPUNIT_TEST(testfunction_ReadFixedCString);
	CPPUNIT_TEST(testfunction_ReadFixedUCString);

	CPPUNIT_TEST(testfunction_MmGetInt32_I);
	CPPUNIT_TEST(testfunction_MmGetUInt32_I);
	CPPUNIT_TEST(testfunction_MmGetInt16_I);
	CPPUNIT_TEST(testfunction_MmGetUInt16_I);
	CPPUNIT_TEST(testfunction_MmGetInt8);
	CPPUNIT_TEST(testfunction_MmGetUInt8);
	CPPUNIT_TEST(testfunction_MmGetFloat_I);
	CPPUNIT_TEST(testfunction_MmGetInt32_N);
	CPPUNIT_TEST(testfunction_MmGetUInt32_N);
	CPPUNIT_TEST(testfunction_MmGetInt16_N);
	CPPUNIT_TEST(testfunction_MmGetUInt16_N);
	CPPUNIT_TEST(testfunction_MmGetFloat_N);

	CPPUNIT_TEST(testfunction_MmSetUInt32_I);
	CPPUNIT_TEST(testfunction_MmSetUInt16_I);
	CPPUNIT_TEST(testfunction_MmSetUInt8);
	CPPUNIT_TEST(testfunction_MmSetFloat_I);
	CPPUNIT_TEST(testfunction_MmSetUInt32_N);
	CPPUNIT_TEST(testfunction_MmSetUInt16_N);
	CPPUNIT_TEST(testfunction_MmSetFloat_N);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testdefine_SIZES();
	void testfunction_affy_swap();
	void testfunction_ReadInt32_I();
	void testfunction_ReadUInt32_I();
	void testfunction_ReadInt16_I();
	void testfunction_ReadUInt16_I();
	void testfunction_ReadInt8();
	void testfunction_ReadUInt8();
	void testfunction_ReadFloat_I();
	void testfunction_ReadCString_I();
	void testfunction_ReadString_I();
	void testfunction_ReadUIntLenString_I();
	void testfunction_ReadInt32_N();
	void testfunction_ReadUInt32_N();
	void testfunction_ReadInt16_N();
	void testfunction_ReadUInt16_N();
	void testfunction_ReadFloat_N();
	void testfunction_ReadCString_N();
	void testfunction_ReadString_N();
	void testfunction_ReadUIntLenString_N();
	void testfunction_ReadNextLine_dosfile();
	void testfunction_ReadNextLine_unixfile();
	void testfunction_ReadFloatFromOldBPMAP();
	void testfunction_MmGetFloatFromOldBPMAP();
	void testfunction_ReadFixedString();
	void testfunction_ReadFixedCString();
	void testfunction_ReadFixedUCString();

	void testfunction_MmGetInt32_I();
	void testfunction_MmGetUInt32_I();
	void testfunction_MmGetInt16_I();
	void testfunction_MmGetUInt16_I();
	void testfunction_MmGetInt8();
	void testfunction_MmGetUInt8();
	void testfunction_MmGetFloat_I();
	void testfunction_MmGetInt32_N();
	void testfunction_MmGetUInt32_N();
	void testfunction_MmGetInt16_N();
	void testfunction_MmGetUInt16_N();
	void testfunction_MmGetFloat_N();

	void testfunction_MmSetUInt32_I();
	void testfunction_MmSetUInt16_I();
	void testfunction_MmSetUInt8();
	void testfunction_MmSetFloat_I();
	void testfunction_MmSetUInt32_N();
	void testfunction_MmSetUInt16_N();
	void testfunction_MmSetFloat_N();
};


#endif // __FILEIOTEST_H_
