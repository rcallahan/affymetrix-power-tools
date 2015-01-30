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

#ifndef __FILEWRITERTEST_H_
#define __FILEWRITERTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class FileWriterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( FileWriterTest );

	CPPUNIT_TEST(testfunction_WriteInt32_I);
	CPPUNIT_TEST(testfunction_WriteUInt32_I);
	CPPUNIT_TEST(testfunction_WriteInt16_I);
	CPPUNIT_TEST(testfunction_WriteUInt16_I);
	CPPUNIT_TEST(testfunction_WriteInt32_N);
	CPPUNIT_TEST(testfunction_WriteUInt32_N);
	CPPUNIT_TEST(testfunction_WriteInt16_N);
	CPPUNIT_TEST(testfunction_WriteUInt16_N);
	CPPUNIT_TEST(testfunction_WriteInt8);
	CPPUNIT_TEST(testfunction_WriteUInt8);
	CPPUNIT_TEST(testfunction_WriteFloat_I);
	CPPUNIT_TEST(testfunction_WriteFloat_N);
	CPPUNIT_TEST(testfunction_WriteFloatLowPrecision);
	CPPUNIT_TEST(testfunction_WriteCharacterArray);
	CPPUNIT_TEST(testfunction_WriteFixedCString);
	CPPUNIT_TEST(testfunction_WriteCString);
	CPPUNIT_TEST(testfunction_WriteFixedString);
	CPPUNIT_TEST(testfunction_WriteString_I);
	CPPUNIT_TEST(testfunction_WriteString_N);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testfunction_WriteInt32_I();
	void testfunction_WriteUInt32_I();
	void testfunction_WriteInt16_I();
	void testfunction_WriteUInt16_I();
	void testfunction_WriteInt32_N();
	void testfunction_WriteUInt32_N();
	void testfunction_WriteInt16_N();
	void testfunction_WriteUInt16_N();
	void testfunction_WriteInt8();
	void testfunction_WriteUInt8();
	void testfunction_WriteFloat_I();
	void testfunction_WriteFloat_N();
	void testfunction_WriteFloatLowPrecision();
	void testfunction_WriteCharacterArray();
	void testfunction_WriteFixedCString();
	void testfunction_WriteCString();
	void testfunction_WriteFixedString();
	void testfunction_WriteString_I();
	void testfunction_WriteString_N();
};


#endif // __FILEWRITERTEST_H_
