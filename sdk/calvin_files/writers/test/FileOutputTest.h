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

#ifndef __FILEOUTPUTTEST_H_
#define __FILEOUTPUTTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class FileOutputTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( FileOutputTest );

	CPPUNIT_TEST( testdefine_SIZES );

	CPPUNIT_TEST( testmethod_WriteInt8 );
	CPPUNIT_TEST( testmethod_WriteUInt8 );
	CPPUNIT_TEST( testmethod_WriteInt16 );
	CPPUNIT_TEST( testmethod_WriteUInt16 );
	CPPUNIT_TEST( testmethod_WriteInt32 );
	CPPUNIT_TEST( testmethod_WriteUInt32 );
	CPPUNIT_TEST( testmethod_WriteFloat );
	CPPUNIT_TEST( testmethod_WriteString8);
	CPPUNIT_TEST( testmethod_WriteString8_fixedlen);
	CPPUNIT_TEST( testmethod_WriteString16);
	CPPUNIT_TEST( testmethod_WriteString16_fixedlen);
	CPPUNIT_TEST( testmethod_WriteBlob);
	CPPUNIT_TEST( testmethod_WriteBlobWithReserve);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testdefine_SIZES();
	void testmethod_WriteInt8();
	void testmethod_WriteUInt8();
	void testmethod_WriteInt16();
	void testmethod_WriteUInt16();
	void testmethod_WriteInt32();
	void testmethod_WriteUInt32();
	void testmethod_WriteFloat();

	void testmethod_WriteString8();
	void testmethod_WriteString8_fixedlen();
	void testmethod_WriteString16();
	void testmethod_WriteString16_fixedlen();
	void testmethod_WriteBlob();
	void testmethod_WriteBlobWithReserve();
};


#endif // __FILEOUTPUTTEST_H_
