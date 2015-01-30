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

#ifndef __FILEINPUTTEST_H_
#define __FILEINPUTTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class FileInputTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( FileInputTest );

	CPPUNIT_TEST( testdefine_SIZES );

	CPPUNIT_TEST( testmethod_ReadInt8 );
	CPPUNIT_TEST( testmethod_ReadInt16 );
	CPPUNIT_TEST( testmethod_ReadInt32 );

	CPPUNIT_TEST( testmethod_ReadUInt8 );
	CPPUNIT_TEST( testmethod_ReadUInt16 );
	CPPUNIT_TEST( testmethod_ReadUInt32 );

	CPPUNIT_TEST( testmethod_ReadInt8_stream );
	CPPUNIT_TEST( testmethod_ReadInt16_stream );
	CPPUNIT_TEST( testmethod_ReadInt32_stream );

	CPPUNIT_TEST( testmethod_ReadUInt8_stream );
	CPPUNIT_TEST( testmethod_ReadUInt16_stream );
	CPPUNIT_TEST( testmethod_ReadUInt32_stream );

	CPPUNIT_TEST( testmethod_ReadFloat );
	CPPUNIT_TEST( testmethod_ReadFloat_stream );

	CPPUNIT_TEST( testmethod_ReadString8);
	CPPUNIT_TEST( testmethod_ReadString8_stream);

	CPPUNIT_TEST( testmethod_ReadString8_fixedlen);
	CPPUNIT_TEST( testmethod_ReadString8_fixedlen_stream);

	CPPUNIT_TEST( testmethod_ReadString16);
	CPPUNIT_TEST( testmethod_ReadString16_stream);

	CPPUNIT_TEST( testmethod_ReadString16_fixedlen);
	CPPUNIT_TEST( testmethod_ReadString16_fixedlen_stream);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testdefine_SIZES();
	void testmethod_ReadInt8();
	void testmethod_ReadInt16();
	void testmethod_ReadInt32();

	void testmethod_ReadUInt8();
	void testmethod_ReadUInt16();
	void testmethod_ReadUInt32();

	void testmethod_ReadInt8_stream();
	void testmethod_ReadInt16_stream();
	void testmethod_ReadInt32_stream();

	void testmethod_ReadUInt8_stream();
	void testmethod_ReadUInt16_stream();
	void testmethod_ReadUInt32_stream();

	void testmethod_ReadFloat();
	void testmethod_ReadFloat_stream();

	void testmethod_ReadString8();
	void testmethod_ReadString8_stream();

	void testmethod_ReadString8_fixedlen();
	void testmethod_ReadString8_fixedlen_stream();

	void testmethod_ReadString16();
	void testmethod_ReadString16_stream();

	void testmethod_ReadString16_fixedlen();
	void testmethod_ReadString16_fixedlen_stream();

};


#endif // __FILEINPUTTEST_H_
