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

#ifndef __GRCFILEDATATEST_H_
#define __GRCFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CGRCFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CGRCFileDataTest );
	
	CPPUNIT_TEST( testCreation );
	CPPUNIT_TEST(testproperty_FileName);
	CPPUNIT_TEST(testproperty_NumCols);
	CPPUNIT_TEST(testproperty_NumRows);
	CPPUNIT_TEST(testproperty_NumB1);
	CPPUNIT_TEST(testproperty_NumB2);
	CPPUNIT_TEST(testproperty_NumNonSynth);
	CPPUNIT_TEST(testmethod_Exists);
	CPPUNIT_TEST(testmethod_Read);
	CPPUNIT_TEST(testmethod_GetB1);
	CPPUNIT_TEST(testmethod_GetB2);
	CPPUNIT_TEST(testmethod_GetNonSynth);
	CPPUNIT_TEST(testmethod_ExistsWhenFileNotExists);
	CPPUNIT_TEST(testmethod_ReadWhenFileNotExists);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();

	void testproperty_FileName();
	void testproperty_NumCols();
	void testproperty_NumRows();
	void testproperty_NumB1();
	void testproperty_NumB2();
	void testproperty_NumNonSynth();

	void testmethod_Exists();
	void testmethod_Read();
	void testmethod_GetB1();
	void testmethod_GetB2();
	void testmethod_GetNonSynth();
	void testmethod_ExistsWhenFileNotExists();
	void testmethod_ReadWhenFileNotExists();
};


#endif // __GRCFILEDATATEST_H_
