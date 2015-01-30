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

#ifndef __EXPFILEDATATEST_H_
#define __EXPFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CEXPFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CEXPFileDataTest );

	// Test interface object creation
	CPPUNIT_TEST( testCreation );

	// Test interface properties
	CPPUNIT_TEST( testproperty_FileName );
	CPPUNIT_TEST( testproperty_ArrayType );

	// Test interface methods
	CPPUNIT_TEST( testmethod_Exists );
	CPPUNIT_TEST( testmethod_Read_full );
	CPPUNIT_TEST( testmethod_Read_bare );
	CPPUNIT_TEST( testmethod_Read_hyb_abort );
	CPPUNIT_TEST( testmethod_ExistsWhenFileNotExists );
	CPPUNIT_TEST( testmethod_ReadWhenFileNotExists );
	CPPUNIT_TEST( testmethod_Clear );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testproperty_FileName();
	void testproperty_ArrayType();
	void testmethod_Exists();
	void testmethod_Read_full();
	void testmethod_Read_bare();
	void testmethod_Read_hyb_abort();
	void testmethod_ExistsWhenFileNotExists();
	void testmethod_ReadWhenFileNotExists();
	void testmethod_Clear();
};


#endif // __EXPFILEDATATEST_H_
