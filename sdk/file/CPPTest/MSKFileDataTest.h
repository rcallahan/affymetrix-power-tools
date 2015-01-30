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

#ifndef __MSKFILEDATATEST_H_
#define __MSKFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CMSKFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CMSKFileDataTest );
	
	// Test interface object creation
	CPPUNIT_TEST( testCreation );

	// Test interface properties
	CPPUNIT_TEST( testproperty_FileName );
	CPPUNIT_TEST( testproperty_ArrayType );

	// Test interface methods
	CPPUNIT_TEST( testmethod_ProbeSetIndiciesListCount );
	CPPUNIT_TEST( testmethod_Exists );
	CPPUNIT_TEST( testmethod_Read );
	CPPUNIT_TEST( testmethod_GetProbeSetIndiciesIterators );
	CPPUNIT_TEST( testmethod_Clear );
	CPPUNIT_TEST( testmethod_GetProbeSetListCount) ;
	CPPUNIT_TEST( testmethod_GetProbeSetListIterators ) ;
	CPPUNIT_TEST( testmethod_ExistsWhenFileNotExists );
	CPPUNIT_TEST( testmethod_ReadWhenFileNotExists );
	CPPUNIT_TEST( testmethod_ReadWhenNoDataIsInTheFile );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testproperty_FileName();
	void testproperty_ArrayType();
	void testmethod_ProbeSetIndiciesListCount();
	void testmethod_Exists();
	void testmethod_Read();
	void testmethod_GetProbeSetIndiciesIterators();
	void testmethod_GetProbeSetListCount();
	void testmethod_GetProbeSetListIterators();
	void testmethod_Clear();
	void testmethod_ExistsWhenFileNotExists();
	void testmethod_ReadWhenFileNotExists();
	void testmethod_ReadWhenNoDataIsInTheFile();
};


#endif // __MSKFILEDATATEST_H_
