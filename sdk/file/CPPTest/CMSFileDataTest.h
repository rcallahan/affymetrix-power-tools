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

#ifndef __CMSFILEDATATEST_H_
#define __CMSFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CCMSFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CCMSFileDataTest );
	
	// Test interface object creation
	CPPUNIT_TEST( testCreation );
	CPPUNIT_TEST( testCreationHeader );

	// Test header interface properties
	CPPUNIT_TEST( testproperty_HeaderVersion );
	CPPUNIT_TEST( testproperty_HeaderArrayCount );
	CPPUNIT_TEST( testproperty_HeaderSNPCount );
	CPPUNIT_TEST( testproperty_HeaderAssay );

	// Test file interface properties
	CPPUNIT_TEST( testproperty_FileName );

	// Test file interface methods
	CPPUNIT_TEST( testmethod_Exists );
	CPPUNIT_TEST( testmethod_ExistsWhenFileNotExists );
	CPPUNIT_TEST( testmethod_Read );
	CPPUNIT_TEST( testmethod_ReadHeader );
	CPPUNIT_TEST( testmethod_ReadWhenFileNotExists );
	CPPUNIT_TEST( testmethod_ReadWhenNoDataIsInTheFile );
	CPPUNIT_TEST( testmethod_ArrayTypeInformation );
	CPPUNIT_TEST( testmethod_SNPInformation );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	// Test interface object creation
	void testCreation();
	void testCreationHeader();

	// Test header interface properties
	void testproperty_HeaderVersion();
	void testproperty_HeaderArrayCount();
	void testproperty_HeaderSNPCount();
	void testproperty_HeaderAssay();

	// Test file interface properties
	void testproperty_Header();
	void testproperty_FileName();

	// Test file interface methods
	void testmethod_Exists();
	void testmethod_ExistsWhenFileNotExists();
	void testmethod_Read();
	void testmethod_ReadHeader();
	void testmethod_ReadWhenFileNotExists();
	void testmethod_ReadWhenNoDataIsInTheFile();
	void testmethod_ArrayTypeInformation();
	void testmethod_SNPInformation();
};


#endif // __CMSFILEDATATEST_H_
