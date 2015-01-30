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

#ifndef __MDLFILEDATATEST_H_
#define __MDLFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CMDLFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CMDLFileDataTest );
	
	// Test interface object creation
	CPPUNIT_TEST( testCreation );
	CPPUNIT_TEST( testCreationHeader );
	CPPUNIT_TEST( testCreationData );

	// Test header interface properties
	CPPUNIT_TEST( testproperty_HeaderVersion );
	CPPUNIT_TEST( testproperty_HeaderNumProbes );
	CPPUNIT_TEST( testproperty_HeaderAnalysisType );
	CPPUNIT_TEST( testproperty_HeaderMasked );
	CPPUNIT_TEST( testproperty_HeaderChipType );
	CPPUNIT_TEST( testproperty_HeaderModifiedDate );
	CPPUNIT_TEST( testproperty_HeaderEqualOperator );

	// Test data interface properties
	CPPUNIT_TEST( testproperty_DataAffinity );
	CPPUNIT_TEST( testproperty_DataRelativeBkg );
	CPPUNIT_TEST( testproperty_DataSaturation );
	CPPUNIT_TEST( testproperty_DataOffset );

	// Test file interface properties
	CPPUNIT_TEST( testproperty_Header );
	CPPUNIT_TEST( testproperty_FileName );

	// Test file interface methods
	CPPUNIT_TEST( testmethod_Exists );
	CPPUNIT_TEST( testmethod_ExistsWhenFileNotExists );
	CPPUNIT_TEST( testmethod_Read );
	CPPUNIT_TEST( testmethod_ReadHeader );
	CPPUNIT_TEST( testmethod_ReadWhenFileNotExists );
	CPPUNIT_TEST( testmethod_ReadWhenNoDataIsInTheFile );
	CPPUNIT_TEST( testmethod_InitializeForWriting );
	CPPUNIT_TEST( testmethod_SetProbeSetData );
	CPPUNIT_TEST( testmethod_GetProbeSetData );
	CPPUNIT_TEST( testmethod_Write );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	// Test the header size
	void testSize();

	// Test interface object creation
	void testCreation();
	void testCreationHeader();
	void testCreationData();

	// Test header interface properties
	void testproperty_HeaderVersion();
	void testproperty_HeaderNumProbes();
	void testproperty_HeaderAnalysisType();
	void testproperty_HeaderMasked();
	void testproperty_HeaderChipType();
	void testproperty_HeaderModifiedDate();
	void testproperty_HeaderEqualOperator();

	// Test data interface properties
	void testproperty_DataAffinity();
	void testproperty_DataRelativeBkg();
	void testproperty_DataSaturation();
	void testproperty_DataOffset();

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
	void testmethod_InitializeForWriting();
	void testmethod_SetProbeSetData();
	void testmethod_GetProbeSetData();
	void testmethod_Write();
};


#endif // __MDLFILEDATATEST_H_
