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

#ifndef __CHPFILEDATATEST_H_
#define __CHPFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CCHPFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CCHPFileDataTest );
	
	CPPUNIT_TEST( testproperty_Header_Dimensions );
	CPPUNIT_TEST( testproperty_Header_NumProbeSets );
	CPPUNIT_TEST( testproperty_Header_AssayType );
	CPPUNIT_TEST( testproperty_Header_AlgName );
	CPPUNIT_TEST( testproperty_Header_AlgVersion );
	CPPUNIT_TEST( testproperty_Header_AlgParams );
	CPPUNIT_TEST( testproperty_Header_SummaryParams );
	CPPUNIT_TEST( testproperty_Header_ChipType );
	CPPUNIT_TEST( testproperty_Header_ParentCellFile );
	CPPUNIT_TEST( testproperty_Header_ProgID );
	CPPUNIT_TEST( testproperty_Header_BackgroundZoneInfo );
	CPPUNIT_TEST( testmethod_Exists);
	CPPUNIT_TEST( testmethod_IsXDACompatibleFile);
	CPPUNIT_TEST( testmethod_ReadHeader_fail);
	CPPUNIT_TEST( testmethod_Read_fail);
	CPPUNIT_TEST( testproperty_Paths);
	CPPUNIT_TEST( testmethod_Read_Exp_Comp);
	CPPUNIT_TEST( testmethod_Read_Exp_Abs);
	CPPUNIT_TEST( testmethod_ReadHeader_Exp_Comp);
	CPPUNIT_TEST( testmethod_ReadHeader_Exp_Abs);
	CPPUNIT_TEST( testmethod_ReadHeader_TagV11);
	CPPUNIT_TEST( testmethod_Read_TagV11);
	CPPUNIT_TEST( testmethod_ReadHeader_TagXDA);
	CPPUNIT_TEST( testmethod_Read_TagXDA);
	CPPUNIT_TEST( testmethod_Read_ReseqXDA_v1);
	CPPUNIT_TEST( testmethod_Read_ReseqXDA_v2);
	CPPUNIT_TEST( testmethod_Read_Reseq_old_file);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testproperty_Header_Dimensions();
	void testproperty_Header_NumProbeSets();
	void testproperty_Header_AssayType();
	void testproperty_Header_AlgName();
	void testproperty_Header_AlgVersion();
	void testproperty_Header_AlgParams();
	void testproperty_Header_SummaryParams();
	void testproperty_Header_ChipType();
	void testproperty_Header_ParentCellFile();
	void testproperty_Header_ProgID();
	void testproperty_Header_BackgroundZoneInfo();
	void testmethod_Exists();
	void testmethod_IsXDACompatibleFile();
	void testmethod_ReadHeader_fail();
	void testmethod_Read_fail();
	void testproperty_Paths();
	void testmethod_Read_Exp_Comp();
	void testmethod_Read_Exp_Abs();
	void testmethod_ReadHeader_Exp_Comp();
	void testmethod_ReadHeader_Exp_Abs();
	void testmethod_ReadHeader_TagV11();
	void testmethod_Read_TagV11();
	void testmethod_ReadHeader_TagXDA();
	void testmethod_Read_TagXDA();
	void testmethod_Read_ReseqXDA_v1();
	void testmethod_Read_ReseqXDA_v2();
	void testmethod_Read_Reseq_old_file();
};


#endif // __CHPFILEDATATEST_H_
