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
#pragma once

#include "calvin_files/fusion/src/FusionCDFData.h"
//
#include "file/CDFFileData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class FusionCdfFileTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( FusionCdfFileTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testExpressionV3 );
	CPPUNIT_TEST ( testExpressionXDA );
	CPPUNIT_TEST ( testMissingFile );
	CPPUNIT_TEST ( testMappingV3 );
	CPPUNIT_TEST ( testMappingXDA );
    CPPUNIT_TEST( test_DMET3ASCII );
    CPPUNIT_TEST( test_DMET3XDA );	
	CPPUNIT_TEST( test_MultichannelASCII );
    CPPUNIT_TEST( test_MultichannelXDA );	
	CPPUNIT_TEST( test_AxiomV6ASCII );
    CPPUNIT_TEST( test_AxiomV4XDA );	
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp(); 
	void tearDown();
	void testRegObjects(affxcdf::CCDFFileData &gcosCdf, affymetrix_fusion_io::FusionCDFData &fusionCdf);
	void testQCObjects(affxcdf::CCDFFileData &gcosCdf, affymetrix_fusion_io::FusionCDFData &fusionCdf);

	void testExpressionV3();
	void testExpressionXDA();
	void testMissingFile();
	void testMappingV3();
	void testMappingXDA();
    void test_DMET3XDA();
	void test_DMET3ASCII();
	void test_MultichannelXDA();
	void test_MultichannelASCII();
	void test_AxiomV4XDA();
	void test_AxiomV6ASCII();
	void testCreation();
};
