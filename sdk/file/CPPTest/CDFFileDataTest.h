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

#ifndef __CDFFILEDATATEST_H_
#define __CDFFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CCDFFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CCDFFileDataTest );
	
	CPPUNIT_TEST( test_CCDFFileHeader );
	CPPUNIT_TEST( test_CCDFProbeInformation );
	CPPUNIT_TEST( test_CCDFProbeGroupInformation );
	CPPUNIT_TEST( test_CCDFProbeSetInformation );
	CPPUNIT_TEST( test_CCDFQCProbeInformation );
	CPPUNIT_TEST( test_CCDFQCProbeSetInformation );
	CPPUNIT_TEST( testmethod_IsXDACompatibleFile );
	CPPUNIT_TEST( testmethod_Exists );
	CPPUNIT_TEST( testmethod_ReadHeader_with_ASCII );
	CPPUNIT_TEST( testmethod_ReadHeader_with_XDA );
	CPPUNIT_TEST( test_ExpressionXDA );
	CPPUNIT_TEST( test_ExpressionASCII );
	CPPUNIT_TEST( test_GenotypingXDA );
	CPPUNIT_TEST( test_GenotypingASCII );
	CPPUNIT_TEST( test_ReseqXDA );
	CPPUNIT_TEST( test_ReseqASCII );
	CPPUNIT_TEST( test_TagXDA );
	CPPUNIT_TEST( test_TagASCII );
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

	void test_CCDFFileHeader();
	void test_CCDFProbeInformation();
	void test_CCDFProbeGroupInformation();
	void test_CCDFProbeSetInformation();
	void test_CCDFQCProbeInformation();
	void test_CCDFQCProbeSetInformation();
	void testmethod_IsXDACompatibleFile();
	void testmethod_Exists();
	void testmethod_ReadHeader_with_ASCII();
	void testmethod_ReadHeader_with_XDA();
	void test_ExpressionXDA();
	void test_ExpressionASCII();
	void test_GenotypingXDA();
	void test_GenotypingASCII();
	void test_ReseqXDA();
	void test_ReseqASCII();
	void test_TagXDA();
	void test_TagASCII();
	void test_DMET3XDA();
	void test_DMET3ASCII();
	void test_MultichannelXDA();
	void test_MultichannelASCII();
	void test_AxiomV4XDA();
	void test_AxiomV6ASCII();
	
};


#endif // __CDFFILEDATATEST_H_
