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

#ifndef __CHPFILEWRITERTEST_H_
#define __CHPFILEWRITERTEST_H_

#include "file/CHPFileWriter.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affxchpwriter;

class CCHPFileWriterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CCHPFileWriterTest );
	
	CPPUNIT_TEST( testproperty_Header_AlgName );
	CPPUNIT_TEST( testproperty_Header_AlgVersion );
	CPPUNIT_TEST( testproperty_Header_AlgParams );
	CPPUNIT_TEST( testproperty_Header_SummaryParams );
	CPPUNIT_TEST( testproperty_Header_ParentCellFile );
	CPPUNIT_TEST( testproperty_Header_ProgID );
	CPPUNIT_TEST( testproperty_Header_BackgroundZoneInfo );
	CPPUNIT_TEST( testfunction_Write_Exp_Comp);
	CPPUNIT_TEST( testfunction_Write_Exp_Abs);
	CPPUNIT_TEST( testfunction_Write_Genotyping);
	CPPUNIT_TEST( testfunction_Write_Universal);
	CPPUNIT_TEST( testfunction_Write_Resequencing);
	CPPUNIT_TEST( testassignment_CResequencingResults );

	CPPUNIT_TEST( testfunction_Write_Exp_Abs_NoBuffer);
	CPPUNIT_TEST( testfunction_Write_Genotyping_NoBuffer);
	CPPUNIT_TEST( testfunction_Write_Universal_NoBuffer);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testproperty_Header_AlgName();
	void testproperty_Header_AlgVersion();
	void testproperty_Header_AlgParams();
	void testproperty_Header_SummaryParams();
	void testproperty_Header_ParentCellFile();
	void testproperty_Header_ProgID();
	void testproperty_Header_BackgroundZoneInfo();
	void testfunction_Write_Exp_Comp();
	void testfunction_Write_Exp_Abs();
	void testassignment_CResequencingResults();
	void testfunction_Write_Genotyping();
	void testfunction_Write_Universal();
	void testfunction_Write_Resequencing();
	void testfunction_Write_Exp_Abs_NoBuffer();
	void testfunction_Write_Genotyping_NoBuffer();
	void testfunction_Write_Universal_NoBuffer();

private:
	void setData(CCHPFileWriter& chp, bool bComparison);
};


#endif // __CHPFILEWRITERTEST_H_
