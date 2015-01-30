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

#ifndef __CHPTILINGDATATEST_H_
#define __CHPTILINGDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CHPTilingDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CHPTilingDataTest );

	CPPUNIT_TEST ( test_FileName );
	CPPUNIT_TEST ( test_NumberSequences );
	CPPUNIT_TEST ( test_AlgName );
	CPPUNIT_TEST ( test_AlgVersion );
	CPPUNIT_TEST ( test_AlgParams );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void test_FileName();
	void test_NumberSequences();
	void test_AlgName();
	void test_AlgVersion();
	void test_AlgParams();
};

#endif // __CHPTILINGDATATEST_H_
