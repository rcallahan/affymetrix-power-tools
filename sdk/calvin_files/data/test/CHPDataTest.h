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

#ifndef __CHPDATATEST_H_
#define __CHPDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CHPDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CHPDataTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( FilenameTest );
	CPPUNIT_TEST ( SetEntryCountTest );
	CPPUNIT_TEST ( ArrayTypeTest );
	CPPUNIT_TEST ( AlgParamTest );
	CPPUNIT_TEST ( ChipSumTest );
	CPPUNIT_TEST ( AlgVersionTest );
	CPPUNIT_TEST ( ProgIdTest );
	CPPUNIT_TEST ( AlgNameTest );
	CPPUNIT_TEST ( RowsTest );
	CPPUNIT_TEST ( ColsTest );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void FilenameTest();
	void SetEntryCountTest();
	void ArrayTypeTest();
	void AlgParamTest();
	void ChipSumTest();
	void AlgVersionTest();
	void ProgIdTest();
	void AlgNameTest();
	void RowsTest();
	void ColsTest();
};

#endif // __CHPDATATEST_H_
