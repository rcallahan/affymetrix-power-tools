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

#ifndef _GRIDCONTROLDATATEST_H_
#define _GRIDCONTROLDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class GridControlDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( GridControlDataTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testproperty_Rows );
	CPPUNIT_TEST ( testproperty_Columns );
	CPPUNIT_TEST ( testproperty_B1 );
	CPPUNIT_TEST ( testproperty_B2 );
	CPPUNIT_TEST ( testproperty_NS );
	CPPUNIT_TEST ( testmethod_Clear );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testproperty_Rows();
	void testproperty_Columns();
	void testproperty_B1();
	void testproperty_B2();
	void testproperty_NS();
	void testmethod_Clear();
};

#endif
