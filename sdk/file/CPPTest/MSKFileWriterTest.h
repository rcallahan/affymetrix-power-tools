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

#ifndef __MSKFILEWRITERTEST_H_
#define __MSKFILEWRITERTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CMSKFileWriterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CMSKFileWriterTest );
	
	CPPUNIT_TEST( testmethod_Save_with_call );
	CPPUNIT_TEST( testmethod_Save_with_comp );
	CPPUNIT_TEST( testmethod_Save_with_both );
	CPPUNIT_TEST( testmethod_Save_empty );
	CPPUNIT_TEST( testmethod_Save_no_array_type );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testmethod_Save_with_call();
	void testmethod_Save_with_comp();
	void testmethod_Save_with_both();
	void testmethod_Save_empty();
	void testmethod_Save_no_array_type();
};


#endif // __MSKFILEWRITERTEST_H_
