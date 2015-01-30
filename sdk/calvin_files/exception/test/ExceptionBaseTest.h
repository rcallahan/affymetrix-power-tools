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

#ifndef __EXCEPTIONBASETEST_H_
#define __EXCEPTIONBASETEST_H_

#include <cppunit/extensions/HelperMacros.h>

class ExceptionTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( ExceptionTest );

	CPPUNIT_TEST( testCreate );
	CPPUNIT_TEST( testToString );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreate();
	void testToString();
};


#endif // __EXCEPTIONBASETEST_H_
