////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#ifndef __SATUINT16TEST_H_
#define __SATUINT16TEST_H_

#include <cppunit/extensions/HelperMacros.h>

class Satuint16Test : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( Satuint16Test );
	
	CPPUNIT_TEST( test_Satuint16 );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void test_Satuint16();
};

#endif // __SATUINT16TEST_H_
