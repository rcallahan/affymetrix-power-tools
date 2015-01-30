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

#ifndef __CELFILEWRITERTEST_H_
#define __CELFILEWRITERTEST_H_

#include "file/CELFileWriter.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affxcel;

class CCELFileWriterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CCELFileWriterTest );
	
	CPPUNIT_TEST(testWrite_v3);
	CPPUNIT_TEST(testWrite_v4);
	CPPUNIT_TEST(testWrite_bcel);
	CPPUNIT_TEST(testWrite_ccel);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testWrite_v3();
	void testWrite_v4();
	void testWrite_bcel();
	void testWrite_ccel();

private:
	void setData(CCELFileWriter& cel);
};

#endif // __CELFILEWRITERTEST_H_
