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
#ifndef __BEDFILEDATATEST_H_
#define __BEDFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class BEDFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( BEDFileDataTest );

	CPPUNIT_TEST ( testIntervalGroupClear );
	CPPUNIT_TEST ( testFormatTrack );
	CPPUNIT_TEST ( testFormatBrowser );
	CPPUNIT_TEST ( testRead );
	CPPUNIT_TEST ( testRead_with_params );
	CPPUNIT_TEST ( testRead_with_no_params );
	CPPUNIT_TEST ( testRead_no_file );
	CPPUNIT_TEST ( testClear );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testIntervalGroupClear();
	void testFormatTrack();
	void testFormatBrowser();
	void testRead();
	void testRead_with_params();
	void testRead_with_no_params();
	void testRead_no_file();
	void testClear();
};

#endif // __BEDFILEDATATEST_H_
