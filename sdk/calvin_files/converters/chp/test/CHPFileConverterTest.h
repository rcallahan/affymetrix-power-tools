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

#ifndef __CHPFILECONVERTERTEST_H_
#define __CHPFILECONVERTERTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class CHPFileConverterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CHPFileConverterTest );

	CPPUNIT_TEST ( test_ConvertFile_to_Calvin );
	CPPUNIT_TEST ( test_Convert_Mas5_File_to_Calvin );
	CPPUNIT_TEST ( test_ConvertFile_to_XDA );
	CPPUNIT_TEST ( test_ConvertFile_fail_to_open_file );
	CPPUNIT_TEST ( test_ConvertFile_already_in_format );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void test_ConvertFile_to_Calvin();
	void test_Convert_Mas5_File_to_Calvin();
	void test_ConvertFile_to_XDA();
	void test_ConvertFile_fail_to_open_file();
	void test_ConvertFile_already_in_format();
};

#endif // __CHPFILECONVERTERTEST_H_
