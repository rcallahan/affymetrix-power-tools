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


#pragma once

#include "calvin_files/fusion/src/FusionArrayFileReader.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_fusion_io;

/*! Test Fixture for testing the Array File Reader. */
class FusionArrayFileReaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( FusionArrayFileReaderTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testmethod_Read_invalid_file );
	CPPUNIT_TEST ( testmethod_Read_missing_file );
	CPPUNIT_TEST ( testmethod_Read_calvin_array_file );
	CPPUNIT_TEST ( testmethod_Read_dtt_array_file );
	CPPUNIT_TEST ( testmethod_Read_mas_array_file );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	/*! Test the ability to create the ArrayFileReader class. */
	void testCreation();
	/*! Test the handling of attempting to read an invalid array file. */
	void testmethod_Read_invalid_file();
	/*! Test the handling of attempting to read an array file that does not exist. */
	void testmethod_Read_missing_file();
	/*! Test the ability to read a valid Calvin array file. */
	void testmethod_Read_calvin_array_file();
	/*! Test the ability to read an Data Transfer Tool formated array file. */
	void testmethod_Read_dtt_array_file();
	/*! Test the ability to read an Micro Array Suit formated array file (EXP).  */
	void testmethod_Read_mas_array_file();
};
