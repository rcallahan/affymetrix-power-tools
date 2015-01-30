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

#include <cppunit/extensions/HelperMacros.h>

class TemplateFileReaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( TemplateFileReaderTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testmethod_Read );
	CPPUNIT_TEST ( testmethod_Read_header_data_only );
	CPPUNIT_TEST ( testmethod_Read_when_file_does_not_exist );
	CPPUNIT_TEST ( testmethod_Read_when_file_is_not_valid );
	CPPUNIT_TEST ( testmethod_IsFileType );
	CPPUNIT_TEST ( testmethod_DataTypeIdentifierStatic );

	CPPUNIT_TEST_SUITE_END();

	void CheckTemplateData(bool headerOnly);

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testmethod_Read();
	void testmethod_Read_header_data_only();
	void testmethod_Read_when_file_does_not_exist();
	void testmethod_Read_when_file_is_not_valid();
	void testmethod_IsFileType();
	void testmethod_DataTypeIdentifierStatic();
};
