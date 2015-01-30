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
#ifndef __DTTARRAYFILEREADERTEST_H_
#define __DTTARRAYFILEREADERTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class DttArrayFileReaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( DttArrayFileReaderTest );

	CPPUNIT_TEST ( testproperty_ArrayType);
	CPPUNIT_TEST ( testproperty_ExpName);
	CPPUNIT_TEST ( testproperty_Attributes);
	CPPUNIT_TEST ( testmethod_Clear);
	CPPUNIT_TEST ( testmethod_Exists_pass );
	CPPUNIT_TEST ( testmethod_Exists_fail );
	CPPUNIT_TEST ( testmethod_Read );
	CPPUNIT_TEST ( testmethod_Read_when_file_does_not_exist );
	CPPUNIT_TEST ( testmethod_Read_when_file_is_not_valid );
	CPPUNIT_TEST ( testmethod_Read_from_gdac_exporter_sample);
	CPPUNIT_TEST ( testmethod_Read_from_gdac_exporter_exp);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testproperty_ArrayType();
	void testproperty_ExpName();
	void testproperty_Attributes();
	void testmethod_Clear();

	void testmethod_Exists_pass();
	void testmethod_Exists_fail();
	void testmethod_Read();
	void testmethod_Read_when_file_does_not_exist();
	void testmethod_Read_when_file_is_not_valid();
	void testmethod_Read_from_gdac_exporter_sample();
	void testmethod_Read_from_gdac_exporter_exp();
};

#endif // __DTTARRAYFILEREADERTEST_H_
