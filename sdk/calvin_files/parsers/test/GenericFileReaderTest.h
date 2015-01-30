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

#ifndef __GENERICFILEREADERTEST_H_
#define __GENERICFILEREADERTEST_H_

#include "calvin_files/data/src/GenericData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class GenericFileReaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( GenericFileReaderTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testproperty_FileName );
	CPPUNIT_TEST ( testmethod_ReadHeaderMinDP_when_file_exists );
	CPPUNIT_TEST ( testmethod_ReadHeaderMinDP_when_file_does_not_exist );
	CPPUNIT_TEST ( testmethod_ReadHeaderMinDP_when_file_is_not_valid );
	CPPUNIT_TEST ( testmethod_ReadHeaderMinDP_when_file_has_wrong_version );
	CPPUNIT_TEST ( testmethod_ReadHeader_when_file_exists );
	CPPUNIT_TEST ( testmethod_ReadHeader_when_file_does_not_exist );
	CPPUNIT_TEST ( testmethod_ReadHeader_when_file_is_not_valid );
	CPPUNIT_TEST ( testmethod_ReadHeader_when_file_has_wrong_version );
	CPPUNIT_TEST ( OpenTest );
	CPPUNIT_TEST ( GetDataGroupCntTest );
	CPPUNIT_TEST ( GetDataGroupReaderByIndexTest );
	CPPUNIT_TEST ( GetDataGroupReaderByNameTest );
	CPPUNIT_TEST ( CloseNoOpenTest );
	CPPUNIT_TEST ( ReadHeaderNoDataGroupHeaderTest );
	CPPUNIT_TEST ( ReadHeaderOfAFileWithMultipleDataGroups );
	CPPUNIT_TEST ( ReadHeaderOfFileWithReserveStringParameters );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testproperty_FileName();
	void testmethod_ReadHeaderMinDP_when_file_exists();
	void testmethod_ReadHeaderMinDP_when_file_does_not_exist();
	void testmethod_ReadHeaderMinDP_when_file_is_not_valid();
	void testmethod_ReadHeaderMinDP_when_file_has_wrong_version();
	void testmethod_ReadHeader_when_file_exists();
	void testmethod_ReadHeader_when_file_does_not_exist();
	void testmethod_ReadHeader_when_file_is_not_valid();
	void testmethod_ReadHeader_when_file_has_wrong_version();
	void OpenTest();
	void GetDataGroupCntTest();
	void GetDataGroupReaderByIndexTest();
	void GetDataGroupReaderByNameTest();
	void CloseNoOpenTest();
	void ReadHeaderNoDataGroupHeaderTest();
	void ReadHeaderOfAFileWithMultipleDataGroups();
	void ReadHeaderOfFileWithReserveStringParameters();

public:
	void Check_TEST_DATA_DAT_FILE_GenericHeader(affymetrix_calvin_io::GenericData& data);
};

#endif // __GENERICFILEREADERTEST_H_
