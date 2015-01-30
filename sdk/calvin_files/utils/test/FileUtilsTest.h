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

#ifndef __FILEUTILSTEST_H_
#define __FILEUTILSTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class FileUtilsTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( FileUtilsTest );

	CPPUNIT_TEST ( testmethod_Exists_when_file_exists );
	CPPUNIT_TEST ( testmethod_Exists_when_file_does_not_exist );
	CPPUNIT_TEST ( testmethod_LockFile_success );
	CPPUNIT_TEST ( testmethod_LockFile_double_lock_failure );
	CPPUNIT_TEST ( testmethod_LockFile_file_does_not_exist );
	CPPUNIT_TEST ( testmethod_UnlockFile_file_does_not_exist );
	CPPUNIT_TEST ( testmethod_ListFiles );
	CPPUNIT_TEST ( testmethod_ListFiles_path_does_not_exist );
	CPPUNIT_TEST ( testmethod_RemoveFile );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testmethod_Exists_when_file_exists();
	void testmethod_Exists_when_file_does_not_exist();
	void testmethod_LockFile_success();
	void testmethod_LockFile_double_lock_failure();
	void testmethod_LockFile_file_does_not_exist();
	void testmethod_UnlockFile_file_does_not_exist();
	void testmethod_ListFiles();
	void testmethod_ListFiles_path_does_not_exist();
	void testmethod_RemoveFile();
};

#endif // __FILEUTILSTEST_H_
