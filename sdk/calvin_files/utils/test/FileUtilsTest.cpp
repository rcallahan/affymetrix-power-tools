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

//
#include "calvin_files/utils/test/FileUtilsTest.h"
//
#include "calvin_files/utils/src/FileUtils.h"
//
#include <cstdio>
#include <fstream>
#include <map>
//

using namespace std;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( FileUtilsTest );

#define TEST_NON_EXISTANT_FILE "./no_file"
#define TEST_FILE "./FileUtilsTest.cpp"

void FileUtilsTest::setUp()
{
}

void FileUtilsTest::tearDown()
{
}

void FileUtilsTest::testmethod_Exists_when_file_exists()
{
	CPPUNIT_ASSERT( FileUtils::Exists(TEST_FILE) == true );
}

void FileUtilsTest::testmethod_Exists_when_file_does_not_exist()
{
	CPPUNIT_ASSERT( FileUtils::Exists(TEST_NON_EXISTANT_FILE) == false );
}

void FileUtilsTest::testmethod_LockFile_success()
{
	const char *fileName = TEST_FILE;
	bool status = FileUtils::LockFile(fileName);
	CPPUNIT_ASSERT( status == true );
	status = FileUtils::UnlockFile(fileName);
	CPPUNIT_ASSERT( status == true );
}

void FileUtilsTest::testmethod_LockFile_double_lock_failure()
{
	const char *fileName = TEST_FILE;
	bool status = FileUtils::LockFile(fileName);
	CPPUNIT_ASSERT( status == true );
	status = FileUtils::LockFile(fileName);
	CPPUNIT_ASSERT( status == false );
	status = FileUtils::UnlockFile(fileName);
	CPPUNIT_ASSERT( status == true );
}

void FileUtilsTest::testmethod_LockFile_file_does_not_exist()
{
	const char *fileName = TEST_NON_EXISTANT_FILE;
	bool status = FileUtils::LockFile(fileName);
	CPPUNIT_ASSERT( status == false );
}

void FileUtilsTest::testmethod_UnlockFile_file_does_not_exist()
{
	const char *fileName = TEST_NON_EXISTANT_FILE;
	bool status = FileUtils::UnlockFile(fileName);
	CPPUNIT_ASSERT( status == true );
}

// rely on string pasting to set the path to the directory of test files.
#define LISTFILES_DIR(x) "../data/listfiles/" x

void FileUtilsTest::testmethod_ListFiles()
{
	const char *path = LISTFILES_DIR("");
	list<string> files;
	list<string>::iterator it;
	string file;
    string filesarray[3];

	files = FileUtils::ListFiles(path, "a");
	CPPUNIT_ASSERT( files.size() == 3 );
    filesarray[0] = "";
    filesarray[1] = "";
    filesarray[2] = "";
	for (it = files.begin(); it!= files.end(); ++it)
    {
        file = *it;
        if (file == LISTFILES_DIR("file1.a"))
            filesarray[0] = file;

        else if (file == LISTFILES_DIR("file2.a"))
            filesarray[1] = file;

        else if (file == LISTFILES_DIR("file3.a"))
            filesarray[2] = file;
    }
	CPPUNIT_ASSERT( filesarray[0] == LISTFILES_DIR("file1.a") );
	CPPUNIT_ASSERT( filesarray[1] == LISTFILES_DIR("file2.a") );
	CPPUNIT_ASSERT( filesarray[2] == LISTFILES_DIR("file3.a") );


	files = FileUtils::ListFiles(path, "b");
	CPPUNIT_ASSERT( files.size() == 2 );    
    filesarray[0] = "";
    filesarray[1] = "";
	for (it = files.begin(); it!= files.end(); ++it)
    {
        file = *it;
        if (file == LISTFILES_DIR("file1.b") )
            filesarray[0] = file;

        else if (file == LISTFILES_DIR("file2.b") )
            filesarray[1] = file;
    }
	CPPUNIT_ASSERT( filesarray[0] == LISTFILES_DIR("file1.b") );
	CPPUNIT_ASSERT( filesarray[1] == LISTFILES_DIR("file2.b") );


	files = FileUtils::ListFiles(path, "c");
	CPPUNIT_ASSERT( files.size() == 1 );
	it = files.begin();
	file = *it;
  CPPUNIT_ASSERT( file == LISTFILES_DIR("file1.c") );
	++it;
	CPPUNIT_ASSERT( it == files.end() );

	files = FileUtils::ListFiles(path, "d");
	CPPUNIT_ASSERT( files.size() == 0 );

	files = FileUtils::ListFiles(path, "");
	CPPUNIT_ASSERT( files.size() == 6 );
}

void FileUtilsTest::testmethod_ListFiles_path_does_not_exist()
{
	const char *path = "../path_does_not_exist";
	list<string> files = FileUtils::ListFiles(path, "a");
	CPPUNIT_ASSERT( files.size() == 0 );
}

void FileUtilsTest::testmethod_RemoveFile()
{
	const char *fname = "./remove_file";
	ofstream o(fname, ios::out);
	o << "test" << endl;
	o.close();

	CPPUNIT_ASSERT(FileUtils::Exists(fname) == true);
	CPPUNIT_ASSERT(FileUtils::RemoveFile(fname) == true);
	CPPUNIT_ASSERT(FileUtils::Exists(fname) == false);
}
