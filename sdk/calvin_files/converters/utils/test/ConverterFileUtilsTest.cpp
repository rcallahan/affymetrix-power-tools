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
#include "calvin_files/converters/utils/test/ConverterFileUtilsTest.h"
//
#include "calvin_files/converters/utils/src/ConverterFileUtils.h"
//
#include "util/Fs.h"
//
#include <cstdio>
#include <fstream>
//

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( ConverterFileUtilsTest );


void ConverterFileUtilsTest::setUp()
{
	CPPUNIT_ASSERT(ConverterRemoveFile("busy") == true);
	CPPUNIT_ASSERT(ConverterRemoveFile("busy2") == true);
}

void ConverterFileUtilsTest::tearDown()
{
	CPPUNIT_ASSERT(ConverterRemoveFile("busy") == true);
	CPPUNIT_ASSERT(ConverterRemoveFile("busy2") == true);
}

void ConverterFileUtilsTest::testMove()
{
	ofstream str("busy", ios::out);
	str << "busy" << endl;
	CPPUNIT_ASSERT(ConverterFileExists("busy") == true);
#ifdef WIN32
	CPPUNIT_ASSERT(ConverterMoveFile("busy", "busy2", false) == false);
#endif
	str.close();
	CPPUNIT_ASSERT(ConverterMoveFile("busy", "busy2", false) == true);
	CPPUNIT_ASSERT(ConverterFileExists("busy2") == true);
	CPPUNIT_ASSERT(ConverterFileExists("busy") == false);

        Fs::aptOpen(str, "busy", ios::out);
	str << "busy" << endl;
	str.close();
	CPPUNIT_ASSERT(ConverterMoveFile("busy2", "busy", false) == false);
	CPPUNIT_ASSERT(ConverterMoveFile("busy2", "busy", true) == true);
}

void ConverterFileUtilsTest::testRemove()
{
	ofstream str("busy", ios::out);
	str << "busy" << endl;
	CPPUNIT_ASSERT(ConverterFileExists("busy") == true);
#ifdef WIN32
	CPPUNIT_ASSERT(ConverterRemoveFile("busy") == false);
#endif
	str.close();
	CPPUNIT_ASSERT(ConverterRemoveFile("busy") == true);
	CPPUNIT_ASSERT(ConverterFileExists("busy") == false);
}

void ConverterFileUtilsTest::testExists()
{
	CPPUNIT_ASSERT(ConverterFileExists("ConverterFileUtilsTest.cpp") == true);
	CPPUNIT_ASSERT(ConverterFileExists("ConverterFileUtilsTest.no_file") == false);
}

