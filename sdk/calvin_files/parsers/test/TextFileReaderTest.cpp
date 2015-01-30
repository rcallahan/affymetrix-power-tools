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
#include "calvin_files/parsers/test/TextFileReaderTest.h"
//
#include "calvin_files/parsers/src/TextFileReader.h"
//
#include <cstring>
#include <fstream>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;

CPPUNIT_TEST_SUITE_REGISTRATION( TextFileReaderTest );

void TextFileReaderTest::setUp()
{
}

void TextFileReaderTest::tearDown()
{
}

void TextFileReaderTest::testCreation()
{
	TextFileReader array;
	CPPUNIT_ASSERT(1);
}

void TextFileReaderTest::testmethod_ReadFile_file_does_not_exist()
{
	TextFileReader reader;
	map<string, string> params;
	std::string name = "file_does_not_exist";
	try
	{
		reader.ReadFile(name, params);
	}
	catch (FileNotFoundException)
	{
		CPPUNIT_ASSERT(1);
	}
	catch (...)
	{
		CPPUNIT_ASSERT(0);
	}
}

void TextFileReaderTest::testmethod_ReadFile()
{
	TextFileReader reader;
	map<string, string> params;
	std::string name = "../data/test.file.text";
	reader.ReadFile(name, params);

	string value;

	CPPUNIT_ASSERT(params.size()==3);

	value = params["a"];
	CPPUNIT_ASSERT( value == "b" );
	value = params["b"];
	CPPUNIT_ASSERT( value == "c" );
	value = params["c"];
	CPPUNIT_ASSERT( value == "d" );
}

void TextFileReaderTest::testmethod_ReadFile_stream()
{
	TextFileReader reader;
	map<string, string> params;
	std::string name = "../data/test.file.text";
	std::ifstream stream(name.c_str(), std::ios::in);
	reader.ReadFile(stream, params);
	stream.close();

	string value;

	CPPUNIT_ASSERT(params.size()==3);

	value = params["a"];
	CPPUNIT_ASSERT( value == "b" );
	value = params["b"];
	CPPUNIT_ASSERT( value == "c" );
	value = params["c"];
	CPPUNIT_ASSERT( value == "d" );
}
