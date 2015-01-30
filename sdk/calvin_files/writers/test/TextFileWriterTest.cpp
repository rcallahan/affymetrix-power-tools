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
#include "calvin_files/writers/test/TextFileWriterTest.h"
//
#include "calvin_files/writers/src/TextFileWriter.h"
//
#include <cstring>
#include <fstream>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;

CPPUNIT_TEST_SUITE_REGISTRATION( TextFileWriterTest );

void TextFileWriterTest::setUp()
{
}

void TextFileWriterTest::tearDown()
{
}

void TextFileWriterTest::testCreation()
{
	TextFileWriter writer;
	CPPUNIT_ASSERT(1);
}

void TextFileWriterTest::testmethod_WriteFile()
{
	TextFileWriter writer;
	map<string, string> params;
	std::string name = "test.file.text";
	params["a"] = "b";
	params["b"] = "c";
	params["c"] = "d";
	writer.WriteFile(name, params);

	const int linelen=32;
	char line[linelen];
	std::ifstream instr(name.c_str());
	instr.getline(line, linelen);
	CPPUNIT_ASSERT( strncmp(line, "a=b", 3) == 0 );
	instr.getline(line, linelen);
	CPPUNIT_ASSERT( strncmp(line, "b=c", 3) == 0 );
	instr.getline(line, linelen);
	CPPUNIT_ASSERT( strncmp(line, "c=d", 3) == 0 );
}

