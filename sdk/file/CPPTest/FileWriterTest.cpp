////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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
#include "file/CPPTest/FileWriterTest.h"
//
#include "file/FileWriter.h"
//
#include <cmath>
#include <cstring>
#include <string.h>
#include <string>
//
#ifdef _MSC_VER
#include <direct.h>
#endif

using namespace std;

#define TEST_FILE "./data/write.test"

CPPUNIT_TEST_SUITE_REGISTRATION( FileWriterTest );

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

#define TEST_DATA_FILE "./data/test.file"

void FileWriterTest::setUp()
{
}

void FileWriterTest::tearDown()
{
}

void FileWriterTest::testfunction_WriteInt32_I()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	int32_t val1=123;
	int32_t val2=0;

	WriteInt32_I(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadInt32_I(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteUInt32_I()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	uint32_t val1=123;
	uint32_t val2=0;

	WriteUInt32_I(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadUInt32_I(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteInt16_I()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	int16_t val1=123;
	int16_t val2=0;

	WriteInt16_I(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadInt16_I(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteUInt16_I()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	uint16_t val1=123;
	uint16_t val2=0;

	WriteUInt16_I(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadUInt16_I(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteInt32_N()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	int32_t val1=123;
	int32_t val2=0;

	WriteInt32_N(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadInt32_N(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteUInt32_N()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	uint32_t val1=123;
	uint32_t val2=0;

	WriteUInt32_N(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadUInt32_N(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteInt16_N()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	int16_t val1=123;
	int16_t val2=0;

	WriteInt16_N(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadInt16_N(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteUInt16_N()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	uint16_t val1=123;
	uint16_t val2=0;

	WriteUInt16_N(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadUInt16_N(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteInt8()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	int8_t val1=123;
	int8_t val2=0;

	WriteInt8(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadInt8(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteUInt8()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	uint8_t val1=123;
	uint8_t val2=0;

	WriteUInt8(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadUInt8(in, val2);

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteFloat_I()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	float val1=123.123f;
	float val2=0;

	WriteFloat_I(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadFloat_I(in, val2);

	CPPUNIT_ASSERT(CompareFloats(val1, val2) == true);
}

void FileWriterTest::testfunction_WriteFloat_N()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	float val1=123.123f;
	float val2=0;

	WriteFloat_N(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadFloat_N(in, val2);

	CPPUNIT_ASSERT(CompareFloats(val1, val2) == true);
}

void FileWriterTest::testfunction_WriteFloatLowPrecision()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	float val1=123.123f;
	float val2=0;

	WriteFloatLowPrecision(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadFloat_I(in, val2);

	CPPUNIT_ASSERT(CompareFloats(val2, 123.1f) == true);
}

void FileWriterTest::testfunction_WriteFixedCString()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	const char *val1 = "test";
	char val2[64]="";

	WriteFixedCString(out, val1, strlen(val1));
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadFixedCString(in, val2, (int32_t)strlen(val1));

	CPPUNIT_ASSERT(strcmp(val1, val2) == 0);
}

void FileWriterTest::testfunction_WriteCString()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	const char *val1;
	char *val2=NULL;

	val1 = "test";
	WriteCString(out, val1);
	val1 = "";
	WriteCString(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadCString_I(in, val2);
	val1 = "test";
	CPPUNIT_ASSERT(strcmp(val1, val2) == 0);
	delete[] val2;

	ReadCString_I(in, val2);
	val1 = "";
	CPPUNIT_ASSERT(strcmp(val1, val2) == 0);
	delete[] val2;
}

void FileWriterTest::testfunction_WriteCharacterArray()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	char val1[] = "test";
	char *pval1 = val1;
	val1[2] = '\0';  // check for an embedded null...
	char val2[] = "xxxx";
	char *pval2 = val2;
 
	WriteCharacterArray(out, pval1, 4);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadCharacterArray(in, pval2, 4);

	CPPUNIT_ASSERT(memcmp(val1, val2, 4)==0);
}


void FileWriterTest::testfunction_WriteFixedString()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	string val1 = "test";
	string val2;
 
	WriteFixedString(out, val1, val1.length());
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadFixedString(in, val2, (int32_t)val1.length());

	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteString_I()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	string val1;
	string val2;

	val1 = "test";
	WriteString_I(out, val1);
	val1 = "";
	WriteString_I(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadString_I(in, val2);
	val1 = "test";
	CPPUNIT_ASSERT(val1 == val2);
	ReadString_I(in, val2);
	val1 = "";
	CPPUNIT_ASSERT(val1 == val2);
}

void FileWriterTest::testfunction_WriteString_N()
{
	ofstream out(TEST_FILE, ios::out | ios::binary);
	CPPUNIT_ASSERT( out != NULL);
	string val1;
	string val2;

	val1 = "test";
	WriteString_N(out, val1);
	val1 = "";
	WriteString_N(out, val1);
	out.close();

	ifstream in(TEST_FILE, ios::in | ios::binary);
	ReadString_N(in, val2);
	val1 = "test";
	CPPUNIT_ASSERT(val1 == val2);
	ReadString_N(in, val2);
	val1 = "";
	CPPUNIT_ASSERT(val1 == val2);
}

