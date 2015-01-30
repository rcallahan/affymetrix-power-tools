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

#ifdef _MSC_VER
#include "windows.h"
#endif
//
#include "calvin_files/writers/test/FileOutputTest.h"
//
#include "calvin_files/parsers/src/FileInput.h"
#include "calvin_files/writers/src/FileOutput.h"
//
#include "util/Fs.h"
//
#include <cmath>
#include <cstring>
#include <istream>
#include <string.h>
#include <string>
//

#ifdef _MSC_VER
#pragma warning(disable: 4996) // ignore deprecated functions warning
#else
#include <sys/types.h>
#include <netinet/in.h>
#include <inttypes.h>
#endif

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( FileOutputTest );

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

#define TEST_DATA_FILE "./data.write.file"

void FileOutputTest::setUp()
{
}

void FileOutputTest::tearDown()
{
}

void FileOutputTest::testdefine_SIZES()
{
	CPPUNIT_ASSERT( sizeof(int8_t) == 1 );
	CPPUNIT_ASSERT( sizeof(int16_t) == 2 );
	CPPUNIT_ASSERT( sizeof(int32_t) == 4 );
	CPPUNIT_ASSERT( sizeof(int64_t) == 8 );
	CPPUNIT_ASSERT( sizeof(u_int8_t) == 1 );
	CPPUNIT_ASSERT( sizeof(u_int16_t) == 2 );
	CPPUNIT_ASSERT( sizeof(u_int32_t) == 4 );
	CPPUNIT_ASSERT( sizeof(u_int64_t) == 8 );
}

void FileOutputTest::testmethod_WriteInt8()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	int8_t val;

	val=11;
	FileOutput::WriteInt8(outstr, val);
	val=-22;
	FileOutput::WriteInt8(outstr, val);
	val=33;
	FileOutput::WriteInt8(outstr, val);
	outstr.close();

	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );
	instr.read((char *) &val, sizeof(val));
	CPPUNIT_ASSERT( val == 11 );
	instr.read((char *) &val, sizeof(val));
	CPPUNIT_ASSERT( val == -22 );
	instr.read((char *) &val, sizeof(val));
	CPPUNIT_ASSERT( val == 33 );
}

void FileOutputTest::testmethod_WriteUInt8()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	u_int8_t val;

	val=11;
	FileOutput::WriteUInt8(outstr, val);
	val=22;
	FileOutput::WriteUInt8(outstr, val);
	val=33;
	FileOutput::WriteUInt8(outstr, val);
	outstr.close();

	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );
	instr.read((char *) &val, sizeof(val));
	CPPUNIT_ASSERT( val == 11 );
	instr.read((char *) &val, sizeof(val));
	CPPUNIT_ASSERT( val == 22 );
	instr.read((char *) &val, sizeof(val));
	CPPUNIT_ASSERT( val == 33 );
}

void FileOutputTest::testmethod_WriteInt16()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	int16_t val;

	val=101;
	FileOutput::WriteInt16(outstr, val);
	val=-202;
	FileOutput::WriteInt16(outstr, val);
	val=303;
	FileOutput::WriteInt16(outstr, val);
	outstr.close();

	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );
	instr.read((char *) &val, sizeof(val));
	val = ntohs(val);
	CPPUNIT_ASSERT( val == 101 );
	instr.read((char *) &val, sizeof(val));
	val = ntohs(val);
	CPPUNIT_ASSERT( val == -202 );
	instr.read((char *) &val, sizeof(val));
	val = ntohs(val);
	CPPUNIT_ASSERT( val == 303 );
}

void FileOutputTest::testmethod_WriteUInt16()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	u_int16_t val;

	val=101;
	FileOutput::WriteUInt16(outstr, val);
	val=202;
	FileOutput::WriteUInt16(outstr, val);
	val=303;
	FileOutput::WriteUInt16(outstr, val);
	outstr.close();

	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );
	instr.read((char *) &val, sizeof(val));
	val = ntohs(val);
	CPPUNIT_ASSERT( val == 101 );
	instr.read((char *) &val, sizeof(val));
	val = ntohs(val);
	CPPUNIT_ASSERT( val == 202 );
	instr.read((char *) &val, sizeof(val));
	val = ntohs(val);
	CPPUNIT_ASSERT( val == 303 );
}

void FileOutputTest::testmethod_WriteInt32()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	int32_t val;

	val=101;
	FileOutput::WriteInt32(outstr, val);
	val=-202;
	FileOutput::WriteInt32(outstr, val);
	val=303;
	FileOutput::WriteInt32(outstr, val);
	outstr.close();

	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );
	instr.read((char *) &val, sizeof(val));
	val = ntohl(val);
	CPPUNIT_ASSERT( val == 101 );
	instr.read((char *) &val, sizeof(val));
	val = ntohl(val);
	CPPUNIT_ASSERT( val == -202 );
	instr.read((char *) &val, sizeof(val));
	val = ntohl(val);
	CPPUNIT_ASSERT( val == 303 );
}

void FileOutputTest::testmethod_WriteUInt32()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	u_int32_t val;

	val=101;
	FileOutput::WriteUInt32(outstr, val);
	val=202;
	FileOutput::WriteUInt32(outstr, val);
	val=303;
	FileOutput::WriteUInt32(outstr, val);
	outstr.close();

	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );
	instr.read((char *) &val, sizeof(val));
	val = ntohl(val);
	CPPUNIT_ASSERT( val == 101 );
	instr.read((char *) &val, sizeof(val));
	val = ntohl(val);
	CPPUNIT_ASSERT( val == 202 );
	instr.read((char *) &val, sizeof(val));
	val = ntohl(val);
	CPPUNIT_ASSERT( val == 303 );
}

void FileOutputTest::testmethod_WriteFloat()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	float val;

	val=1.1f;
	FileOutput::WriteFloat(outstr, val);
	val=-2.2f;
	FileOutput::WriteFloat(outstr, val);
	val=3.3f;
	FileOutput::WriteFloat(outstr, val);
	outstr.close();

	int32_t ival;
	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );
	instr.read((char *) &ival, sizeof(ival));
	ival = ntohl(ival);
    memcpy(&val, &ival, sizeof(ival));
	CPPUNIT_ASSERT( CompareFloats(val, 1.1f) );
	instr.read((char *) &ival, sizeof(ival));
	ival = ntohl(ival);
    memcpy(&val, &ival, sizeof(ival));
	CPPUNIT_ASSERT( CompareFloats(val, -2.2f) );
	instr.read((char *) &ival, sizeof(ival));
	ival = ntohl(ival);
    memcpy(&val, &ival, sizeof(ival));
	CPPUNIT_ASSERT( CompareFloats(val, 3.3f) );
}

void FileOutputTest::testmethod_WriteString8()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	string val;

	val = "string";
	FileOutput::WriteString8(outstr, val);
	val = "test";
	FileOutput::WriteString8(outstr, val);
	val = "case";
	FileOutput::WriteString8(outstr, val);
	outstr.close();

	int32_t ival;
	char *buf=NULL;
	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );

	instr.read((char *) &ival, sizeof(ival));
	ival = ntohl(ival);
	buf = new char[ival+1];
	instr.read((char *) buf, ival);
	buf[ival] = 0;
	CPPUNIT_ASSERT ( strcmp("string", buf) == 0);
	delete[] buf;

	instr.read((char *) &ival, sizeof(ival));
	ival = ntohl(ival);
	buf = new char[ival+1];
	instr.read((char *) buf, ival);
	buf[ival] = 0;
	CPPUNIT_ASSERT ( strcmp("test", buf) == 0);
	delete[] buf;

	instr.read((char *) &ival, sizeof(ival));
	ival = ntohl(ival);
	buf = new char[ival+1];
	instr.read((char *) buf, ival);
	buf[ival] = 0;
	CPPUNIT_ASSERT ( strcmp("case", buf) == 0);
	delete[] buf;
}

void FileOutputTest::testmethod_WriteString8_fixedlen()
{
	ofstream outstr;
	string file = TEST_DATA_FILE + string("string8");
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	string val;

	val = "string";
	FileOutput::WriteString8(outstr, val.c_str(), (int32_t) val.length() * 2);
	val = "test";
	FileOutput::WriteString8(outstr, val.c_str(), (int32_t) val.length() * 2);
	val = "case";
	FileOutput::WriteString8(outstr, val.c_str(), (int32_t) val.length() * 2);
	outstr.close();

	int32_t ival;
	char buf[128];
	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );

	ival = (int32_t) strlen("string");
	instr.read((char *) buf, ival * 2);
	buf[ival]=0;
	CPPUNIT_ASSERT ( strcmp("string", buf) == 0);

	ival = (int32_t) strlen("test");
	instr.read((char *) buf, ival * 2);
	buf[ival]=0;
	CPPUNIT_ASSERT ( strcmp("test", buf) == 0);

	ival = (int32_t) strlen("case");
	instr.read((char *) buf, ival * 2);
	buf[ival]=0;
	CPPUNIT_ASSERT ( strcmp("case", buf) == 0);

}

void FileOutputTest::testmethod_WriteString16()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	wstring val;

	val = L"string";
	FileOutput::WriteString16(outstr, val);
	val = L"test";
	FileOutput::WriteString16(outstr, val);
	val = L"case";
	FileOutput::WriteString16(outstr, val);
	outstr.close();

	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );

	val = FileInput::ReadString16(instr);
	CPPUNIT_ASSERT ( val == L"string");

	val = FileInput::ReadString16(instr);
	CPPUNIT_ASSERT ( val == L"test");

	val = FileInput::ReadString16(instr);
	CPPUNIT_ASSERT ( val == L"case");
}

void FileOutputTest::testmethod_WriteString16_fixedlen()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	wstring val;

	val = L"string";
	FileOutput::WriteString16(outstr, val.c_str(), (int32_t) val.length());
	val = L"test";
	FileOutput::WriteString16(outstr, val.c_str(), (int32_t) val.length());
	val = L"case";
	FileOutput::WriteString16(outstr, val.c_str(), (int32_t) val.length());
	outstr.close();


	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );

	val = FileInput::ReadString16(instr, (int32_t) strlen("string"));
	CPPUNIT_ASSERT ( val == L"string");

	val = FileInput::ReadString16(instr, (int32_t) strlen("test"));
	CPPUNIT_ASSERT ( val == L"test");

	val = FileInput::ReadString16(instr, (int32_t) strlen("case"));
	CPPUNIT_ASSERT ( val == L"case");
}

void FileOutputTest::testmethod_WriteBlob()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	wstring val;

	val = L"this is some random string";
	FileOutput::WriteBlob(outstr, (void*)val.c_str(), sizeof(wchar_t)*val.length());

	val = L"The Sharks are on the verge of eliminating the Predators";
	FileOutput::WriteBlob(outstr, (void*)val.c_str(), sizeof(wchar_t)*val.length());

	outstr.close();

	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );

	wstring expectedValue;
	const void* value = 0;
	int32_t size = FileInput::ReadBlob(instr, value);
	expectedValue = L"this is some random string";
	CPPUNIT_ASSERT(size == expectedValue.length()*sizeof(wchar_t));
	CPPUNIT_ASSERT(memcmp(value, (void*)expectedValue.c_str(), size) == 0);

	size = FileInput::ReadBlob(instr, value);
	expectedValue = L"The Sharks are on the verge of eliminating the Predators";
	CPPUNIT_ASSERT(size == expectedValue.length()*sizeof(wchar_t));
	CPPUNIT_ASSERT(memcmp(value, (void*)expectedValue.c_str(), size) == 0);
}

void FileOutputTest::testmethod_WriteBlobWithReserve()
{
	ofstream outstr;
	string file = TEST_DATA_FILE;
	Fs::aptOpen(outstr,file, ios::out);
	CPPUNIT_ASSERT( outstr.is_open() );
	wstring val;

	val = L"no sugar tonight - The Guess Who";
	FileOutput::WriteBlob(outstr, (void*)val.c_str(), sizeof(wchar_t)*val.length(), sizeof(wchar_t)*val.length()+10);

	val = L"pinball wizard - The Who";
	FileOutput::WriteBlob(outstr, (void*)val.c_str(), sizeof(wchar_t)*val.length(), sizeof(wchar_t)*val.length()+10);

	outstr.close();

	ifstream instr;
	Fs::aptOpen(instr, file, ios::in);
	CPPUNIT_ASSERT( instr.is_open() );

	wstring expectedValue;
	const void* value = 0;
	int32_t size = FileInput::ReadBlob(instr, value);
	expectedValue = L"no sugar tonight - The Guess Who";
	CPPUNIT_ASSERT(size == expectedValue.length()*sizeof(wchar_t)+10);
	CPPUNIT_ASSERT(memcmp(value, (void*)expectedValue.c_str(), expectedValue.length()*sizeof(wchar_t)) == 0);

	size = FileInput::ReadBlob(instr, value);
	expectedValue = L"pinball wizard - The Who";
	CPPUNIT_ASSERT(size == expectedValue.length()*sizeof(wchar_t)+10);
	CPPUNIT_ASSERT(memcmp(value, (void*)expectedValue.c_str(), expectedValue.length()*sizeof(wchar_t)) == 0);
}
