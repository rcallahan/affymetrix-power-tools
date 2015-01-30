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
#include "file/CPPTest/FileIOTest.h"
//
#include "file/FileIO.h"
//
#include <cmath>
#include <cstring>
#include <istream>
#include <string.h>
#include <string>
//
#ifdef _MSC_VER
#include <direct.h>
#endif


using namespace std;

#ifdef _MSC_VER
#pragma warning(disable: 4996) // don't show deprecated warnings.
#endif

CPPUNIT_TEST_SUITE_REGISTRATION( CFileIOTest );

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

#define TEST_DATA_FILE "./data/test.file"

void CFileIOTest::setUp()
{
}

void CFileIOTest::tearDown()
{
}

void CFileIOTest::testdefine_SIZES()
{
	CPPUNIT_ASSERT( INT32_SIZE == 4 );
	CPPUNIT_ASSERT( INT16_SIZE == 2 );
	CPPUNIT_ASSERT( INT8_SIZE == 1 );
	CPPUNIT_ASSERT(	UINT32_SIZE == 4 );
	CPPUNIT_ASSERT( UINT16_SIZE == 2 );
	CPPUNIT_ASSERT( UINT8_SIZE == 1 );
	CPPUNIT_ASSERT( FLOAT_SIZE == 4 );
	CPPUNIT_ASSERT( INT32_SIZE  == 4 );
	CPPUNIT_ASSERT( LONG_SIZE  == 4 );
	CPPUNIT_ASSERT( ULONG_SIZE  == 4 );
	CPPUNIT_ASSERT( INT16_SIZE  == 2 );
	CPPUNIT_ASSERT( UINT16_SIZE  == 2 );
	CPPUNIT_ASSERT( INT8_SIZE  == 1 );
	CPPUNIT_ASSERT( UINT8_SIZE  == 1 );
}

void CFileIOTest::testfunction_affy_swap()
{
	const int istartSwapped = -91003647;
	int istartValue = 23434234;
	int isrc = istartValue;
	int idest = 0;
	idest = affy_swap32(isrc);
	CPPUNIT_ASSERT( idest == istartSwapped );
	isrc = affy_swap32(idest);
	CPPUNIT_ASSERT( isrc == istartValue );

	const uint16_t sstartSwapped = 255;
	uint16_t sstartValue = 65280;
	uint16_t ssrc = sstartValue;
	uint16_t sdest = 0;
	sdest = affy_swap16(ssrc);
	CPPUNIT_ASSERT( sdest == sstartSwapped );
	ssrc = affy_swap16(sdest);
	CPPUNIT_ASSERT( ssrc == sstartValue );

}

void CFileIOTest::testfunction_ReadInt32_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".int");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	int32_t ival=0;
	ReadInt32_I(instr, ival);
	CPPUNIT_ASSERT(ival == 10 );
	ReadInt32_I(instr, ival);
	CPPUNIT_ASSERT(ival == 20 );
	ReadInt32_I(instr, ival);
	CPPUNIT_ASSERT(ival == 30 );
	instr.close();
}

void CFileIOTest::testfunction_ReadUInt32_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".uint");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint32_t ival=0;
	ReadUInt32_I(instr, ival);
	CPPUNIT_ASSERT(ival == 101 );
	ReadUInt32_I(instr, ival);
	CPPUNIT_ASSERT(ival == 202 );
	ReadUInt32_I(instr, ival);
	CPPUNIT_ASSERT(ival == 303 );
	instr.close();
}

void CFileIOTest::testfunction_ReadInt16_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".short");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	short val=0;
	ReadInt16_I(instr, val);
	CPPUNIT_ASSERT(val == -21 );
	ReadInt16_I(instr, val);
	CPPUNIT_ASSERT(val == 21 );
	ReadInt16_I(instr, val);
	CPPUNIT_ASSERT(val == -31 );
	instr.close();
}

void CFileIOTest::testfunction_ReadUInt16_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".ushort");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint16_t val=0;
	ReadUInt16_I(instr, val);
	CPPUNIT_ASSERT(val == 45 );
	ReadUInt16_I(instr, val);
	CPPUNIT_ASSERT(val == 54 );
	ReadUInt16_I(instr, val);
	CPPUNIT_ASSERT(val == 66 );
	instr.close();
}

void CFileIOTest::testfunction_ReadInt8()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".char");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	int8_t val=0;
	ReadInt8(instr, val);
	CPPUNIT_ASSERT(val == -11 );
	ReadInt8(instr, val);
	CPPUNIT_ASSERT(val == -22 );
	ReadInt8(instr, val);
	CPPUNIT_ASSERT(val == -33 );
	instr.close();
}

void CFileIOTest::testfunction_ReadUInt8()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".uchar");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint8_t val=0;
	ReadUInt8(instr, val);
	CPPUNIT_ASSERT(val == 11 );
	ReadUInt8(instr, val);
	CPPUNIT_ASSERT(val == 22 );
	ReadUInt8(instr, val);
	CPPUNIT_ASSERT(val == 33 );
	instr.close();
}

void CFileIOTest::testfunction_ReadFloat_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".float");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	float val=0.0f;
	ReadFloat_I(instr, val);
	CPPUNIT_ASSERT( CompareFloats(val, 321.123f) );
	ReadFloat_I(instr, val);
	CPPUNIT_ASSERT( CompareFloats(val, 123.123f) );
	ReadFloat_I(instr, val);
	CPPUNIT_ASSERT( CompareFloats(val, 543.56f) );
	instr.close();
}

void CFileIOTest::testfunction_ReadCString_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	char *val = NULL;
	ReadCString_I(instr, val);
	CPPUNIT_ASSERT( val == string("string") );
	delete[] val;
	val = NULL;
	ReadCString_I(instr, val);
	CPPUNIT_ASSERT( val == string("test") );
	delete[] val;
	val = NULL;
	ReadCString_I(instr, val);
	CPPUNIT_ASSERT( val == string("case") );
	delete[] val;
	val = NULL;
	instr.close();
}

void CFileIOTest::testfunction_ReadString_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	string val;
	ReadString_I(instr, val);
	CPPUNIT_ASSERT( val == string("string") );
	ReadString_I(instr, val);
	CPPUNIT_ASSERT( val == string("test") );
	ReadString_I(instr, val);
	CPPUNIT_ASSERT( val == string("case") );
	instr.close();
}

void CFileIOTest::testfunction_ReadUIntLenString_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	string val;
	ReadString_I(instr, val);
	CPPUNIT_ASSERT( val == string("string") );
	ReadString_I(instr, val);
	CPPUNIT_ASSERT( val == string("test") );
	ReadString_I(instr, val);
	CPPUNIT_ASSERT( val == string("case") );
	instr.close();
}

void CFileIOTest::testfunction_ReadInt32_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_int");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	int32_t ival=0;
	ReadInt32_N(instr, ival);
	CPPUNIT_ASSERT(ival == 10 );
	ReadInt32_N(instr, ival);
	CPPUNIT_ASSERT(ival == 20 );
	ReadInt32_N(instr, ival);
	CPPUNIT_ASSERT(ival == 30 );
	instr.close();
}

void CFileIOTest::testfunction_ReadUInt32_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_uint");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint32_t ival=0;
	ReadUInt32_N(instr, ival);
	CPPUNIT_ASSERT(ival == 101 );
	ReadUInt32_N(instr, ival);
	CPPUNIT_ASSERT(ival == 202 );
	ReadUInt32_N(instr, ival);
	CPPUNIT_ASSERT(ival == 303 );
	instr.close();
}

void CFileIOTest::testfunction_ReadInt16_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_short");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	int16_t val=0;
	ReadInt16_N(instr, val);
	CPPUNIT_ASSERT(val == -21 );
	ReadInt16_N(instr, val);
	CPPUNIT_ASSERT(val == 21 );
	ReadInt16_N(instr, val);
	CPPUNIT_ASSERT(val == -31 );
	instr.close();
}

void CFileIOTest::testfunction_ReadUInt16_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_ushort");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint16_t val=0;
	ReadUInt16_N(instr, val);
	CPPUNIT_ASSERT(val == 45 );
	ReadUInt16_N(instr, val);
	CPPUNIT_ASSERT(val == 54 );
	ReadUInt16_N(instr, val);
	CPPUNIT_ASSERT(val == 66 );
	instr.close();
}

void CFileIOTest::testfunction_ReadFloat_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_float");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	float val=0.0f;
	ReadFloat_N(instr, val);
	CPPUNIT_ASSERT( CompareFloats(val, 321.123f) );
	ReadFloat_N(instr, val);
	CPPUNIT_ASSERT( CompareFloats(val, 123.123f) );
	ReadFloat_N(instr, val);
	CPPUNIT_ASSERT( CompareFloats(val, 543.56f) );
	instr.close();
}

void CFileIOTest::testfunction_ReadCString_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_string");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	char *val = NULL;
	ReadCString_N(instr, val);
	CPPUNIT_ASSERT( val == string("string") );
	delete[] val;
	val = NULL;
	ReadCString_N(instr, val);
	CPPUNIT_ASSERT( val == string("test") );
	delete[] val;
	val = NULL;
	ReadCString_N(instr, val);
	CPPUNIT_ASSERT( val == string("case") );
	delete[] val;
	val = NULL;
	instr.close();
}

void CFileIOTest::testfunction_ReadString_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_string");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	string val;
	ReadString_N(instr, val);
	CPPUNIT_ASSERT( val == string("string") );
	ReadString_N(instr, val);
	CPPUNIT_ASSERT( val == string("test") );
	ReadString_N(instr, val);
	CPPUNIT_ASSERT( val == string("case") );
	instr.close();
}

void CFileIOTest::testfunction_ReadNextLine_dosfile()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".dos");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	const int len = 100;
	char val[len];
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("line 1") );
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("line 2") );
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("line 3") );
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("") );
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("") );
	instr.close();
}

void CFileIOTest::testfunction_ReadNextLine_unixfile()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".unix");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	const int len = 100;
	char val[len];
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("line 1") );
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("line 2") );
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("line 3") );
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("") );
	ReadNextLine(instr, val, len);
	CPPUNIT_ASSERT( val == string("") );
	instr.close();
}

void CFileIOTest::testfunction_ReadFloatFromOldBPMAP()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".bpmap_float");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	float val=0.0f;
	ReadFloatFromOldBPMAP_N(instr, val);
	CPPUNIT_ASSERT( CompareFloats(val, 1.0f) );
	instr.close();
}


void CFileIOTest::testfunction_MmGetFloatFromOldBPMAP()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".bpmap_float");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	const float origVal=1.0f;
	float val=0.0f;
	float res=0.0f;
	instr.read((char *)&val, FLOAT_SIZE);
	res = MmGetFloatFromOldBPMAP_N((float *)&val);
	CPPUNIT_ASSERT( CompareFloats(res, origVal) );
}

void CFileIOTest::testfunction_ReadFixedCString()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".fixed_string");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	char val[7];
	ReadFixedCString(instr, val, 5);
	CPPUNIT_ASSERT( strncmp(val, "fixed", strlen("fixed")) == 0 );
	ReadFixedCString(instr, val, 6);
	CPPUNIT_ASSERT( strncmp(val, "string", strlen("string")) == 0 );
	instr.close();
}

void CFileIOTest::testfunction_ReadFixedUCString()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".fixed_string");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	unsigned char val[7];
	ReadFixedUCString(instr, val, 5);
	CPPUNIT_ASSERT( strncmp((char *)val, "fixed", strlen("fixed")) == 0 );
	ReadFixedUCString(instr, val, 6);
	CPPUNIT_ASSERT( strncmp((char *)val, "string", strlen("string")) == 0 );
	instr.close();
}

void CFileIOTest::testfunction_ReadFixedString()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".fixed_string");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	string val;
	ReadFixedString(instr, val, 5);
	CPPUNIT_ASSERT( strncmp(val.c_str(), "fixed", strlen("fixed")) == 0 );
	ReadFixedString(instr, val, 6);
	CPPUNIT_ASSERT( strncmp(val.c_str(), "string", strlen("string")) == 0 );
	instr.close();
}

void CFileIOTest::testfunction_ReadUIntLenString_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_string");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	string val;
	ReadUIntLenString_N(instr, val);
	CPPUNIT_ASSERT( val == string("string") );
	ReadUIntLenString_N(instr, val);
	CPPUNIT_ASSERT( val == string("test") );
	ReadUIntLenString_N(instr, val);
	CPPUNIT_ASSERT( val == string("case") );
	instr.close();
}

void CFileIOTest::testfunction_MmGetInt32_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".int");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	int32_t vals[3];
	instr.read((char *)vals, 3*INT32_SIZE);
	int32_t val=0;
	val = MmGetInt32_I(&vals[0]);
	CPPUNIT_ASSERT(val == 10 );
	val = MmGetInt32_I(&vals[1]);
	CPPUNIT_ASSERT(val == 20 );
	val = MmGetInt32_I(&vals[2]);
	CPPUNIT_ASSERT(val == 30 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetUInt32_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".uint");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint32_t vals[3];
	instr.read((char *)vals, 3*INT32_SIZE);
	uint32_t val=0;
	val = MmGetUInt32_I(&vals[0]);
	CPPUNIT_ASSERT(val == 101 );
	val = MmGetUInt32_I(&vals[1]);
	CPPUNIT_ASSERT(val == 202 );
	val = MmGetUInt32_I(&vals[2]);
	CPPUNIT_ASSERT(val == 303 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetInt16_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".short");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	int16_t vals[3];
	instr.read((char *)vals, 3*INT16_SIZE);
	int16_t val=0;
	val = MmGetInt16_I(&vals[0]);
	CPPUNIT_ASSERT(val == -21 );
	val = MmGetInt16_I(&vals[1]);
	CPPUNIT_ASSERT(val == 21 );
	val = MmGetInt16_I(&vals[2]);
	CPPUNIT_ASSERT(val == -31 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetUInt16_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".ushort");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint16_t vals[3];
	instr.read((char *)vals, 3*INT16_SIZE);
	uint16_t val=0;
	val = MmGetUInt16_I(&vals[0]);
	CPPUNIT_ASSERT(val == 45 );
	val = MmGetUInt16_I(&vals[1]);
	CPPUNIT_ASSERT(val == 54 );
	val = MmGetUInt16_I(&vals[2]);
	CPPUNIT_ASSERT(val == 66 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetInt8()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".char");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	int8_t vals[3];
	instr.read((char *)vals, 3*INT8_SIZE);
	int8_t val=0;
	val = MmGetInt8(&vals[0]);
	CPPUNIT_ASSERT(val == -11 );
	val = MmGetInt8(&vals[1]);
	CPPUNIT_ASSERT(val == -22 );
	val = MmGetInt8(&vals[2]);
	CPPUNIT_ASSERT(val == -33 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetUInt8()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".uchar");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint8_t vals[3];
	instr.read((char *)vals, 3*INT8_SIZE);
	uint8_t val=0;
	val = MmGetUInt8(&vals[0]);
	CPPUNIT_ASSERT(val == 11 );
	val = MmGetUInt8(&vals[1]);
	CPPUNIT_ASSERT(val == 22 );
	val = MmGetUInt8(&vals[2]);
	CPPUNIT_ASSERT(val == 33 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetFloat_I()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".float");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	float vals[3];
	instr.read((char *)vals, 3*FLOAT_SIZE);
	float val=0.0f;
	val = MmGetFloat_I(&vals[0]);
	CPPUNIT_ASSERT( CompareFloats(val, 321.123f) );
	val = MmGetFloat_I(&vals[1]);
	CPPUNIT_ASSERT( CompareFloats(val, 123.123f) );
	val = MmGetFloat_I(&vals[2]);
	CPPUNIT_ASSERT( CompareFloats(val, 543.56f) );
	instr.close();
}

void CFileIOTest::testfunction_MmGetInt32_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_int");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	int32_t vals[3];
	instr.read((char *)vals, 3*INT32_SIZE);
	int32_t val=0;
	val = MmGetInt32_N(&vals[0]);
	CPPUNIT_ASSERT(val == 10 );
	val = MmGetInt32_N(&vals[1]);
	CPPUNIT_ASSERT(val == 20 );
	val = MmGetInt32_N(&vals[2]);
	CPPUNIT_ASSERT(val == 30 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetUInt32_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_uint");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint32_t vals[3];
	instr.read((char *)vals, 3*INT32_SIZE);
	uint32_t val=0;
	val = MmGetUInt32_N(&vals[0]);
	CPPUNIT_ASSERT(val == 101 );
	val = MmGetUInt32_N(&vals[1]);
	CPPUNIT_ASSERT(val == 202 );
	val = MmGetUInt32_N(&vals[2]);
	CPPUNIT_ASSERT(val == 303 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetInt16_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_short");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	int16_t vals[3];
	instr.read((char *)vals, 3*INT16_SIZE);
	int16_t val=0;
	val = MmGetInt16_N(&vals[0]);
	CPPUNIT_ASSERT(val == -21 );
	val = MmGetInt16_N(&vals[1]);
	CPPUNIT_ASSERT(val == 21 );
	val = MmGetInt16_N(&vals[2]);
	CPPUNIT_ASSERT(val == -31 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetUInt16_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_ushort");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	uint16_t vals[3];
	instr.read((char *)vals, 3*INT16_SIZE);
	uint16_t val=0;
	val = MmGetUInt16_N(&vals[0]);
	CPPUNIT_ASSERT(val == 45 );
	val = MmGetUInt16_N(&vals[1]);
	CPPUNIT_ASSERT(val == 54 );
	val = MmGetUInt16_N(&vals[2]);
	CPPUNIT_ASSERT(val == 66 );
	instr.close();
}

void CFileIOTest::testfunction_MmGetFloat_N()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".big_float");
	instr.open(file.c_str(), ios::in);
	CPPUNIT_ASSERT( instr );
	float vals[3];
	instr.read((char *)vals, 3*FLOAT_SIZE);
	float val=0.0f;
	val = MmGetFloat_N(&vals[0]);
	CPPUNIT_ASSERT( CompareFloats(val, 321.123f) );
	val = MmGetFloat_N(&vals[1]);
	CPPUNIT_ASSERT( CompareFloats(val, 123.123f) );
	val = MmGetFloat_N(&vals[2]);
	CPPUNIT_ASSERT( CompareFloats(val, 543.56f) );
	instr.close();
}

void CFileIOTest::testfunction_MmSetUInt32_I()
{
	uint32_t val1;
	uint32_t val2;
	val1 = 123;
	MmSetUInt32_I(&val2, val1);
	CPPUNIT_ASSERT(MmGetUInt32_I(&val2) == val1);
}

void CFileIOTest::testfunction_MmSetUInt16_I()
{
	uint16_t val1;
	uint16_t val2;
	val1 = 123;
	MmSetUInt16_I(&val2, val1);
	CPPUNIT_ASSERT(MmGetUInt16_I(&val2) == val1);
}

void CFileIOTest::testfunction_MmSetUInt8()
{
	uint8_t val1;
	uint8_t val2;
	val1 = 123;
	MmSetUInt8(&val2, val1);
	CPPUNIT_ASSERT(MmGetUInt8(&val2) == val1);
}

void CFileIOTest::testfunction_MmSetFloat_I()
{
	float val1;
	float val2;
	val1 = 123.0f;
	MmSetFloat_I(&val2, val1);
	CPPUNIT_ASSERT(CompareFloats(MmGetFloat_I(&val2), val1));
}

void CFileIOTest::testfunction_MmSetUInt32_N()
{
	uint32_t val1;
	uint32_t val2;
	val1 = 123;
	MmSetUInt32_I(&val2, val1);
	CPPUNIT_ASSERT(MmGetUInt32_I(&val2) == val1);
}

void CFileIOTest::testfunction_MmSetUInt16_N()
{
	uint16_t val1;
	uint16_t val2;
	val1 = 123;
	MmSetUInt16_I(&val2, val1);
	CPPUNIT_ASSERT(MmGetUInt16_I(&val2) == val1);
}

void CFileIOTest::testfunction_MmSetFloat_N()
{
	float val1;
	float val2;
	val1 = 123.0f;
	MmSetFloat_N(&val2, val1);
	CPPUNIT_ASSERT(CompareFloats(MmGetFloat_N(&val2), val1));
}
