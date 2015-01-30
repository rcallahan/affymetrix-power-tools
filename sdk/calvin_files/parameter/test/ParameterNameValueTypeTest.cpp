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
#include "calvin_files/parameter/test/ParameterNameValueTypeTest.h"
//
#include "calvin_files/parameter/src/ParameterException.h"
//
#include <cstring>
#include <string.h>
#include <string>
//

#ifndef _MSC_VER
#include <sys/types.h>
#include <netinet/in.h>
#include <inttypes.h>
#endif

const u_int32_t NUMBER_BUFFER_LEN = 16;

// Copied from Parameter.cpp
const wchar_t Int8MIMEType[] = L"text/x-calvin-integer-8";
const wchar_t UInt8MIMEType[] = L"text/x-calvin-unsigned-integer-8";
const wchar_t Int16MIMEType[] = L"text/x-calvin-integer-16";
const wchar_t UInt16MIMEType[] = L"text/x-calvin-unsigned-integer-16";
const wchar_t Int32MIMEType[] = L"text/x-calvin-integer-32";
const wchar_t UInt32MIMEType[] = L"text/x-calvin-unsigned-integer-32";
const wchar_t FloatMIMEType[] = L"text/x-calvin-float";
const wchar_t TextMIMEType[] = L"text/plain";
const wchar_t AsciiMIMEType[] = L"text/ascii";


CPPUNIT_TEST_SUITE_REGISTRATION( ParameterNameValueTypeTest );

using namespace affymetrix_calvin_parameter;

void ParameterNameValueTypeTest::setUp()
{
}

void ParameterNameValueTypeTest::tearDown()
{
}

void ParameterNameValueTypeTest::CreationTest()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT(1);
}

void ParameterNameValueTypeTest::ConstructorOfInt8TypeTest()
{
	MIMEValue v;
	FillMIMEValue(-5, v);
	ParameterNameValueType nvt(L"fred", v, Int8MIMEType);
	CPPUNIT_ASSERT(nvt.GetName() == L"fred");
	CPPUNIT_ASSERT(nvt.GetValueInt8() == -5);
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::Int8Type);
	CPPUNIT_ASSERT_THROW(nvt.GetValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);
}

void ParameterNameValueTypeTest::ConstructorOfUInt8TypeTest()
{
	MIMEValue v;
	FillMIMEValue(5, v);
	ParameterNameValueType nvt(L"wilma", v, UInt8MIMEType);
	CPPUNIT_ASSERT(nvt.GetName() == L"wilma");
	CPPUNIT_ASSERT(nvt.GetValueUInt8() == 5);
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::UInt8Type);
	CPPUNIT_ASSERT_THROW(nvt.GetValueInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);
}

void ParameterNameValueTypeTest::ConstructorOfInt16TypeTest()
{
	MIMEValue v;
	FillMIMEValue(-23245, v);
	ParameterNameValueType nvt(L"beckham", v, Int16MIMEType);
	CPPUNIT_ASSERT(nvt.GetName() == L"beckham");
	CPPUNIT_ASSERT(nvt.GetValueInt16() == -23245);
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::Int16Type);
	CPPUNIT_ASSERT_THROW(nvt.GetValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);
}

void ParameterNameValueTypeTest::ConstructorOfUInt16TypeTest()
{
	MIMEValue v;
	FillMIMEValue(54444, v);
	ParameterNameValueType nvt(L"Ronaldinho", v, UInt16MIMEType);
	CPPUNIT_ASSERT(nvt.GetName() == L"Ronaldinho");
	CPPUNIT_ASSERT(nvt.GetValueUInt16() == 54444);
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::UInt16Type);
	CPPUNIT_ASSERT_THROW(nvt.GetValueInt16(), affymetrix_calvin_exceptions::ParameterMismatchException);
}

void ParameterNameValueTypeTest::ConstructorOfInt32TypeTest()
{
	MIMEValue v;
	FillMIMEValue(-23245345, v);
	ParameterNameValueType nvt(L"raul", v, Int32MIMEType);
	CPPUNIT_ASSERT(nvt.GetName() == L"raul");
	CPPUNIT_ASSERT(nvt.GetValueInt32() == -23245345);
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::Int32Type);
	CPPUNIT_ASSERT_THROW(nvt.GetValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);
}

void ParameterNameValueTypeTest::ConstructorOfUInt32TypeTest()
{
	MIMEValue v;
	FillMIMEValue(4000000000, v);
	ParameterNameValueType nvt(L"ZIzou", v, UInt32MIMEType);
	CPPUNIT_ASSERT(nvt.GetName() == L"ZIzou");
	CPPUNIT_ASSERT(nvt.GetValueUInt32() == 4000000000);
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::UInt32Type);
	CPPUNIT_ASSERT_THROW(nvt.GetValueInt16(), affymetrix_calvin_exceptions::ParameterMismatchException);
}

void ParameterNameValueTypeTest::ConstructorOfFloatTypeTest()
{
	MIMEValue v;
	float f = 2.34567f;
	type_punned pun;
	pun.v_float=f;
	FillMIMEValue(pun.v_uint32, v);
	ParameterNameValueType nvt(L"Henry", v, FloatMIMEType);
	CPPUNIT_ASSERT(nvt.GetName() == L"Henry");
	CPPUNIT_ASSERT(nvt.GetValueFloat() == 2.34567f);
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::FloatType);
	CPPUNIT_ASSERT_THROW(nvt.GetValueInt16(), affymetrix_calvin_exceptions::ParameterMismatchException);
}

void ParameterNameValueTypeTest::ConstructorOfTextTypeTest()
{
	MIMEValue v;
	ConvertWStringToMIMEText(v, L"Arsenal Football Club");
	ParameterNameValueType nvt(L"Viera", v, TextMIMEType);
	CPPUNIT_ASSERT(nvt.GetName() == L"Viera");
	CPPUNIT_ASSERT(nvt.GetValueText() == L"Arsenal Football Club");
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::TextType);
	CPPUNIT_ASSERT_THROW(nvt.GetValueInt16(), affymetrix_calvin_exceptions::ParameterMismatchException);
}

void ParameterNameValueTypeTest::ConstructorOfAsciiTypeTest()
{
	MIMEValue v;
	ConvertStringToMIMEText(v, "Arsenal Football Club");
	ParameterNameValueType nvt(L"Viera", v, AsciiMIMEType);
	CPPUNIT_ASSERT(nvt.GetName() == L"Viera");
	CPPUNIT_ASSERT(nvt.GetValueAscii() == "Arsenal Football Club");
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::AsciiType);
	CPPUNIT_ASSERT_THROW(nvt.GetValueInt16(), affymetrix_calvin_exceptions::ParameterMismatchException);
}

void ParameterNameValueTypeTest::CopyConstructorTest()
{
	MIMEValue v;
	float f = 2.34567f;
	type_punned pun;
	pun.v_float=f;
	FillMIMEValue(pun.v_uint32, v);

	ParameterNameValueType nvt(L"Henry", v, FloatMIMEType);
	ParameterNameValueType nvtCopy(nvt);
	CPPUNIT_ASSERT(nvt.GetName() == L"Henry");
	CPPUNIT_ASSERT(nvt.GetValueFloat() == 2.34567f);
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::FloatType);
}

void ParameterNameValueTypeTest::AssignmentOperatorTest()
{
	MIMEValue v;
	ConvertWStringToMIMEText(v, L"Arsenal Football Club");
	ParameterNameValueType nvt(L"Viera", v, TextMIMEType);
	ParameterNameValueType nvtAssign;
	nvtAssign = nvt;

	CPPUNIT_ASSERT(nvt.GetName() == L"Viera");
	CPPUNIT_ASSERT(nvt.GetValueText() == L"Arsenal Football Club");
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::TextType);
}

void ParameterNameValueTypeTest::EqualityOperatorTest()
{
	MIMEValue v;
	ConvertWStringToMIMEText(v, L"Arsenal Football Club");
	ParameterNameValueType nvt(L"Viera", v, TextMIMEType);

	MIMEValue v2;
	FillMIMEValue(6, v2);
	ParameterNameValueType nvtSameName(L"Viera", v2, Int32MIMEType);
	CPPUNIT_ASSERT(nvt == nvtSameName);
	CPPUNIT_ASSERT(nvt == L"Viera");
}

void ParameterNameValueTypeTest::InequalityOperatorTest()
{
	MIMEValue v;
	ConvertWStringToMIMEText(v, L"Arsenal Football Club");
	ParameterNameValueType nvt(L"Viera", v, TextMIMEType);
	ParameterNameValueType nvtDifferentName(L"Robben", v, TextMIMEType);
	CPPUNIT_ASSERT(nvt != nvtDifferentName);
	CPPUNIT_ASSERT(nvt != L"Lampard");
}

void ParameterNameValueTypeTest::LessThanOperatorTest()
{
	MIMEValue v;
	ConvertWStringToMIMEText(v, L"Arsenal Football Club");
	ParameterNameValueType nvt(L"Viera", v, TextMIMEType);

	MIMEValue v2;
	FillMIMEValue(6, v2);
	ParameterNameValueType nvtSameName(L"Viera", v2, Int32MIMEType);
	ParameterNameValueType nvtDifferentName(L"Robben", v, TextMIMEType);
	CPPUNIT_ASSERT_ASSERTION_FAIL(CPPUNIT_ASSERT(nvt < nvtSameName));
	CPPUNIT_ASSERT(nvtDifferentName < nvt);
}

void ParameterNameValueTypeTest::GreaterThanOperatorTest()
{
	MIMEValue v;
	ConvertWStringToMIMEText(v, L"Arsenal Football Club");
	ParameterNameValueType nvt(L"Viera", v, TextMIMEType);

	MIMEValue v2;
	FillMIMEValue(6, v2);
	ParameterNameValueType nvtSameName(L"Viera", v2, Int32MIMEType);
	ParameterNameValueType nvtDifferentName(L"Robben", v, TextMIMEType);
	CPPUNIT_ASSERT_ASSERTION_FAIL(CPPUNIT_ASSERT(nvt > nvtSameName));
	CPPUNIT_ASSERT(nvt > nvtDifferentName);
}

void ParameterNameValueTypeTest::NameTest()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetName(L"PS3"));
	CPPUNIT_ASSERT(nvt.GetName() == L"PS3");
}

void ParameterNameValueTypeTest::GetTypeUnknownTest()
{
	MIMEValue v;
	ConvertWStringToMIMEText(v, L"Liverpool Football Club");
	ParameterNameValueType nvt(L"Gerard", v, L"Whatever");
	CPPUNIT_ASSERT(nvt.GetParameterType() == ParameterNameValueType::UnknownType);
}

void ParameterNameValueTypeTest::ValueInt8Test()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueInt8(127));
	CPPUNIT_ASSERT(nvt.GetValueInt8() == 127);
	CPPUNIT_ASSERT_THROW(nvt.GetValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(127, v);
	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ValueUInt8Test()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueUInt8(255));
	CPPUNIT_ASSERT(nvt.GetValueUInt8() == 255);
	CPPUNIT_ASSERT_THROW(nvt.GetValueUInt16(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(255, v);
	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ValueInt16Test()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueInt16(-30002));
	CPPUNIT_ASSERT(nvt.GetValueInt16() == -30002);
	CPPUNIT_ASSERT_THROW(nvt.GetValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(-30002, v);
	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ValueUInt16Test()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueUInt16(55213));
	CPPUNIT_ASSERT(nvt.GetValueUInt16() == 55213);
	CPPUNIT_ASSERT_THROW(nvt.GetValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(55213, v);
	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ValueInt32Test()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueInt32(-521423654));
	CPPUNIT_ASSERT(nvt.GetValueInt32() == -521423654);
	CPPUNIT_ASSERT_THROW(nvt.GetValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(-521423654, v);
	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ValueUInt32Test()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueUInt32(3000123456));
	CPPUNIT_ASSERT(nvt.GetValueUInt32() == 3000123456);
	CPPUNIT_ASSERT_THROW(nvt.GetValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(3000123456, v);
	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ValueFloatTest()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueFloat(874065.7845f));
	CPPUNIT_ASSERT(nvt.GetValueFloat() == 874065.7845f);
	CPPUNIT_ASSERT_THROW(nvt.GetValueText(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	float f = 874065.7845f;
	type_punned pun;
	pun.v_float=f;
	FillMIMEValue(pun.v_uint32, v);

	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ValueTextTest()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueText(L"The Los Angeles Angels of Anaheim"));
	CPPUNIT_ASSERT(nvt.GetValueText() == L"The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT_THROW(nvt.GetValueFloat(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	ConvertWStringToMIMEText(v, L"The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ValueAsciiTest()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueAscii("The Los Angeles Angels of Anaheim"));
	CPPUNIT_ASSERT(nvt.GetValueAscii() == "The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT_THROW(nvt.GetValueFloat(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	ConvertStringToMIMEText(v, "The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ValueTextWithReserveTest()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueText(L"The Los Angeles Angels of Anaheim", 100));
	CPPUNIT_ASSERT(nvt.GetValueText() == L"The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT(nvt.GetMIMEValue().Size() == 100*sizeof(u_int16_t));
}

void ParameterNameValueTypeTest::ValueAsciiWithReserveTest()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetValueAscii("The Los Angeles Angels of Anaheim", 100));
	CPPUNIT_ASSERT(nvt.GetValueAscii() == "The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT(nvt.GetMIMEValue().Size() == 100);
}

void ParameterNameValueTypeTest::MIMETypeTest()
{
	ParameterNameValueType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetMIMEType(L"text/xml"));
	CPPUNIT_ASSERT(nvt.GetMIMEType() == L"text/xml");
}

void ParameterNameValueTypeTest::MIMEValueTest()
{
	ParameterNameValueType nvt;
	char buf[4];
	buf[0] = 4;
	buf[1] = 127;
	buf[2] = 88;
	buf[3] = 1;
	MIMEValue v(buf, 4);

	CPPUNIT_ASSERT_NO_THROW(nvt.SetMIMEValue(v));
	CPPUNIT_ASSERT(nvt.GetMIMEValue() == v);
}

void ParameterNameValueTypeTest::ToStringTest()
{
	ParameterNameValueType nvt;
	nvt.SetName(L"fred");
	nvt.SetValueInt8(1);
	CPPUNIT_ASSERT(nvt.ToString() == L"1");
	nvt.SetValueUInt8(2);
	CPPUNIT_ASSERT(nvt.ToString() == L"2");
	nvt.SetValueInt16(3);
	CPPUNIT_ASSERT(nvt.ToString() == L"3");
	nvt.SetValueUInt16(4);
	CPPUNIT_ASSERT(nvt.ToString() == L"4");
	nvt.SetValueInt32(5);
	CPPUNIT_ASSERT(nvt.ToString() == L"5");
	nvt.SetValueUInt32(6);
	CPPUNIT_ASSERT(nvt.ToString() == L"6");
	nvt.SetValueFloat(1015.0f);
	CPPUNIT_ASSERT(nvt.ToString() == L"1015.000000");
	nvt.SetValueText(L"test text 16");
	CPPUNIT_ASSERT(nvt.ToString() == L"test text 16");
	nvt.SetValueAscii("test text 8");
	CPPUNIT_ASSERT(nvt.ToString() == L"test text 8");

	MIMEValue mimeValue;
	mimeValue.SetValue(0, 0);
	nvt.SetMIMEType(L"foo");
	nvt.SetMIMEValue(mimeValue);
	CPPUNIT_ASSERT(nvt.ToString() == L"");
}

void ParameterNameValueTypeTest::ConvertWStringToMIMEText(MIMEValue& value, const std::wstring& s)
{
	u_int32_t len = (u_int32_t) s.length();
	u_int16_t* buf = new u_int16_t[len];
	for (u_int32_t i=0; i<len; ++i)
	{
		buf[i] = (u_int16_t) s[i];
		buf[i] = htons(buf[i]);
	}
	value.SetValue(buf, len*sizeof(u_int16_t));
	delete[] buf;
}

void ParameterNameValueTypeTest::ConvertStringToMIMEText(MIMEValue& value, const std::string& s)
{
	u_int32_t len = (u_int32_t) s.length();
	value.SetValue(s.c_str(), len);
}

void ParameterNameValueTypeTest::FillMIMEValue(u_int32_t value, MIMEValue& mime)
{
	char buf[NUMBER_BUFFER_LEN];
	memset(buf, 0, NUMBER_BUFFER_LEN);	// for aesthetics
	value = htonl(value);
	memcpy(buf, &value, sizeof(u_int32_t));
	mime.SetValue(buf, NUMBER_BUFFER_LEN);
}

void ParameterNameValueTypeTest::testproperties_ParameterNameValueDefaultRequiredType()
{
	ParameterNameValueDefaultRequiredType param;
	MIMEValue v;
	MIMEValue d;
	FillMIMEValue(-5, v);
	FillMIMEValue(-10, d);

	param.SetMIMEValue(v);
	param.DefaultMIMEValue() = d;
	CPPUNIT_ASSERT(param.DefaultMIMEValue() == d);
	CPPUNIT_ASSERT(param.GetMIMEValue() == v);

	param.RequiredFlag() = true;
	CPPUNIT_ASSERT(param.RequiredFlag() == true);

	param.RequiredFlag() = false;
	CPPUNIT_ASSERT(param.RequiredFlag() == false);
}

void ParameterNameValueTypeTest::testconstructors_ParameterNameValueDefaultRequiredType()
{
	int32_t v1=-5;
	int32_t d1=-10;
	MIMEValue v;
	MIMEValue d;
	v.SetValue(&v1, sizeof(v1));
	d.SetValue(&d1, sizeof(d1));
	ParameterNameValueDefaultRequiredType p1(L"fred", v, Int32MIMEType, d, true);
	CPPUNIT_ASSERT(p1.DefaultMIMEValue() == d);
	CPPUNIT_ASSERT(p1.GetMIMEValue() == v);
	CPPUNIT_ASSERT(p1.RequiredFlag() == true);

	ParameterNameValueDefaultRequiredType p2(L"fred", &v1, sizeof(int32_t), Int32MIMEType, &d1, sizeof(int32_t), false);
	CPPUNIT_ASSERT(p2.DefaultMIMEValue() == d);
	CPPUNIT_ASSERT(p2.GetMIMEValue() == v);
	CPPUNIT_ASSERT(p2.RequiredFlag() == false);
}

void ParameterNameValueTypeTest::testassignment_ParameterNameValueDefaultRequiredType()
{
	MIMEValue v;
	MIMEValue d;
	FillMIMEValue(-5, v);
	FillMIMEValue(-10, d);
	ParameterNameValueDefaultRequiredType p1(L"fred", v, Int8MIMEType, d, true);

	ParameterNameValueDefaultRequiredType p2 = p1;
	ParameterNameValueDefaultRequiredType p3;
	p3 = p2;

	CPPUNIT_ASSERT(p3.GetMIMEValue() == p1.GetMIMEValue());
	CPPUNIT_ASSERT(p3.DefaultMIMEValue() == p1.DefaultMIMEValue());
	CPPUNIT_ASSERT(p3.GetMIMEValue() == p1.GetMIMEValue());
	CPPUNIT_ASSERT(p3.RequiredFlag() == p1.RequiredFlag());
}

void ParameterNameValueTypeTest::testToString_ParameterNameValueDefaultRequiredType()
{
	ParameterNameValueDefaultRequiredType nvt;

	nvt.SetName(L"fred");
	nvt.SetValueInt8(10);
	nvt.SetDefaultValueInt8(1);
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"1");
	nvt.SetValueUInt8(12);
	nvt.SetDefaultValueUInt8(2);
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"2");
	nvt.SetValueInt16(13);
	nvt.SetDefaultValueInt16(3);
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"3");
	nvt.SetValueUInt16(14);
	nvt.SetDefaultValueUInt16(4);
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"4");
	nvt.SetValueInt32(15);
	nvt.SetDefaultValueInt32(5);
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"5");
	nvt.SetValueUInt32(16);
	nvt.SetDefaultValueUInt32(6);
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"6");
	nvt.SetValueFloat(15.0f);
	nvt.SetDefaultValueFloat(1015.0f);
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"1015.000000");
	nvt.SetValueText(L"test text 16b");
	nvt.SetDefaultValueText(L"test text 16");
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"test text 16");
	nvt.SetValueAscii("test text 8b");
	nvt.SetDefaultValueAscii("test text 8");
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"test text 8");

	MIMEValue mimeValue;
	mimeValue.SetValue(0, 0);
	nvt.SetMIMEType(L"foo");
	nvt.DefaultMIMEValue() = mimeValue;
	CPPUNIT_ASSERT(nvt.DefaultToString() == L"");
}

void ParameterNameValueTypeTest::DefaultValueInt8Test()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueInt8(127));
	CPPUNIT_ASSERT(nvt.GetDefaultValueInt8() == 127);
	CPPUNIT_ASSERT_THROW(nvt.GetDefaultValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(127, v);
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue() == v);
}

void ParameterNameValueTypeTest::DefaultValueUInt8Test()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueUInt8(255));
	CPPUNIT_ASSERT(nvt.GetDefaultValueUInt8() == 255);
	CPPUNIT_ASSERT_THROW(nvt.GetDefaultValueUInt16(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(255, v);
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue() == v);
}

void ParameterNameValueTypeTest::DefaultValueInt16Test()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueInt16(-30002));
	CPPUNIT_ASSERT(nvt.GetDefaultValueInt16() == -30002);
	CPPUNIT_ASSERT_THROW(nvt.GetDefaultValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(-30002, v);
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue() == v);
}

void ParameterNameValueTypeTest::DefaultValueUInt16Test()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueUInt16(55213));
	CPPUNIT_ASSERT(nvt.GetDefaultValueUInt16() == 55213);
	CPPUNIT_ASSERT_THROW(nvt.GetDefaultValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(55213, v);
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue() == v);
}

void ParameterNameValueTypeTest::DefaultValueInt32Test()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueInt32(-521423654));
	CPPUNIT_ASSERT(nvt.GetDefaultValueInt32() == -521423654);
	CPPUNIT_ASSERT_THROW(nvt.GetDefaultValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(-521423654, v);
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue() == v);
}

void ParameterNameValueTypeTest::DefaultValueUInt32Test()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueUInt32(3000123456));
	CPPUNIT_ASSERT(nvt.GetDefaultValueUInt32() == 3000123456);
	CPPUNIT_ASSERT_THROW(nvt.GetDefaultValueUInt8(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	FillMIMEValue(3000123456, v);
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue() == v);
}

void ParameterNameValueTypeTest::DefaultValueFloatTest()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueFloat(874065.7845f));
	CPPUNIT_ASSERT(nvt.GetDefaultValueFloat() == 874065.7845f);
	CPPUNIT_ASSERT_THROW(nvt.GetDefaultValueText(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	float f = 874065.7845f;
	type_punned pun;
	pun.v_float=f;
	FillMIMEValue(pun.v_uint32, v);

	CPPUNIT_ASSERT(nvt.DefaultMIMEValue() == v);
}

void ParameterNameValueTypeTest::DefaultValueTextTest()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueText(L"The Los Angeles Angels of Anaheim"));
	CPPUNIT_ASSERT(nvt.GetDefaultValueText() == L"The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT_THROW(nvt.GetDefaultValueFloat(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	ConvertWStringToMIMEText(v, L"The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue() == v);
}

void ParameterNameValueTypeTest::DefaultValueAsciiTest()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueAscii("The Los Angeles Angels of Anaheim"));
	CPPUNIT_ASSERT(nvt.GetDefaultValueAscii() == "The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT_THROW(nvt.GetDefaultValueFloat(), affymetrix_calvin_exceptions::ParameterMismatchException);

	MIMEValue v;
	ConvertStringToMIMEText(v, "The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue() == v);
}

void ParameterNameValueTypeTest::DefaultValueTextWithReserveTest()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueText(L"The Los Angeles Angels of Anaheim", 100));
	CPPUNIT_ASSERT(nvt.GetDefaultValueText() == L"The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue().Size() == 100*sizeof(u_int16_t));
}

void ParameterNameValueTypeTest::DefaultValueAsciiWithReserveTest()
{
	ParameterNameValueDefaultRequiredType nvt;
	CPPUNIT_ASSERT_NO_THROW(nvt.SetDefaultValueAscii("The Los Angeles Angels of Anaheim", 100));
	CPPUNIT_ASSERT(nvt.GetDefaultValueAscii() == "The Los Angeles Angels of Anaheim");
	CPPUNIT_ASSERT(nvt.DefaultMIMEValue().Size() == 100);
}

/*
 * Convert the type to a string.
 */
void ParameterNameValueTypeTest::ParameterValueTypeToStringTest()
{
	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeToString(
		ParameterNameValueDefaultRequiredType::NoParameterType) == L"");

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeToString(
		ParameterNameValueDefaultRequiredType::IntegerParameterType) == L"Int");
		
	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeToString(
		ParameterNameValueDefaultRequiredType::FloatParameterType) == L"Float");

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeToString(
		ParameterNameValueDefaultRequiredType::TextParameterType) == L"String");

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeToString(
		ParameterNameValueDefaultRequiredType::DateParameterType) == L"Date");

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeToString(
		ParameterNameValueDefaultRequiredType::TimeParameterType) == L"Time");

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeToString(
		ParameterNameValueDefaultRequiredType::DateTimeParameterType) == L"DateTime");

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeToString(
		ParameterNameValueDefaultRequiredType::ControlSingleParameterType) == L"SingleControl");
		
	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeToString(
		ParameterNameValueDefaultRequiredType::ControlMultiParameterType) == L"MultiControl");
}

/*
 * Convert the string to a type.
 */
void ParameterNameValueTypeTest::ParameterValueTypeFromStringTest()
{
	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeFromString(
		L"Int") == ParameterNameValueDefaultRequiredType::IntegerParameterType);
		
	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeFromString(
		L"Float") == ParameterNameValueDefaultRequiredType::FloatParameterType);

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeFromString(
		L"String") == ParameterNameValueDefaultRequiredType::TextParameterType);

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeFromString(
		L"Date") == ParameterNameValueDefaultRequiredType::DateParameterType);

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeFromString(
		L"Time") == ParameterNameValueDefaultRequiredType::TimeParameterType);

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeFromString(
		L"DateTime") == ParameterNameValueDefaultRequiredType::DateTimeParameterType);

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeFromString(
		L"SingleControl") == ParameterNameValueDefaultRequiredType::ControlSingleParameterType);
		
	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeFromString(
		L"MultiControl") == ParameterNameValueDefaultRequiredType::ControlMultiParameterType);

	CPPUNIT_ASSERT(ParameterNameValueDefaultRequiredType::ParameterValueTypeFromString(
		L"") == ParameterNameValueDefaultRequiredType::NoParameterType);
}

