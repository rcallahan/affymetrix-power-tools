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

#ifndef __PARAMETERNAMEVALUETYPETEST_H_
#define __PARAMETERNAMEVALUETYPETEST_H_

#include "calvin_files/parameter/src/ParameterNameValueType.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class ParameterNameValueTypeTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE (ParameterNameValueTypeTest);

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (ConstructorOfInt8TypeTest);
	CPPUNIT_TEST (ConstructorOfUInt8TypeTest);
	CPPUNIT_TEST (ConstructorOfInt16TypeTest);
	CPPUNIT_TEST (ConstructorOfUInt16TypeTest);
	CPPUNIT_TEST (ConstructorOfInt32TypeTest);
	CPPUNIT_TEST (ConstructorOfUInt32TypeTest);
	CPPUNIT_TEST (ConstructorOfFloatTypeTest);
	CPPUNIT_TEST (ConstructorOfTextTypeTest);
	CPPUNIT_TEST (ConstructorOfAsciiTypeTest);
	CPPUNIT_TEST (CopyConstructorTest);
	CPPUNIT_TEST (AssignmentOperatorTest);
	CPPUNIT_TEST (EqualityOperatorTest);
	CPPUNIT_TEST (InequalityOperatorTest);
	CPPUNIT_TEST (LessThanOperatorTest);
	CPPUNIT_TEST (GreaterThanOperatorTest);
	CPPUNIT_TEST (NameTest);
	CPPUNIT_TEST (GetTypeUnknownTest);
	CPPUNIT_TEST (ValueInt8Test);
	CPPUNIT_TEST (ValueUInt8Test);
	CPPUNIT_TEST (ValueInt16Test);
	CPPUNIT_TEST (ValueUInt16Test);
	CPPUNIT_TEST (ValueInt32Test);
	CPPUNIT_TEST (ValueUInt32Test);
	CPPUNIT_TEST (ValueFloatTest);
	CPPUNIT_TEST (ValueTextTest);
	CPPUNIT_TEST (ValueAsciiTest);
	CPPUNIT_TEST (ValueTextWithReserveTest);
	CPPUNIT_TEST (ValueAsciiWithReserveTest);
	CPPUNIT_TEST (MIMETypeTest);
	CPPUNIT_TEST (MIMEValueTest);
	CPPUNIT_TEST (ToStringTest);

	CPPUNIT_TEST (testproperties_ParameterNameValueDefaultRequiredType);
	CPPUNIT_TEST (testconstructors_ParameterNameValueDefaultRequiredType);
	CPPUNIT_TEST (testassignment_ParameterNameValueDefaultRequiredType);
	CPPUNIT_TEST (testToString_ParameterNameValueDefaultRequiredType);
	CPPUNIT_TEST (DefaultValueInt8Test);
	CPPUNIT_TEST (DefaultValueUInt8Test);
	CPPUNIT_TEST (DefaultValueInt16Test);
	CPPUNIT_TEST (DefaultValueUInt16Test);
	CPPUNIT_TEST (DefaultValueInt32Test);
	CPPUNIT_TEST (DefaultValueUInt32Test);
	CPPUNIT_TEST (DefaultValueFloatTest);
	CPPUNIT_TEST (DefaultValueTextTest);
	CPPUNIT_TEST (DefaultValueAsciiTest);
	CPPUNIT_TEST (DefaultValueTextWithReserveTest);
	CPPUNIT_TEST (DefaultValueAsciiWithReserveTest);
	CPPUNIT_TEST (ParameterValueTypeFromStringTest);
	CPPUNIT_TEST (ParameterValueTypeToStringTest);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void CreationTest();
	void ConstructorOfInt8TypeTest();
	void ConstructorOfUInt8TypeTest();
	void ConstructorOfInt16TypeTest();
	void ConstructorOfUInt16TypeTest();
	void ConstructorOfInt32TypeTest();
	void ConstructorOfUInt32TypeTest();
	void ConstructorOfFloatTypeTest();
	void ConstructorOfTextTypeTest();
	void ConstructorOfAsciiTypeTest();
	void CopyConstructorTest();
	void AssignmentOperatorTest();
	void EqualityOperatorTest();
	void InequalityOperatorTest();
	void LessThanOperatorTest();
	void GreaterThanOperatorTest();
	void NameTest();
	void GetTypeUnknownTest();
	void ValueInt8Test();
	void ValueUInt8Test();
	void ValueInt16Test();
	void ValueUInt16Test();
	void ValueInt32Test();
	void ValueUInt32Test();
	void ValueFloatTest();
	void ValueTextTest();
	void ValueAsciiTest();
	void ValueTextWithReserveTest();
	void ValueAsciiWithReserveTest();
	void MIMETypeTest();
	void MIMEValueTest();
	void ToStringTest();
	void testproperties_ParameterNameValueDefaultRequiredType();
	void testconstructors_ParameterNameValueDefaultRequiredType();
	void testassignment_ParameterNameValueDefaultRequiredType();
	void testToString_ParameterNameValueDefaultRequiredType();
	void DefaultValueInt8Test();
	void DefaultValueUInt8Test();
	void DefaultValueInt16Test();
	void DefaultValueUInt16Test();
	void DefaultValueInt32Test();
	void DefaultValueUInt32Test();
	void DefaultValueFloatTest();
	void DefaultValueTextTest();
	void DefaultValueAsciiTest();
	void DefaultValueTextWithReserveTest();
	void DefaultValueAsciiWithReserveTest();
	void ParameterValueTypeToStringTest();
	void ParameterValueTypeFromStringTest();

// support methods
	void ConvertWStringToMIMEText(affymetrix_calvin_parameter::MIMEValue& value, const std::wstring& s);
	void ConvertStringToMIMEText(affymetrix_calvin_parameter::MIMEValue& value, const std::string& s);
	void FillMIMEValue(u_int32_t value, affymetrix_calvin_parameter::MIMEValue& mime);
};

#endif // __PARAMETERNAMEVALUETYPETEST_H_
