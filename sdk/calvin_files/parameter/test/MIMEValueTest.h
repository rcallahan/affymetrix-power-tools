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

#ifndef _MIMEVALUETEST_H_
#define _MIMEVALUETEST_H_

//
#include "calvin_files/parameter/src/ParameterNameValueType.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class MIMEValueTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( MIMEValueTest );

	CPPUNIT_TEST( CreationTest );
	CPPUNIT_TEST( GetValueTest );
	CPPUNIT_TEST( CopyConstructorTest );
	CPPUNIT_TEST( AssignmentOperatorTest );
	CPPUNIT_TEST( SizeTest );
	CPPUNIT_TEST( OperatorEqualsTest );
	CPPUNIT_TEST( OperatorNotEqualsTest );
	CPPUNIT_TEST( SetValueTest );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void CreationTest();
	void GetValueTest();
	void CopyConstructorTest();
	void AssignmentOperatorTest();
	void SizeTest();
	void OperatorEqualsTest();
	void OperatorNotEqualsTest();
	void SetValueTest();

private:		
	char arr[10];
	affymetrix_calvin_parameter::MIMEValue* mimeValue;
};

#endif // _MIMEVALUETEST_H_
