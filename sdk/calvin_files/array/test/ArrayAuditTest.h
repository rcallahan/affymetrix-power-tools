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
#pragma once

#include <cppunit/extensions/HelperMacros.h>

class ArrayAuditEntryTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( ArrayAuditEntryTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testproperty_UserName );
	CPPUNIT_TEST ( testproperty_DateTime );
	CPPUNIT_TEST ( testproperty_ArrayGuid );
	CPPUNIT_TEST ( testproperty_InputFileGuids ); 	 
	CPPUNIT_TEST ( testproperty_OutputFileGuids );
	CPPUNIT_TEST ( testproperty_ActionParameters );
	CPPUNIT_TEST ( testproperty_ActionType );
	CPPUNIT_TEST ( testmethod_Clear );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testproperty_UserName();
	void testproperty_DateTime();
	void testproperty_ArrayGuid();
	void testproperty_InputFileGuids(); 	 
	void testproperty_OutputFileGuids();
	void testproperty_ActionParameters();
	void testproperty_ActionType();
	void testmethod_Clear();
};
