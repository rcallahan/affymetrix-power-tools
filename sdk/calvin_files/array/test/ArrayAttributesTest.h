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

class ArrayAttributesTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( ArrayAttributesTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testproperty_Identifier);
	CPPUNIT_TEST ( testproperty_Attributes);
	CPPUNIT_TEST ( testmethod_Clear );
	CPPUNIT_TEST ( testproperty_ArrayBarcode );
	CPPUNIT_TEST ( testproperty_ArrayName );
	CPPUNIT_TEST ( testproperty_Comment );
	CPPUNIT_TEST ( testproperty_Media );
	CPPUNIT_TEST ( testproperty_CreatedBy );
	CPPUNIT_TEST ( testproperty_CreationDateTime );
	CPPUNIT_TEST ( testproperty_MediaFileName );
	CPPUNIT_TEST ( testproperty_MediaFileGUID );
	CPPUNIT_TEST ( testproperty_LibraryPackageName );
	CPPUNIT_TEST ( testproperty_MasterFile );
	CPPUNIT_TEST ( testproperty_MediaCol );
	CPPUNIT_TEST ( testproperty_MediaRow );
	CPPUNIT_TEST ( testproperty_PatAssignment );
	CPPUNIT_TEST ( testproperty_MasterFileId );
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testproperty_Identifier();
	void testproperty_Attributes();
	void testproperty_ArrayBarcode();
	void testproperty_ArrayName();
	void testproperty_Comment();
	void testproperty_Media();
	void testproperty_CreatedBy();
	void testproperty_CreationDateTime();
	void testproperty_MasterFile();
	void testproperty_MediaCol();
	void testproperty_MediaRow();
	void testproperty_PatAssignment();
	void testproperty_MasterFileId();
	void testproperty_MediaFileName();
	void testproperty_MediaFileGUID();
	void testproperty_LibraryPackageName();

	void testmethod_Clear();
};
