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

#ifndef __GENERICDATATEST_FILEINDEPENDENT_H_
#define __GENERICDATATEST_FILEINDEPENDENT_H_

#include "calvin_files/data/src/GenericData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

/*
 * These are the GenericData tests that do not rely on reading a file.
 * The setUp method prepares the GenericData object for testing.
 * File-based GenericData tests are elsewhere.
 */
class GenericDataTest_FileIndependent : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( GenericDataTest_FileIndependent );

	CPPUNIT_TEST (testCreation);
	CPPUNIT_TEST (testproperty_Header);
	CPPUNIT_TEST (testmethod_DataSet_ByIndex);
	CPPUNIT_TEST (testmethod_DataSet_ByName);
	CPPUNIT_TEST (testmethod_FindDataSetHeader_ByName);
	CPPUNIT_TEST (testmethod_FindDataSetHeader_ByIndex);
	CPPUNIT_TEST (testmethod_DataGroupNames);
	CPPUNIT_TEST (testproperty_DataGroupCnt);
	CPPUNIT_TEST (testproperty_ArrayFileIdentifier);
	CPPUNIT_TEST (testproperty_FileIdentifier);
	CPPUNIT_TEST (testproperty_ArrayIdentifier);
	CPPUNIT_TEST (testmethod_Clear);
	CPPUNIT_TEST (testmethod_DataSetCnt);
	CPPUNIT_TEST (testmethod_DataSetNames);
	CPPUNIT_TEST (testmethod_FindDataGroupHeader_ByName);
	CPPUNIT_TEST (testmethod_FindDataGroupHeader_ByIndex);

	CPPUNIT_TEST_SUITE_END();

public:
	GenericDataTest_FileIndependent();
	~GenericDataTest_FileIndependent();

	void setUp();
	void tearDown();

	void testCreation();
	void testproperty_Header();
	void testmethod_DataSet_ByIndex();
	void testmethod_DataSet_ByName();
	void testmethod_FindDataSetHeader_ByName();
	void testmethod_FindDataSetHeader_ByIndex();
	void testmethod_DataGroupNames();
	void testproperty_DataGroupCnt();
	void testproperty_ArrayFileIdentifier();
	void testproperty_FileIdentifier();
	void testproperty_ArrayIdentifier();
	void testmethod_Clear();
	void testmethod_DataSetCnt();
	void testmethod_DataSetNames();
	void testmethod_FindDataGroupHeader_ByName();
	void testmethod_FindDataGroupHeader_ByIndex();

private:
	affymetrix_calvin_io::GenericData* data;
	affymetrix_calvin_io::GenericDataHeader* header;
	affymetrix_calvin_io::GenericDataHeader* parent;
	affymetrix_calvin_io::DataGroupHeader* dch;
	affymetrix_calvin_io::DataSetHeader* dphPI;
	affymetrix_calvin_io::DataSetHeader* dphGrid;

};

// tests to write:
// IsDPHPartiallyRead
// ReadFullDataSetHeader
// OpenFStream

#endif // __GENERICDATATEST_FILEINDEPENDENT_H_
