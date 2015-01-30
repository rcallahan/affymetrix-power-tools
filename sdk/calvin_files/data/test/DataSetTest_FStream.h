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

#ifndef __DATASETTEST_FSTREAM_H_
#define __DATASETTEST_FSTREAM_H_

#include "calvin_files/data/src/GenericData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

/*
 * This tests the basic functions of the DataSet class using fstream to read data.
 * The fixture creates DataSets that access data using a small memory footprint.
 * Other tests exercise features to handle different
 * data types, column arrangements, memory footprints
 * and multiple DataSets referencing the same file.
 */
class DataSetTest_FStream : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE(DataSetTest_FStream);

	CPPUNIT_TEST (testCreation);
	CPPUNIT_TEST (testmethod_Delete);
	CPPUNIT_TEST (testmethod_Open);
	CPPUNIT_TEST (testmethod_Close);
	CPPUNIT_TEST (testproperty_BytesPerRow);
	CPPUNIT_TEST (testproperty_Rows);
	CPPUNIT_TEST (testproperty_Cols);
	CPPUNIT_TEST (testproperty_Header);
	CPPUNIT_TEST (testmethod_ReadUInt16_Single);
	CPPUNIT_TEST (testmethod_ReadUInt16_Vector);
	CPPUNIT_TEST (testmethod_ReadUInt16_Vector_TooManyRows);
	CPPUNIT_TEST (testmethod_ReadUInt16_Vector_InChunks);
	CPPUNIT_TEST (testmethod_ReadUInt16_Raw);
	CPPUNIT_TEST (testmethod_ReadUInt16_Raw_TooManyRows);
	CPPUNIT_TEST (testmethod_ReadUInt16_Raw_InChunks);
	CPPUNIT_TEST (testmethod_CheckRowColumnAndType_DataSetNotOpenException);
	CPPUNIT_TEST (testmethod_CheckRowColumnAndType_ColumnOutOfBoundsException);
	CPPUNIT_TEST (testmethod_CheckRowColumnAndType_RowOutofBoundsException);
	CPPUNIT_TEST (testmethod_CheckRowColumnAndType_UnexpectedColumnTypeException);
	CPPUNIT_TEST (testmethod_ReadUInt16_Single_DataSetNotOpenException);
	CPPUNIT_TEST (testmethod_ReadUInt16_Vector_DataSetNotOpenException);
	CPPUNIT_TEST (testmethod_ReadUInt16_Raw_DataSetNotOpenException);

	CPPUNIT_TEST_SUITE_END();

public:
	DataSetTest_FStream();
	~DataSetTest_FStream();

	void setUp();
	void tearDown();

	void testCreation();
	void testmethod_Delete();
	void testmethod_Open();
	void testmethod_Close();
	void testproperty_BytesPerRow();
	void testproperty_Rows();
	void testproperty_Cols();
	void testproperty_Header();
	void testmethod_ReadUInt16_Single();
	void testmethod_ReadUInt16_Vector();
	void testmethod_ReadUInt16_Vector_TooManyRows();
	void testmethod_ReadUInt16_Vector_InChunks();
	void testmethod_ReadUInt16_Raw();
	void testmethod_ReadUInt16_Raw_TooManyRows();
	void testmethod_ReadUInt16_Raw_InChunks();
	void testmethod_CheckRowColumnAndType_DataSetNotOpenException();
	void testmethod_CheckRowColumnAndType_ColumnOutOfBoundsException();
	void testmethod_CheckRowColumnAndType_RowOutofBoundsException();
	void testmethod_CheckRowColumnAndType_UnexpectedColumnTypeException();
	void testmethod_ReadUInt16_Single_DataSetNotOpenException();
	void testmethod_ReadUInt16_Vector_DataSetNotOpenException();
	void testmethod_ReadUInt16_Raw_DataSetNotOpenException();

private:
	affymetrix_calvin_io::GenericData* data;
	affymetrix_calvin_io::DataSet* dataPlane;

	std::ifstream ifs;
};

#endif // __DATASETTEST_FSTREAM_H_
