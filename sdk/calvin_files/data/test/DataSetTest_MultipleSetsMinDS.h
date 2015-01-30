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

#ifndef __DATASETTEST_MULTIPLESETSMINDS_H_
#define __DATASETTEST_MULTIPLESETSMINDS_H_

#include "calvin_files/data/src/GenericData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

/*
 * This tests access by multiple planes of the same file concurrently
 * that that have been read with ReadHeaderMinDP method
 */
class DataSetTest_MultipleSetsMinDS : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE(DataSetTest_MultipleSetsMinDS);

	CPPUNIT_TEST (testCreation);
	CPPUNIT_TEST (testCreateSecondDataSet);
	CPPUNIT_TEST (testOpenTwoDataSets);
	CPPUNIT_TEST (testAccessDataFromTwoDifferentDataSets);
	CPPUNIT_TEST (testAccessDataFromTwoInstancesOfSamePlane);

	CPPUNIT_TEST_SUITE_END();

public:
	DataSetTest_MultipleSetsMinDS();
	~DataSetTest_MultipleSetsMinDS();

	void setUp();
	void tearDown();

	void testCreation();
	void testCreateSecondDataSet();
	void testOpenTwoDataSets();
	void testAccessDataFromTwoDifferentDataSets();
	void testAccessDataFromTwoInstancesOfSamePlane();


private:
	affymetrix_calvin_io::GenericData* data;
	affymetrix_calvin_io::DataSet* dataPlane;
};

#endif // __DATASETTEST_MULTIPLESETSMINDS_H_
