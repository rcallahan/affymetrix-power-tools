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


#ifndef __DATAGROUPREADERTEST_H_
#define __DATAGROUPREADERTEST_H_

#include "calvin_files/data/src/FileHeader.h"
#include "calvin_files/parsers/src/DataGroupReader.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

using namespace affymetrix_calvin_io;

class DataGroupReaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( DataGroupReaderTest );

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (GetDataGroupNameTest);
	CPPUNIT_TEST (GetDataSetReaderByIndexTest);
	CPPUNIT_TEST (GetDataSetReaderByNameTest);

	CPPUNIT_TEST_SUITE_END();

public:

	void setUp();
	void tearDown();

	void CreationTest();
	void GetDataGroupNameTest();
	void GetDataSetReaderByIndexTest();
	void GetDataSetReaderByNameTest();

private:

	std::ifstream is;
	FileHeader fh;

};

#endif // __DATAGROUPREADERTEST_H_
