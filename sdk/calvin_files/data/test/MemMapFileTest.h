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

#include "calvin_files/data/src/MemMapFile.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

// MemMapFile has more capabilities than are currently tested.  These
// are basic tests intended at exercising the functionality of this
// class that is currently used.  The MemMapFile class is based
// on the DATIO7 in terms of its use of memory mapping and so MemMapFile
// inherited an richer set of functionality than it needed.

using namespace affymetrix_calvin_io;

class MemMapFileTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE (MemMapFileTest);

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (FilenameTest);
	CPPUNIT_TEST (CreateTest);
	CPPUNIT_TEST (OpenTest);
	CPPUNIT_TEST (MapDataTest);
	CPPUNIT_TEST (CloseTest);
	CPPUNIT_TEST (GetDataPtrTest);
	CPPUNIT_TEST (GetFirstMappedBytePosTest);
	CPPUNIT_TEST (GetBytesMappedTest);
	CPPUNIT_TEST (GetLastErrorTest);
	CPPUNIT_TEST (IsViewMappedTest);
	CPPUNIT_TEST (CreateNoSetFilenameTest);
	CPPUNIT_TEST (OpenNoCreateTest);
	CPPUNIT_TEST (WriteValuesTest);

	CPPUNIT_TEST_SUITE_END();

public:

	void setUp();
	void tearDown();
	void CreationTest();
	void FilenameTest();
	void CreateTest();
	void OpenTest();
	void MapDataTest();
	void CloseTest();
	void GetDataPtrTest();
	void GetFirstMappedBytePosTest();
	void GetBytesMappedTest();
	void GetLastErrorTest();
	void IsViewMappedTest();
	void CreateNoSetFilenameTest();
	void OpenNoCreateTest();
	void WriteValuesTest();

};
