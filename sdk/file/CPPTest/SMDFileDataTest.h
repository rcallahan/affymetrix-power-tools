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

#ifndef __SMDFILEDATATEST_H_
#define __SMDFILEDATATEST_H_

#include <cppunit/extensions/HelperMacros.h>

class SMDFileDataTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( SMDFileDataTest );
	
	CPPUNIT_TEST(CreationTest);
	CPPUNIT_TEST(FileNameTest);
	CPPUNIT_TEST(ExistsTest);
	CPPUNIT_TEST(ReadTest);
	CPPUNIT_TEST(ReadFailedTest);
	CPPUNIT_TEST(FrameRowsMemberTest);
	CPPUNIT_TEST(FrameColsMemberTest);
	CPPUNIT_TEST(NumFramesTest);
	CPPUNIT_TEST(GetFrameTest);
	CPPUNIT_TEST(GetFrameOutOfBoundsTest);
	CPPUNIT_TEST(FailedReadClearsTest);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void CreationTest();
	void FileNameTest();
	void ExistsTest();
	void ReadTest();
	void ReadFailedTest();
	void FrameRowsMemberTest();
	void FrameColsMemberTest();
	void NumFramesTest();
	void GetFrameTest();
	void GetFrameOutOfBoundsTest();
	void FailedReadClearsTest();

};
#endif // __SMDFILEDATATEST_H_
