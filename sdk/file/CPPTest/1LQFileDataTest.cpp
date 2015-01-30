////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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


#include "file/CPPTest/1LQFileDataTest.h"
//
#include "file/1LQFileData.h"
//

CPPUNIT_TEST_SUITE_REGISTRATION( C1LQFileDataTest );

#define TEST_FILE "./data/test3x2.1lq"
#define TEST_IGNORE_LINES_FILE "./data/test3x2-ignorelines.1lq"

using namespace affx1lq;

void C1LQFileDataTest::setUp()
{
}

void C1LQFileDataTest::tearDown()
{
}

void C1LQFileDataTest::testCreation()
{
	C1LQFileData data;
	CPPUNIT_ASSERT( 1 );
}

void C1LQFileDataTest::testproperty_FileName()
{
	C1LQFileData data;
	std::string path = TEST_FILE;
	data.SetFileName(path.c_str());
	CPPUNIT_ASSERT( path == data.GetFileName() );
}

void C1LQFileDataTest::testmethod_Exists()
{
	C1LQFileData data;
	std::string path = TEST_FILE;
	data.SetFileName(path.c_str());
	CPPUNIT_ASSERT( data.Exists() == true );
}

void C1LQFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	C1LQFileData data;
	std::string path = "no_file";
	data.SetFileName(path.c_str());
	CPPUNIT_ASSERT( data.Exists() == false );
}

void C1LQFileDataTest::testmethod_Read_ignoring_comment_lines()
{
	C1LQFileData data;
	std::string path = TEST_IGNORE_LINES_FILE;
	data.SetFileName(path.c_str());
	CPPUNIT_ASSERT( data.Read() );

	CPPUNIT_ASSERT( data.GetNumberRows() == 1 );
	CPPUNIT_ASSERT( data.GetNumberColumns() == 2 );
	CPPUNIT_ASSERT( data.GetEntries().size() == 2);

	const std::list<DataEntry> &entries = data.GetEntries();
	std::list<DataEntry>::const_iterator it = entries.begin();
	DataEntry entry;

	entry = *it;

	CPPUNIT_ASSERT( entry.x == 0 );
	CPPUNIT_ASSERT( entry.y == 0 );
	CPPUNIT_ASSERT( entry.probe == "AAA" );
	CPPUNIT_ASSERT( entry.destype == 3 );
	CPPUNIT_ASSERT( entry.feature == "c1" );
	CPPUNIT_ASSERT( entry.qualifier == "q1" );
	CPPUNIT_ASSERT( entry.expos == 11 );
	CPPUNIT_ASSERT( entry.plength == 15 );
	CPPUNIT_ASSERT( entry.position == 9 );
	CPPUNIT_ASSERT( entry.cbase == 'a' );
	CPPUNIT_ASSERT( entry.pbase == 'b' );
	CPPUNIT_ASSERT( entry.tbase == 'c' );
	CPPUNIT_ASSERT( entry.ipbase == 'd' );
	CPPUNIT_ASSERT( entry.unit == 0 );
	CPPUNIT_ASSERT( entry.block == 1 );
	CPPUNIT_ASSERT( entry.atom == 2 );

	++it;
	entry = *it;

	CPPUNIT_ASSERT( entry.x == 1 );
	CPPUNIT_ASSERT( entry.y == 0 );
	CPPUNIT_ASSERT( entry.probe == "GGG" );
	CPPUNIT_ASSERT( entry.destype == -3 );
	CPPUNIT_ASSERT( entry.feature == "c2" );
	CPPUNIT_ASSERT( entry.qualifier == "q2" );
	CPPUNIT_ASSERT( entry.expos == 22 );
	CPPUNIT_ASSERT( entry.plength == 25 );
	CPPUNIT_ASSERT( entry.position == 8 );
	CPPUNIT_ASSERT( entry.cbase == 'e' );
	CPPUNIT_ASSERT( entry.pbase == 'f' );
	CPPUNIT_ASSERT( entry.tbase == 'g' );
	CPPUNIT_ASSERT( entry.ipbase == 'h' );
	CPPUNIT_ASSERT( entry.unit == 1 );
	CPPUNIT_ASSERT( entry.block == 2 );
	CPPUNIT_ASSERT( entry.atom == 3 );

	++it;

	CPPUNIT_ASSERT( it == entries.end() );
}

void C1LQFileDataTest::testmethod_Read()
{
	C1LQFileData data;
	std::string path = TEST_FILE;
	data.SetFileName(path.c_str());
	CPPUNIT_ASSERT( data.Read() );

	CPPUNIT_ASSERT( data.GetNumberRows() == 2 );
	CPPUNIT_ASSERT( data.GetNumberColumns() == 3 );
	CPPUNIT_ASSERT( data.GetEntries().size() == 6);

	const std::list<DataEntry> &entries = data.GetEntries();
	std::list<DataEntry>::const_iterator it = entries.begin();
	DataEntry entry;

	entry = *it;

	CPPUNIT_ASSERT( entry.x == 0 );
	CPPUNIT_ASSERT( entry.y == 0 );
	CPPUNIT_ASSERT( entry.probe == "AAA" );
	CPPUNIT_ASSERT( entry.destype == 3 );
	CPPUNIT_ASSERT( entry.feature == "c1" );
	CPPUNIT_ASSERT( entry.qualifier == "q1" );
	CPPUNIT_ASSERT( entry.expos == 11 );
	CPPUNIT_ASSERT( entry.plength == 15 );
	CPPUNIT_ASSERT( entry.position == 9 );
	CPPUNIT_ASSERT( entry.cbase == 'a' );
	CPPUNIT_ASSERT( entry.pbase == 'b' );
	CPPUNIT_ASSERT( entry.tbase == 'c' );
	CPPUNIT_ASSERT( entry.ipbase == 'd' );
	CPPUNIT_ASSERT( entry.unit == 0 );
	CPPUNIT_ASSERT( entry.block == 1 );
	CPPUNIT_ASSERT( entry.atom == 2 );

	++it;
	entry = *it;

	CPPUNIT_ASSERT( entry.x == 1 );
	CPPUNIT_ASSERT( entry.y == 0 );
	CPPUNIT_ASSERT( entry.probe == "GGG" );
	CPPUNIT_ASSERT( entry.destype == -3 );
	CPPUNIT_ASSERT( entry.feature == "c2" );
	CPPUNIT_ASSERT( entry.qualifier == "q2" );
	CPPUNIT_ASSERT( entry.expos == 22 );
	CPPUNIT_ASSERT( entry.plength == 25 );
	CPPUNIT_ASSERT( entry.position == 8 );
	CPPUNIT_ASSERT( entry.cbase == 'e' );
	CPPUNIT_ASSERT( entry.pbase == 'f' );
	CPPUNIT_ASSERT( entry.tbase == 'g' );
	CPPUNIT_ASSERT( entry.ipbase == 'h' );
	CPPUNIT_ASSERT( entry.unit == 1 );
	CPPUNIT_ASSERT( entry.block == 2 );
	CPPUNIT_ASSERT( entry.atom == 3 );

	++it;
	++it;
	++it;
	++it;
	++it;

	CPPUNIT_ASSERT( it == entries.end() );

	data.Rotate();

	// After rotate, 
	//    y' becomes MaxX - x - 1
	//    x' becomes y
	// MaxX becomes MaxY and vice versa.

	CPPUNIT_ASSERT( data.GetNumberRows() == 3 );
	CPPUNIT_ASSERT( data.GetNumberColumns() == 2 );

	// start over
	it = entries.begin();

	entry = *it;

	CPPUNIT_ASSERT( entry.x == 0 );
	CPPUNIT_ASSERT( entry.y == 2 );

	++it;
	entry = *it;

	CPPUNIT_ASSERT( entry.x == 0 );
	CPPUNIT_ASSERT( entry.y == 1 );

	++it;
	entry = *it;

	CPPUNIT_ASSERT( entry.x == 0 );
	CPPUNIT_ASSERT( entry.y == 0 );

	++it;	
	entry = *it;

	CPPUNIT_ASSERT( entry.x == 1 );
	CPPUNIT_ASSERT( entry.y == 2 );


	++it;
	entry = *it;

	CPPUNIT_ASSERT( entry.x == 1 );
	CPPUNIT_ASSERT( entry.y == 1 );

	++it;
	entry = *it;

	CPPUNIT_ASSERT( entry.x == 1 );
	CPPUNIT_ASSERT( entry.y == 0 );

	++it;

	CPPUNIT_ASSERT( it == entries.end() );

}

void C1LQFileDataTest::testmethod_Clear()
{
	C1LQFileData data;
	std::string path = TEST_FILE;
	data.SetFileName(path.c_str());
	data.Read();
	data.Clear();
	CPPUNIT_ASSERT( data.GetNumberRows() == 0 );
	CPPUNIT_ASSERT( data.GetNumberColumns() == 0 );
	CPPUNIT_ASSERT( data.GetEntries().size() == 0);
}
