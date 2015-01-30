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


//
#include "file/CPPTest/MSKFileWriterTest.h"
//
#include "file/MSKFileData.h"
#include "file/MSKFileWriter.h"
//
#include <cstring>
#include <string.h>
//

using namespace affxmsk;
using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( CMSKFileWriterTest );

void CMSKFileWriterTest::setUp()
{
}

void CMSKFileWriterTest::tearDown()
{
}

void CMSKFileWriterTest::testmethod_Save_with_call()
{
	const char *file = "./data/new.msk";
	CMSKFileWriter msk;
	msk.Clear();
	msk.SetFileName(file);
	msk.SetArrayType("Hu6800");

	ProbeSetIndiciesType indicies;
	indicies.probeSetName = "A28102_at";
	indicies.indicies.push_back(17);
	indicies.indicies.push_back(18);
	indicies.indicies.push_back(19);
	indicies.indicies.push_back(20);
	msk.AddProbeSetIndicies(indicies);

	indicies.indicies.clear();
	indicies.probeSetName = "AB000381_s_at";
	indicies.indicies.push_back(7);
	msk.AddProbeSetIndicies(indicies);

	CPPUNIT_ASSERT(msk.Save() == true);

	////////////////////////////

	msk.Clear();
	msk.SetFileName(file);
	CPPUNIT_ASSERT(msk.Read() == true);

	ProbeSetIndiciesType ind;
	list<int>::iterator index_iter;
	ProbeSetIndiciesListConstIt iter;
	ProbeSetIndiciesListConstIt begin;
	ProbeSetIndiciesListConstIt end;

	msk.GetProbeSetIndiciesIterators(begin, end);
	iter = begin;

	CPPUNIT_ASSERT( msk.GetProbeSetIndiciesListCount() == 2);

	ind = *iter;
	CPPUNIT_ASSERT( ind.probeSetName == "A28102_at" );
	index_iter = ind.indicies.begin();
	CPPUNIT_ASSERT( *index_iter == 17 );
	++index_iter;
	CPPUNIT_ASSERT( *index_iter == 18 );
	++index_iter;
	CPPUNIT_ASSERT( *index_iter == 19 );
	++index_iter;
	CPPUNIT_ASSERT( *index_iter == 20 );
	++index_iter;
	CPPUNIT_ASSERT( index_iter == ind.indicies.end() );

	++iter;
	ind = *iter;
	CPPUNIT_ASSERT( ind.probeSetName == "AB000381_s_at" );
	index_iter = ind.indicies.begin();
	CPPUNIT_ASSERT( *index_iter == 7 );
	++index_iter;
	CPPUNIT_ASSERT( index_iter == ind.indicies.end() );


	CPPUNIT_ASSERT( msk.GetProbeSetListCount() == 0);

	CPPUNIT_ASSERT( strcmp(msk.GetArrayType(), "Hu6800") == 0 );
}

void CMSKFileWriterTest::testmethod_Save_with_comp()
{
	const char *file = "./data/new.msk";
	CMSKFileWriter msk;
	msk.Clear();
	msk.SetFileName(file);
	msk.SetArrayType("Hu6800");
	msk.AddProbeSet("ABC_at");
	msk.AddProbeSet("XYZ_at");

	CPPUNIT_ASSERT(msk.Save() == true);

	////////////////////////////

	msk.Clear();
	msk.SetFileName(file);
	CPPUNIT_ASSERT(msk.Read() == true);

	CPPUNIT_ASSERT( msk.GetProbeSetIndiciesListCount() == 0);

	CPPUNIT_ASSERT( msk.GetProbeSetListCount() == 2);

	ProbeSetListConstIt siter;
	ProbeSetListConstIt sbegin;
	ProbeSetListConstIt send;

	msk.GetProbeSetIterators(sbegin, send);

	siter = sbegin;
	CPPUNIT_ASSERT( *siter == "ABC_at" );

	++siter;
	CPPUNIT_ASSERT( *siter == "XYZ_at" );

	++siter;
	CPPUNIT_ASSERT( siter == send );

	CPPUNIT_ASSERT( strcmp(msk.GetArrayType(), "Hu6800") == 0 );
}

void CMSKFileWriterTest::testmethod_Save_with_both()
{
	const char *file = "./data/new.msk";
	CMSKFileWriter msk;
	msk.Clear();
	msk.SetFileName(file);
	msk.SetArrayType("Hu6800");
	msk.AddProbeSet("ABC_at");
	msk.AddProbeSet("XYZ_at");

	ProbeSetIndiciesType indicies;
	indicies.probeSetName = "A28102_at";
	indicies.indicies.push_back(17);
	indicies.indicies.push_back(18);
	indicies.indicies.push_back(19);
	indicies.indicies.push_back(20);
	msk.AddProbeSetIndicies(indicies);

	indicies.indicies.clear();
	indicies.probeSetName = "AB000381_s_at";
	indicies.indicies.push_back(7);
	msk.AddProbeSetIndicies(indicies);

	CPPUNIT_ASSERT(msk.Save() == true);

	////////////////////////////

	msk.Clear();
	msk.SetFileName(file);
	CPPUNIT_ASSERT(msk.Read() == true);

	ProbeSetIndiciesType ind;
	list<int>::iterator index_iter;
	ProbeSetIndiciesListConstIt iter;
	ProbeSetIndiciesListConstIt begin;
	ProbeSetIndiciesListConstIt end;

	msk.GetProbeSetIndiciesIterators(begin, end);
	iter = begin;

	CPPUNIT_ASSERT( msk.GetProbeSetIndiciesListCount() == 2);

	ind = *iter;
	CPPUNIT_ASSERT( ind.probeSetName == "A28102_at" );
	index_iter = ind.indicies.begin();
	CPPUNIT_ASSERT( *index_iter == 17 );
	++index_iter;
	CPPUNIT_ASSERT( *index_iter == 18 );
	++index_iter;
	CPPUNIT_ASSERT( *index_iter == 19 );
	++index_iter;
	CPPUNIT_ASSERT( *index_iter == 20 );
	++index_iter;
	CPPUNIT_ASSERT( index_iter == ind.indicies.end() );

	++iter;
	ind = *iter;
	CPPUNIT_ASSERT( ind.probeSetName == "AB000381_s_at" );
	index_iter = ind.indicies.begin();
	CPPUNIT_ASSERT( *index_iter == 7 );
	++index_iter;
	CPPUNIT_ASSERT( index_iter == ind.indicies.end() );




	CPPUNIT_ASSERT( msk.GetProbeSetListCount() == 2);

	ProbeSetListConstIt siter;
	ProbeSetListConstIt sbegin;
	ProbeSetListConstIt send;

	msk.GetProbeSetIterators(sbegin, send);

	siter = sbegin;
	CPPUNIT_ASSERT( *siter == "ABC_at" );

	++siter;
	CPPUNIT_ASSERT( *siter == "XYZ_at" );

	++siter;
	CPPUNIT_ASSERT( siter == send );

	CPPUNIT_ASSERT( strcmp(msk.GetArrayType(), "Hu6800") == 0 );
}

void CMSKFileWriterTest::testmethod_Save_empty()
{
	const char *file = "./data/new.msk";
	CMSKFileWriter msk;
	msk.Clear();
	msk.SetFileName(file);
	msk.SetArrayType("Test3");
	CPPUNIT_ASSERT(msk.Save() == true);

	msk.Clear();
	msk.SetFileName(file);
	CPPUNIT_ASSERT(msk.Read() == true);
	CPPUNIT_ASSERT(strcmp(msk.GetArrayType(), "Test3") == 0);
	CPPUNIT_ASSERT( msk.GetProbeSetListCount() == 0);
	CPPUNIT_ASSERT( msk.GetProbeSetIndiciesListCount() == 0);
}

void CMSKFileWriterTest::testmethod_Save_no_array_type()
{
	const char *file = "./data/new.msk";
	CMSKFileWriter msk;
	msk.Clear();
	msk.SetFileName(file);
	CPPUNIT_ASSERT(msk.Save() == false);
}
