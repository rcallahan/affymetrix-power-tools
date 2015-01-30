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

//
#include "calvin_files/parsers/test/CHPTilingFileReaderTest.h"
//
#include "calvin_files/data/src/CHPTilingEntry.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/data/src/TilingResultData.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/parsers/src/CHPTilingFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_parameter;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPTilingFileReaderTest );

void CHPTilingFileReaderTest::setUp()
{
}

void CHPTilingFileReaderTest::tearDown()
{
}

void CHPTilingFileReaderTest::testCreation()
{
	CHPTilingFileReader reader;
	CPPUNIT_ASSERT(1);
}

void CHPTilingFileReaderTest::testRead()
{
	CHPTilingData data;
	CHPTilingFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename("../data/CHP_tiling_file"));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == "../data/CHP_tiling_file");
	CPPUNIT_ASSERT(data.GetNumberSequences() == 2);
	CPPUNIT_ASSERT(data.GetAlgName() == L"tile");
	CPPUNIT_ASSERT(data.GetAlgVersion() == L"1.0");

	ParameterNameValueTypeList params = data.GetAlgParams();
	CPPUNIT_ASSERT(params.size() == 1);
	ParameterNameValueTypeList::iterator it=params.begin();
	ParameterNameValueType &param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"p1");
	CPPUNIT_ASSERT(param.GetValueText() == L"v1");

	const double eps=0.00001;
	CHPTilingEntry e;
	TilingSequenceData seq;

	data.OpenTilingSequenceDataSet(0);
	seq = data.GetTilingSequenceData();
	
	CPPUNIT_ASSERT(seq.name == L"n1");
	CPPUNIT_ASSERT(seq.groupName == L"g1");
	CPPUNIT_ASSERT(seq.version == L"v1");
	CPPUNIT_ASSERT(seq.parameters.size() == 1);
	it = seq.parameters.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"seq1_p1");
	CPPUNIT_ASSERT(param.GetValueText() == L"seq1_v1");

	CPPUNIT_ASSERT(data.GetTilingSequenceEntryCount(0) == 2);
	data.GetTilingSequenceEntry(0, e);
	CPPUNIT_ASSERT(e.position == 10);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 10.0f, eps);

	data.GetTilingSequenceEntry(1, e);
	CPPUNIT_ASSERT(e.position == 20);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 20.0f, eps);


	data.OpenTilingSequenceDataSet(1);
	seq = data.GetTilingSequenceData();
	
	CPPUNIT_ASSERT(seq.name == L"n2");
	CPPUNIT_ASSERT(seq.groupName == L"g2");
	CPPUNIT_ASSERT(seq.version == L"v2");
	CPPUNIT_ASSERT(seq.parameters.size() == 1);
	it = seq.parameters.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"seq2_p1");
	CPPUNIT_ASSERT(param.GetValueText() == L"seq2_v1");

	CPPUNIT_ASSERT(data.GetTilingSequenceEntryCount(1) == 3);
	data.GetTilingSequenceEntry(0, e);
	CPPUNIT_ASSERT(e.position == 11);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 11.0f, eps);

	data.GetTilingSequenceEntry(1, e);
	CPPUNIT_ASSERT(e.position == 21);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 21.0f, eps);

	data.GetTilingSequenceEntry(2, e);
	CPPUNIT_ASSERT(e.position == 31);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 31.0f, eps);
}
