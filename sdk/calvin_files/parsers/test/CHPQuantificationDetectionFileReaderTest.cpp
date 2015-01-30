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
#include "calvin_files/parsers/test/CHPQuantificationDetectionFileReaderTest.h"
//
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/parsers/src/CHPQuantificationDetectionFileReader.h"
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

CPPUNIT_TEST_SUITE_REGISTRATION( CHPQuantificationDetectionFileReaderTest );

void CHPQuantificationDetectionFileReaderTest::setUp()
{
}

void CHPQuantificationDetectionFileReaderTest::tearDown()
{
}

void CHPQuantificationDetectionFileReaderTest::testCreation()
{
	CHPQuantificationDetectionFileReader reader;
	CPPUNIT_ASSERT(1);
}

void CHPQuantificationDetectionFileReaderTest::testRead()
{
	CHPQuantificationDetectionData data;
	CHPQuantificationDetectionFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename("../data/CHP_quantification_detection_file"));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == "../data/CHP_quantification_detection_file");
	CPPUNIT_ASSERT(data.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data.GetEntryCount() == 2);
	
	ParameterNameValueTypeList params = data.GetAlgParams();
	CPPUNIT_ASSERT(params.size() == 1);
	ParameterNameValueTypeList::iterator it=params.begin();
	ParameterNameValueType param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	params = data.GetSummaryParams();
	CPPUNIT_ASSERT(params.size() == 1);
	it=params.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");

	ProbeSetQuantificationDetectionData e;
	data.GetQuantificationDetectionEntry(0, e);
	CPPUNIT_ASSERT(e.name == "abc");
	CPPUNIT_ASSERT(e.quantification == 10.0f);
	CPPUNIT_ASSERT(e.pvalue == 0.1f);
	CPPUNIT_ASSERT(e.id == -1);
	data.GetQuantificationDetectionEntry(1, e);
	CPPUNIT_ASSERT(e.name == "xyz");
	CPPUNIT_ASSERT(e.quantification == 20.0f);
	CPPUNIT_ASSERT(e.pvalue == 0.2f);
	CPPUNIT_ASSERT(e.id == -1);
}

void CHPQuantificationDetectionFileReaderTest::testReadId()
{
	CHPQuantificationDetectionData data;
	CHPQuantificationDetectionFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename("../data/CHP_quantification_detection_file_id"));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == "../data/CHP_quantification_detection_file_id");
	CPPUNIT_ASSERT(data.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data.GetEntryCount() == 2);
	
	ParameterNameValueTypeList params = data.GetAlgParams();
	CPPUNIT_ASSERT(params.size() == 1);
	ParameterNameValueTypeList::iterator it=params.begin();
	ParameterNameValueType param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	params = data.GetSummaryParams();
	CPPUNIT_ASSERT(params.size() == 1);
	it=params.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");

	ProbeSetQuantificationDetectionData e;
	data.GetQuantificationDetectionEntry(0, e);
	CPPUNIT_ASSERT(e.id == 10);
	CPPUNIT_ASSERT(e.quantification == 10.0f);
	CPPUNIT_ASSERT(e.pvalue == 0.1f);
	CPPUNIT_ASSERT(e.name == "");
	data.GetQuantificationDetectionEntry(1, e);
	CPPUNIT_ASSERT(e.id == 20);
	CPPUNIT_ASSERT(e.quantification == 20.0f);
	CPPUNIT_ASSERT(e.pvalue == 0.2f);
	CPPUNIT_ASSERT(e.name == "");
}
