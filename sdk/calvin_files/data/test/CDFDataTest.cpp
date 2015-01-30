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
#include "calvin_files/data/test/CDFDataTest.h"
//
#include "calvin_files/data/src/CDFData.h"
#include "calvin_files/utils/src/AffyStlCollectionTypes.h"
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;

CPPUNIT_TEST_SUITE_REGISTRATION( CDFDataTest );

void CDFDataTest::setUp()
{
	data = new CDFData("CDFMetaData");
}

void CDFDataTest::tearDown()
{
	delete data;
}

void CDFDataTest::testCreation()
{
	CDFData hdr("CDFMetaData");
	CPPUNIT_ASSERT(1);
}

void CDFDataTest::FilenameTest()
{
	std::string p1 = "CDFMetaData";
	data->SetFilename(p1);
	std::string p2 = data->GetFilename();
	CPPUNIT_ASSERT(p1 == p2);
}

void CDFDataTest::DataTypeIdTest()
{
	data->SetProbeSetCnt(10, Expression);
	CPPUNIT_ASSERT(data->GetDataTypeId() == AFFY_EXPR_PS);

	data->SetProbeSetCnt(10, Genotyping);
	CPPUNIT_ASSERT(data->GetDataTypeId() == AFFY_GENO_PS);

	data->SetProbeSetCnt(10, Tag);
	CPPUNIT_ASSERT(data->GetDataTypeId() == AFFY_TAG_PS);

	data->SetProbeSetCnt(10, Resequencing);
	CPPUNIT_ASSERT(data->GetDataTypeId() == AFFY_RESEQ_PS);

	data->SetProbeSetCnt(10, Control);
	CPPUNIT_ASSERT(data->GetDataTypeId() == AFFY_CNTRL_PS);
}

void CDFDataTest::ProbeSetCntTest()
{
	data->SetProbeSetCnt(100, Expression);
	CPPUNIT_ASSERT(data->GetProbeSetCnt() == 100);
	data->SetProbeSetCnt(10, Genotyping);
	CPPUNIT_ASSERT(data->GetProbeSetCnt() == 10);
}

void CDFDataTest::ArrayRowsTest()
{
	data->SetArrayRows(1001);
	CPPUNIT_ASSERT(data->GetArrayRows() == 1001);
}

void CDFDataTest::ArrayColsTest()
{
	data->SetArrayCols(2001);
	CPPUNIT_ASSERT(data->GetArrayCols() == 2001);
}

void CDFDataTest::RefSequenceTest()
{
	std::string seq("this really could be a std::string, too bad GenericDataHeader parameter list does not support it.");
	data->SetRefSequence(seq);
	CPPUNIT_ASSERT(data->GetRefSequence() == seq);
}

void CDFDataTest::FileHeaderTest()
{
	// Put something into the header.
	data->SetArrayCols(45);

	// Get the header.
	FileHeader* fh = data->GetFileHeader();
	CPPUNIT_ASSERT(fh);

	// Try to read the thing from the header.
	ParameterNameValueType pvt;
	CPPUNIT_ASSERT(fh->GetGenericDataHdr()->FindNameValParam(CDF_COLS_PARAM, pvt));
	CPPUNIT_ASSERT(pvt.GetValueUInt32() == 45);
}

void CDFDataTest::GetGenericDataTest()
{
	// Put something into the header.
	data->SetArrayCols(33);

	// Get the generic data.
	GenericData& gd = data->GetGenericData();

	// Try to read the thing from the header.
	ParameterNameValueType pvt;
	CPPUNIT_ASSERT(gd.Header().GetGenericDataHdr()->FindNameValParam(CDF_COLS_PARAM, pvt));
	CPPUNIT_ASSERT(pvt.GetValueUInt32() == 33);
}

