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

#include "file/CPPTest/BARFileWriterTest.h"
//
#include "file/BARFileWriter.h"
//
#include <cmath>
//

#define TEST_BAR_FILE_NAME "./data/test_new.bar"

CPPUNIT_TEST_SUITE_REGISTRATION( CBARFileWriterTest );

using namespace affxbarwriter;
using namespace affxbar;

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

void CBARFileWriterTest::setUp()
{
}

void CBARFileWriterTest::tearDown()
{
}

void CBARFileWriterTest::testCreation()
{
	CBARFileWriter bar;
	CPPUNIT_ASSERT( 1 );
}

void CBARFileWriterTest::testmethod_CreateNewFile_and_Save()
{
	CBARFileWriter barWrite;

	barWrite.SetFileName(TEST_BAR_FILE_NAME);
	CPPUNIT_ASSERT(barWrite.CreateNewFile() == true);

	TagValuePairType param;
	param.Tag="Param1";
	param.Value="NoVal";
	barWrite.AddAlgorithmParameter(param.Tag.c_str(), param.Value.c_str());

	// Affy BAR files currently expect 1 int for the position and 1 or more floats.
	barWrite.AddColumn(BAR_DATA_INTEGER);
	barWrite.AddColumn(BAR_DATA_FLOAT);

	barWrite.SetNumberSequences(1);
	BarSequenceResultData data;
	for (int k=0; k<barWrite.GetNumberSequences(); k++)
	{
		CGDACSequenceResultItem *p = barWrite.GetResultsPtr(k);
		p->SetNumberDataPoints(5);
		p->SetName("Seq1");
		p->SetGroupName("Group1");
		p->SetVersion("1.0");
		p->AddParameter("SeqParam1", "NoVal");
		for (int l=0; l<p->GetNumberDataPoints(); l++)
		{
			data.iValue = l;
			p->SetDataPoint(l, 0, data);
			data.fValue = (float)(2*l)+0.5f;
			p->SetDataPoint(l, 1, data);
		}
	}

	CPPUNIT_ASSERT(barWrite.Save() == true);
	barWrite.Close();


	CBARFileWriter barRead;
	barRead.SetFileName(TEST_BAR_FILE_NAME);
	CPPUNIT_ASSERT(barRead.Read() == true);

	CPPUNIT_ASSERT( CompareFloats(barRead.GetVersion(), 2.0f) );
	CPPUNIT_ASSERT( barRead.GetNumberSequences() == 1);
	CPPUNIT_ASSERT( barRead.GetNumberColumns() == 2);
	CPPUNIT_ASSERT( barRead.GetNumberParameters() == 1);

	CPPUNIT_ASSERT( barRead.GetParameter(0).Tag == std::string("Param1") );
	CPPUNIT_ASSERT( barRead.GetParameter(0).Value == std::string("NoVal") );
	CPPUNIT_ASSERT( barRead.GetColumnTypes(0) == affxbar::BAR_DATA_INTEGER );
	CPPUNIT_ASSERT( barRead.GetColumnTypes(1) == affxbar::BAR_DATA_FLOAT );

	const double eps = 1e-5;
	CGDACSequenceResultItem seq;
	barRead.GetResults(0, seq);
	CPPUNIT_ASSERT(seq.GetColumnType(0) == 	affxbar::BAR_DATA_INTEGER);
	CPPUNIT_ASSERT(seq.GetColumnType(1) == affxbar::BAR_DATA_FLOAT);
	CPPUNIT_ASSERT(seq.GetNumberColumns() == 2);
	CPPUNIT_ASSERT(seq.GetName() == std::string("Seq1") );
	CPPUNIT_ASSERT(seq.GetVersion() == std::string("1.0") );
	CPPUNIT_ASSERT(seq.GetGroupName() == std::string("Group1") );
	CPPUNIT_ASSERT(seq.GetNumberDataPoints() == 5);
	CPPUNIT_ASSERT(seq.GetNumberParameters() == 1);
	CPPUNIT_ASSERT(seq.GetParameter(0).Tag == std::string("SeqParam1") );
	CPPUNIT_ASSERT(seq.GetParameter(0).Value == std::string("NoVal") );
	seq.GetData(0, 0, data);
	CPPUNIT_ASSERT(data.iValue == 0);
	seq.GetData(0, 1, data);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.fValue, .5, eps);
	seq.GetData(1, 0, data);
	CPPUNIT_ASSERT(data.iValue == 1);
	seq.GetData(1, 1, data);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.fValue, 2.5, eps);
	seq.GetData(2, 0, data);
	CPPUNIT_ASSERT(data.iValue == 2);
	seq.GetData(2, 1, data);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.fValue, 4.5, eps);
	seq.GetData(3, 0, data);
	CPPUNIT_ASSERT(data.iValue == 3);
	seq.GetData(3, 1, data);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.fValue, 6.5, eps);
	seq.GetData(4, 0, data);
	CPPUNIT_ASSERT(data.iValue == 4);
	seq.GetData(4, 1, data);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.fValue, 8.5, eps);
	barRead.Close();
}
