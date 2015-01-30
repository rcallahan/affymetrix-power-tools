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
#include "calvin_files/writers/test/CDFFileWriterTest.h"
//
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/writers/src/CDFFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CDFFileWriterTest );

void CDFFileWriterTest::setUp()
{
	CDFData data("CDF_file");
	data.SetProbeSetCnt(10, Expression);
	writer = new CDFFileWriter(data);
}

void CDFFileWriterTest::tearDown()
{
	delete writer;
}

void CDFFileWriterTest::testCreation()
{
	CDFData c("CDF_file");
	c.SetProbeSetCnt(10, Expression);
	CDFFileWriter* w = new CDFFileWriter(c);
	CPPUNIT_ASSERT(1);
	delete w;
}

void CDFFileWriterTest::WriteTest()
{
	writer->OpenDataGroup(L"affy probe set 1", 1);
	CDFProbeSetWriter* probeWriter = writer->CreateProbeSetWriter(L"xda block name",100,100,33,17,3009,100);
	probeWriter->WriteHeader();
	probeWriter->Write(10,10,100,45,11,11);
	probeWriter->Close();
	delete probeWriter;
	writer->CloseDataGroup();

	writer->OpenDataGroup(L"affy probe set 2", 1);
	probeWriter = writer->CreateProbeSetWriter(L"xda block name",100,100,33,17,3009,100);
	probeWriter->WriteHeader();
	probeWriter->Write(10,10,100,45,11,11);
	probeWriter->Close();
	delete probeWriter;
	writer->CloseDataGroup();

	writer->OpenDataGroup(L"affy probe set 3", 1);
	probeWriter = writer->CreateProbeSetWriter(L"xda block name",100,100,33,17,3009,100);
	probeWriter->WriteHeader();
	probeWriter->Write(10,10,100,45,11,11);
	probeWriter->Close();
	delete probeWriter;
	writer->CloseDataGroup();

	writer->OpenDataGroup(L"affy probe set 4", 1);
	probeWriter = writer->CreateProbeSetWriter(L"xda block name",100,100,33,17,3009,100);
	probeWriter->WriteHeader();
	probeWriter->Write(10,10,100,45,11,11);
	probeWriter->Close();
	delete probeWriter;
	writer->CloseDataGroup();

	CPPUNIT_ASSERT(1);
}
