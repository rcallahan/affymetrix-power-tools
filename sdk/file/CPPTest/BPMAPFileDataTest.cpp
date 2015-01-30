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

#include "file/CPPTest/BPMAPFileDataTest.h"
//
#include "file/BPMAPFileData.h"
//
#include <cmath>
#include <iostream>
//

CPPUNIT_TEST_SUITE_REGISTRATION( CBPMAPFileDataTest );

#ifdef _MSC_VER
#define TEST_BPMAP_V1_FILE_NAME "TilingArrays\\size_test.sary-full.2004-09-27-v1.bpmap"
#define TEST_BPMAP_V2_FILE_NAME "TilingArrays\\size_test.sary-full.2004-09-27-v2.bpmap"
#define TEST_BPMAP_V3_FILE_NAME "TilingArrays\\test.v3.bpmap"
#define TEST_BPMAP_INCORRECT_HEADER_FILE_NAME "TilingArrays\\incorrect_header.bpmap"
#else
#define TEST_BPMAP_V1_FILE_NAME "TilingArrays/size_test.sary-full.2004-09-27-v1.bpmap"
#define TEST_BPMAP_V2_FILE_NAME "TilingArrays/size_test.sary-full.2004-09-27-v2.bpmap"
#define TEST_BPMAP_V2_FILE_NAME "TilingArrays/size_test.sary-full.2004-09-27-v2.bpmap"
#define TEST_BPMAP_V3_FILE_NAME "TilingArrays/test.v3.bpmap"
#define TEST_BPMAP_INCORRECT_HEADER_FILE_NAME "TilingArrays/incorrect_header.bpmap"
#endif

#define PM_MM 0
#define PM_ONLY 1

static std::string GetDataPath(const char *relPath)
{
	extern std::string externalDataPath; // The path to the data.
//	std::string path = externalDataPath;
//#ifdef _MSC_VER
//	path += "\\";
//#else
//	path += "/";
//#endif
//	path += relPath;
//	return path;
  return Fs::join(externalDataPath,relPath);
}

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

void CBPMAPFileDataTest::setUp()
{
}

void CBPMAPFileDataTest::tearDown()
{
}

void CBPMAPFileDataTest::testCreation()
{
	affxbpmap::CBPMAPFileData bpmap;
	CPPUNIT_ASSERT( 1 );
}

void CBPMAPFileDataTest::testproperty_FileName()
{
	affxbpmap::CBPMAPFileData bpmap;
	std::string path = "test";
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT( path == bpmap.GetFileName() );
}

void CBPMAPFileDataTest::testmethod_Exists()
{
	affxbpmap::CBPMAPFileData bpmap;
	std::string path = GetDataPath(TEST_BPMAP_V1_FILE_NAME);
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bpmap.Exists() == true );
}

void CBPMAPFileDataTest::testmethod_ExistsWhenFileNotExists()
{
	affxbpmap::CBPMAPFileData bpmap;
	std::string path = "test";
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bpmap.Exists() == false );
}

void CBPMAPFileDataTest::testmethod_Read()
{
	affxbpmap::CBPMAPFileData bpmap;
	affxbpmap::CGDACSequenceItem seq;
	affxbpmap::GDACSequenceHitItemType hit;

	std::string path = GetDataPath(TEST_BPMAP_V1_FILE_NAME);
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", bpmap.Exists() == true );
	CPPUNIT_ASSERT( bpmap.Read() );
	CPPUNIT_ASSERT( bpmap.GetNumberSequences() == 68);
	CPPUNIT_ASSERT( CompareFloats(bpmap.GetVersion(), 1.0f) );
	bpmap.GetSequenceItem(0, seq);
	CPPUNIT_ASSERT( seq.GetName() == std::string("at:Jun_2003;chr1") );
	CPPUNIT_ASSERT( seq.GroupName() == std::string("") );
	CPPUNIT_ASSERT( seq.GetSeqVersion() == std::string("") );
	CPPUNIT_ASSERT( seq.GetNumber() == 0 );
	CPPUNIT_ASSERT( seq.GetNumberHits() == 5204 );
	CPPUNIT_ASSERT( seq.GetNumberParameters() == 0 );
	seq.GetHitItem(0, hit, true);
	CPPUNIT_ASSERT( hit.PMX == 394 );
	CPPUNIT_ASSERT( hit.PMY == 637 );
	CPPUNIT_ASSERT( hit.MMX == 394 );
	CPPUNIT_ASSERT( hit.MMY == 638 );
	CPPUNIT_ASSERT( CompareFloats(hit.MatchScore, 1.0f) );
	CPPUNIT_ASSERT( hit.PMProbe == std::string("TATAATATTTCACCGCCACGAATTA") );
	CPPUNIT_ASSERT( hit.Position == 3754875 );
	CPPUNIT_ASSERT( hit.getStartPosition() == 3754875 );
	CPPUNIT_ASSERT( hit.getCenterPosition() == 3754875 + (hit.ProbeLength-1)/2 );
	CPPUNIT_ASSERT( hit.ProbeLength == 25 );
	CPPUNIT_ASSERT( hit.TopStrand == 0 );

	seq.GetHitItem(seq.GetNumberHits()-1, hit, true);
	CPPUNIT_ASSERT( hit.PMX == 488 );
	CPPUNIT_ASSERT( hit.PMY == 43 );
	CPPUNIT_ASSERT( hit.MMX == 488 );
	CPPUNIT_ASSERT( hit.MMY == 44 );
	CPPUNIT_ASSERT( CompareFloats(hit.MatchScore, 1.0f) );
	CPPUNIT_ASSERT( hit.PMProbe == std::string("ATGATCTGAGGTAATGGCGAAAACG") );
	CPPUNIT_ASSERT( hit.Position == 21015718 );
	CPPUNIT_ASSERT( hit.getStartPosition() == 21015718 );
	CPPUNIT_ASSERT( hit.getCenterPosition() == 21015718 + (hit.ProbeLength-1)/2 );
	CPPUNIT_ASSERT( hit.ProbeLength == 25 );
	CPPUNIT_ASSERT( hit.TopStrand == 1 );
	bpmap.Close();
	
	path = GetDataPath(TEST_BPMAP_V2_FILE_NAME);
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", bpmap.Exists() == true );
	CPPUNIT_ASSERT( bpmap.Read() );
	CPPUNIT_ASSERT( bpmap.GetNumberSequences() == 68);
	CPPUNIT_ASSERT( CompareFloats(bpmap.GetVersion(), 2.0f) );
	bpmap.GetSequenceItem(0, seq);
	CPPUNIT_ASSERT( seq.GetName() == std::string("chr1") );
	CPPUNIT_ASSERT( seq.GroupName() == std::string("At") );
	CPPUNIT_ASSERT( seq.GetSeqVersion() == std::string("Jun_2003") );
	CPPUNIT_ASSERT( seq.GetNumber() == 0 );
	CPPUNIT_ASSERT( seq.GetNumberHits() == 5204 );
	CPPUNIT_ASSERT( seq.GetNumberParameters() == 0 );
	seq.GetHitItem(0, hit, true);
	CPPUNIT_ASSERT( hit.PMX == 394 );
	CPPUNIT_ASSERT( hit.PMY == 637 );
	CPPUNIT_ASSERT( hit.MMX == 394 );
	CPPUNIT_ASSERT( hit.MMY == 638 );
	CPPUNIT_ASSERT( CompareFloats(hit.MatchScore, 1.0f) );
	CPPUNIT_ASSERT( hit.PMProbe == std::string("TATAATATTTCACCGCCACGAATTA") );
	CPPUNIT_ASSERT( hit.Position == 3754875 );
	CPPUNIT_ASSERT( hit.getStartPosition() == 3754875 );
	CPPUNIT_ASSERT( hit.getCenterPosition() == 3754875 + (hit.ProbeLength-1)/2 );
	CPPUNIT_ASSERT( hit.ProbeLength == 25 );
	CPPUNIT_ASSERT( hit.TopStrand == 0 );

	seq.GetHitItem(seq.GetNumberHits()-1, hit, true);
	CPPUNIT_ASSERT( hit.PMX == 488 );
	CPPUNIT_ASSERT( hit.PMY == 43 );
	CPPUNIT_ASSERT( hit.MMX == 488 );
	CPPUNIT_ASSERT( hit.MMY == 44 );
	CPPUNIT_ASSERT( CompareFloats(hit.MatchScore, 1.0f) );
	CPPUNIT_ASSERT( hit.PMProbe == std::string("ATGATCTGAGGTAATGGCGAAAACG") );
	CPPUNIT_ASSERT( hit.Position == 21015718 );
	CPPUNIT_ASSERT( hit.getStartPosition() == 21015718 );
	CPPUNIT_ASSERT( hit.getCenterPosition() == 21015718 + (hit.ProbeLength-1)/2 );
	CPPUNIT_ASSERT( hit.ProbeLength == 25 );
	CPPUNIT_ASSERT( hit.TopStrand == 1 );
	bpmap.Close();
}

void CBPMAPFileDataTest::testmethod_ReadV3File()
{
	affxbpmap::CBPMAPFileData bpmap;
	affxbpmap::CGDACSequenceItem seq;
	affxbpmap::GDACSequenceHitItemType hit;
	TagValuePairType tv;

	// note that TEST_BPMAP_V3_FILE_NAME is rewritten in testing
	// BPMAPFileWriterTest.cpp and this occurs in the test suite
	// after this test, a failure in the BPMAPWriter test that 
	// breaks this test won't be seen until the next time the
	// test suite is run. The test suite should be fixed to run the
	// BPMAPWriter test first... not sure how to do this.

	std::string path = GetDataPath(TEST_BPMAP_V3_FILE_NAME);
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", bpmap.Exists() == true );
	CPPUNIT_ASSERT( bpmap.Read() );
	CPPUNIT_ASSERT( bpmap.GetNumberSequences() == 4);
	CPPUNIT_ASSERT( CompareFloats(bpmap.GetVersion(), 3.0f) );

	bpmap.GetSequenceItem(0, seq);
	CPPUNIT_ASSERT( seq.GetName() == std::string("chr1") );
	CPPUNIT_ASSERT( seq.GroupName() == std::string("B_subtilis") );
	CPPUNIT_ASSERT( seq.GetSeqVersion() == std::string("1/1/00:11am") );
	CPPUNIT_ASSERT( seq.GetNumber() == 0 );
	CPPUNIT_ASSERT( seq.GetNumberHits() == 4 );
	CPPUNIT_ASSERT( seq.GetNumberParameters()== 1 );
	CPPUNIT_ASSERT( seq.GetProbeMapping() == PM_MM );

	tv = seq.GetParameter(0);
	CPPUNIT_ASSERT( tv.Tag == "tag1" );
	CPPUNIT_ASSERT( tv.Value == "controls");

	seq.GetHitItem(0, hit, true);
	CPPUNIT_ASSERT( hit.PMX == 1069 );
	CPPUNIT_ASSERT( hit.PMY == 1943 );
	CPPUNIT_ASSERT( hit.MMX == 3069 );
	CPPUNIT_ASSERT( hit.MMY == 1944 );
	CPPUNIT_ASSERT( CompareFloats(hit.MatchScore, 0.8f) );
	CPPUNIT_ASSERT( hit.PMProbe == std::string("CACACCCTAACACTACCCTAACACT") );
	CPPUNIT_ASSERT( hit.Position == 356 );
	CPPUNIT_ASSERT( hit.getStartPosition() == 356 );
	CPPUNIT_ASSERT( hit.getCenterPosition() == 356 + (hit.ProbeLength-1)/2 );
	CPPUNIT_ASSERT( hit.ProbeLength == 25 );
	CPPUNIT_ASSERT( hit.TopStrand == 0 );

	seq.GetHitItem(seq.GetNumberHits()-1, hit, true);
	CPPUNIT_ASSERT( hit.PMX == 1996 );
	CPPUNIT_ASSERT( hit.PMY == 2209 );
	CPPUNIT_ASSERT( hit.MMX == 1996 );
	CPPUNIT_ASSERT( hit.MMY == 2210 );
	CPPUNIT_ASSERT( CompareFloats(hit.MatchScore, 1.0f) );
	CPPUNIT_ASSERT( hit.PMProbe == std::string("AAAAATAGTGTTAGGGTAGTGTTAG") );
	CPPUNIT_ASSERT( hit.Position == 2012000111 );
	CPPUNIT_ASSERT( hit.getStartPosition() == 2012000111 );
	CPPUNIT_ASSERT( hit.getCenterPosition() == 2012000111 + (hit.ProbeLength-1)/2 );
	CPPUNIT_ASSERT( hit.ProbeLength == 25 );
	CPPUNIT_ASSERT( hit.TopStrand == 1 );

	bpmap.GetSequenceItem(1, seq);
	CPPUNIT_ASSERT( seq.GetName() == std::string("chr1") );
	CPPUNIT_ASSERT( seq.GroupName() == std::string("HS1") );
	CPPUNIT_ASSERT( seq.GetSeqVersion() == std::string("11_Nov_2005") );
	CPPUNIT_ASSERT( seq.GetNumber() == 1 );
	CPPUNIT_ASSERT( seq.GetNumberHits() == 8 );
	CPPUNIT_ASSERT( seq.GetNumberParameters() == 2 );
	CPPUNIT_ASSERT( seq.GetProbeMapping() == PM_MM );

	tv = seq.GetParameter(0);
	CPPUNIT_ASSERT( tv.Tag == "tag1" );
	CPPUNIT_ASSERT( tv.Value == "value1");
	tv = seq.GetParameter(1);
	CPPUNIT_ASSERT( tv.Tag == "tag2" );
	CPPUNIT_ASSERT( tv.Value == "value2");

	seq.GetHitItem(0, hit, true);
	CPPUNIT_ASSERT( hit.PMX == 1745 );
	CPPUNIT_ASSERT( hit.PMY == 1095 );
	CPPUNIT_ASSERT( hit.MMX == 1745 );
	CPPUNIT_ASSERT( hit.MMY == 1096 );
	CPPUNIT_ASSERT( CompareFloats(hit.MatchScore, 1.0f) );
	CPPUNIT_ASSERT( hit.PMProbe == std::string("TAGGGCTGTGTTAGGGTAGTGTTAG") );
	CPPUNIT_ASSERT( hit.Position == 64 );
	CPPUNIT_ASSERT( hit.ProbeLength == 25 );
	CPPUNIT_ASSERT( hit.TopStrand == 1 );

	seq.GetHitItem(4, hit, true);
	CPPUNIT_ASSERT( hit.PMX == 354 );
	CPPUNIT_ASSERT( hit.PMY == 1437 );
	CPPUNIT_ASSERT( hit.MMX == 354 );
	CPPUNIT_ASSERT( hit.MMY == 1438 );
	CPPUNIT_ASSERT( CompareFloats(hit.MatchScore, 1.0f) );
	CPPUNIT_ASSERT( hit.PMProbe == std::string("GTCTCTCAACTTACCCTCCATTACC") );
	CPPUNIT_ASSERT( hit.Position == 108 );
	CPPUNIT_ASSERT( hit.ProbeLength == 25 );
	CPPUNIT_ASSERT( hit.TopStrand == 0 );

	bpmap.GetSequenceItem(2, seq);
	CPPUNIT_ASSERT( seq.GetName() == std::string("chr2") );
	CPPUNIT_ASSERT( seq.GroupName() == std::string("HS1") );
	CPPUNIT_ASSERT( seq.GetSeqVersion() == std::string("11_Nov_2005") );
	CPPUNIT_ASSERT( seq.GetNumber() == 2 );
	CPPUNIT_ASSERT( seq.GetNumberHits() == 5 );
	CPPUNIT_ASSERT( seq.GetNumberParameters() == 2 );
	CPPUNIT_ASSERT( seq.GetProbeMapping() == PM_MM );

	tv = seq.GetParameter(1);
	CPPUNIT_ASSERT( tv.Tag == "tag2" );
	CPPUNIT_ASSERT( tv.Value == "value2");

	bpmap.GetSequenceItem(3, seq);
	CPPUNIT_ASSERT( seq.GetName() == std::string("chr3") );
	CPPUNIT_ASSERT( seq.GroupName() == std::string("HS1") );
	CPPUNIT_ASSERT( seq.GetSeqVersion() == std::string("11_Nov_2005") );
	CPPUNIT_ASSERT( seq.GetNumber() == 3 );
	CPPUNIT_ASSERT( seq.GetNumberHits() == 11 );
	CPPUNIT_ASSERT( seq.GetNumberParameters() == 2 );
	CPPUNIT_ASSERT( seq.GetProbeMapping() == PM_ONLY );

	seq.GetHitItem(1, hit, true);
	CPPUNIT_ASSERT( hit.PMX == 1297 );
	CPPUNIT_ASSERT( hit.PMY == 981 );
	CPPUNIT_ASSERT( CompareFloats(hit.MatchScore, 0.9f) );
	CPPUNIT_ASSERT( hit.PMProbe == std::string("TAAGTAGAGAGATGGATGGTGGTTG") );
	CPPUNIT_ASSERT( hit.Position == 477 );
	CPPUNIT_ASSERT( hit.ProbeLength == 25 );
	CPPUNIT_ASSERT( hit.TopStrand == 1 );

	bpmap.Close();
}
void CBPMAPFileDataTest::testmethod_ReadWhenFileNotExists()
{
	affxbpmap::CBPMAPFileData bpmap;
	std::string path = "test";
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bpmap.Read() == false );
}

void CBPMAPFileDataTest::testmethod_ReadHeader()
{
	affxbpmap::CBPMAPFileData bpmap;
	
	std::string path = GetDataPath(TEST_BPMAP_V1_FILE_NAME);
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", bpmap.Exists() == true );
	CPPUNIT_ASSERT( bpmap.ReadHeader() );
	CPPUNIT_ASSERT( bpmap.GetNumberSequences() == 68);
	CPPUNIT_ASSERT( CompareFloats(bpmap.GetVersion(), 1.0f) );
	bpmap.Close();

	path = GetDataPath(TEST_BPMAP_V2_FILE_NAME);
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", bpmap.Exists() == true );
	CPPUNIT_ASSERT( bpmap.ReadHeader() );
	CPPUNIT_ASSERT( bpmap.GetNumberSequences() == 68);
	CPPUNIT_ASSERT( CompareFloats(bpmap.GetVersion(), 2.0f) );
	bpmap.Close();

	path = GetDataPath(TEST_BPMAP_INCORRECT_HEADER_FILE_NAME);
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", bpmap.Exists() == true );
	CPPUNIT_ASSERT( !bpmap.ReadHeader() );
	bpmap.Close();
}


void CBPMAPFileDataTest::testmethod_ReadHeaderWhenFileNotExists()
{
	affxbpmap::CBPMAPFileData bpmap;
	std::string path = "test";
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bpmap.ReadHeader() == false );
}

void CBPMAPFileDataTest::testmethod_GetErrorFromReadError()
{
	affxbpmap::CBPMAPFileData bpmap;
	std::string path = "test";
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bpmap.Read() == false );
	CPPUNIT_ASSERT( bpmap.GetError() == std::string("Unable to open the file.") );
}

void CBPMAPFileDataTest::testmethod_GetErrorFromReadHeaderError()
{
	affxbpmap::CBPMAPFileData bpmap;
	std::string path = "test";
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT( bpmap.ReadHeader() == false );
	CPPUNIT_ASSERT( bpmap.GetError() == std::string("Unable to open the file.") );
}

void CBPMAPFileDataTest::testmethod_FullName()
{
	affxbpmap::CBPMAPFileData bpmap;
	affxbpmap::CGDACSequenceItem seq;
	affxbpmap::GDACSequenceHitItemType hit;
	TagValuePairType tv;

	// note that TEST_BPMAP_V3_FILE_NAME is rewritten in testing
	// BPMAPFileWriterTest.cpp and this occurs in the test suite
	// after this test, a failure in the BPMAPWriter test that 
	// breaks this test won't be seen until the next time the
	// test suite is run. The test suite should be fixed to run the
	// BPMAPWriter test first... not sure how to do this.

	std::string path = GetDataPath(TEST_BPMAP_V3_FILE_NAME);
	bpmap.SetFileName(path.c_str());
	CPPUNIT_ASSERT_MESSAGE( "file does not exist", bpmap.Exists() == true );
	CPPUNIT_ASSERT( bpmap.Read() );
	CPPUNIT_ASSERT( bpmap.GetNumberSequences() == 4);
	CPPUNIT_ASSERT( CompareFloats(bpmap.GetVersion(), 3.0f) );

	bpmap.GetSequenceItem(0, seq);
	CPPUNIT_ASSERT( seq.GetName() == std::string("chr1") );
	CPPUNIT_ASSERT( seq.GroupName() == std::string("B_subtilis") );
	CPPUNIT_ASSERT( seq.GetSeqVersion() == std::string("1/1/00:11am") );
	CPPUNIT_ASSERT(seq.FullName() == "B_subtilis:1/1/00:11am;chr1");

}
