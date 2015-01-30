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
#include "file/CPPTest/CDFFileDataTest.h" 
//
#include "file/CDFFileData.h"
#include "util/Fs.h"
//
#include <iostream>

CPPUNIT_TEST_SUITE_REGISTRATION( CCDFFileDataTest );

using namespace affxcdf;
using namespace std;

#ifdef _MSC_VER
#define EXP_XDA_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\Test3-xda.CDF"
#define EXP_ASCII_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\Test3-ascii.CDF"
#define NO_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\NoFile.CDF"
#define MAP_XDA_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\Mapping10K_Xba131-xda.CDF"
#define MAP_ASCII_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\Mapping10K_Xba131-ascii.CDF"
#define DMET3_ASCII_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\DMET3_Plus_test_ascii.cdf"
#define DMET3_XDA_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\DMET3_Plus_test_xda.cdf"
#define MULTICHANNEL_ASCII_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\Multichannel_test_ascii.cdf"
#define MULTICHANNEL_XDA_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\Multichannel_test_xda.cdf"
#define AXIOM_V6_ASCII_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\Axiom_v6_test_ascii.r1.cdf"
#define AXIOM_V4_XDA_FILE "..\\..\\calvin_files\\fusion\\data\\CDF\\Axiom_v4_test_xda.r1.cdf"
#else
#define EXP_XDA_FILE "../../calvin_files/fusion/data/CDF/Test3-xda.CDF"
#define EXP_ASCII_FILE "../../calvin_files/fusion/data/CDF/Test3-ascii.CDF"
#define NO_FILE ".../../calvin_files/fusion/data/CDF/NoFile.CDF"
#define MAP_XDA_FILE "../../calvin_files/fusion/data/CDF/Mapping10K_Xba131-xda.CDF"
#define MAP_ASCII_FILE "../../calvin_files/fusion/data/CDF/Mapping10K_Xba131-ascii.CDF"
#define DMET3_ASCII_FILE "../../calvin_files/fusion/data/CDF/DMET3_Plus_test_ascii.cdf"
#define DMET3_XDA_FILE "../../calvin_files/fusion/data/CDF/DMET3_Plus_test_xda.cdf"
#define MULTICHANNEL_ASCII_FILE "../../calvin_files/fusion/data/CDF/Multichannel_test_ascii.cdf"
#define MULTICHANNEL_XDA_FILE "../../calvin_files/fusion/data/CDF/Multichannel_test_xda.cdf"
#define AXIOM_V6_ASCII_FILE "../../calvin_files/fusion/data/CDF/Axiom_v6_test_ascii.r1.cdf"
#define AXIOM_V4_XDA_FILE "../../calvin_files/fusion/data/CDF/Axiom_v4_test_xda.r1.cdf"
#endif

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

void CCDFFileDataTest::setUp()
{
}

void CCDFFileDataTest::tearDown()
{
}

void CCDFFileDataTest::test_CCDFFileHeader()
{
	CCDFFileHeader header;
	CPPUNIT_ASSERT(header.GetCols() == 0);
	CPPUNIT_ASSERT(header.GetRows() == 0);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 0);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 0);
	CPPUNIT_ASSERT(header.GetReference() == "");
}

void CCDFFileDataTest::test_CCDFProbeInformation()
{
	CCDFProbeInformation probe;
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 0);
	CPPUNIT_ASSERT(probe.GetY() == 0);
	CPPUNIT_ASSERT(probe.GetPBase() == ' ');
	CPPUNIT_ASSERT(probe.GetTBase() == ' ');
}

void CCDFFileDataTest::test_CCDFProbeGroupInformation()
{
	CCDFProbeGroupInformation group;
	CPPUNIT_ASSERT(group.GetDirection() == NoDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 0);
	CPPUNIT_ASSERT(group.GetNumCells() == 0);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 0);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 0);
	CPPUNIT_ASSERT(group.GetName() == "");
}

void CCDFFileDataTest::test_CCDFProbeSetInformation()
{
	CCDFProbeSetInformation set;
	CPPUNIT_ASSERT( set.GetProbeSetType() == UnknownProbeSetType );
	CPPUNIT_ASSERT( set.GetDirection() == NoDirection);
	CPPUNIT_ASSERT( set.GetNumLists() == 0 );
	CPPUNIT_ASSERT( set.GetNumGroups() == 0 );
	CPPUNIT_ASSERT( set.GetNumCells() == 0 );
	CPPUNIT_ASSERT( set.GetNumCellsPerList() == 0 );
	CPPUNIT_ASSERT( set.GetProbeSetNumber() == 0 );
}

void CCDFFileDataTest::test_CCDFQCProbeInformation()
{
	CCDFQCProbeInformation probe;
	CPPUNIT_ASSERT(probe.GetX() == 0);
	CPPUNIT_ASSERT(probe.GetY() == 0);
	CPPUNIT_ASSERT(probe.GetPLen() == 0);
	CPPUNIT_ASSERT(probe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(probe.IsBackgroundProbe() == false );
}

void CCDFFileDataTest::test_CCDFQCProbeSetInformation()
{
	CCDFQCProbeSetInformation set;
	CPPUNIT_ASSERT(set.GetNumCells() == 0);
	CPPUNIT_ASSERT(set.GetQCProbeSetType() == UnknownQCProbeSetType);
}

void CCDFFileDataTest::testmethod_IsXDACompatibleFile()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(EXP_XDA_FILE).c_str());
	CPPUNIT_ASSERT(cdf.IsXDACompatibleFile() == true );
	cdf.SetFileName(GetDataPath(EXP_ASCII_FILE).c_str());
	CPPUNIT_ASSERT(cdf.IsXDACompatibleFile() == false );
}

void CCDFFileDataTest::testmethod_Exists()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(EXP_XDA_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Exists() == true );
	cdf.SetFileName(GetDataPath(NO_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Exists() == false );
}

void CCDFFileDataTest::testmethod_ReadHeader_with_ASCII()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(EXP_ASCII_FILE).c_str());
	CPPUNIT_ASSERT(cdf.ReadHeader() == true );
	CPPUNIT_ASSERT(cdf.GetChipType() == "Test3-ascii");
	CCDFFileHeader &header = cdf.GetHeader();
	CPPUNIT_ASSERT(header.GetCols() == 126);
	CPPUNIT_ASSERT(header.GetRows() == 126);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 345);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 13);
	CPPUNIT_ASSERT(header.GetReference() == "");
}

void CCDFFileDataTest::testmethod_ReadHeader_with_XDA()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(EXP_XDA_FILE).c_str());
	CPPUNIT_ASSERT(cdf.ReadHeader() == true );
	CPPUNIT_ASSERT(cdf.GetChipType() == "Test3-xda");
	CCDFFileHeader &header = cdf.GetHeader();
	CPPUNIT_ASSERT(header.GetCols() == 126);
	CPPUNIT_ASSERT(header.GetRows() == 126);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 345);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 13);
	CPPUNIT_ASSERT(header.GetReference() == "");
}

void CCDFFileDataTest::test_ExpressionXDA()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(EXP_XDA_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );

	CPPUNIT_ASSERT(cdf.GetChipType() == "Test3-xda");
	CCDFFileHeader &header = cdf.GetHeader();
	CPPUNIT_ASSERT(header.GetCols() == 126);
	CPPUNIT_ASSERT(header.GetRows() == 126);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 345);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 13);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "Pae_16SrRNA_s_at");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "Pae_23SrRNA_s_at");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "PA1178_oprH_at");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "AFFX_ratb2/X14115_at");

	for (int i=0; i<header.GetNumProbeSets(); i++)
	{
		CPPUNIT_ASSERT(cdf.GetProbeSetType(i) == ExpressionProbeSetType);
	}

	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;

	cdf.GetProbeSetInformation(0, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 16 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 32 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 1000 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 16);
	CPPUNIT_ASSERT(group.GetNumCells() == 32);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "Pae_16SrRNA_s_at");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 111);
	CPPUNIT_ASSERT(probe.GetY() == 79);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 111);
	CPPUNIT_ASSERT(probe.GetY() == 80);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 93);
	CPPUNIT_ASSERT(probe.GetY() == 82);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');



	cdf.GetProbeSetInformation(1, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 16 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 32 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 1001 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 16);
	CPPUNIT_ASSERT(group.GetNumCells() == 32);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "Pae_23SrRNA_s_at");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 124);
	CPPUNIT_ASSERT(probe.GetY() == 95);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 124);
	CPPUNIT_ASSERT(probe.GetY() == 96);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 36);
	CPPUNIT_ASSERT(probe.GetY() == 8);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');



	cdf.GetProbeSetInformation(header.GetNumProbeSets()-1, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 20 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 40 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 3101 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 20);
	CPPUNIT_ASSERT(group.GetNumCells() == 40);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "AFFX_ratb2/X14115_at");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 20);
	CPPUNIT_ASSERT(probe.GetY() == 113);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 20);
	CPPUNIT_ASSERT(probe.GetY() == 114);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 119);
	CPPUNIT_ASSERT(probe.GetY() == 58);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');


	CCDFQCProbeSetInformation qcset;
	CCDFQCProbeInformation qcprobe;

	cdf.GetQCProbeSetInformation(0, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 300);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == GeneExpNegativeQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 77);
	CPPUNIT_ASSERT(qcprobe.GetY() == 82);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 20);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 77);
	CPPUNIT_ASSERT(qcprobe.GetY() == 83);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 20);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == true);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 15);
	CPPUNIT_ASSERT(qcprobe.GetY() == 86);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 1);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == true );


	cdf.GetQCProbeSetInformation(header.GetNumQCProbeSets()-1, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 9);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == CentralCrossNegativeQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 60);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 61);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 64);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );


	cdf.GetQCProbeSetInformation(CentralCrossNegativeQCProbeSetType, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 9);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == CentralCrossNegativeQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 60);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 61);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 64);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );
}

void CCDFFileDataTest::test_ExpressionASCII()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(EXP_ASCII_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );

	CPPUNIT_ASSERT(cdf.GetChipType() == "Test3-ascii");
	CCDFFileHeader &header = cdf.GetHeader();
	CPPUNIT_ASSERT(header.GetCols() == 126);
	CPPUNIT_ASSERT(header.GetRows() == 126);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 345);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 13);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "Pae_16SrRNA_s_at");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "Pae_23SrRNA_s_at");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "PA1178_oprH_at");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "AFFX_ratb2/X14115_at");

	for (int i=0; i<header.GetNumProbeSets(); i++)
	{
		CPPUNIT_ASSERT(cdf.GetProbeSetType(i) == ExpressionProbeSetType);
	}

	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;

	cdf.GetProbeSetInformation(0, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 16 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 32 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 1000 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 16);
	CPPUNIT_ASSERT(group.GetNumCells() == 32);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "Pae_16SrRNA_s_at");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 111);
	CPPUNIT_ASSERT(probe.GetY() == 79);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 111);
	CPPUNIT_ASSERT(probe.GetY() == 80);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 93);
	CPPUNIT_ASSERT(probe.GetY() == 82);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');



	cdf.GetProbeSetInformation(1, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 16 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 32 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 1001 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 16);
	CPPUNIT_ASSERT(group.GetNumCells() == 32);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "Pae_23SrRNA_s_at");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 124);
	CPPUNIT_ASSERT(probe.GetY() == 95);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 124);
	CPPUNIT_ASSERT(probe.GetY() == 96);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 36);
	CPPUNIT_ASSERT(probe.GetY() == 8);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');



	cdf.GetProbeSetInformation(header.GetNumProbeSets()-1, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 20 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 40 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 3101 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 20);
	CPPUNIT_ASSERT(group.GetNumCells() == 40);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "AFFX_ratb2/X14115_at");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 20);
	CPPUNIT_ASSERT(probe.GetY() == 113);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 20);
	CPPUNIT_ASSERT(probe.GetY() == 114);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 119);
	CPPUNIT_ASSERT(probe.GetY() == 58);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');


	CCDFQCProbeSetInformation qcset;
	CCDFQCProbeInformation qcprobe;

	cdf.GetQCProbeSetInformation(0, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 300);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == GeneExpNegativeQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 77);
	CPPUNIT_ASSERT(qcprobe.GetY() == 82);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 20);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 77);
	CPPUNIT_ASSERT(qcprobe.GetY() == 83);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 20);

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 15);
	CPPUNIT_ASSERT(qcprobe.GetY() == 86);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 1);

	cdf.GetQCProbeSetInformation(header.GetNumQCProbeSets()-1, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 9);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == CentralCrossNegativeQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 60);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 61);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 64);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	cdf.GetQCProbeSetInformation(CentralCrossNegativeQCProbeSetType, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 9);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == CentralCrossNegativeQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 60);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 61);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 62);
	CPPUNIT_ASSERT(qcprobe.GetY() == 64);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);
}

void CCDFFileDataTest::test_GenotypingXDA()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(MAP_XDA_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );

	CPPUNIT_ASSERT(cdf.GetChipType() == "Mapping10K_Xba131-xda");
	CCDFFileHeader &header = cdf.GetHeader();
	CPPUNIT_ASSERT(header.GetCols() == 712);
	CPPUNIT_ASSERT(header.GetRows() == 712);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 11564);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 9);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "AFFX-5Q-123");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "AFFX-5Q-456");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "AFFX-5Q-789");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "SNP_A-1508078");

	for (int i=0; i<header.GetNumProbeSets(); i++)
	{
		CPPUNIT_ASSERT(cdf.GetProbeSetType(i) == GenotypingProbeSetType);
	}

	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;

	cdf.GetProbeSetInformation(0, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == GenotypingProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 30 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 60 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 41 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 30);
	CPPUNIT_ASSERT(group.GetNumCells() == 60);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "AFFX-5Q-123");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 12);
	CPPUNIT_ASSERT(probe.GetX() == 323);
	CPPUNIT_ASSERT(probe.GetY() == 386);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 12);
	CPPUNIT_ASSERT(probe.GetX() == 323);
	CPPUNIT_ASSERT(probe.GetY() == 387);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 46);
	CPPUNIT_ASSERT(probe.GetX() == 337);
	CPPUNIT_ASSERT(probe.GetY() == 389);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');


	cdf.GetProbeSetInformation(header.GetNumProbeSets()-1, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == GenotypingProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 20 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 4 );
	CPPUNIT_ASSERT(set.GetNumCells() == 40 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 12548 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 5);
	CPPUNIT_ASSERT(group.GetNumCells() == 10);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 4);
	CPPUNIT_ASSERT(group.GetName() == "SNP_A-1508078A");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 642);
	CPPUNIT_ASSERT(probe.GetY() == 665);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 642);
	CPPUNIT_ASSERT(probe.GetY() == 666);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 21);
	CPPUNIT_ASSERT(probe.GetX() == 172);
	CPPUNIT_ASSERT(probe.GetY() == 699);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');


	set.GetGroupInformation(3, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 5);
	CPPUNIT_ASSERT(group.GetNumCells() == 10);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 15);
	CPPUNIT_ASSERT(group.GetStop() == 19);
	CPPUNIT_ASSERT(group.GetName() == "SNP_A-1508078G");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 15);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 648);
	CPPUNIT_ASSERT(probe.GetY() == 604);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 15);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 648);
	CPPUNIT_ASSERT(probe.GetY() == 603);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 19);
	CPPUNIT_ASSERT(probe.GetExpos() == 21);
	CPPUNIT_ASSERT(probe.GetX() == 220);
	CPPUNIT_ASSERT(probe.GetY() == 686);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');


	CCDFQCProbeSetInformation qcset;
	CCDFQCProbeInformation qcprobe;

	cdf.GetQCProbeSetInformation(0, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 24);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == HybPositiveQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 473);
	CPPUNIT_ASSERT(qcprobe.GetY() == 2);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);
	CPPUNIT_ASSERT(qcprobe.IsPerfectMatchProbe() == false);
	CPPUNIT_ASSERT(qcprobe.IsBackgroundProbe() == false );

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 236);
	CPPUNIT_ASSERT(qcprobe.GetY() == 2);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 711);
	CPPUNIT_ASSERT(qcprobe.GetY() == 472);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

	cdf.GetQCProbeSetInformation(header.GetNumQCProbeSets()-1, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 9);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == CentralCrossNegativeQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 356);
	CPPUNIT_ASSERT(qcprobe.GetY() == 354);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 356);
	CPPUNIT_ASSERT(qcprobe.GetY() == 355);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 356);
	CPPUNIT_ASSERT(qcprobe.GetY() == 358);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	cdf.GetQCProbeSetInformation(HybPositiveQCProbeSetType, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 24);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == HybPositiveQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 473);
	CPPUNIT_ASSERT(qcprobe.GetY() == 2);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 236);
	CPPUNIT_ASSERT(qcprobe.GetY() == 2);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 711);
	CPPUNIT_ASSERT(qcprobe.GetY() == 472);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

}

void CCDFFileDataTest::test_GenotypingASCII()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(MAP_ASCII_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );

	CPPUNIT_ASSERT(cdf.GetChipType() == "Mapping10K_Xba131-ascii");
	CCDFFileHeader &header = cdf.GetHeader();
	CPPUNIT_ASSERT(header.GetCols() == 712);
	CPPUNIT_ASSERT(header.GetRows() == 712);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 11564);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 9);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "AFFX-5Q-123");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "AFFX-5Q-456");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "AFFX-5Q-789");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "SNP_A-1508078");

	for (int i=0; i<header.GetNumProbeSets(); i++)
	{
		CPPUNIT_ASSERT(cdf.GetProbeSetType(i) == GenotypingProbeSetType);
	}

	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;

	cdf.GetProbeSetInformation(0, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == GenotypingProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 30 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 60 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 41 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 30);
	CPPUNIT_ASSERT(group.GetNumCells() == 60);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "AFFX-5Q-123");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 12);
	CPPUNIT_ASSERT(probe.GetX() == 323);
	CPPUNIT_ASSERT(probe.GetY() == 386);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 12);
	CPPUNIT_ASSERT(probe.GetX() == 323);
	CPPUNIT_ASSERT(probe.GetY() == 387);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 46);
	CPPUNIT_ASSERT(probe.GetX() == 337);
	CPPUNIT_ASSERT(probe.GetY() == 389);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');


	cdf.GetProbeSetInformation(header.GetNumProbeSets()-1, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == GenotypingProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 20 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 4 );
	CPPUNIT_ASSERT(set.GetNumCells() == 40 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 12548 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 5);
	CPPUNIT_ASSERT(group.GetNumCells() == 10);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 4);
	CPPUNIT_ASSERT(group.GetName() == "SNP_A-1508078A");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 642);
	CPPUNIT_ASSERT(probe.GetY() == 665);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 642);
	CPPUNIT_ASSERT(probe.GetY() == 666);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 21);
	CPPUNIT_ASSERT(probe.GetX() == 172);
	CPPUNIT_ASSERT(probe.GetY() == 699);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');


	set.GetGroupInformation(3, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 5);
	CPPUNIT_ASSERT(group.GetNumCells() == 10);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 15);
	CPPUNIT_ASSERT(group.GetStop() == 19);
	CPPUNIT_ASSERT(group.GetName() == "SNP_A-1508078G");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 15);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 648);
	CPPUNIT_ASSERT(probe.GetY() == 604);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 15);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 648);
	CPPUNIT_ASSERT(probe.GetY() == 603);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 19);
	CPPUNIT_ASSERT(probe.GetExpos() == 21);
	CPPUNIT_ASSERT(probe.GetX() == 220);
	CPPUNIT_ASSERT(probe.GetY() == 686);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');


	CCDFQCProbeSetInformation qcset;
	CCDFQCProbeInformation qcprobe;

	cdf.GetQCProbeSetInformation(0, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 24);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == HybPositiveQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 473);
	CPPUNIT_ASSERT(qcprobe.GetY() == 2);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 236);
	CPPUNIT_ASSERT(qcprobe.GetY() == 2);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 711);
	CPPUNIT_ASSERT(qcprobe.GetY() == 472);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

	cdf.GetQCProbeSetInformation(header.GetNumQCProbeSets()-1, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 9);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == CentralCrossNegativeQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 356);
	CPPUNIT_ASSERT(qcprobe.GetY() == 354);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 356);
	CPPUNIT_ASSERT(qcprobe.GetY() == 355);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 356);
	CPPUNIT_ASSERT(qcprobe.GetY() == 358);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 25);

	cdf.GetQCProbeSetInformation(HybPositiveQCProbeSetType, qcset);
	CPPUNIT_ASSERT(qcset.GetNumCells() == 24);
	CPPUNIT_ASSERT(qcset.GetQCProbeSetType() == HybPositiveQCProbeSetType);

	qcset.GetProbeInformation(0, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 473);
	CPPUNIT_ASSERT(qcprobe.GetY() == 2);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

	qcset.GetProbeInformation(1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 236);
	CPPUNIT_ASSERT(qcprobe.GetY() == 2);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);

	qcset.GetProbeInformation(qcset.GetNumCells()-1, qcprobe);
	CPPUNIT_ASSERT(qcprobe.GetX() == 711);
	CPPUNIT_ASSERT(qcprobe.GetY() == 472);
	CPPUNIT_ASSERT(qcprobe.GetPLen() == 16);
}

void CCDFFileDataTest::test_ReseqXDA()
{
}

void CCDFFileDataTest::test_ReseqASCII()
{
}

void CCDFFileDataTest::test_TagXDA()
{
}

void CCDFFileDataTest::test_TagASCII()
{
}

void CCDFFileDataTest::test_DMET3ASCII()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(DMET3_ASCII_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );

	CPPUNIT_ASSERT(cdf.GetChipType() == "DMET3_Plus_test_ascii");
	CCDFFileHeader &header = cdf.GetHeader();
   	CPPUNIT_ASSERT(header.GetCols() == 1050);
	CPPUNIT_ASSERT(header.GetRows() == 1050);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 21);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 4);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "AFFX-5Q-123");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "AFFX-RandomGC23");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "chr3.70290");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "AFR_fake_no_opt-08");

	
	CPPUNIT_ASSERT(cdf.GetProbeSetType(0) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(1) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(2) == CopyNumberProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(3) == CopyNumberProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(4) == CopyNumberProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(5) == CopyNumberProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(6) == GenotypingProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(7) == GenotypingProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(8) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(9) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(10) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(11) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(12) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(13) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(14) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(15) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(16) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(17) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(18) == GenotypingProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(19) == GenotypingProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(20) == ExpressionProbeSetType);
	
	 
	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;
    
	//ExpressionProbeSetType
	cdf.GetProbeSetInformation(0, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 30 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 60 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 41 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 30);
	CPPUNIT_ASSERT(group.GetNumCells() == 60);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "AFFX-5Q-123");

	group.GetCell(0, probe);//PM
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 590);
	CPPUNIT_ASSERT(probe.GetY() == 532);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
    CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);//MM
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 591);
	CPPUNIT_ASSERT(probe.GetY() == 532);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

	group.GetCell(group.GetNumCells()-2, probe);//PM
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 592); 
	CPPUNIT_ASSERT(probe.GetY() == 518);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
    CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

	group.GetCell(group.GetNumCells()-1, probe);//MM
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 593); 
	CPPUNIT_ASSERT(probe.GetY() == 518);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
    CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

    
	//CopyNumberProbeSetType
	cdf.GetProbeSetInformation(2, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == CopyNumberProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 24 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 24 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 2084 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 24);
	CPPUNIT_ASSERT(group.GetNumCells() == 24);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "chr3.70290");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 67);
	CPPUNIT_ASSERT(probe.GetX() == 560);
	CPPUNIT_ASSERT(probe.GetY() == 606);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 14);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 67);
	CPPUNIT_ASSERT(probe.GetX() == 993);
	CPPUNIT_ASSERT(probe.GetY() == 658);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 16);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 241);
	CPPUNIT_ASSERT(probe.GetX() == 925); 
	CPPUNIT_ASSERT(probe.GetY() == 517);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 11);


	//GenotypingProbeSetType
	cdf.GetProbeSetInformation(6, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == GenotypingProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == EitherDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 540 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 4 );
	CPPUNIT_ASSERT(set.GetNumCells() == 540 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 5011 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 134);
	CPPUNIT_ASSERT(group.GetName() == "C");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 375);
	CPPUNIT_ASSERT(probe.GetY() == 8);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 2);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 43);
	CPPUNIT_ASSERT(probe.GetY() == 9);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 3);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 9);
	CPPUNIT_ASSERT(probe.GetY() == 550);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 50);


	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 135);
	CPPUNIT_ASSERT(group.GetStop() == 269);
	CPPUNIT_ASSERT(group.GetName() == "T");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 135);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 374);
	CPPUNIT_ASSERT(probe.GetY() == 8);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 2);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 136);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 42);
	CPPUNIT_ASSERT(probe.GetY() == 9);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 3);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 269);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 8);
	CPPUNIT_ASSERT(probe.GetY() == 550);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 50);


	set.GetGroupInformation(2, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 270);
	CPPUNIT_ASSERT(group.GetStop() == 404);
	CPPUNIT_ASSERT(group.GetName() == "C");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 270);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 411);
	CPPUNIT_ASSERT(probe.GetY() == 227);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 22);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 271);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 900);
	CPPUNIT_ASSERT(probe.GetY() == 250);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 30);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 404);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 774);
	CPPUNIT_ASSERT(probe.GetY() == 992);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 97);


	set.GetGroupInformation(3, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 405);
	CPPUNIT_ASSERT(group.GetStop() == 539);
	CPPUNIT_ASSERT(group.GetName() == "T");

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 405);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 410);
	CPPUNIT_ASSERT(probe.GetY() == 227);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 22);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 406);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 899);
	CPPUNIT_ASSERT(probe.GetY() == 250);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 30);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 539);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 773);
	CPPUNIT_ASSERT(probe.GetY() == 992);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 97);

    //MarkerProbeSetType//CYP2A6_27_G
    cdf.GetProbeSetInformation(8, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == EitherDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 270 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 270 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 5346 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 134);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 11);
	CPPUNIT_ASSERT(probe.GetX() == 971);
	CPPUNIT_ASSERT(probe.GetY() == 99);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 31);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 11);
	CPPUNIT_ASSERT(probe.GetX() == 682);
	CPPUNIT_ASSERT(probe.GetY() == 709);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 135);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 25);
	CPPUNIT_ASSERT(probe.GetX() == 220);
	CPPUNIT_ASSERT(probe.GetY() == 294);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 66); 
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 135);
	CPPUNIT_ASSERT(group.GetStop() == 269);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 135);
	CPPUNIT_ASSERT(probe.GetExpos() == 11);
	CPPUNIT_ASSERT(probe.GetX() == 709);
	CPPUNIT_ASSERT(probe.GetY() == 73);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 18);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 136);
	CPPUNIT_ASSERT(probe.GetExpos() == 11);
	CPPUNIT_ASSERT(probe.GetX() == 1040);
	CPPUNIT_ASSERT(probe.GetY() == 103);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 33);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 269);
	CPPUNIT_ASSERT(probe.GetExpos() == 25);
	CPPUNIT_ASSERT(probe.GetX() == 1025);
	CPPUNIT_ASSERT(probe.GetY() == 1029);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 187);
  
	//MarkerProbeSetType//DMET3B10768
    cdf.GetProbeSetInformation(12, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == EitherDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6480 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 48 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6480 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 7845 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 134);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 1);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 207);
	CPPUNIT_ASSERT(probe.GetY() == 79);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 40);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 543);
	CPPUNIT_ASSERT(probe.GetY() == 239);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 143);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 107); 
	CPPUNIT_ASSERT(probe.GetY() == 761);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 475); 
    
	set.GetGroupInformation(47, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 6345);
	CPPUNIT_ASSERT(group.GetStop() == 6479);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 9);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6345);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 388);
	CPPUNIT_ASSERT(probe.GetY() == 323);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 202);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6346);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 915);
	CPPUNIT_ASSERT(probe.GetY() == 502);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 313);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6479);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 755);
	CPPUNIT_ASSERT(probe.GetY() == 905);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 557);
}

void CCDFFileDataTest::test_DMET3XDA()
{
	CCDFFileData cdf; 
	cdf.SetFileName(GetDataPath(DMET3_XDA_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );
    CPPUNIT_ASSERT(cdf.IsXDACompatibleFile() == true );
	CPPUNIT_ASSERT(cdf.GetChipType() == "DMET3_Plus_test_xda");
	CCDFFileHeader &header = cdf.GetHeader();
   	CPPUNIT_ASSERT(header.GetCols() == 1050);
	CPPUNIT_ASSERT(header.GetRows() == 1050);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 21);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 4);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "AFFX-5Q-123");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "AFFX-RandomGC23");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "chr3.70290");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "AFR_fake_no_opt-08");
    	
	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;

	CPPUNIT_ASSERT(cdf.GetProbeSetType(0) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(1) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(2) == CopyNumberProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(3) == CopyNumberProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(4) == CopyNumberProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(5) == CopyNumberProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(6) == GenotypingProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(7) == GenotypingProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(8) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(9) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(10) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(11) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(12) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(13) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(14) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(15) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(16) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(17) == MarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(18) == GenotypingProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(19) == GenotypingProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(20) == ExpressionProbeSetType);
    
	//ExpressionProbeSetType
	cdf.GetProbeSetInformation(0, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 30 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 60 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 2 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 41 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 30);
	CPPUNIT_ASSERT(group.GetNumCells() == 60);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 2);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "AFFX-5Q-123");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
    
	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 591);
	CPPUNIT_ASSERT(probe.GetY() == 532);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 590);
	CPPUNIT_ASSERT(probe.GetY() == 532);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 592);
	CPPUNIT_ASSERT(probe.GetY() == 518);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
    CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

    //CopyNumberProbeSetType
	cdf.GetProbeSetInformation(2, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == CopyNumberProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 24 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 24 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 2084 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 24);
	CPPUNIT_ASSERT(group.GetNumCells() == 24);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "chr3.70290");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 67);
	CPPUNIT_ASSERT(probe.GetX() == 560);
	CPPUNIT_ASSERT(probe.GetY() == 606);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 14);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 67);
	CPPUNIT_ASSERT(probe.GetX() == 993);
	CPPUNIT_ASSERT(probe.GetY() == 658);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 16);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 241);
	CPPUNIT_ASSERT(probe.GetX() == 925); 
	CPPUNIT_ASSERT(probe.GetY() == 517);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 11);
 
	//GenotypingProbeSetType
	cdf.GetProbeSetInformation(6, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == GenotypingProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == EitherDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 540 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 4 );
	CPPUNIT_ASSERT(set.GetNumCells() == 540 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 5011 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 134);
	CPPUNIT_ASSERT(group.GetName() == "C");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 375);
	CPPUNIT_ASSERT(probe.GetY() == 8);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 2);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 43);
	CPPUNIT_ASSERT(probe.GetY() == 9);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 3);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 9);
	CPPUNIT_ASSERT(probe.GetY() == 550);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 50);


	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 135);
	CPPUNIT_ASSERT(group.GetStop() == 269);
	CPPUNIT_ASSERT(group.GetName() == "T");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 135);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 374);
	CPPUNIT_ASSERT(probe.GetY() == 8);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 2);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 136);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 42);
	CPPUNIT_ASSERT(probe.GetY() == 9);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 3);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 269);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 8);
	CPPUNIT_ASSERT(probe.GetY() == 550);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 50);


	set.GetGroupInformation(2, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 270);
	CPPUNIT_ASSERT(group.GetStop() == 404);
	CPPUNIT_ASSERT(group.GetName() == "C");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 270);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 411);
	CPPUNIT_ASSERT(probe.GetY() == 227);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 22);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 271);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 900);
	CPPUNIT_ASSERT(probe.GetY() == 250);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 30);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 404);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 774);
	CPPUNIT_ASSERT(probe.GetY() == 992);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 97);


	set.GetGroupInformation(3, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 405);
	CPPUNIT_ASSERT(group.GetStop() == 539);
	CPPUNIT_ASSERT(group.GetName() == "T");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 405);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 410);
	CPPUNIT_ASSERT(probe.GetY() == 227);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 22);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 406);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 899);
	CPPUNIT_ASSERT(probe.GetY() == 250);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 30);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 539);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 773);
	CPPUNIT_ASSERT(probe.GetY() == 992);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 97);

	//MarkerProbeSetType//CYP2A6_27_G
    cdf.GetProbeSetInformation(8, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == EitherDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 270 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 270 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 5346 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 134);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 11);
	CPPUNIT_ASSERT(probe.GetX() == 971);
	CPPUNIT_ASSERT(probe.GetY() == 99);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 31);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 11);
	CPPUNIT_ASSERT(probe.GetX() == 682);
	CPPUNIT_ASSERT(probe.GetY() == 709);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 135);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 25);
	CPPUNIT_ASSERT(probe.GetX() == 220);
	CPPUNIT_ASSERT(probe.GetY() == 294);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 66); 
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 135);
	CPPUNIT_ASSERT(group.GetStop() == 269);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 135);
	CPPUNIT_ASSERT(probe.GetExpos() == 11);
	CPPUNIT_ASSERT(probe.GetX() == 709);
	CPPUNIT_ASSERT(probe.GetY() == 73);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 18);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 136);
	CPPUNIT_ASSERT(probe.GetExpos() == 11);
	CPPUNIT_ASSERT(probe.GetX() == 1040);
	CPPUNIT_ASSERT(probe.GetY() == 103);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 33);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 269);
	CPPUNIT_ASSERT(probe.GetExpos() == 25);
	CPPUNIT_ASSERT(probe.GetX() == 1025);
	CPPUNIT_ASSERT(probe.GetY() == 1029);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 187);
  
	//MarkerProbeSetType DMET3B10768
    cdf.GetProbeSetInformation(12, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == EitherDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6480 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 48 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6480 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 7845 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 134);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 1);
	CPPUNIT_ASSERT((char)group.GetAlleleCode() == 'A');

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 207);
	CPPUNIT_ASSERT(probe.GetY() == 79);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 40);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 543);
	CPPUNIT_ASSERT(probe.GetY() == 239);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 143);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 107); 
	CPPUNIT_ASSERT(probe.GetY() == 761);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 475); 
    
	set.GetGroupInformation(47, group);
	CPPUNIT_ASSERT(group.GetDirection() == AntiSenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 135);
	CPPUNIT_ASSERT(group.GetNumCells() == 135);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 6345);
	CPPUNIT_ASSERT(group.GetStop() == 6479);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 9);
	CPPUNIT_ASSERT((char)group.GetAlleleCode() == 'B');

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6345);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 388);
	CPPUNIT_ASSERT(probe.GetY() == 323);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 202);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6346);
	CPPUNIT_ASSERT(probe.GetExpos() == 15);
	CPPUNIT_ASSERT(probe.GetX() == 915);
	CPPUNIT_ASSERT(probe.GetY() == 502);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 21);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 313);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6479);
	CPPUNIT_ASSERT(probe.GetExpos() == 23);
	CPPUNIT_ASSERT(probe.GetX() == 755);
	CPPUNIT_ASSERT(probe.GetY() == 905);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 29);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 557);
}

void CCDFFileDataTest::test_MultichannelXDA()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(MULTICHANNEL_XDA_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );

	CPPUNIT_ASSERT(cdf.GetChipType() == "Multichannel_test_xda");
	CCDFFileHeader &header = cdf.GetHeader();
   	CPPUNIT_ASSERT(header.GetCols() == 1462);
	CPPUNIT_ASSERT(header.GetRows() == 1487);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 25);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 0);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "10151091AC_r");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "10151091GT_f");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "10500817AG_f");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "snp_rs1022296_0w_r");
	
	CPPUNIT_ASSERT(cdf.GetProbeSetType(0) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(1) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(2) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(3) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(4) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(5) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(6) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(7) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(8) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(9) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(header.GetNumProbeSets()-1) == MultichannelMarkerProbeSetType);	
	 
	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;
    
	//ExpressionProbeSetType
	cdf.GetProbeSetInformation(4, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 200 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 200 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 4 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 200);
	CPPUNIT_ASSERT(group.GetNumCells() == 200);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "AFFX-BkGr-GC11");

	group.GetCell(0, probe);//PM
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 692);
	CPPUNIT_ASSERT(probe.GetY() == 1167);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
    CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

	group.GetCell(group.GetNumCells()-1, probe);//PM
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 489); 
	CPPUNIT_ASSERT(probe.GetY() == 1059);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
    CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);    

    //MultichannelProbeSetType //indel_rs10707934_1wa0_r
    cdf.GetProbeSetInformation(13, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 531269 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1456);
	CPPUNIT_ASSERT(probe.GetY() == 1229);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 220);
	CPPUNIT_ASSERT(probe.GetY() == 1230);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 471);
	CPPUNIT_ASSERT(probe.GetY() == 1230);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0); 
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 6);
	CPPUNIT_ASSERT(group.GetStop() == 8);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1456);
	CPPUNIT_ASSERT(probe.GetY() == 1229);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 7);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 220);
	CPPUNIT_ASSERT(probe.GetY() == 1230);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 8);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 471);
	CPPUNIT_ASSERT(probe.GetY() == 1230);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0); 
  
	//MultichannelProbeSetType //snp_rs1022296_0w_f
    cdf.GetProbeSetInformation(23, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 77939 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 579);
	CPPUNIT_ASSERT(probe.GetY() == 97);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 224);
	CPPUNIT_ASSERT(probe.GetY() == 43);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1148); 
	CPPUNIT_ASSERT(probe.GetY() == 111);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0); 
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 3);
	CPPUNIT_ASSERT(group.GetStop() == 5);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 3);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 579);
	CPPUNIT_ASSERT(probe.GetY() == 96);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 4);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 224);
	CPPUNIT_ASSERT(probe.GetY() == 42);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 5);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1148); 
	CPPUNIT_ASSERT(probe.GetY() == 110);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0); 

	//MultichannelProbeSetType //AX_as_test1_f
    cdf.GetProbeSetInformation(8, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 95001 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == MixedRepType);

	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 292);
	CPPUNIT_ASSERT(probe.GetY() == 81);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 3);
	CPPUNIT_ASSERT(group.GetStop() == 5);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == MixedRepType);

	group.GetCell(2, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 5);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 700);
	CPPUNIT_ASSERT(probe.GetY() == 156);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

	//MultichannelProbeSetType //AX_as_test2_f
    cdf.GetProbeSetInformation(9, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 95002 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == DifferentRepType);

	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 292);
	CPPUNIT_ASSERT(probe.GetY() == 81);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 3);
	CPPUNIT_ASSERT(group.GetStop() == 5);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == DifferentRepType);

	group.GetCell(2, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 5);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 700);
	CPPUNIT_ASSERT(probe.GetY() == 156);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

	//MultichannelProbeSetType //AX_nas_test1_r
    cdf.GetProbeSetInformation(10, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 95003 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == DifferentRepType);

	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 18);
	CPPUNIT_ASSERT(probe.GetX() == 799);
	CPPUNIT_ASSERT(probe.GetY() == 1191);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 35);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 6);
	CPPUNIT_ASSERT(group.GetStop() == 8);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == DifferentRepType);

	group.GetCell(2, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 8);
	CPPUNIT_ASSERT(probe.GetExpos() == 18);
	CPPUNIT_ASSERT(probe.GetX() == 991);
	CPPUNIT_ASSERT(probe.GetY() == 1191);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 35);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
}

void CCDFFileDataTest::test_MultichannelASCII()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(MULTICHANNEL_ASCII_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );

	CPPUNIT_ASSERT(cdf.GetChipType() == "Multichannel_test_ascii");
	CCDFFileHeader &header = cdf.GetHeader();
   	CPPUNIT_ASSERT(header.GetCols() == 1462);
	CPPUNIT_ASSERT(header.GetRows() == 1487);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 25);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 0);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "10151091AC_r");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "10151091GT_f");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "10500817AG_f");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "snp_rs1022296_0w_r");
	
	CPPUNIT_ASSERT(cdf.GetProbeSetType(0) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(1) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(2) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(3) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(4) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(5) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(6) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(7) == ExpressionProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(8) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(9) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(header.GetNumProbeSets()-1) == MultichannelMarkerProbeSetType);	
	 
	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;
    
	//ExpressionProbeSetType
	cdf.GetProbeSetInformation(4, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == ExpressionProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 200 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 1 );
	CPPUNIT_ASSERT(set.GetNumCells() == 200 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 4 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 200);
	CPPUNIT_ASSERT(group.GetNumCells() == 200);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(group.GetName() == "AFFX-BkGr-GC11");

	group.GetCell(0, probe);//PM
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 0);
	CPPUNIT_ASSERT(probe.GetX() == 692);
	CPPUNIT_ASSERT(probe.GetY() == 1167);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
    CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

	group.GetCell(group.GetNumCells()-1, probe);//PM
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetX() == 489); 
	CPPUNIT_ASSERT(probe.GetY() == 1059);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 't');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'a');
    CPPUNIT_ASSERT(probe.GetProbeLength() == 25);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);    

    //MultichannelProbeSetType //indel_rs10707934_1wa0_r
    cdf.GetProbeSetInformation(13, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 531269 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 0);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1456);
	CPPUNIT_ASSERT(probe.GetY() == 1229);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 220);
	CPPUNIT_ASSERT(probe.GetY() == 1230);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 471);
	CPPUNIT_ASSERT(probe.GetY() == 1230);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0); 
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 6);
	CPPUNIT_ASSERT(group.GetStop() == 8);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1456);
	CPPUNIT_ASSERT(probe.GetY() == 1229);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 7);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 220);
	CPPUNIT_ASSERT(probe.GetY() == 1230);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 8);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 471);
	CPPUNIT_ASSERT(probe.GetY() == 1230);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0); 
  
	//MultichannelProbeSetType //snp_rs1022296_0w_f
    cdf.GetProbeSetInformation(23, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 77939 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 579);
	CPPUNIT_ASSERT(probe.GetY() == 97);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 224);
	CPPUNIT_ASSERT(probe.GetY() == 43);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1148); 
	CPPUNIT_ASSERT(probe.GetY() == 111);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0); 
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 3);
	CPPUNIT_ASSERT(group.GetStop() == 5);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 1);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 3);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 579);
	CPPUNIT_ASSERT(probe.GetY() == 96);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 4);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 224);
	CPPUNIT_ASSERT(probe.GetY() == 42);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 5);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1148); 
	CPPUNIT_ASSERT(probe.GetY() == 110);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0); 

	//MultichannelProbeSetType //AX_as_test1_f
    cdf.GetProbeSetInformation(8, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 95001 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == MixedRepType);

	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 292);
	CPPUNIT_ASSERT(probe.GetY() == 81);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 3);
	CPPUNIT_ASSERT(group.GetStop() == 5);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == MixedRepType);

	group.GetCell(2, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 5);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 700);
	CPPUNIT_ASSERT(probe.GetY() == 156);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

	//MultichannelProbeSetType //AX_as_test2_f
    cdf.GetProbeSetInformation(9, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 95002 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == DifferentRepType);

	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 292);
	CPPUNIT_ASSERT(probe.GetY() == 81);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 3);
	CPPUNIT_ASSERT(group.GetStop() == 5);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == DifferentRepType);

	group.GetCell(2, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 5);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 700);
	CPPUNIT_ASSERT(probe.GetY() == 156);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);

	//MultichannelProbeSetType //AX_nas_test1_r
    cdf.GetProbeSetInformation(10, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 95003 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == DifferentRepType);

	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 18);
	CPPUNIT_ASSERT(probe.GetX() == 799);
	CPPUNIT_ASSERT(probe.GetY() == 1191);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 35);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 6);
	CPPUNIT_ASSERT(group.GetStop() == 8);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == DifferentRepType);

	group.GetCell(2, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 8);
	CPPUNIT_ASSERT(probe.GetExpos() == 18);
	CPPUNIT_ASSERT(probe.GetX() == 991);
	CPPUNIT_ASSERT(probe.GetY() == 1191);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 35);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
}

void CCDFFileDataTest::test_AxiomV4XDA()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(AXIOM_V4_XDA_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );

	CPPUNIT_ASSERT(cdf.GetGUID() == "000030a3-4695-44df-5795-002ce7002203");
	CPPUNIT_ASSERT(cdf.GetIntegrityMd5() == "");
	CPPUNIT_ASSERT(cdf.GetChipType() == "Test_Axiom_GW_SNP");
	std::vector<std::string> chiptypes = cdf.GetChipTypes();
	CPPUNIT_ASSERT(chiptypes.at(0) == "Test_Axiom_GW_SNP");
	CPPUNIT_ASSERT(chiptypes.at(1) == "Test_Axiom_GW_SNP.v1");
	CPPUNIT_ASSERT(chiptypes.at(2) == "Test_Axiom_GW_SNP.r1");

	CCDFFileHeader &header = cdf.GetHeader();
   	CPPUNIT_ASSERT(header.GetCols() == 1190);
	CPPUNIT_ASSERT(header.GetRows() == 1190);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 11);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 3);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "AFFX-LCP-11715001");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "AFFX-LCP-11715002");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "AFFX-NP-11083590");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "AX-12510582");
	
	CPPUNIT_ASSERT(cdf.GetProbeSetType(0) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(1) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(2) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(3) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(header.GetNumProbeSets()-1) == MultichannelMarkerProbeSetType);	
	 
	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;

    //MultichannelProbeSetType //AX-11219884
    cdf.GetProbeSetInformation(5, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 774492 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 339);
	CPPUNIT_ASSERT(probe.GetY() == 584);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 997);
	CPPUNIT_ASSERT(probe.GetY() == 11);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 1);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 5);
	CPPUNIT_ASSERT(probe.GetY() == 1);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 2); 
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 4);
	CPPUNIT_ASSERT(group.GetStop() == 6);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 4);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 339);
	CPPUNIT_ASSERT(probe.GetY() == 583);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 4);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 5);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 997);
	CPPUNIT_ASSERT(probe.GetY() == 10);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 5);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 5);
	CPPUNIT_ASSERT(probe.GetY() == 0);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 6); 
  
	//MultichannelProbeSetType //AX-12398238
    cdf.GetProbeSetInformation(9, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 2 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 2 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 299927 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 1);
	CPPUNIT_ASSERT(group.GetNumCells() == 1);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 1);
	CPPUNIT_ASSERT(group.GetStop() == 1);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1056);
	CPPUNIT_ASSERT(probe.GetY() == 565);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 1);
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 1);
	CPPUNIT_ASSERT(group.GetNumCells() == 1);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 3);
	CPPUNIT_ASSERT(group.GetStop() == 3);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 3);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1056);
	CPPUNIT_ASSERT(probe.GetY() == 565);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 1);
}

void CCDFFileDataTest::test_AxiomV6ASCII()
{
	CCDFFileData cdf;
	cdf.SetFileName(GetDataPath(AXIOM_V6_ASCII_FILE).c_str());
	CPPUNIT_ASSERT(cdf.Read() == true );

	CPPUNIT_ASSERT(cdf.GetGUID() == "000030a3-4695-44df-5795-002ce7002203");
	CPPUNIT_ASSERT(cdf.GetIntegrityMd5() == "");
	CPPUNIT_ASSERT(cdf.GetChipType() == "Test_Axiom_GW_SNP");
	std::vector<std::string> chiptypes = cdf.GetChipTypes();
	CPPUNIT_ASSERT(chiptypes.at(0) == "Test_Axiom_GW_SNP");
	CPPUNIT_ASSERT(chiptypes.at(1) == "Test_Axiom_GW_SNP.v1");
	CPPUNIT_ASSERT(chiptypes.at(2) == "Test_Axiom_GW_SNP.r1");

	CCDFFileHeader &header = cdf.GetHeader();
   	CPPUNIT_ASSERT(header.GetCols() == 1190);
	CPPUNIT_ASSERT(header.GetRows() == 1190);
	CPPUNIT_ASSERT(header.GetNumProbeSets() == 11);
	CPPUNIT_ASSERT(header.GetNumQCProbeSets() == 3);
	CPPUNIT_ASSERT(header.GetReference() == "");

	CPPUNIT_ASSERT(cdf.GetProbeSetName(0) == "AFFX-LCP-11715001");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(1) == "AFFX-LCP-11715002");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(2) == "AFFX-NP-11083590");
	CPPUNIT_ASSERT(cdf.GetProbeSetName(header.GetNumProbeSets()-1) == "AX-12510582");
	
	CPPUNIT_ASSERT(cdf.GetProbeSetType(0) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(1) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(2) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(3) == MultichannelMarkerProbeSetType);
	CPPUNIT_ASSERT(cdf.GetProbeSetType(header.GetNumProbeSets()-1) == MultichannelMarkerProbeSetType);	
	 
	CCDFProbeSetInformation set;
	CCDFProbeGroupInformation group;
	CCDFProbeInformation probe;

    //MultichannelProbeSetType //AX-11219884
    cdf.GetProbeSetInformation(5, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 6 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 6 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 774492 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 0);
	CPPUNIT_ASSERT(group.GetStop() == 2);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 0);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 339);
	CPPUNIT_ASSERT(probe.GetY() == 584);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 0);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 997);
	CPPUNIT_ASSERT(probe.GetY() == 11);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 1);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == group.GetNumLists()-1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 5);
	CPPUNIT_ASSERT(probe.GetY() == 1);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'g');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'c');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 2); 
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 3);
	CPPUNIT_ASSERT(group.GetNumCells() == 3);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 4);
	CPPUNIT_ASSERT(group.GetStop() == 6);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
    CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 4);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 339);
	CPPUNIT_ASSERT(probe.GetY() == 583);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 4);
	group.GetCell(1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 5);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 997);
	CPPUNIT_ASSERT(probe.GetY() == 10);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 5);
	group.GetCell(group.GetNumCells()-1, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 6);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 5);
	CPPUNIT_ASSERT(probe.GetY() == 0);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'c');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 'g');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 6); 
  
	//MultichannelProbeSetType //AX-12398238
    cdf.GetProbeSetInformation(9, set);
	CPPUNIT_ASSERT(set.GetProbeSetType() == MultichannelMarkerProbeSetType );
	CPPUNIT_ASSERT(set.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(set.GetNumLists() == 2 );
	CPPUNIT_ASSERT(set.GetNumGroups() == 2 );
	CPPUNIT_ASSERT(set.GetNumCells() == 2 );
	CPPUNIT_ASSERT(set.GetNumCellsPerList() == 1 );
	CPPUNIT_ASSERT(set.GetProbeSetNumber() == 299927 );

	set.GetGroupInformation(0, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 1);
	CPPUNIT_ASSERT(group.GetNumCells() == 1);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 1);
	CPPUNIT_ASSERT(group.GetStop() == 1);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 0);
	CPPUNIT_ASSERT(group.GetChannel() == 1);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 1);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1056);
	CPPUNIT_ASSERT(probe.GetY() == 565);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 1);
    
	set.GetGroupInformation(1, group);
	CPPUNIT_ASSERT(group.GetDirection() == SenseDirection);
	CPPUNIT_ASSERT(group.GetNumLists() == 1);
	CPPUNIT_ASSERT(group.GetNumCells() == 1);
	CPPUNIT_ASSERT(group.GetNumCellsPerList() == 1);
	CPPUNIT_ASSERT(group.GetStart() == 3);
	CPPUNIT_ASSERT(group.GetStop() == 3);
	CPPUNIT_ASSERT(group.GetName() == "NONE");
	CPPUNIT_ASSERT(group.GetWobbleSituation() == 0);
	CPPUNIT_ASSERT(group.GetAlleleCode() == 1);
	CPPUNIT_ASSERT(group.GetChannel() == 0);
	CPPUNIT_ASSERT(group.GetRepType() == IdenticalRepType);

	group.GetCell(0, probe);
	CPPUNIT_ASSERT(probe.GetListIndex() == 3);
	CPPUNIT_ASSERT(probe.GetExpos() == 16);
	CPPUNIT_ASSERT(probe.GetX() == 1056);
	CPPUNIT_ASSERT(probe.GetY() == 565);
	CPPUNIT_ASSERT(tolower(probe.GetPBase()) == 'a');
	CPPUNIT_ASSERT(tolower(probe.GetTBase()) == 't');
	CPPUNIT_ASSERT(probe.GetProbeLength() == 30);
	CPPUNIT_ASSERT(probe.GetProbeGrouping() == 1);
}
