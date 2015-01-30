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
#include "calvin_files/fusion/test/FusionCdfFileTest.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affxcdf;
using namespace affymetrix_fusion_io;

const string EXP_XDA_FILE		= "../data/CDF/Test3-xda.CDF";
const string EXP_ASCII_FILE		= "../data/CDF/Test3-ascii.CDF";
const string NO_FILE			= "../data/CDF/NoFile.CDF";
const string MAP_XDA_FILE		= "../data/CDF/Mapping10K_Xba131-xda.CDF";
const string MAP_ASCII_FILE		= "../data/CDF/Mapping10K_Xba131-ascii.CDF";
const string DMET3_XDA_FILE     = "../data/CDF/DMET3_Plus_test_xda.cdf";
const string DMET3_ASCII_FILE	= "../data/CDF/DMET3_Plus_test_ascii.cdf";
const string MULTICHANNEL_XDA_FILE     = "../data/CDF/Multichannel_test_xda.cdf";
const string MULTICHANNEL_ASCII_FILE	= "../data/CDF/Multichannel_test_ascii.cdf";
const string AXIOM_V4_XDA_FILE     = "../data/CDF/Axiom_v4_test_xda.r1.cdf";
const string AXIOM_V6_ASCII_FILE	= "../data/CDF/Axiom_v6_test_ascii.r1.cdf";




CPPUNIT_TEST_SUITE_REGISTRATION( FusionCdfFileTest );

void FusionCdfFileTest::setUp()
{}

void FusionCdfFileTest::tearDown()
{}

void FusionCdfFileTest::testCreation()
{
	FusionCDFData *gcosCdf = new FusionCDFData;
	CPPUNIT_ASSERT(gcosCdf != NULL);
	delete gcosCdf;
}

void FusionCdfFileTest::testQCObjects(affxcdf::CCDFFileData &gcosCdf, affymetrix_fusion_io::FusionCDFData &fusionCdf)
{
	affxcdf::CCDFQCProbeSetInformation gcosSet;
	affxcdf::CCDFQCProbeInformation gcosCel;

	affymetrix_fusion_io::FusionCDFQCProbeSetInformation fusionSet;
	affymetrix_fusion_io::FusionCDFQCProbeInformation fusionCel;

	CPPUNIT_ASSERT(gcosCdf.GetHeader().GetCols() == fusionCdf.GetHeader().GetCols());
	CPPUNIT_ASSERT(gcosCdf.GetHeader().GetRows() == fusionCdf.GetHeader().GetRows());
	CPPUNIT_ASSERT(gcosCdf.GetHeader().GetNumQCProbeSets() == fusionCdf.GetHeader().GetNumQCProbeSets());
	CPPUNIT_ASSERT(gcosCdf.GetChipType() == fusionCdf.GetChipType());
	CPPUNIT_ASSERT(gcosCdf.GetFileName() == fusionCdf.GetFileName());

	int nps = gcosCdf.GetHeader().GetNumQCProbeSets();
	for (int ips=0; ips<nps; ips++)
	{
		gcosCdf.GetQCProbeSetInformation(ips, gcosSet);
		fusionCdf.GetQCProbeSetInformation(ips, fusionSet);

		CPPUNIT_ASSERT(gcosSet.GetQCProbeSetType() == fusionSet.GetQCProbeSetType());
		CPPUNIT_ASSERT(gcosSet.GetNumCells() == fusionSet.GetNumCells());

		int nc = gcosSet.GetNumCells();
		for (int ic=0; ic<nc; ic++)
		{
			gcosSet.GetProbeInformation(ic, gcosCel);
			fusionSet.GetProbeInformation(ic, fusionCel);

			CPPUNIT_ASSERT(gcosCel.GetX() == fusionCel.GetX());
			CPPUNIT_ASSERT(gcosCel.GetY() == fusionCel.GetY());
			CPPUNIT_ASSERT(gcosCel.GetPLen() == fusionCel.GetPLen());
			CPPUNIT_ASSERT(gcosCel.IsPerfectMatchProbe() == fusionCel.IsPerfectMatchProbe());
			CPPUNIT_ASSERT(gcosCel.IsBackgroundProbe() == fusionCel.IsBackgroundProbe());
		}
	}
}
void FusionCdfFileTest::testRegObjects(affxcdf::CCDFFileData &gcosCdf, affymetrix_fusion_io::FusionCDFData &fusionCdf)
{
	affxcdf::CCDFFileHeader &gcosHeader = gcosCdf.GetHeader();
	affxcdf::CCDFProbeSetInformation gcosSet;
	affxcdf::CCDFProbeGroupInformation gcosGroup;
	affxcdf::CCDFProbeInformation gcosCel;

	affymetrix_fusion_io::FusionCDFFileHeader &fusionHeader = fusionCdf.GetHeader();
	affymetrix_fusion_io::FusionCDFProbeSetInformation fusionSet;
	affymetrix_fusion_io::FusionCDFProbeGroupInformation fusionGroup;
	affymetrix_fusion_io::FusionCDFProbeInformation fusionCel;

	CPPUNIT_ASSERT(gcosHeader.GetCols() == fusionHeader.GetCols());
	CPPUNIT_ASSERT(gcosHeader.GetRows() == fusionHeader.GetRows());
	CPPUNIT_ASSERT(gcosHeader.GetNumProbeSets() == fusionHeader.GetNumProbeSets());
	CPPUNIT_ASSERT(gcosHeader.GetNumQCProbeSets() == fusionHeader.GetNumQCProbeSets());
	CPPUNIT_ASSERT(gcosHeader.GetReference() == fusionHeader.GetReference());
	CPPUNIT_ASSERT(gcosCdf.GetChipType() == fusionCdf.GetChipType());
	CPPUNIT_ASSERT(gcosCdf.GetFileName() == fusionCdf.GetFileName());

	int nps = gcosHeader.GetNumProbeSets();
	for (int ips=0; ips<nps; ips++)
	{
		CPPUNIT_ASSERT(gcosCdf.GetProbeSetName(ips) == fusionCdf.GetProbeSetName(ips));
		CPPUNIT_ASSERT(gcosCdf.GetProbeSetType(ips) == fusionCdf.GetProbeSetType(ips));
	}

	for (int ips=0; ips<nps; ips++)
	{
		gcosCdf.GetProbeSetInformation(ips, gcosSet);
		fusionCdf.GetProbeSetInformation(ips, fusionSet);

		CPPUNIT_ASSERT(gcosSet.GetNumLists() == fusionSet.GetNumLists());
		CPPUNIT_ASSERT(gcosSet.GetNumGroups() == fusionSet.GetNumGroups());
		CPPUNIT_ASSERT(gcosSet.GetNumCells() == fusionSet.GetNumCells());
		CPPUNIT_ASSERT(gcosSet.GetNumCellsPerList() == fusionSet.GetNumCellsPerList());
		CPPUNIT_ASSERT(gcosSet.GetProbeSetNumber() == fusionSet.GetProbeSetNumber());
		CPPUNIT_ASSERT(gcosSet.GetProbeSetType() == fusionSet.GetProbeSetType());
		CPPUNIT_ASSERT(gcosSet.GetDirection() == fusionSet.GetDirection());

		int ng = gcosSet.GetNumGroups();
		for (int ig=0; ig<ng; ig++)
		{
			gcosSet.GetGroupInformation(ig, gcosGroup);
			fusionSet.GetGroupInformation(ig, fusionGroup);

			CPPUNIT_ASSERT(gcosGroup.GetDirection() == fusionGroup.GetDirection());
			CPPUNIT_ASSERT(gcosGroup.GetNumLists() == fusionGroup.GetNumLists());
			CPPUNIT_ASSERT(gcosGroup.GetNumCells() == fusionGroup.GetNumCells());
			CPPUNIT_ASSERT(gcosGroup.GetNumCellsPerList() == fusionGroup.GetNumCellsPerList());
			CPPUNIT_ASSERT(gcosGroup.GetStart() == fusionGroup.GetStart());
			CPPUNIT_ASSERT(gcosGroup.GetStop() == fusionGroup.GetStop());
			CPPUNIT_ASSERT(gcosGroup.GetName() == fusionGroup.GetName());
			CPPUNIT_ASSERT(gcosGroup.GetAlleleCode() == fusionGroup.GetAlleleCode());
			CPPUNIT_ASSERT(gcosGroup.GetWobbleSituation() == fusionGroup.GetWobbleSituation());
			CPPUNIT_ASSERT(gcosGroup.GetChannel() == fusionGroup.GetChannel());
			CPPUNIT_ASSERT(gcosGroup.GetRepType() == fusionGroup.GetRepType());
			
			int nc=gcosGroup.GetNumCells();
			for (int ic=0; ic<nc; ic++)
			{
				gcosGroup.GetCell(ic, gcosCel);
				fusionGroup.GetCell(ic, fusionCel);

				CPPUNIT_ASSERT(gcosCel.GetX() == fusionCel.GetX());
				CPPUNIT_ASSERT(gcosCel.GetY() == fusionCel.GetY());
				CPPUNIT_ASSERT(gcosCel.GetListIndex() == fusionCel.GetListIndex());
				CPPUNIT_ASSERT(gcosCel.GetExpos() == fusionCel.GetExpos());
				CPPUNIT_ASSERT(gcosCel.GetPBase() == fusionCel.GetPBase());
				CPPUNIT_ASSERT(gcosCel.GetTBase() == fusionCel.GetTBase());
				CPPUNIT_ASSERT(gcosCel.GetProbeLength() == fusionCel.GetProbeLength());
				CPPUNIT_ASSERT(gcosCel.GetProbeGrouping() == fusionCel.GetProbeGrouping());
			}
		}
	}
}

void FusionCdfFileTest::testExpressionV3()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(EXP_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(EXP_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
}

void FusionCdfFileTest::testExpressionXDA()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(EXP_XDA_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(EXP_XDA_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
}

void FusionCdfFileTest::testMissingFile()
{
	affymetrix_fusion_io::FusionCDFData fusionCdf;
	fusionCdf.SetFileName(NO_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == false);
}

void FusionCdfFileTest::testMappingV3()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(MAP_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(MAP_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
}

void FusionCdfFileTest::testMappingXDA()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(MAP_XDA_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(MAP_XDA_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
}

void FusionCdfFileTest::test_DMET3ASCII()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(DMET3_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(DMET3_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
	
}

void FusionCdfFileTest::test_DMET3XDA()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(DMET3_XDA_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(DMET3_XDA_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
}


void FusionCdfFileTest::test_MultichannelASCII()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(MULTICHANNEL_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(MULTICHANNEL_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
	
}

void FusionCdfFileTest::test_MultichannelXDA()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(MULTICHANNEL_XDA_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(MULTICHANNEL_XDA_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
}

void FusionCdfFileTest::test_AxiomV6ASCII()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(AXIOM_V6_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(AXIOM_V6_ASCII_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
	
}

void FusionCdfFileTest::test_AxiomV4XDA()
{
	affxcdf::CCDFFileData gcosCdf;
	affymetrix_fusion_io::FusionCDFData fusionCdf;

	gcosCdf.SetFileName(AXIOM_V4_XDA_FILE.c_str());
	CPPUNIT_ASSERT(gcosCdf.Read() == true);
	fusionCdf.SetFileName(AXIOM_V4_XDA_FILE.c_str());
	CPPUNIT_ASSERT(fusionCdf.Read() == true);

	testQCObjects(gcosCdf, fusionCdf);
	testRegObjects(gcosCdf, fusionCdf);
}

