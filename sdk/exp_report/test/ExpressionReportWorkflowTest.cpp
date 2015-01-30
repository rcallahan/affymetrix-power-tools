////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#include "exp_report/workflow/ExpressionReportWorkflow.h"
//
#include "calvin_files/parsers/src/CHPFileReader.h"
#include "calvin_files/parsers/src/CHPQuantificationFileReader.h"
#include "file/CHPFileData.h"
//
#include <cppunit/extensions/HelperMacros.h>
//
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
//

using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;
using namespace ExpressionReport;
using namespace affxchp;

#ifdef _MSC_VER
#pragma warning(disable: 4996)
#endif

class ExpressionReportWorkflowTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( ExpressionReportWorkflowTest );
    ///@todo exon data missing from CVS
	//CPPUNIT_TEST( testExonChpFile );
	CPPUNIT_TEST( testQuantificationChpFile );
	CPPUNIT_TEST( testCCLegChpFile );
	CPPUNIT_TEST( testXDAChpFile );
	CPPUNIT_TEST( testFail );
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
	void testQuantificationChpFile();
    void testCCLegChpFile();
    void testXDAChpFile();
    void testExonChpFile();
    void testFail();
};


CPPUNIT_TEST_SUITE_REGISTRATION( ExpressionReportWorkflowTest );

void ExpressionReportWorkflowTest::setUp()
{
}

void ExpressionReportWorkflowTest::tearDown()
{
}

void copy(const char *f1, const char *f2)
{
    FILE *infp = fopen(f1, "rb");
    FILE *outfp = fopen(f2, "wb");
    int c;
    while ((c = getc(infp)) != EOF)
    {
        putc(c, outfp);
    }
    fclose(outfp);
    fclose(infp);
}

void ExpressionReportWorkflowTest::testFail()
{
	ExpressionReportWorkflow wrk;
    copy("../data/xda-mas5.chp", "../data/test.chp");
    affymetrix_calvin_parameter::ParameterNameValueTypeList extra;
    CPPUNIT_ASSERT(wrk.Run("../data/nofile.chp", "../data/nofile.cel", "../data", "../data/Hu6800.default.report_controls", NULL, extra, false) == false);
}

void ExpressionReportWorkflowTest::testXDAChpFile()
{
	ExpressionReportWorkflow wrk;

    copy("../data/xda-mas5.chp", "../data/test.chp");
    affymetrix_calvin_parameter::ParameterNameValueTypeList extra;
    CPPUNIT_ASSERT(wrk.Run("../data/test.chp", "../data/test.CEL", "../data", "../data/Hu6800.default.report_controls", NULL, extra, false) == true);


    CCHPFileData chp1;
    CCHPFileData chp2;

    chp1.SetFileName("../data/test.chp");
    chp2.SetFileName("../data/xda-mas5-ref.chp");
    CPPUNIT_ASSERT(chp1.Read() == true);
    CPPUNIT_ASSERT(chp2.Read() == true);

    CCHPFileHeader &h1 = chp1.GetHeader();
    CCHPFileHeader &h2 = chp2.GetHeader();

    CPPUNIT_ASSERT(h1.GetAlgName() == h2.GetAlgName());
    CPPUNIT_ASSERT(h1.GetAlgVersion() == h2.GetAlgVersion());
    CPPUNIT_ASSERT(h1.GetParentCellFile() == h2.GetParentCellFile());
    CPPUNIT_ASSERT(h1.GetProgID() == h2.GetProgID());
    CPPUNIT_ASSERT(h1.GetAssayType() == h1.GetAssayType());
    CPPUNIT_ASSERT(h1.GetChipType() == h2.GetChipType());
    CPPUNIT_ASSERT(h1.GetCols() == h2.GetCols());
    CPPUNIT_ASSERT(h1.GetRows() == h2.GetRows());
    CPPUNIT_ASSERT(h1.GetNumProbeSets() == h2.GetNumProbeSets());

    TagValuePairTypeList::iterator it1=h1.AlgorithmParameters().begin();
    TagValuePairTypeList::iterator it2=h2.AlgorithmParameters().begin();
    CPPUNIT_ASSERT(h1.AlgorithmParameters().size() == h2.AlgorithmParameters().size());
    while (it1 != h1.AlgorithmParameters().end())
    {
        CPPUNIT_ASSERT(it1->Tag == it2->Tag);
        CPPUNIT_ASSERT(it1->Value == it2->Value);
        ++it1;
        ++it2;
    }

    it1=h1.SummaryParameters().begin();
    it2=h2.SummaryParameters().begin();
    CPPUNIT_ASSERT(h1.SummaryParameters().size() == h2.SummaryParameters().size());

    std::map<std::string, std::string> h1sum;
    while (it1 != h1.SummaryParameters().end())
    {
        h1sum[it1->Tag] = it1->Value;
        ++it1;
    }
    std::map<std::string, std::string> h2sum;
    while (it2 != h2.SummaryParameters().end())
    {
        h1sum[it2->Tag] = it2->Value;
        ++it2;
    }

    std::map<std::string, std::string>::iterator mit = h1sum.begin();
    while (mit != h1sum.end())
    {
        CPPUNIT_ASSERT( h2sum[ mit->second ] == h1sum[ mit->second ]);
        ++mit;
    }

    /*
    while (it1 != h1.SummaryParameters().end())
    {
        CPPUNIT_ASSERT(it1->Tag == it2->Tag);
        CPPUNIT_ASSERT(it1->Value == it2->Value);
        ++it1;
        ++it2;
    }
    */


    int nsets = h1.GetNumProbeSets();
    CExpressionProbeSetResults *e1;
    CExpressionProbeSetResults *e2;
	for (int iset=0; iset<nsets; iset++)
	{
        e1 = chp1.GetExpressionResults(iset);
        e2 = chp2.GetExpressionResults(iset);
        CPPUNIT_ASSERT(e1->Detection == e2->Detection);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->Signal, e2->Signal, 0.0001f);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->DetectionPValue, e2->DetectionPValue, 0.0001f);
        CPPUNIT_ASSERT(e1->NumPairs == e2->NumPairs);
        CPPUNIT_ASSERT(e1->NumUsedPairs== e2->NumUsedPairs);
        CPPUNIT_ASSERT(e1->m_HasCompResults == e2->m_HasCompResults);
    }

    BackgroundZoneInfo &b1 = h1.GetBackgroundZoneInfo();
    BackgroundZoneInfo &b2 = h2.GetBackgroundZoneInfo();
    CPPUNIT_ASSERT(b1.number_zones == b2.number_zones);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(b1.smooth_factor, b2.smooth_factor, 0.0001f);

    BackgroundZoneTypeList::iterator z1it = b1.zones.begin();
    BackgroundZoneTypeList::iterator z2it = b2.zones.begin();

    while (z1it != b1.zones.end())
    {
        CPPUNIT_ASSERT(z1it->centerx == z2it->centerx);
        CPPUNIT_ASSERT(z1it->centery == z2it->centery);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(z1it->background, z2it->background, 0.0001f);
        ++z1it;
        ++z2it;
    }
   

    remove("../data/test.chp");
}

void ExpressionReportWorkflowTest::testCCLegChpFile()
{
	ExpressionReportWorkflow wrk;

    copy("../data/ec-mas5.chp", "../data/test.chp");
    affymetrix_calvin_parameter::ParameterNameValueTypeList extra;
    CPPUNIT_ASSERT(wrk.Run("../data/test.chp", "../data/test.CEL", "../data", "../data/Hu6800.default.report_controls", NULL, extra, false) == true);

    CHPData *chp1 = new CHPData;
    CHPData *chp2 = new CHPData;
    //try
    //{
        // Read the input file.
        CHPFileReader reader;
        reader.SetFilename("../data/test.chp");
        reader.Read(*chp1);
        reader.SetFilename("../data/ec-mas5-ref.chp");
        reader.Read(*chp2);

        // Check the results.
        CPPUNIT_ASSERT(
            chp1->GetGenericData().Header().GetGenericDataHdr()->GetParentCnt() ==
            chp2->GetGenericData().Header().GetGenericDataHdr()->GetParentCnt());
        CPPUNIT_ASSERT(chp1->GetEntryCount() == chp2->GetEntryCount());
        CPPUNIT_ASSERT(chp1->GetAlgName() == chp2->GetAlgName());
        CPPUNIT_ASSERT(chp1->GetAlgVersion() == chp2->GetAlgVersion());
        CPPUNIT_ASSERT(chp1->GetArrayType() == chp2->GetArrayType());

        ParameterNameValueTypeVector l1 = chp1->GetAlgParams();
        ParameterNameValueTypeVector l2 = chp2->GetAlgParams();
        CPPUNIT_ASSERT(l1.size() == l2.size());
        ParameterNameValueTypeVector::iterator it1 = l1.begin();
        ParameterNameValueTypeVector::iterator it2 = l2.begin();
        while (it1 != l1.end() && it2 != l2.end())
        {
            CPPUNIT_ASSERT(it1->GetName() == it2->GetName());
            CPPUNIT_ASSERT(it1->ToString() == it2->ToString());
            ++it1;
            ++it2;
        }

        l1 = chp1->GetChipSums();
        l2 = chp2->GetChipSums();
        CPPUNIT_ASSERT(l1.size() == l2.size());

        it1 = l1.begin();
        std::map<std::wstring, std::wstring> h1sum;
        while (it1 != l1.end())
        {
            h1sum[it1->GetName()] = it1->ToString();
            ++it1;
        }
        it2 = l2.begin();
        std::map<std::wstring, std::wstring> h2sum;
        while (it2 != l2.end())
        {
            h2sum[it2->GetName()] = it2->ToString();
            ++it2;
        }

        std::map<std::wstring, std::wstring>::iterator mit = h2sum.begin();
        while (mit != h2sum.end())
        {
            CPPUNIT_ASSERT( h2sum[ mit->second ] == h1sum[ mit->second ]);
            ++mit;
        }

       
        int nsets = chp1->GetEntryCount();
        CHPExpressionEntry e1;
        CHPExpressionEntry e2;
		for (int iset=0; iset<nsets; iset++)
		{
            chp1->GetEntry(iset, e1);
            chp2->GetEntry(iset, e2);
            CPPUNIT_ASSERT(e1.GetDetection() == e2.GetDetection());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(e1.GetSignal(), e2.GetSignal(), 0.0001f);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(e1.GetDetectionPValue(), e2.GetDetectionPValue(), 0.0001f);
            CPPUNIT_ASSERT(e1.GetNumPairs() == e2.GetNumPairs());
            CPPUNIT_ASSERT(e1.GetNumPairsUsed() == e2.GetNumPairsUsed());
            CPPUNIT_ASSERT(e1.GetHasComparisonData() == e2.GetHasComparisonData());
        }

        int nz = chp1->GetBackgroundZoneCnt();
        CHPBackgroundZone z1;
        CHPBackgroundZone z2;
        for (int iz=0; iz<nz; iz++)
        {
            chp1->GetBackgroundZone(iz, z1);
            chp1->GetBackgroundZone(iz, z2);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(z1.GetBackground(), z2.GetBackground(), 0.0001f);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(z1.GetCenterX(), z2.GetCenterX(), 0.0001f);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(z1.GetCenterY(), z2.GetCenterY(), 0.0001f);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(z1.GetSmoothFactor(), z2.GetSmoothFactor(), 0.0001f);
        }
    //}
    //catch (...)
    //{
    //}

    chp2->Clear();
    delete chp2;
    chp2 = NULL;
            
    chp1->Clear();
    delete chp1;
    chp1 = NULL;


    remove("../data/test.chp");
}

void ExpressionReportWorkflowTest::testExonChpFile()
{
	ExpressionReportWorkflow wrk;

    copy("../data/exon.chp", "../data/test.chp");
    affymetrix_calvin_parameter::ParameterNameValueTypeList extra;
    CPPUNIT_ASSERT(wrk.Run("../data/test.chp", "../data/test.CEL", "../data", NULL, "../data/exon_report_probesets", extra, false) == true);

}

void ExpressionReportWorkflowTest::testQuantificationChpFile()
{
	ExpressionReportWorkflow wrk;

    copy("../data/ec-sig.chp", "../data/test.chp");
    affymetrix_calvin_parameter::ParameterNameValueTypeList extra;
    CPPUNIT_ASSERT(wrk.Run("../data/test.chp", "../data/test.CEL", "../data", "../data/Hu6800.default.report_controls", NULL, extra, false) == true);

    CHPQuantificationData *chp1 = new CHPQuantificationData("../data/test.chp");
    CHPQuantificationData *chp2 = new CHPQuantificationData("../data/ec-sig-ref.chp");
    //try
    //{
        // Read the input file.
        CHPQuantificationFileReader reader;
        reader.Read(*chp1);
        reader.Read(*chp2);

        // Check the results.
        CPPUNIT_ASSERT(
            chp1->GetGenericData().Header().GetGenericDataHdr()->GetParentCnt() ==
            chp2->GetGenericData().Header().GetGenericDataHdr()->GetParentCnt());
        CPPUNIT_ASSERT(chp1->GetEntryCount() == chp2->GetEntryCount());
        CPPUNIT_ASSERT(chp1->GetAlgName() == chp2->GetAlgName());
        CPPUNIT_ASSERT(chp1->GetAlgVersion() == chp2->GetAlgVersion());
        CPPUNIT_ASSERT(chp1->GetArrayType() == chp2->GetArrayType());

        ParameterNameValueTypeList l1 = chp1->GetAlgParams();
        ParameterNameValueTypeList l2 = chp2->GetAlgParams();
        CPPUNIT_ASSERT(l1.size() == l2.size());
        ParameterNameValueTypeList::iterator it1 = l1.begin();
        ParameterNameValueTypeList::iterator it2 = l2.begin();
        while (it1 != l1.end() && it2 != l2.end())
        {
            CPPUNIT_ASSERT(it1->GetName() == it2->GetName());
            CPPUNIT_ASSERT(it1->ToString() == it2->ToString());
            ++it1;
            ++it2;
        }

        l1 = chp1->GetSummaryParams();
        l2 = chp2->GetSummaryParams();
        CPPUNIT_ASSERT(l1.size() == l2.size());
        it1 = l1.begin();
        it2 = l2.begin();
        while (it1 != l1.end() && it2 != l2.end())
        {
            CPPUNIT_ASSERT(it1->GetName() == it2->GetName());
            CPPUNIT_ASSERT(it1->ToString() == it2->ToString());
            ++it1;
            ++it2;
        }
        
        int nsets = chp1->GetEntryCount();
        ProbeSetQuantificationData e1;
        ProbeSetQuantificationData e2;
		for (int iset=0; iset<nsets; iset++)
		{
            chp1->GetQuantificationEntry(iset, e1);
            chp2->GetQuantificationEntry(iset, e2);
            CPPUNIT_ASSERT(e1.name == e2.name);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(e1.quantification, e2.quantification, 0.0001f);

        }
    //}
    //catch (...)
    //{
    //}

    chp2->Clear();
    delete chp2;
    chp2 = NULL;
            
    chp1->Clear();
    delete chp1;
    chp1 = NULL;

    remove("../data/test.chp");

}
