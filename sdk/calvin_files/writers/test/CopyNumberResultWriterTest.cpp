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

#define _CRT_SECURE_NO_WARNINGS

#include "calvin_files/writers/test/CopyNumberResultWriterTest.h"
//
#include "calvin_files/fusion/src/FusionCHPData.h"
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
#include "calvin_files/writers/src/CopyNumberResultWriter.h"
//

#ifdef _MSC_VER
#pragma warning(disable: 4996) // don't show deprecated warnings.
#ifdef HAVE_SNPRINTF // If not using visual c++'s _snprintf include snprintf.
extern "C" {
#include "snprintf.h"
} 
#else // otherwise use _snprintf where normally use snprintf.
#define snprintf _snprintf
#endif // HAVE_SNPRINTF
#endif // _MSC_VER

using namespace std;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_exceptions;

CPPUNIT_TEST_SUITE_REGISTRATION( CopyNumberResultWriterTest );

void CopyNumberResultWriterTest::setUp()
{
}

void CopyNumberResultWriterTest::tearDown()
{
}

void CopyNumberResultWriterTest::testWrite()
{
	CopyNumberResultWriter writer;
    const char *fileNames[] = {"test1.cn", "test2.cn"};

    list<ParameterNameValueType> algParams;
    list<ParameterNameValueType> sumParams;
    vector<ColumnInfo> cols;

    writer.MaximumProbeSetNameLength(12);
    writer.MaximumCytoRegionNameLength(12);
	writer.MaximumGenotypeProbeSetNameLength(12);
    writer.AlgName() = "MYALG";
    writer.AlgVersion() = "1.0";
    writer.NumberProbeSets() = 10;
    writer.NumberCytoRegions() = 10;
	writer.NumberGenotypeProbeSets() = 10;
    writer.Columns() = cols;
    writer.AlgParams() = algParams;
    writer.SetChromosomeProbeSetIndexInformation(X_CHR, 0, writer.NumberProbeSets());
    for (int i=0; i<2; i++)
    {
        FusionCELData cel;
        cel.SetFileName("../data/small_cel_file");
        cel.Read(false);

        writer.SummaryParams() = sumParams;
        writer.CreateResultFile(cel, fileNames[i]);
        cel.Close();

        ProbeSetMultiDataCopyNumberData entry;
        char buf[64];
        for (int j=0; j<writer.NumberProbeSets(); j++)
        {
            entry.chr = X_CHR;
            entry.position = j+i;
            snprintf(buf, 64, "%d", j+i);
            entry.name = buf;
            writer.WriteProbeSetResult(entry);
        }
        ProbeSetMultiDataCytoRegionData cy;
        for (int j=0; j<writer.NumberCytoRegions(); j++)
        {
            cy.call = 1;
            cy.confidenceScore = (float)(j+i);
            snprintf(buf, 64, "%d", j+i);
            cy.name = buf;
            writer.WriteCytoRegionResult(cy);
        }
        ProbeSetMultiDataGenotypeData gt;
        for (int j=0; j<writer.NumberGenotypeProbeSets(); j++)
        {
            gt.call = 1;
			gt.confidence = (float)(j+i);
            snprintf(buf, 64, "%d", j+i);
            gt.name = buf;
			writer.WriteGenotypeProbeSetResult(gt);
        }
        writer.CloseResultsFile();
    }

    for (int i=0; i<2; i++)
    {
        FusionCHPData *chp = FusionCHPDataReg::Read(fileNames[i]);
	    CPPUNIT_ASSERT(chp != NULL);
	    FusionCHPMultiDataData *genoChp = FusionCHPMultiDataData::FromBase(chp); 
	    CPPUNIT_ASSERT(genoChp != NULL);

	    CPPUNIT_ASSERT(genoChp->GetAlgName() == L"MYALG");
	    CPPUNIT_ASSERT(genoChp->GetAlgVersion() == L"1.0");
	    CPPUNIT_ASSERT(genoChp->GetArrayType() == L"Hg-small");
	    CPPUNIT_ASSERT(genoChp->GetEntryCount(ExpressionMultiDataType) == 0);
	    CPPUNIT_ASSERT(genoChp->GetEntryCount(GenotypeMultiDataType) == 10);
	    CPPUNIT_ASSERT(genoChp->GetEntryCount(CopyNumberMultiDataType) == 10);
	    CPPUNIT_ASSERT(genoChp->GetEntryCount(CytoMultiDataType) == 10);

        DataSetHeader *dsh = genoChp->GetDataSetHeader(CopyNumberMultiDataType);
        CPPUNIT_ASSERT(dsh->GetNameValParamCnt() == 3);
        ParameterNameValueTypeConstIt begin;
        ParameterNameValueTypeConstIt end;
        ParameterNameValueTypeConstIt it;
        dsh->GetNameValIterators(begin, end);
        it = begin;
        CPPUNIT_ASSERT(it->GetValueInt32() == 0);
        ++it;
        CPPUNIT_ASSERT(it->GetValueInt32() == 10);
        ++it;
        CPPUNIT_ASSERT(it->GetValueAscii() == "X");
        
        ProbeSetMultiDataCopyNumberData entry;
        char buf[64];
        for (int j=0; j<10; j++)
        {
            genoChp->GetCopyNumberEntry(CopyNumberMultiDataType, j, entry);
            CPPUNIT_ASSERT(entry.chr == X_CHR);
            CPPUNIT_ASSERT(entry.position == j+i);
            snprintf(buf, 64, "%d", j+i);
            CPPUNIT_ASSERT(entry.name.compare(buf) == 0);
        }

        ProbeSetMultiDataCytoRegionData cy;
        for (int j=0; j<10; j++)
        {
            genoChp->GetCytoRegionEntry(CytoMultiDataType, j, cy);
            CPPUNIT_ASSERT(cy.call == 1);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(cy.confidenceScore, j+i, 0.0001f);
            snprintf(buf, 64, "%d", j+i);
            CPPUNIT_ASSERT(cy.name.compare(buf) == 0);
        }

        ProbeSetMultiDataGenotypeData gt;
        for (int j=0; j<10; j++)
        {
            genoChp->GetGenotypeEntry(GenotypeMultiDataType, j, gt);
            CPPUNIT_ASSERT(gt.call == 1);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(gt.confidence, j+i, 0.0001f);
            snprintf(buf, 64, "%d", j+i);
            CPPUNIT_ASSERT(gt.name.compare(buf) == 0);
        }

        delete chp;
    }
}
