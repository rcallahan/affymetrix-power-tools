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
#include "calvin_files/fusion/test/FusionCDFQCProbeSetNamesTest.h"
//
#include "calvin_files/fusion/src/FusionCDFQCProbeSetNames.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affxcdf;
using namespace affymetrix_fusion_io;

CPPUNIT_TEST_SUITE_REGISTRATION( FusionCDFQCProbeSetNamesTest );

void FusionCDFQCProbeSetNamesTest::setUp()
{}

void FusionCDFQCProbeSetNamesTest::tearDown()
{}

void FusionCDFQCProbeSetNamesTest::testGetStaticCDFQCProbeSetName()
{
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(UnknownQCProbeSetType) == L"UnknownQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CheckerboardNegativeQCProbeSetType) == L"CheckerboardNegativeQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CheckerboardPositiveQCProbeSetType) == L"CheckerboardPositiveQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(HybNegativeQCProbeSetType) == L"HybNegativeQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(HybPositiveQCProbeSetType) == L"HybPositiveQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(TextFeaturesNegativeQCProbeSetType) == L"TextFeaturesNegativeQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(TextFeaturesPositiveQCProbeSetType) == L"TextFeaturesPositiveQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CentralNegativeQCProbeSetType) == L"CentralNegativeQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CentralPositiveQCProbeSetType) == L"CentralPositiveQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(GeneExpNegativeQCProbeSetType) == L"GeneExpNegativeQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(GeneExpPositiveQCProbeSetType) == L"GeneExpPositiveQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CycleFidelityNegativeQCProbeSetType) == L"CycleFidelityNegativeQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CycleFidelityPositiveQCProbeSetType) == L"CycleFidelityPositiveQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CentralCrossNegativeQCProbeSetType) == L"CentralCrossNegativeQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CentralCrossPositiveQCProbeSetType) == L"CentralCrossPositiveQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CrossHybNegativeQCProbeSetType) == L"CrossHybNegativeQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(CrossHybPositiveQCProbeSetType) == L"CrossHybPositiveQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(SpatialNormalizationNegativeQCProbeSetType) == L"SpatialNormalizationNegativeQC");
	CPPUNIT_ASSERT(FusionCDFQCProbeSetNames::GetStaticCDFQCProbeSetName(SpatialNormalizationPositiveQCProbeSetType) == L"SpatialNormalizationPositiveQC");
}
