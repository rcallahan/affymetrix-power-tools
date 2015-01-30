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

/**
 * @file   CytogeneticstrioAnalysisTest.cpp
 *
 * @brief  Testing the util functions.
 *
 */

//
#include "copynumber/CytogeneticsTrioAnalysis.h"
#include "file/TsvFile/TsvFile.h"
#include "util/AffxByteArray.h"
#include "util/RegressionCheck.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
//
#include "util/CPPTest/Setup.h"

/**
 * @class CytogeneticsTrioAnalysisTest
 * @brief cppunit class for testing functions.
 */
class CytogeneticsTrioAnalysisTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( CytogeneticsTrioAnalysisTest );
  CPPUNIT_TEST( LOD_PaternityTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void LOD_PaternityTest();
};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( CytogeneticsTrioAnalysisTest );

void CytogeneticsTrioAnalysisTest::LOD_PaternityTest() {
	
	Verbose::out(1, "");
	Verbose::out(1, "CytogeneticsTrioAnalysisTest::LOD_PaternityTest()");
	AffxString strFileName = "input/LOD_PaternityTest.txt";
        affx::TsvFile tsv;
        tsv.m_optAutoTrim = true;
	unsigned int uiCount = 0;
	if (tsv.open(strFileName) == affx::TSV_OK) {
          while(tsv.nextLevel(0) == affx::TSV_OK) {
            uiCount++;
          }
          tsv.clear();
	} else {Err::errAbort("Cannot open file: " + strFileName);}
	std::vector<char> vSampleGenotypeCalls(uiCount);
	std::vector<char> vKnownParentGenotypeCalls(uiCount);
	std::vector<char> vAllegedParentGenotypeCalls(uiCount);
	std::vector<float> vProbabilityAAlleles(uiCount);
	std::vector<float> vProbabilityBAlleles(uiCount);
	uiCount = 0;
	if (tsv.open(strFileName) == affx::TSV_OK) {
          while(tsv.nextLevel(0) == affx::TSV_OK) {
            int colIdx = 1;
            std::string col;
            tsv.get(0,colIdx++, col);
            vSampleGenotypeCalls[uiCount] = (char)AffxByteArray(col).parseInt();
            tsv.get(0,colIdx++, col);
            vKnownParentGenotypeCalls[uiCount] = (char)AffxByteArray(col).parseInt();
            tsv.get(0,colIdx++, col);
            vAllegedParentGenotypeCalls[uiCount] = (char)AffxByteArray(col).parseInt();
            tsv.get(0,colIdx++, col);
            vProbabilityAAlleles[uiCount] = (float)AffxByteArray(col).parseDouble();
            tsv.get(0,colIdx++, col);
            vProbabilityBAlleles[uiCount] = (float)AffxByteArray(col).parseDouble();
            uiCount++;
          }
          tsv.clear();
	} else {Err::errAbort("Cannot open file: " + strFileName);}
	
	time_t startTime = time(NULL);	
	CytogeneticsTrioAnalysis obj;
	double dLOD = obj.calculateLODPaternity(vSampleGenotypeCalls, vKnownParentGenotypeCalls, vAllegedParentGenotypeCalls, vProbabilityAAlleles, vProbabilityBAlleles, 0.01);
	Verbose::out(1, "LODPaternity = " + ::getDouble(dLOD) + "\tExpected: 1602.1812311429");
	CPPUNIT_ASSERT(fabs(1602.1812311429 - dLOD) <= 0.01);
	
	dLOD = obj.calculateLODPaternity(vSampleGenotypeCalls, vKnownParentGenotypeCalls, vProbabilityAAlleles, vProbabilityBAlleles, 0.01);
	Verbose::out(1, "LODPaternity(Single Known) = " + ::getDouble(dLOD) + "\tExpected: 935.65432393134");
	CPPUNIT_ASSERT(fabs(935.65432393134 - dLOD) <= 0.01);
	
	dLOD = obj.calculateLODPaternity(vSampleGenotypeCalls, vAllegedParentGenotypeCalls, vProbabilityAAlleles, vProbabilityBAlleles, 0.01);
	Verbose::out(1, "LODPaternity(Single Alleged) = " + ::getDouble(dLOD) + "\tExpected: 953.69509137911");;
	CPPUNIT_ASSERT(fabs(953.69509137911 - dLOD) <= 0.01);

	time_t endTime = time(NULL);
	int t = int(  (float)(endTime - startTime) * 100); // convert to minutes
	Verbose::out(1, ToStr("Run took approximately: ") + ToStr((float)t/100) + ToStr(" seconds."));
	Verbose::out(1, "");
}
