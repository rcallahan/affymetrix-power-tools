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


#include "calvin_files/fusion/test/FusionCHPDataGCOSTest.h"
//
#include "calvin_files/exception/src/DevelopmentException.h"
#include "calvin_files/fusion/src/FusionCHPLegacyData.h"
#include "calvin_files/parsers/src/FileException.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include "file/CHPFileData.h"	
//
#include <cstring>
#include <map>
#include <string>
//

using namespace std;
using namespace affxchp;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_utilities;

// Underlying CCelFileData class expects backslashes rather than forward slashes on Windows hmmmm
const string testFiles[] = 
{
	"doesnt_exist.CHP",														// file does not exist - for testing when a file does not exist
	"../../../TestDataFiles/CHP/GCOS104-U133a.CHP",
	"../../../TestDataFiles/CHP/pcwrs041a_pd01.CHP",
	"../../../TestDataFiles/CHP/NA02705_00243_01_3001197.CHP",
	"../../../TestDataFiles/CHP/10k_1.CHP",
	"../../../TestDataFiles/CHP/T1_r1_stat.CHP",
	"../../../TestDataFiles/CHP/T1_r1_v_T2_r1_stat.CHP",
	"../../../TestDataFiles/CHP/DCN_01_xda.CHP",
	"../../../TestDataFiles/CHP/GCOS104-U133a-test_rel.CHP",
	"../../../TestDataFiles/CHP/StatData_baseline.CHP",
	"../../../TestDataFiles/CHP/StatData.CHP",
	"../../../TestDataFiles/CHP/PLIERdata1.CHP"
};
const int FILE_NAME_LIST_COUNT = 11; // does not count the does not exist file case because we do not explicitly test it

const int DOES_NOT_EXIST = 0;	// index into the file name array above of the file does not exist name

CPPUNIT_TEST_SUITE_REGISTRATION(FusionCHPDataGCOSTest);

// has to be called before each test
void FusionCHPDataGCOSTest::setUp()
{
}

// has to be called after each test
void FusionCHPDataGCOSTest::tearDown()
{
}

void FusionCHPDataGCOSTest::TestFileId()
{
	if (FileUtils::Exists(testFiles[5].c_str()) == true)
	{
		FusionCHPData *chp = FusionCHPDataReg::ReadHeader(testFiles[5]);
		FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
		CPPUNIT_ASSERT(fusionChp->FileId() == "");
		CPPUNIT_ASSERT(fusionChp->GetGenericData() == NULL);
	}
}

// test to see if we can even creat the class
void FusionCHPDataGCOSTest::CreationTest(string fileName)
{
	//FusionCHPLegacyData * fusionChp = new FusionCHPLegacyData;
	//CPPUNIT_ASSERT_MESSAGE("Failed: Create fusion object!",fusionChp);
	//delete fusionChp;
}

// main test entry method, calls the rest.
// had to do it this way, because CPPUINT will not accept test method parameters
void FusionCHPDataGCOSTest::TestAll()
{
	// loop over each chp file to test
	for(int i = 1; i < FILE_NAME_LIST_COUNT; ++i)
	{
		string fileName = testFiles[i];
		if (FileUtils::Exists(fileName.c_str()) == false)
			continue;
		CreationTest(fileName);
		tearDown();
		setUp();
		ReadHeaderTest(fileName);
		tearDown();
		setUp();
		RowsColsTest(fileName);
		tearDown();
		setUp();
		VersionTest(fileName);
		tearDown();
		setUp();
		AlgorithmNameAndVersionTest(fileName);
		tearDown();
		setUp();
		ChipTypeTest(fileName);
		tearDown();
		setUp();
		ParentCelFileName(fileName);
		tearDown();
		setUp();
		ProgIDTest(fileName);
		tearDown();
		setUp();
		AlgParamsTest(fileName);
		tearDown();
		setUp();
		SummaryParamsTest(fileName);
		tearDown();
		setUp();
		NumberOfSetsTest(fileName);
		tearDown();
		setUp();
		ReadHeaderAndVerifyAlgParamsTest(fileName);
		tearDown();
		setUp();
		ReadHeaderAndVerifySummaryParamsTest(fileName);
		tearDown();
		setUp();
		ReadTest(fileName);
		tearDown();
		setUp();
		CallMethodsWhenObjectIsNotReadyTest(fileName);
		tearDown();
		setUp();
		DataCompareTest(fileName);
		tearDown();
		setUp();
		BackgroundDataTest(fileName);
	}
}

// test that we can read the chp header information
void FusionCHPDataGCOSTest::ReadHeaderTest(string fileName)
{
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);
	CCHPFileData gcosChp;

	gcosChp.SetFileName(fileName.c_str());
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	delete fusionChp;
}

// test the rows and colums
void FusionCHPDataGCOSTest::RowsColsTest(string fileName)
{
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # rows!",fusionChp->GetHeader().GetRows() == gcosChp.GetHeader().GetRows());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # cols!",fusionChp->GetHeader().GetCols() == gcosChp.GetHeader().GetCols());
	delete fusionChp;
}

// test the chp file version
void FusionCHPDataGCOSTest::VersionTest(string fileName)
{
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Version",fusionChp->GetHeader().GetVersion() == gcosChp.GetHeader().GetVersionNumber());
	delete fusionChp;
}

// test the algorithm and alg. verion info
void FusionCHPDataGCOSTest::AlgorithmNameAndVersionTest(string fileName)
{
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg name!", fusionChp->GetHeader().GetAlgName() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetAlgName()));
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg version!", fusionChp->GetHeader().GetAlgVersion() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetAlgVersion()));
	delete fusionChp;
}

// test the chiptype name
void FusionCHPDataGCOSTest::ChipTypeTest(string fileName)
{
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg version!", fusionChp->GetHeader().GetChipType() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetChipType()));
	delete fusionChp;
}

// test the parent cell file name
void FusionCHPDataGCOSTest::ParentCelFileName(string fileName)
{
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - parent cell file name!", fusionChp->GetHeader().GetParentCellFile() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetParentCellFile()));
	delete fusionChp;
}

// test the prog id
void FusionCHPDataGCOSTest::ProgIDTest(string fileName)
{
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - ProgID!", fusionChp->GetHeader().GetProgID() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetProgID()));
	delete fusionChp;
}

// test the algorithm params
void FusionCHPDataGCOSTest::AlgParamsTest(string fileName)
{
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

    FusionTagValuePairTypeList algParamsFusion;
    fusionChp->GetHeader().AlgorithmParameters(algParamsFusion);
    FusionTagValuePairTypeList::iterator beginFusion = algParamsFusion.begin();
    FusionTagValuePairTypeList::iterator endFusion = algParamsFusion.end();

	TagValuePairTypeList& algParamsGCOS = gcosChp.GetHeader().AlgorithmParameters();
	TagValuePairTypeList::iterator beginGCOS = algParamsGCOS.begin(), endGCOS = algParamsGCOS.end();

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg count",fusionChp->GetHeader().AlgorithmParameterCount() == gcosChp.GetHeader().AlgorithmParameters().size());

	typedef map<wstring,wstring,less<wstring> > theFusionMap;
	theFusionMap map;

	// need the map to have an easy lookup for the key to see if the fusion has the key
	for (; beginFusion != endFusion; beginFusion++)
	{
		map.insert(theFusionMap::value_type(beginFusion->Tag,beginFusion->Value));
	}

	theFusionMap::const_iterator findIt;
	for(; beginGCOS != endGCOS; beginGCOS++)
	{
		// test single alg value call
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare Algorithm Parameter!", fusionChp->GetHeader().GetAlgorithmParameter(StringUtils::ConvertMBSToWCS(beginGCOS->Tag).c_str()) == StringUtils::ConvertMBSToWCS(beginGCOS->Value));

		// test to see if the keys and values are the same between what was retrieved from fusion get algorithm params and
		// with what was retrieved by the gcos api
		findIt = map.find(StringUtils::ConvertMBSToWCS(beginGCOS->Tag));
		CPPUNIT_ASSERT_MESSAGE("Faile: Algorithm parameter not found in Fusion list", findIt != map.end());
		CPPUNIT_ASSERT_MESSAGE("Faile: Algorithm parameter not found in Fusion list", findIt->second == StringUtils::ConvertMBSToWCS(beginGCOS->Value));
	}
	delete fusionChp;
}

// tests the summary params
void FusionCHPDataGCOSTest::SummaryParamsTest(string fileName)
{
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

    FusionTagValuePairTypeList sumParamsFusion;
    fusionChp->GetHeader().SummaryParameters(sumParamsFusion);
    FusionTagValuePairTypeList::iterator beginFusion = sumParamsFusion.begin();
    FusionTagValuePairTypeList::iterator endFusion = sumParamsFusion.end();

	TagValuePairTypeList& sumParamsGCOS = gcosChp.GetHeader().SummaryParameters();
	TagValuePairTypeList::iterator beginGCOS = sumParamsGCOS.begin(), endGCOS = sumParamsGCOS.end();

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - sum count",fusionChp->GetHeader().SummaryParameterCount() == gcosChp.GetHeader().SummaryParameters().size());

	typedef map<wstring,wstring,less<wstring> > theFusionMap;
	theFusionMap map;

	// need the map to have an easy lookup for the key to see if the fusion has the key
	for (; beginFusion != endFusion; beginFusion++)
	{
		map.insert(theFusionMap::value_type(beginFusion->Tag,beginFusion->Value));
	}

	theFusionMap::const_iterator findIt;
	for(; beginGCOS != endGCOS; beginGCOS++)
	{
		// test single alg value call
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare Summary Parameter!", fusionChp->GetHeader().GetSummaryParameter(StringUtils::ConvertMBSToWCS(beginGCOS->Tag).c_str()) == StringUtils::ConvertMBSToWCS(beginGCOS->Value));

		// test to see if the keys and values are the same between what was retrieved from fusion get algorithm params and
		// with what was retrieved by the gcos api
		findIt = map.find(StringUtils::ConvertMBSToWCS(beginGCOS->Tag));
		CPPUNIT_ASSERT_MESSAGE("Faile: Algorithm parameter not found in Fusion list", findIt != map.end());
		CPPUNIT_ASSERT_MESSAGE("Faile: Algorithm parameter not found in Fusion list", findIt->second == StringUtils::ConvertMBSToWCS(beginGCOS->Value));
	}
	delete fusionChp;
}

// tests the number of probesets
void FusionCHPDataGCOSTest::NumberOfSetsTest(string fileName)
{
	// test the number of sets
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # probesets",fusionChp->GetHeader().GetNumProbeSets() == gcosChp.GetHeader().GetNumProbeSets());
	delete fusionChp;
}

// tests the algorithm parameters
void FusionCHPDataGCOSTest::ReadHeaderAndVerifyAlgParamsTest(string fileName)
{
	// read the header with the fusion api
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	// read the header with the gcos api
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);

	// get the gcos list
	TagValuePairTypeList listGCOS = gcosChp.GetHeader().AlgorithmParameters();
	FusionTagValuePairTypeList listFusion;
	// get the fusion list
	fusionChp->GetHeader().AlgorithmParameters(listFusion);

	// test that the number of parameters are the same between a pure gcos file sdk call and the fusdion api call
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg param count",listFusion.size() == listGCOS.size());

	// this section will test that the parameters read by both the pure gcos api and the fusion are the same
	// first we get the interators for the loop
	FusionTagValuePairTypeList::iterator begin = listFusion.begin();
	FusionTagValuePairTypeList::iterator end = listFusion.end();

	// loop over every parameter in the fusion parameter list and use it to get the value thru the getalgorithmparameter call
	// we could of used the value property of the of the iterator that is from the fusion algorithm list retrieved
	// early, but this is to exerise the get individual alg params api for both.
	// However, I will test the interator Value property from the fusion list and compare it to the other values returned from the
	// the get alg parameter calls
	for (; begin != end; begin++)
	{
		// get the fusion value for the parameter name
		wstring tempFusion = fusionChp->GetHeader().GetAlgorithmParameter(begin->Tag.c_str());
		// get the gcos value for the parameter name
		// have to convert the string to wstring first though
		wstring tempGCOS = StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetAlgorithmParameter(StringUtils::ConvertWCSToMBS(begin->Tag).c_str()));

		// double check that the value from the fusion list is the same as the other values retieved from the get alg params calls
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params",tempFusion == tempGCOS && begin->Value == tempGCOS);
	}
	delete fusionChp;
}

// tests the summary parameters
void FusionCHPDataGCOSTest::ReadHeaderAndVerifySummaryParamsTest(string fileName)
{
	// read the header with the fusion api
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	// read the header with the gcos api
	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);

	// get the gcos list
	TagValuePairTypeList listGCOS = gcosChp.GetHeader().SummaryParameters();
	FusionTagValuePairTypeList listFusion;
	// get the fusion list
	fusionChp->GetHeader().SummaryParameters(listFusion);

	// test that the number of parameters are the same between a pure gcos file sdk call and the fusdion api call
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - param list count",listFusion.size() == listGCOS.size());

	// this section will test that the parameters read by both the pure gcos api and the fusion are the same
	// first we get the interators for the loop
	FusionTagValuePairTypeList::iterator begin = listFusion.begin();
	FusionTagValuePairTypeList::iterator end = listFusion.end();

	// loop over every parameter in the fusion parameter list and use it to get the value thru the getsummaryparameter call
	// we could of used the value property of the of the iterator that is from the fusion summary list retrieved
	// early, but this is to exerise the get individual summary params api for both.
	// However, I will test the interator Value property from the fusion list and compare it to the other values returned from the
	// the get summary parameter calls
	for (; begin != end; begin++)
	{
		// get the fusion value for the parameter name
		wstring tempFusion = fusionChp->GetHeader().GetSummaryParameter(begin->Tag.c_str());
		// get the gcos value for the parameter name
		// have to convert the string to wstring first though
		wstring tempGCOS = StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetSummaryParameter(StringUtils::ConvertWCSToMBS(begin->Tag).c_str()));

		// double check that the value from the fusion list is the same as the other values retieved from the get summary params calls
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params",tempFusion == tempGCOS && begin->Value == tempGCOS);
	}
	delete fusionChp;
}

// test the read method of fusion. This differs only from previous readheader tests in that instead of readheader being called
// read is called.
void FusionCHPDataGCOSTest::ReadTest(string fileName)
{
	FusionCHPData *chp = FusionCHPDataReg::Read(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);

	// Check some header values
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - version",fusionChp->GetHeader().GetVersion() == gcosChp.GetHeader().GetVersionNumber());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - cols",fusionChp->GetHeader().GetCols() == gcosChp.GetHeader().GetCols());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - rows",fusionChp->GetHeader().GetRows() == gcosChp.GetHeader().GetRows());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # probesets",fusionChp->GetHeader().GetNumProbeSets() == gcosChp.GetHeader().GetNumProbeSets());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg name",fusionChp->GetHeader().GetAlgName() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetAlgName()));
    CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg version",fusionChp->GetHeader().GetAlgVersion() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetAlgVersion()));
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - chip type",fusionChp->GetHeader().GetChipType() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetChipType()));
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - parent cell name",fusionChp->GetHeader().GetParentCellFile() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetParentCellFile()));
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - prog id",fusionChp->GetHeader().GetProgID() == StringUtils::ConvertMBSToWCS(gcosChp.GetHeader().GetProgID()));
	delete fusionChp;
}

void FusionCHPDataGCOSTest::CallMethodsWhenObjectIsNotReadyTest(string fileName)
{
   // the idea here is that without performing the set file name and without calling the readxxx methods
    // any of the following calls shall throw an exception.
    
 //////////   CPPUNIT_ASSERT_THROW(fusionChp->ReadHeader(), FileNotOpenException);

	//////////CPPUNIT_ASSERT_THROW(fusionChp->GetHeader(), FileNotOpenException);
	//////////CPPUNIT_ASSERT_THROW(fusionChp->GetHeader().GetRows(), FileNotOpenException);
	//////////CPPUNIT_ASSERT_THROW(fusionChp->GetHeader().GetCols(), FileNotOpenException);
	//////////CPPUNIT_ASSERT_THROW(fusionChp->GetHeader().GetRows(), FileNotOpenException);
}

// tests the main data of the chp file
// currently supports expression, genotyping, universal, and soon resequencing
void FusionCHPDataGCOSTest::DataCompareTest(string fileName)
{
	// here we compare the data contents with both the gcos reader and the fusion reader and they better be the same
	FusionCHPData *chp = FusionCHPDataReg::Read(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	CPPUNIT_ASSERT_MESSAGE("Failed: Read!", gcosChp.Read());

	CExpressionProbeSetResults resultsGCOSExpression;
	FusionExpressionProbeSetResults resultsFusionExpression;

	CGenotypeProbeSetResults resultsGCOSGeno;
	FusionGenotypeProbeSetResults resultsFusionGeno;

	CUniversalProbeSetResults resultsGCOSUni;
	FusionUniversalProbeSetResults resultsFusionUni;

	u_int32_t count = fusionChp->GetHeader().GetNumProbeSets(); 
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # probesets",count == (u_int32_t)gcosChp.GetHeader().GetNumProbeSets());


	if(fusionChp->GetHeader().GetAssayType() == FusionExpression)
	{
		// for each type supported, get the probeset results
		// loop over the number of probesets
		for(u_int32_t i = 0; i < count; ++i)
		{
			CPPUNIT_ASSERT(fusionChp->GetExpressionResults(i, resultsFusionExpression));
			resultsGCOSExpression = *gcosChp.GetExpressionResults((int)i);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - DetectionPValue",resultsFusionExpression.GetDetectionPValue() == resultsGCOSExpression.DetectionPValue);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Signal",resultsFusionExpression.GetSignal() == resultsGCOSExpression.Signal);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Signal",resultsFusionExpression.GetNumPairs() == resultsGCOSExpression.NumPairs);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - NumPairs",resultsFusionExpression.GetNumUsedPairs() == resultsGCOSExpression.NumUsedPairs);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Detection",resultsFusionExpression.GetDetection() == resultsGCOSExpression.Detection);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - m_HasCompResults",resultsFusionExpression.HasCompResults() == resultsGCOSExpression.m_HasCompResults);
			if (resultsFusionExpression.HasCompResults())
			{
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - ChangePValue",resultsFusionExpression.GetChangePValue() == resultsGCOSExpression.ChangePValue);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - SignalLogRatio",resultsFusionExpression.GetSignalLogRatio() == resultsGCOSExpression.SignalLogRatio);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - SignalLogRatioLow",resultsFusionExpression.GetSignalLogRatioLow() == resultsGCOSExpression.SignalLogRatioLow);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - SignalLogRatioHigh",resultsFusionExpression.GetSignalLogRatioHigh() == resultsGCOSExpression.SignalLogRatioHigh);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - NumCommonPairs",resultsFusionExpression.GetNumCommonPairs() == resultsGCOSExpression.NumCommonPairs);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Change",resultsFusionExpression.GetChange() == resultsGCOSExpression.Change);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - GetDetectionString",resultsFusionExpression.GetDetectionString() == resultsGCOSExpression.GetDetectionString());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - GetChangeString",resultsFusionExpression.GetChangeString() == resultsGCOSExpression.GetChangeString());
			}
		}
	}
	else if(fusionChp->GetHeader().GetAssayType() == FusionGenotyping)
	{
		// for genotyping, we have a special case
		// for 10k we only check curtain parameters are checked and in 100K other are checked.
		// because we are not guaranteed to have properly initialized data from the files, we 
		// have to make checks against the algorithm names.
		// there are currently two.
		// if there are any algs added in the future, we will need to update this else if to include it.
		bool bWholeGenenome = false, bDynamicModel = false;

		if(fusionChp->GetHeader().GetAlgName() == L"WholeGenome")
		{
			bWholeGenenome = true;
		}
		else if(fusionChp->GetHeader().GetAlgName() == L"DynamicModel")
		{
			bDynamicModel = true;
		}

		// for each type supported, get the probeset results
		// loop over the number of probesets
		for(u_int32_t i = 0; i < count; ++i)
		{
			CPPUNIT_ASSERT(fusionChp->GetGenotypingResults(i, resultsFusionGeno));
			resultsGCOSGeno = *gcosChp.GetGenotypingResults((int)i);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - AlleleCall",resultsFusionGeno.GetAlleleCall() == resultsGCOSGeno.AlleleCall);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Confidence",resultsFusionGeno.GetConfidence() == resultsGCOSGeno.Confidence);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - AlleleCallString",resultsFusionGeno.GetAlleleCallString() == resultsGCOSGeno.GetAlleleCallString());

			if(bWholeGenenome)
			{
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - RAS1",resultsFusionGeno.GetRAS1() == resultsGCOSGeno.RAS1);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - RAS2",resultsFusionGeno.GetRAS2() == resultsGCOSGeno.RAS2);
			}
			else if(bDynamicModel)
			{
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - pvalue_AA",resultsFusionGeno.GetPValueAA() == resultsGCOSGeno.pvalue_AA);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - pvalue_AB",resultsFusionGeno.GetPValueAB() == resultsGCOSGeno.pvalue_AB);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - pvalue_BB",resultsFusionGeno.GetPValueBB() == resultsGCOSGeno.pvalue_BB);
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - pvalue_NoCall",resultsFusionGeno.GetPValueNoCall() == resultsGCOSGeno.pvalue_NoCall);
			}
		}
	}
	else if(fusionChp->GetHeader().GetAssayType() == FusionUniversal)
	{
		// for each type supported, get the probeset results
		// loop over the number of probesets
		for(u_int32_t i = 0; i < count; ++i)
		{
			CPPUNIT_ASSERT(fusionChp->GetUniversalResults(i, resultsFusionUni));
			resultsGCOSUni = *gcosChp.GetUniversalResults((int)i);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Background",resultsFusionUni.GetBackground() == resultsGCOSUni.GetBackground());
		}
	}
	else if(fusionChp->GetHeader().GetAssayType() == FusionResequencing)
	{
		FusionResequencingResults fusionResults;
		fusionChp->GetReseqResults(fusionResults);
		CResequencingResults *gcosResults = gcosChp.GetResequencingResults();

		CPPUNIT_ASSERT(fusionResults.GetCalledBasesSize() == gcosResults->GetCalledBasesSize());
		CPPUNIT_ASSERT(fusionResults.GetScoresSize() == gcosResults->GetScoresSize());
		CPPUNIT_ASSERT(fusionResults.GetForceCallsSize() == gcosResults->GetForceCallsSize());
		CPPUNIT_ASSERT(fusionResults.GetOrigCallsSize() == gcosResults->GetOrigCallsSize());

		const Int8Vector &fusionCalls = fusionResults.GetCalledBases();
		const FloatVector &fusionScores = fusionResults.GetScores();
		const std::vector<char> &gcosCalls = gcosResults->GetCalledBases();
		const FloatVector &gcosScores = gcosResults->GetScores();

		int sz = fusionResults.GetCalledBasesSize();
		for (int i=0; i<sz; i++)
		{
			CPPUNIT_ASSERT(fusionResults.GetCalledBase(i) == gcosResults->GetCalledBase(i));
			CPPUNIT_ASSERT(fusionResults.GetScore(i) == gcosResults->GetScore(i));
			CPPUNIT_ASSERT(fusionCalls[i] == gcosCalls[i]);
			CPPUNIT_ASSERT(fusionScores[i] == gcosScores[i]);
		}

		ForceCallType gcosForce;
		FusionForceCallType fusionForce;
		const FusionForceCallVector& fusionForces = fusionResults.GetForceCalls();
		const std::vector<ForceCallType>& gcosForces = gcosResults->GetForceCalls();

		sz = fusionResults.GetForceCallsSize();
		for (int i=0; i<sz; i++)
		{
			fusionForce = fusionResults.GetForceCall(i);
			gcosForce = gcosResults->GetForceCall(i);

			CPPUNIT_ASSERT(fusionForce.GetPosition() == gcosForce.position);
			CPPUNIT_ASSERT(fusionForce.GetCall() == gcosForce.call);
			CPPUNIT_ASSERT(fusionForce.GetReason() == gcosForce.reason);

			CPPUNIT_ASSERT(fusionForces[i].GetPosition() == gcosForces[i].position);
			CPPUNIT_ASSERT(fusionForces[i].GetCall() == gcosForces[i].call);
			CPPUNIT_ASSERT(fusionForces[i].GetReason() == gcosForces[i].reason);
		}

		BaseCallType gcosOrig;
		FusionBaseCallType fusionOrig;
		const FusionBaseCallVector& fusionOrigs = fusionResults.GetOrigCalls();
		const std::vector<BaseCallType>& gcosOrigs = gcosResults->GetOrigCalls();

		sz = fusionResults.GetOrigCallsSize();
		for (int i=0; i<sz; i++)
		{
			fusionOrig = fusionResults.GetOrigCall(i);
			gcosOrig = gcosResults->GetOrigCall(i);

			CPPUNIT_ASSERT(fusionOrig.GetPosition() == gcosOrig.position);
			CPPUNIT_ASSERT(fusionOrig.GetCall() == gcosOrig.call);

			CPPUNIT_ASSERT(fusionOrigs[i].GetPosition() == gcosOrigs[i].position);
			CPPUNIT_ASSERT(fusionOrigs[i].GetCall() == gcosOrigs[i].call);
		}
	}
	delete fusionChp;
}

// tests the background data for the chp file. This tests the backgroundzones, backgroundinfo and backgroundtype for each zone
void FusionCHPDataGCOSTest::BackgroundDataTest(string fileName)
{
	CPPUNIT_ADDITIONALMESSAGE_H
	// here we compare the background data contents with both the gcos reader and the fusion reader and they better be the same
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CCHPFileData gcosChp;
	gcosChp.SetFileName(fileName.c_str());
	CPPUNIT_ASSERT_MESSAGE("Failed: Read header!", gcosChp.ReadHeader() == true);

	BackgroundZoneInfo zoneInfoFusion;
	fusionChp->GetHeader().GetBackgroundZoneInfo(zoneInfoFusion);
	BackgroundZoneInfo zoneInfoGCOS = gcosChp.GetHeader().GetBackgroundZoneInfo();

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - number_zones",zoneInfoFusion.number_zones == zoneInfoGCOS.number_zones);
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - smooth_factor",zoneInfoFusion.smooth_factor == zoneInfoGCOS.smooth_factor);
	BackgroundZoneTypeList zonesFusion = zoneInfoFusion.zones;
	BackgroundZoneTypeList zonesGCOS = zoneInfoGCOS.zones;

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - GetChangeString",zonesFusion.size() == zonesGCOS.size());

	BackgroundZoneTypeList::iterator beginFusion = zonesFusion.begin(), endFusion = zonesFusion.end();
	BackgroundZoneTypeList::iterator beginGCOS = zonesGCOS.begin(), endGCOS = zonesGCOS.end();

	for(; beginFusion != endFusion && beginGCOS != endGCOS; beginFusion++, beginGCOS++)
	{
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - centerx",beginFusion->centerx == beginGCOS->centerx);
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - centery",beginFusion->centery == beginGCOS->centery);
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - background",beginFusion->background == beginGCOS->background);
	}
	delete fusionChp;
}
