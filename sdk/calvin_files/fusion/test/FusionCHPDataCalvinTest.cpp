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
#include "calvin_files/fusion/test/FusionCHPDataCalvinTest.h"
//
#include "calvin_files/exception/src/DevelopmentException.h"
#include "calvin_files/fusion/src/FusionCHPLegacyData.h"
#include "calvin_files/parsers/src/CHPFileReader.h"	
#include "calvin_files/parsers/src/FileException.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cstring>
#include <map>
#include <string>
//

using namespace std;
using namespace affxchp;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_utilities;

// Underlying CCelFileData class expects backslashes rather than forward slashes on Windows hmmmm
const string testFiles[] = 
{
	"doesnt_exist.CHP",														// file does not exist - for testing when a file does not exist
//	//"../../../TestDataFiles/CHP/GCOS104-U133a.CHP",
//	//"../../../TestDataFiles/CHP/pcwrs041a_pd01.CHP",
//	//"../../../TestDataFiles/CHP/NA02705_00243_01_3001197.CHP",
	"../../../TestDataFiles/CHP/10k_1-calvin.CHP",
//	//"../../../TestDataFiles/CHP/T1_r1_stat.CHP",
//	//"../../../TestDataFiles/CHP/T1_r1_v_T2_r1_stat.CHP",
	"../../../TestDataFiles/CHP/DCN_01_xda-calvin.CHP",
//	//"../../../TestDataFiles/CHP/GCOS104-U133a-test_rel.CHP",
//	//"../../../TestDataFiles/CHP/StatData_baseline.CHP",
	"../../../TestDataFiles/CHP/StatData-calvin.CHP",
	"../../../TestDataFiles/CHP/PLIERdata1-calvin.CHP"
};
const int FILE_NAME_LIST_COUNT = 5; // does not count the does not exist file case because we do not explicitly test it

const int DOES_NOT_EXIST = 0;	// index into the file name array above of the file does not exist name

CPPUNIT_TEST_SUITE_REGISTRATION(FusionCHPDataCalvinTest);

// has to be called before each test
void FusionCHPDataCalvinTest::setUp()
{
}

// has to be called after each test
void FusionCHPDataCalvinTest::tearDown()
{
}

void FusionCHPDataCalvinTest::TestFileId()
{
	if (FileUtils::Exists(testFiles[3].c_str()) == true)
	{
		FusionCHPData *chp = FusionCHPDataReg::ReadHeader(testFiles[3]);
		FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
		CPPUNIT_ASSERT(fusionChp->FileId() == "0000065535-1131145034-0000006334-0000018467-0000000041");
		CPPUNIT_ASSERT(fusionChp->GetGenericData()->FileIdentifier() == "0000065535-1131145034-0000006334-0000018467-0000000041");
	}
}

// test to see if we can even creat the class
void FusionCHPDataCalvinTest::CreationTest(string fileName)
{
	//FusionCHPLegacyData * fusionChp = new FusionCHPLegacyData;
	//CPPUNIT_ASSERT_MESSAGE("Failed: Create fusion object!",fusionChp);
	//delete fusionChp;
}

// main test entry method, calls the rest.
// had to do it this way, because CPPUINT will not accept test method parameters
void FusionCHPDataCalvinTest::TestAll()
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
void FusionCHPDataCalvinTest::ReadHeaderTest(string fileName)
{
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	delete fusionChp;
}

// test the rows and colums
void FusionCHPDataCalvinTest::RowsColsTest(string fileName)
{
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # rows!",fusionChp->GetHeader().GetRows() == calvinChp.GetRows());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # cols!",fusionChp->GetHeader().GetCols() == calvinChp.GetCols());
	delete fusionChp;
}

// test the chp file version
void FusionCHPDataCalvinTest::VersionTest(string fileName)
{
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Version",fusionChp->GetHeader().GetVersion() == calvinChp.GetVersion());
	delete fusionChp;
}

// test the algorithm and alg. verion info
void FusionCHPDataCalvinTest::AlgorithmNameAndVersionTest(string fileName)
{
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg name!", fusionChp->GetHeader().GetAlgName() == calvinChp.GetAlgName());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg version!", fusionChp->GetHeader().GetAlgVersion() == calvinChp.GetAlgVersion());
	delete fusionChp;
}

// test the chiptype name
void FusionCHPDataCalvinTest::ChipTypeTest(string fileName)
{
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg version!", fusionChp->GetHeader().GetChipType() == calvinChp.GetArrayType());
	delete fusionChp;
}

// test the parent cell file name
void FusionCHPDataCalvinTest::ParentCelFileName(string fileName)
{
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - parent cell file name!", fusionChp->GetHeader().GetParentCellFile() == calvinChp.GetParentCell());
	delete fusionChp;
}

// test the prog id
void FusionCHPDataCalvinTest::ProgIDTest(string fileName)
{
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	std::wstring progid;
	ParameterNameValueType paramType;
	CPPUNIT_ASSERT_MESSAGE("Failed: ProgID not found in calvin file", calvinChp.GetFileHeader()->GetGenericDataHdr()->FindNameValParam(CHP_PROGID, paramType));
	progid = paramType.GetValueText();
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - ProgID!", fusionChp->GetHeader().GetProgID() == progid);
	delete fusionChp;
}

// test the algorithm params
void FusionCHPDataCalvinTest::AlgParamsTest(string fileName)
{
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

  FusionTagValuePairTypeList algParamsFusion;
  fusionChp->GetHeader().AlgorithmParameters(algParamsFusion);
  FusionTagValuePairTypeList::iterator beginFusion = algParamsFusion.begin();
  FusionTagValuePairTypeList::iterator endFusion = algParamsFusion.end();

	const ParameterNameValueTypeVector& algParamsCalvin = calvinChp.GetAlgParams();
	ParameterNameValueTypeVector::const_iterator beginCalvin = algParamsCalvin.begin(), endCalvin = algParamsCalvin.end();

	CPPUNIT_ASSERT(algParamsFusion.size() == algParamsCalvin.size());

	typedef map<wstring,wstring,less<wstring> > theFusionMap;
	theFusionMap map;

	// need the map to have an easy lookup for the key to see if the fusion has the key
	for (; beginFusion != endFusion; beginFusion++)
	{
		map.insert(theFusionMap::value_type(beginFusion->Tag,beginFusion->Value));
	}

	theFusionMap::const_iterator findIt;
	for(; beginCalvin != endCalvin; beginCalvin++)
	{
		// test single alg value call
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare Algorithm Parameter!", fusionChp->GetHeader().GetAlgorithmParameter(beginCalvin->GetName().c_str()) == beginCalvin->ToString());

		// test to see if the keys and values are the same between what was retrieved from fusion get algorithm params and
		// with what was retrieved by the Calvin api
		findIt = map.find(beginCalvin->GetName());
		CPPUNIT_ASSERT_MESSAGE("Faile: Algorithm parameter not found in Fusion list", findIt != map.end());
		CPPUNIT_ASSERT_MESSAGE("Faile: Algorithm parameter not found in Fusion list", findIt->second == beginCalvin->ToString());
	}
	delete fusionChp;
}

// tests the summary params
void FusionCHPDataCalvinTest::SummaryParamsTest(string fileName)
{
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

  FusionTagValuePairTypeList sumParamsFusion;
  fusionChp->GetHeader().SummaryParameters(sumParamsFusion);
  FusionTagValuePairTypeList::iterator beginFusion = sumParamsFusion.begin();
  FusionTagValuePairTypeList::iterator endFusion = sumParamsFusion.end();

	const ParameterNameValueTypeVector& sumParamsCalvin = calvinChp.GetChipSums();
	ParameterNameValueTypeVector::const_iterator beginCalvin = sumParamsCalvin.begin(), endCalvin = sumParamsCalvin.end();

	typedef map<wstring,wstring,less<wstring> > theFusionMap;
	theFusionMap map;

	// need the map to have an easy lookup for the key to see if the fusion has the key
	for (; beginFusion != endFusion; beginFusion++)
	{
		map.insert(theFusionMap::value_type(beginFusion->Tag,beginFusion->Value));
	}

	theFusionMap::const_iterator findIt;
	for(; beginCalvin != endCalvin; beginCalvin++)
	{
		// test single alg value call
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare Summary Parameter!", fusionChp->GetHeader().GetSummaryParameter(beginCalvin->GetName().c_str()) == beginCalvin->ToString());

		// test to see if the keys and values are the same between what was retrieved from fusion get algorithm params and
		// with what was retrieved by the calvin api
		findIt = map.find(beginCalvin->GetName());
		CPPUNIT_ASSERT_MESSAGE("Faile: Algorithm parameter not found in Fusion list", findIt != map.end());
		CPPUNIT_ASSERT_MESSAGE("Faile: Algorithm parameter not found in Fusion list", findIt->second == beginCalvin->ToString());
	}
	delete fusionChp;
}

// tests the number of probesets
void FusionCHPDataCalvinTest::NumberOfSetsTest(string fileName)
{
	// test the number of sets
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # probesets",fusionChp->GetHeader().GetNumProbeSets() == calvinChp.GetEntryCount());
	delete fusionChp;
}

// tests the algorithm parameters
void FusionCHPDataCalvinTest::ReadHeaderAndVerifyAlgParamsTest(string fileName)
{
	// read the header with the fusion api
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	// read the header with the calvin api
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	// get the calvin list
	ParameterNameValueTypeVector listCalvin = calvinChp.GetAlgParams();
	FusionTagValuePairTypeList listFusion;
	// get the fusion list
	fusionChp->GetHeader().AlgorithmParameters(listFusion);

	// test that the number of parameters are the same between a pure calvin file sdk call and the fusdion api call
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg param count",listFusion.size() == listCalvin.size());

	// this section will test that the parameters read by both the pure calvin api and the fusion are the same
	// first we get the interators for the loop
	FusionTagValuePairTypeList::iterator beginFusion = listFusion.begin();
	FusionTagValuePairTypeList::iterator endFusion = listFusion.end();
	ParameterNameValueTypeVector::iterator beginCalvin = listCalvin.begin();
	ParameterNameValueTypeVector::iterator endCalvin = listCalvin.end();

	// loop over every parameter in the fusion parameter list and use it to get the value thru the getalgorithmparameter call
	// we could of used the value property of the of the iterator that is from the fusion algorithm list retrieved
	// early, but this is to exerise the get individual alg params api for both.
	// However, I will test the interator Value property from the fusion list and compare it to the other values returned from the
	// the get alg parameter calls
	for (; beginFusion != endFusion && beginCalvin != endCalvin ; beginFusion++, beginCalvin++)
	{
		// get the fusion value for the parameter name
		wstring tempFusion = fusionChp->GetHeader().GetAlgorithmParameter(beginFusion->Tag.c_str());
		// get the calvin value for the parameter name
		// have to convert the string to wstring first though
		ParameterNameValueType tempCalvin = calvinChp.GetAlgParam(beginFusion->Tag);

		// double check that the value from the fusion list is the same as the other values retieved from the get alg params calls
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params", tempFusion == tempCalvin.ToString() && beginFusion->Value == tempCalvin.ToString());

		// check the Detail information
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params", beginFusion->DetailedType().GetName() == beginCalvin->GetName());
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params", beginFusion->DetailedType().GetParameterType() == beginCalvin->GetParameterType());
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params", beginFusion->DetailedType().ToString() == beginCalvin->ToString());
	}
	delete fusionChp;
}

// tests the summary parameters
void FusionCHPDataCalvinTest::ReadHeaderAndVerifySummaryParamsTest(string fileName)
{
	// read the header with the fusion api
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	// read the header with the calvin api
	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	// get the calvin list
	ParameterNameValueTypeVector listCalvin = calvinChp.GetChipSums();
	FusionTagValuePairTypeList listFusion;
	// get the fusion list
	fusionChp->GetHeader().SummaryParameters(listFusion);

	// test that the number of parameters are the same between a pure calvin file sdk call and the fusdion api call
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - param list count",listFusion.size() == listCalvin.size());

	// this section will test that the parameters read by both the pure calvin api and the fusion are the same
	// first we get the interators for the loop
	FusionTagValuePairTypeList::iterator beginFusion = listFusion.begin();
	FusionTagValuePairTypeList::iterator endFusion = listFusion.end();
	ParameterNameValueTypeVector::iterator beginCalvin = listCalvin.begin();
	ParameterNameValueTypeVector::iterator endCalvin = listCalvin.end();

	// loop over every parameter in the fusion parameter list and use it to get the value thru the getsummaryparameter call
	// we could of used the value property of the of the iterator that is from the fusion summary list retrieved
	// early, but this is to exerise the get individual summary params api for both.
	// However, I will test the interator Value property from the fusion list and compare it to the other values returned from the
	// the get summary parameter calls
	for (; beginFusion != endFusion && beginCalvin != endCalvin ; beginFusion++, beginCalvin++)
	{
		// get the fusion value for the parameter name
		wstring tempFusion = fusionChp->GetHeader().GetSummaryParameter(beginFusion->Tag.c_str());
		// get the calvin value for the parameter name
		// have to convert the string to wstring first though
		ParameterNameValueType tempCalvin = calvinChp.GetChipSum(beginFusion->Tag);

		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params", tempFusion == tempCalvin.ToString() && beginFusion->Value == tempCalvin.ToString());

		// check the Detail information
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params", beginFusion->DetailedType().GetName() == beginCalvin->GetName());
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params", beginFusion->DetailedType().GetParameterType() == beginCalvin->GetParameterType());
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - params", beginFusion->DetailedType().ToString() == beginCalvin->ToString());
	}
	delete fusionChp;
}

// test the read method of fusion. This differs only from previous readheader tests in that instead of readheader being called
// read is called.
void FusionCHPDataCalvinTest::ReadTest(string fileName)
{
	FusionCHPData *chp = FusionCHPDataReg::Read(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	// Check some header values
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - version",fusionChp->GetHeader().GetVersion() == calvinChp.GetVersion());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - cols",fusionChp->GetHeader().GetCols() == calvinChp.GetCols());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - rows",fusionChp->GetHeader().GetRows() == calvinChp.GetRows());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # probesets",fusionChp->GetHeader().GetNumProbeSets() == calvinChp.GetEntryCount());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg name",fusionChp->GetHeader().GetAlgName() == calvinChp.GetAlgName());
  CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg version",fusionChp->GetHeader().GetAlgVersion() == calvinChp.GetAlgVersion());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - chip type",fusionChp->GetHeader().GetChipType() == calvinChp.GetArrayType());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - parent cell name",fusionChp->GetHeader().GetParentCellFile() == calvinChp.GetParentCell());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - alg count",fusionChp->GetHeader().AlgorithmParameterCount() == calvinChp.GetAlgParams().size());
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - sum count",fusionChp->GetHeader().SummaryParameterCount() == calvinChp.GetChipSums().size());
	//CPPUNIT_ASSERT_MESSAGE("Failed: Compare - prog id",fusionChp->GetHeader().GetProgID() == StringUtils::ConvertMBSToWCS(calvinChp.GetHeader().GetProgID()));
	delete fusionChp;
}

void FusionCHPDataCalvinTest::CallMethodsWhenObjectIsNotReadyTest(string fileName)
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
void FusionCHPDataCalvinTest::DataCompareTest(string fileName)
{
	// here we compare the data contents with both the calvin reader and the fusion reader and they better be the same
	FusionCHPData *chp = FusionCHPDataReg::Read(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	CHPExpressionEntry resultsCalvinExpression;
	FusionExpressionProbeSetResults resultsFusionExpression;

	CHPGenotypeEntry resultsCalvinGeno;
	FusionGenotypeProbeSetResults resultsFusionGeno;

	CHPUniversalEntry resultsCalvinUni;
	FusionUniversalProbeSetResults resultsFusionUni;

	u_int32_t count = fusionChp->GetHeader().GetNumProbeSets(); 
	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - # probesets",count == calvinChp.GetEntryCount());


	if(fusionChp->GetHeader().GetAssayType() == FusionExpression)
	{
		// for each type supported, get the probeset results
		// loop over the number of probesets
		for(u_int32_t i = 0; i < count; ++i)
		{
			CPPUNIT_ASSERT(fusionChp->GetExpressionResults(i, resultsFusionExpression));
			calvinChp.GetEntry(i, resultsCalvinExpression);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - DetectionPValue",resultsFusionExpression.GetDetectionPValue() == resultsCalvinExpression.GetDetectionPValue());
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Signal",resultsFusionExpression.GetSignal() == resultsCalvinExpression.GetSignal());
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Signal",resultsFusionExpression.GetNumPairs() == resultsCalvinExpression.GetNumPairs());
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - NumPairs",resultsFusionExpression.GetNumUsedPairs() == resultsCalvinExpression.GetNumPairsUsed());
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Detection",resultsFusionExpression.GetDetection() == resultsCalvinExpression.GetDetection());
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - m_HasCompResults",resultsFusionExpression.HasCompResults() == resultsCalvinExpression.GetHasComparisonData());
			if (resultsFusionExpression.HasCompResults())
			{
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - ChangePValue",resultsFusionExpression.GetChangePValue() == resultsCalvinExpression.GetChangePValue());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - SignalLogRatio",resultsFusionExpression.GetSignalLogRatio() == resultsCalvinExpression.GetSigLogRatio());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - SignalLogRatioLow",resultsFusionExpression.GetSignalLogRatioLow() == resultsCalvinExpression.GetSigLogRatioLo());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - SignalLogRatioHigh",resultsFusionExpression.GetSignalLogRatioHigh() == resultsCalvinExpression.GetSigLogRatioHi());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - NumCommonPairs",resultsFusionExpression.GetNumCommonPairs() == resultsCalvinExpression.GetCommonPairs());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Change",resultsFusionExpression.GetChange() == resultsCalvinExpression.GetChange());
	//			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - GetDetectionString",resultsFusionExpression.GetDetectionString() == resultsCalvinExpression.GetDetectionString());
	//			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - GetChangeString",resultsFusionExpression.GetChangeString() == resultsCalvinExpression.GetChangeString());
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
			calvinChp.GetEntry(i, resultsCalvinGeno);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - AlleleCall",resultsFusionGeno.GetAlleleCall() == resultsCalvinGeno.GetCall());
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Confidence",resultsFusionGeno.GetConfidence() == resultsCalvinGeno.GetConfidence());
	//		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - AlleleCallString",resultsFusionGeno.GetAlleleCallString() == resultsCalvinGeno.GetAlleleCallString());

			if(bWholeGenenome)
			{
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - RAS1",resultsFusionGeno.GetRAS1() == resultsCalvinGeno.GetRAS1());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - RAS2",resultsFusionGeno.GetRAS2() == resultsCalvinGeno.GetRAS2());
			}
			else if(bDynamicModel)
			{
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - pvalue_AA",resultsFusionGeno.GetPValueAA() == resultsCalvinGeno.GetAACall());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - pvalue_AB",resultsFusionGeno.GetPValueAB() == resultsCalvinGeno.GetABCall());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - pvalue_BB",resultsFusionGeno.GetPValueBB() == resultsCalvinGeno.GetBBCall());
				CPPUNIT_ASSERT_MESSAGE("Failed: Compare - pvalue_NoCall",resultsFusionGeno.GetPValueNoCall() == resultsCalvinGeno.GetNoCall());
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
			calvinChp.GetEntry(i, resultsCalvinUni);
			CPPUNIT_ASSERT_MESSAGE("Failed: Compare - Background",resultsFusionUni.GetBackground() == resultsCalvinUni.GetBackground());
		}
	}
	else if(fusionChp->GetHeader().GetAssayType() == FusionResequencing)
	{
		FusionResequencingResults fusionResults;
		fusionChp->GetReseqResults(fusionResults);

		CPPUNIT_ASSERT(fusionResults.GetCalledBasesSize() == calvinChp.GetEntryCount());
		CPPUNIT_ASSERT(fusionResults.GetForceCallsSize() == calvinChp.GetForceCnt());
		CPPUNIT_ASSERT(fusionResults.GetOrigCallsSize() == calvinChp.GetOrigCnt());

		const Int8Vector &fusionCalls = fusionResults.GetCalledBases();
		const FloatVector &fusionScores = fusionResults.GetScores();

		int sz = fusionResults.GetCalledBasesSize();
		for (int i=0; i<sz; i++)
		{
			CHPReseqEntry calvinResult;
			calvinChp.GetEntry(i, calvinResult);
			CPPUNIT_ASSERT(fusionResults.GetCalledBase(i) == calvinResult.call);
			CPPUNIT_ASSERT(fusionResults.GetScore(i) == calvinResult.score);
			CPPUNIT_ASSERT(fusionCalls[i] == calvinResult.call);
			CPPUNIT_ASSERT(fusionScores[i] == calvinResult.score);
		}

		FusionForceCallType fusionForce;
		const FusionForceCallVector& fusionForces = fusionResults.GetForceCalls();

		sz = fusionResults.GetForceCallsSize();
		for (int i=0; i<sz; i++)
		{
			fusionForce = fusionResults.GetForceCall(i);

			CHPReseqForceCall calvinForce;
			calvinChp.GetForceCall(i, calvinForce);

			CPPUNIT_ASSERT(fusionForce.GetPosition() == calvinForce.position);
			CPPUNIT_ASSERT(fusionForce.GetCall() == calvinForce.call);
			CPPUNIT_ASSERT(fusionForce.GetReason() == calvinForce.reason);

			CPPUNIT_ASSERT(fusionForces[i].GetPosition() == calvinForce.position);
			CPPUNIT_ASSERT(fusionForces[i].GetCall() == calvinForce.call);
			CPPUNIT_ASSERT(fusionForces[i].GetReason() == calvinForce.reason);
		}

		FusionBaseCallType fusionOrig;
		const FusionBaseCallVector& fusionOrigs = fusionResults.GetOrigCalls();

		sz = fusionResults.GetOrigCallsSize();
		for (int i=0; i<sz; i++)
		{
			fusionOrig = fusionResults.GetOrigCall(i);

			CHPReseqOrigCall calvinOrig;
			calvinChp.GetOrigCall(i, calvinOrig);

			CPPUNIT_ASSERT(fusionOrig.GetPosition() == calvinOrig.position);
			CPPUNIT_ASSERT(fusionOrig.GetCall() == calvinOrig.call);

			CPPUNIT_ASSERT(fusionOrigs[i].GetPosition() == calvinOrig.position);
			CPPUNIT_ASSERT(fusionOrigs[i].GetCall() == calvinOrig.call);
		}
	}
	delete fusionChp;
}

// tests the background data for the chp file. This tests the backgroundzones, backgroundinfo and backgroundtype for each zone
void FusionCHPDataCalvinTest::BackgroundDataTest(string fileName)
{
	CPPUNIT_ADDITIONALMESSAGE_H
	// here we compare the background data contents with both the calvin reader and the fusion reader and they better be the same
	FusionCHPData *chp = FusionCHPDataReg::ReadHeader(fileName);
	FusionCHPLegacyData *fusionChp = FusionCHPLegacyData::FromBase(chp);
	CPPUNIT_ASSERT(fusionChp != NULL);

	CHPData calvinChp;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(fileName));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(calvinChp));

	BackgroundZoneInfo zoneInfoFusion;
	fusionChp->GetHeader().GetBackgroundZoneInfo(zoneInfoFusion);

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - number_zones",zoneInfoFusion.number_zones == calvinChp.GetBackgroundZoneCnt());
//	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - smooth_factor",zoneInfoFusion.smooth_factor == zoneInfocalvin.smooth_factor);
	BackgroundZoneTypeList zonesFusion = zoneInfoFusion.zones;

//	CHPBackgroundZoneVector calvinZones;
//	calvinChp.GetBackgroundZones(0, ?, calvinZones);

	CPPUNIT_ASSERT_MESSAGE("Failed: Compare - GetChangeString",zonesFusion.size() == calvinChp.GetBackgroundZoneCnt());

	BackgroundZoneTypeList::iterator beginFusion = zonesFusion.begin(), endFusion = zonesFusion.end();

	for(int32_t izone = 0; beginFusion != endFusion && izone < calvinChp.GetBackgroundZoneCnt(); beginFusion++, izone++)
	{
		CHPBackgroundZone calvinZone;
		calvinChp.GetBackgroundZone(izone, calvinZone);
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - centerx",beginFusion->centerx == calvinZone.GetCenterX());
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - centery",beginFusion->centery == calvinZone.GetCenterY());
		CPPUNIT_ASSERT_MESSAGE("Failed: Compare - background",beginFusion->background == calvinZone.GetBackground());
	}
	delete fusionChp;
}
