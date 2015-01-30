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
#include "calvin_files/parsers/test/CDFFileReaderTest.h"
//
#include "calvin_files/data/src/CDFProbeGroupInformation.h"
#include "calvin_files/data/src/CDFProbeInformation.h"
#include "calvin_files/data/src/CDFQCProbeInformation.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/parsers/src/CDFFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_parameter;

CPPUNIT_TEST_SUITE_REGISTRATION( CDFFileReaderTest );

const std::string SMALL_CDF_FILE = "../data/small_CDF_file";
const std::string SMALL_QCCDF_FILE = "../data/small_QCCDF_file";

void CDFFileReaderTest::setUp()
{
}

void CDFFileReaderTest::tearDown()
{
}

void CDFFileReaderTest::CreationTest()
{
	CDFFileReader reader;
	CPPUNIT_ASSERT(1);
}

void CDFFileReaderTest::ReadSmallCDFFileBasicTest()
{
	CDFData data;
	CDFFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_CDF_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == SMALL_CDF_FILE);

	FileHeader* fh = data.GetFileHeader();
	u_int32_t pos = fh->GetFirstDataGroupFilePos();
	CPPUNIT_ASSERT(pos == 0x1BD);
	CPPUNIT_ASSERT(data.GetProbeSetCnt() == 10);
	CPPUNIT_ASSERT(data.GetArrayRows() == 10*2);
	CPPUNIT_ASSERT(data.GetArrayCols() == 11);
}

void CDFFileReaderTest::ReadSmallCDFFileSeqModeTest()
{
	CDFData data;
	CDFFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_CDF_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data, CDFFileReader::ReadSequential));

	int32_t probeSetCnt = data.GetProbeSetCnt();

	for (int32_t i = 0; i < probeSetCnt; ++i)
	{
		CDFProbeSetInformation info;
		data.GetProbeSetInformation(i, info);

		CheckSmallCDFProbeSetInformation(i, info);
	}

	CheckGetProbeSetName(data);
}

void CDFFileReaderTest::ReadSmallCDFFileProbeSetNumberModeTest()
{
	CDFData data;
	CDFFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_CDF_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data, CDFFileReader::ReadByProbeSetNumber));

	int32_t probeSetCnt = data.GetProbeSetCnt();

	// go backwards
	for (int32_t i = probeSetCnt-1; i >= 0; --i)
	{
		CDFProbeSetInformation info;
		data.GetProbeSetInformation(i, info);

		CheckSmallCDFProbeSetInformation(i, info);
	}

	CheckGetProbeSetName(data);
}

void CDFFileReaderTest::ReadSmallCDFFileProbeSetNameModeTest()
{
	CDFData data;
	CDFFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_CDF_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data, CDFFileReader::ReadByProbeSetName));

	int32_t probeSetCnt = data.GetProbeSetCnt();

	// go backwards
	for (int32_t i = probeSetCnt-1; i >= 0; --i)
	{
		wchar_t name[100];
		FormatString1(name, 100, L"biob_%d", i);

		CDFProbeSetInformation info;
		data.GetProbeSetInformation(name, info);

		CheckSmallCDFProbeSetInformation(i, info);
	}

	CheckGetProbeSetName(data);
}

void CDFFileReaderTest::CheckGetProbeSetName(CDFData& data)
{
	int32_t probeSetCnt = data.GetProbeSetCnt();
	wchar_t name[100];

	for (int32_t i = 0; i < probeSetCnt; ++i)
	{
		FormatString1(name, 100, L"biob_%d", i);

		CPPUNIT_ASSERT(data.GetProbeSetName(i) == name);
	}
}

void CDFFileReaderTest::CheckSmallCDFProbeSetInformation(int32_t index, CDFProbeSetInformation& psi)
{
	// values expected in the file
	// u_int8_t unitType = 3;
	u_int32_t numGroups = 1;
	DirectionType direction = ProbeSenseDirection;
	u_int32_t numLists = 11;
	u_int8_t cellsPerList = 2;
	u_int32_t numCells = numLists*cellsPerList;

	wchar_t name[100];
	FormatString1(name, 100, L"biob_%d", index);

	// Check the CDFProbeSetInformation
	// TODO: do we need to add a method to get unitType from CDFProbeSetInformation.  It goes in but doesn't come out.
	CPPUNIT_ASSERT(psi.GetName() == name);
	CPPUNIT_ASSERT(direction == psi.GetDirection());
	CPPUNIT_ASSERT(numLists == psi.GetNumLists());
	CPPUNIT_ASSERT(numGroups == psi.GetNumGroups());
	CPPUNIT_ASSERT(numCells == psi.GetNumCells());
	CPPUNIT_ASSERT(cellsPerList == psi.GetNumCellsPerList());
	CPPUNIT_ASSERT(index == psi.GetProbeSetNumber());

	// Check the CDFProbeGroupInformation
	CDFProbeGroupInformation pgi;
	CPPUNIT_ASSERT_NO_THROW(psi.GetGroupInformation(0, pgi));
	CPPUNIT_ASSERT(pgi.GetName() == name);
	CPPUNIT_ASSERT(direction == pgi.GetDirection());
	CPPUNIT_ASSERT(numLists == pgi.GetNumLists());
	CPPUNIT_ASSERT(numCells == pgi.GetNumCells());
	CPPUNIT_ASSERT(cellsPerList == pgi.GetNumCellsPerList());

	// Check the CDFProbeInformation
	CDFProbeInformation pi;
	for (u_int32_t list = 0; list < numLists; ++list)
	{
		u_int32_t startCell = list*cellsPerList;
		for (u_int32_t cellInList = 0; cellInList < cellsPerList; ++cellInList)
		{
			CPPUNIT_ASSERT_NO_THROW(pgi.GetCell(startCell+cellInList, pi));
			CPPUNIT_ASSERT(pi.GetX() == list);
			CPPUNIT_ASSERT(pi.GetY() == index*cellsPerList+cellInList);
			CPPUNIT_ASSERT(pi.GetListIndex() == list);
			CPPUNIT_ASSERT(pi.GetExpos() == list);
			CPPUNIT_ASSERT(pi.GetPBase() == 'C');
			CPPUNIT_ASSERT(pi.GetTBase() == 'G');
		}
	}
}

void CDFFileReaderTest::ReadSmallQCCDFFileBasicTest()
{
	CDFData data;
	CDFFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_QCCDF_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == SMALL_QCCDF_FILE);

	FileHeader* fh = data.GetFileHeader();
	u_int32_t pos = fh->GetFirstDataGroupFilePos();
	CPPUNIT_ASSERT(pos == 0x1BA);
	CPPUNIT_ASSERT(data.GetProbeSetCnt() == 10);
	CPPUNIT_ASSERT(data.GetArrayRows() == 10*2);
	CPPUNIT_ASSERT(data.GetArrayCols() == 11);
}

void CDFFileReaderTest::ReadSmallQCCDFFileSeqModeTest()
{
	CDFData data;
	CDFFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_QCCDF_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data, CDFFileReader::ReadSequential));

	int32_t probeSetCnt = data.GetProbeSetCnt();

	for (int32_t i = 0; i < probeSetCnt; ++i)
	{
		CDFQCProbeSetInformation info;
		data.GetQCProbeSetInformation(i, info);

		CheckSmallQCCDFProbeSetInformation(i, info);
	}

	CheckQCGetProbeSetName(data);
}

void CDFFileReaderTest::ReadSmallQCCDFFileProbeSetNumberModeTest()
{
	CDFData data;
	CDFFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_QCCDF_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data, CDFFileReader::ReadByProbeSetNumber));

	int32_t probeSetCnt = data.GetProbeSetCnt();

	// go backwards
	for (int32_t i = probeSetCnt-1; i >= 0; --i)
	{
		CDFQCProbeSetInformation info;
		data.GetQCProbeSetInformation(i, info);

		CheckSmallQCCDFProbeSetInformation(i, info);
	}

	CheckQCGetProbeSetName(data);
}

void CDFFileReaderTest::ReadSmallQCCDFFileProbeSetNameModeTest()
{
	CDFData data;
	CDFFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(SMALL_QCCDF_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data, CDFFileReader::ReadByProbeSetName));

	int32_t probeSetCnt = data.GetProbeSetCnt();

	// go backwards
	for (int32_t i = probeSetCnt-1; i >= 0; --i)
	{
		wchar_t name[100];
		FormatString1(name, 100, L"control_%d", i);

		CDFQCProbeSetInformation info;
		data.GetQCProbeSetInformation(name, info);

		CheckSmallQCCDFProbeSetInformation(i, info);
	}

	CheckQCGetProbeSetName(data);
}

void CDFFileReaderTest::CheckQCGetProbeSetName(CDFData& data)
{
	int32_t probeSetCnt = data.GetProbeSetCnt();
	wchar_t name[100];

	for (int32_t i = 0; i < probeSetCnt; ++i)
	{
		FormatString1(name, 100, L"control_%d", i);

		CPPUNIT_ASSERT(data.GetProbeSetName(i) == name);
	}
}

void CDFFileReaderTest::CheckSmallQCCDFProbeSetInformation(int32_t index, CDFQCProbeSetInformation& psi)
{
	// values expected in the file
	//u_int32_t numGroups = 1;
	u_int32_t numLists = 11;
	u_int8_t cellsPerList = 2;
	u_int32_t numCells = numLists*cellsPerList;

	wchar_t qcType[100];
	FormatString1(qcType, 100, L"control_%d", index);

	// Check the CDFQCProbeSetInformation
	CPPUNIT_ASSERT(psi.GetQCProbeSetType() == qcType);
	CPPUNIT_ASSERT(numCells == psi.GetNumCells());

	// Check the CDFQCProbeInformation
	CDFQCProbeInformation pi;
	for (u_int32_t list = 0; list < numLists; ++list)
	{
		u_int32_t startCell = list*cellsPerList;
		for (u_int32_t cellInList = 0; cellInList < cellsPerList; ++cellInList)
		{
			bool expectedPM = (cellInList % 2)? true : false;
			bool expectedBG = (list % 2) ? true : false;

			CPPUNIT_ASSERT_NO_THROW(psi.GetProbeInformation(startCell+cellInList, pi));
			CPPUNIT_ASSERT(pi.GetX() == list);
			CPPUNIT_ASSERT(pi.GetY() == index*cellsPerList+cellInList);
			CPPUNIT_ASSERT(pi.GetPLen() == 25);
			CPPUNIT_ASSERT(pi.IsPerfectMatchProbe() == expectedPM);
			CPPUNIT_ASSERT(pi.IsBackgroundProbe() == expectedBG);
		}
	}
}

void CDFFileReaderTest::BadFilenameTest()
{
	CDFData data;
	CDFFileReader reader;
	data.SetFilename("CDFMetaDataBad");
	CPPUNIT_ASSERT_THROW(reader.Read(data, CDFFileReader::ReadSequential), affymetrix_calvin_exceptions::FileNotFoundException);
}

void CDFFileReaderTest::UseGetQCProbeSetInformationToReadCDFFileTest()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_CDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetNumber);
	CDFQCProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetQCProbeSetInformation(0, info), ProbeSetNotFoundException);
	CPPUNIT_ASSERT_THROW(data.GetQCProbeSetInformation(L"hello", info), ProbeSetNotFoundException);
}

void CDFFileReaderTest::UseGetProbeSEtInformationTOReadQCCDFFileTest()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_QCCDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetNumber);
	CDFProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetProbeSetInformation(0, info), ProbeSetNotFoundException);
	CPPUNIT_ASSERT_THROW(data.GetProbeSetInformation(L"hello", info), ProbeSetNotFoundException);
}

void CDFFileReaderTest::ReadCDFProbeSetInformationOpenSeqModeInWrongMode()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_CDF_FILE);
	reader.Read(data, CDFFileReader::ReadSequential);

	CDFProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetProbeSetInformation(4, info), CDFAccessNotSupportedByModeException);
	CPPUNIT_ASSERT_THROW(data.GetProbeSetInformation(L"hello", info), CDFAccessNotSupportedByModeException);
}

void CDFFileReaderTest::ReadCDFProbeSetInformationOpenProbeIndexModeInWrongMode()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_CDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetNumber);

	CDFProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetProbeSetInformation(L"hello", info), CDFAccessNotSupportedByModeException);
}

void CDFFileReaderTest::ReadCDFProbeSetInformationOpenProbeNameModeInWrongMode()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_CDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetName);

	CDFProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetProbeSetInformation(3, info), CDFAccessNotSupportedByModeException);
}

void CDFFileReaderTest::ReadQCCDFProbeSetInformationOpenSeqModeInWrongMode()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_QCCDF_FILE);
	reader.Read(data, CDFFileReader::ReadSequential);

	CDFQCProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetQCProbeSetInformation(4, info), CDFAccessNotSupportedByModeException);
	CPPUNIT_ASSERT_THROW(data.GetQCProbeSetInformation(L"hello", info), CDFAccessNotSupportedByModeException);
}

void CDFFileReaderTest::ReadQCCDFProbeSetInformationOpenProbeIndexModeInWrongMode()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_QCCDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetNumber);

	CDFQCProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetQCProbeSetInformation(L"hello", info), CDFAccessNotSupportedByModeException);
}

void CDFFileReaderTest::ReadQCCDFProbeSetInformationOpenProbeNameModeInWrongMode()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_QCCDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetName);

	CDFQCProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetQCProbeSetInformation(3, info), CDFAccessNotSupportedByModeException);
}

void CDFFileReaderTest::GetProbeSetInformationWithProbeSetNumberOutOfBoundsTest()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_CDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetNumber);

	CDFProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetProbeSetInformation(11, info), ProbeSetNotFoundException);
}

void CDFFileReaderTest::GetQCProbeSetInformationWithProbeSetNumberOutOfBoundsTest()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_QCCDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetNumber);

	CDFQCProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetQCProbeSetInformation(11, info), ProbeSetNotFoundException);
}

void CDFFileReaderTest::UnknownProbeSetNameTest()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_CDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetName);

	CDFProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetProbeSetInformation(L"hello", info), ProbeSetNotFoundException);
}

void CDFFileReaderTest::UnknownQCProbeSetNameTest()
{
	CDFData data;
	CDFFileReader reader;
	reader.SetFilename(SMALL_QCCDF_FILE);
	reader.Read(data, CDFFileReader::ReadByProbeSetName);

	CDFQCProbeSetInformation info;
	CPPUNIT_ASSERT_THROW(data.GetQCProbeSetInformation(L"hello", info), ProbeSetNotFoundException);
}

