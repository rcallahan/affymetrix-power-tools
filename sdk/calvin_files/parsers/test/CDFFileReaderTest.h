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

#ifndef __CDFFILEREADERTEST_H_
#define __CDFFILEREADERTEST_H_

//
#include "calvin_files/data/src/CDFData.h"
#include "calvin_files/data/src/CDFProbeSetInformation.h"
#include "calvin_files/data/src/CDFQCProbeSetInformation.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

//class CDFProbeSetInformation;

class CDFFileReaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CDFFileReaderTest );

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (ReadSmallCDFFileBasicTest);
	CPPUNIT_TEST (ReadSmallCDFFileSeqModeTest);
	CPPUNIT_TEST (ReadSmallCDFFileProbeSetNumberModeTest);
	CPPUNIT_TEST (ReadSmallCDFFileProbeSetNameModeTest);
	CPPUNIT_TEST (ReadSmallQCCDFFileBasicTest);
	CPPUNIT_TEST (ReadSmallQCCDFFileSeqModeTest);
	CPPUNIT_TEST (ReadSmallQCCDFFileProbeSetNumberModeTest);
	CPPUNIT_TEST (ReadSmallQCCDFFileProbeSetNameModeTest);
	CPPUNIT_TEST (BadFilenameTest);
	CPPUNIT_TEST (UseGetQCProbeSetInformationToReadCDFFileTest);
	CPPUNIT_TEST (UseGetProbeSEtInformationTOReadQCCDFFileTest);
	CPPUNIT_TEST (ReadCDFProbeSetInformationOpenSeqModeInWrongMode);
	CPPUNIT_TEST (ReadCDFProbeSetInformationOpenProbeIndexModeInWrongMode);
	CPPUNIT_TEST (ReadCDFProbeSetInformationOpenProbeNameModeInWrongMode);
	CPPUNIT_TEST (ReadQCCDFProbeSetInformationOpenSeqModeInWrongMode);
	CPPUNIT_TEST (ReadQCCDFProbeSetInformationOpenProbeIndexModeInWrongMode);
	CPPUNIT_TEST (ReadQCCDFProbeSetInformationOpenProbeNameModeInWrongMode);
	CPPUNIT_TEST (GetProbeSetInformationWithProbeSetNumberOutOfBoundsTest);
	CPPUNIT_TEST (GetQCProbeSetInformationWithProbeSetNumberOutOfBoundsTest);
	CPPUNIT_TEST (UnknownProbeSetNameTest);
	CPPUNIT_TEST (UnknownQCProbeSetNameTest);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
	void CreationTest();
	void ReadSmallCDFFileBasicTest();
	void ReadSmallCDFFileSeqModeTest();
	void ReadSmallCDFFileProbeSetNumberModeTest();
	void ReadSmallCDFFileProbeSetNameModeTest();
	void ReadSmallQCCDFFileBasicTest();
	void ReadSmallQCCDFFileSeqModeTest();
	void ReadSmallQCCDFFileProbeSetNumberModeTest();
	void ReadSmallQCCDFFileProbeSetNameModeTest();
	void BadFilenameTest();
	void UseGetQCProbeSetInformationToReadCDFFileTest();
	void UseGetProbeSEtInformationTOReadQCCDFFileTest();
	void ReadCDFProbeSetInformationOpenSeqModeInWrongMode();
	void ReadCDFProbeSetInformationOpenProbeIndexModeInWrongMode();
	void ReadCDFProbeSetInformationOpenProbeNameModeInWrongMode();
	void ReadQCCDFProbeSetInformationOpenSeqModeInWrongMode();
	void ReadQCCDFProbeSetInformationOpenProbeIndexModeInWrongMode();
	void ReadQCCDFProbeSetInformationOpenProbeNameModeInWrongMode();
	void GetProbeSetInformationWithProbeSetNumberOutOfBoundsTest();
	void GetQCProbeSetInformationWithProbeSetNumberOutOfBoundsTest();
	void UnknownProbeSetNameTest();
	void UnknownQCProbeSetNameTest();

// support methods
	void CheckSmallCDFProbeSetInformation(int32_t index, affymetrix_calvin_io::CDFProbeSetInformation& info);
	void CheckSmallQCCDFProbeSetInformation(int32_t index, affymetrix_calvin_io::CDFQCProbeSetInformation& psi);
	void CheckGetProbeSetName(affymetrix_calvin_io::CDFData& data);
	void CheckQCGetProbeSetName(affymetrix_calvin_io::CDFData& data);

};

#endif // __CDFFILEREADERTEST_H_
