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
#include "calvin_files/parsers/test/ArrayFileReaderTest.h"
//
#include "calvin_files/parsers/src/ArrayFileReader.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_array;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( ArrayFileReaderTest );

#define TEST_NON_EXISTANT_FILE "../data/file_does_not_exist";
#define TEST_DATA_FILE_INVALID "../data/test.file.data_header_only"
#define TEST_DATA_FILE "../data/ArraySetFile.xml"

void ArrayFileReaderTest::setUp()
{
}

void ArrayFileReaderTest::tearDown()
{
}

void ArrayFileReaderTest::testCreation()
{
	ArrayFileReader reader;
	CPPUNIT_ASSERT(1);
}

void ArrayFileReaderTest::CheckArrayData(bool headerOnly)
{
	ArrayFileReader reader;
	ArrayData array;
	std::string name = TEST_DATA_FILE;
	CPPUNIT_ASSERT(reader.Read(name, array, headerOnly) == true);

	CPPUNIT_ASSERT( array.DataTypeIdentifier() == ARRAY_SET_FILE_TYPE_IDENTIFIER );
	CPPUNIT_ASSERT(	array.ArraySetFileIdentifier() == std::string("432-432-432-432"));
	CPPUNIT_ASSERT( array.CreatedStep() == ArrayRegistrationStep);
	CPPUNIT_ASSERT( array.InitialProject() == L"my_project");
	CPPUNIT_ASSERT( array.CreationDateTime() == L"8/12/2005 9:00AM");
	CPPUNIT_ASSERT( array.CreatedBy() == L"ljevon");

	if (headerOnly)
		return;

	CPPUNIT_ASSERT( array.PhysicalArraysAttributes().size() == 1);

	ArrayAttributes &atts = array.PhysicalArraysAttributes()[0];

	CPPUNIT_ASSERT( atts.Identifier() == "123-123-123-123");
	CPPUNIT_ASSERT( atts.ArrayName() == "mychip");
	CPPUNIT_ASSERT( atts.ArrayBarcode() == "@1234567890");
	CPPUNIT_ASSERT( atts.Media() == PlateOrStripMedia);
	CPPUNIT_ASSERT( atts.MediaRow() == 1);
	CPPUNIT_ASSERT( atts.MediaCol() == 12);
	CPPUNIT_ASSERT( atts.MasterFile() == "Test3.master");
	CPPUNIT_ASSERT( atts.MasterFileId() == "123-123-123-123");
	CPPUNIT_ASSERT( atts.PatAssignment() == AffyBarcodeAssignment);
	CPPUNIT_ASSERT( atts.CreationDateTime() == L"8/12/2005 10:00AM");
	CPPUNIT_ASSERT( atts.CreatedBy() == L"ljevon");
	CPPUNIT_ASSERT( atts.CreatedStep() == ScanningStep);
	CPPUNIT_ASSERT( atts.Comment() == L"here is a comment");
	CPPUNIT_ASSERT( atts.LibraryPackageName() == "libpackage");
	CPPUNIT_ASSERT( atts.MediaFileGUID() == "mediaguid");
	CPPUNIT_ASSERT( atts.MediaFileName() == "mediafile");

	ParameterNameValuePairVector &params = atts.Attributes();
	ParameterNameValuePairVector::iterator paramIt = params.begin();
	CPPUNIT_ASSERT( params.size() == 2);
	CPPUNIT_ASSERT( (*paramIt).Name == L"SampleID");
	CPPUNIT_ASSERT( (*paramIt).Value == L"433232");
	++paramIt;
	CPPUNIT_ASSERT( (*paramIt).Name == L"LIMS System");
	CPPUNIT_ASSERT( (*paramIt).Value == L"Nautilis");


	ParameterNameValueDefaultRequiredTypeList &userParams = array.UserAttributes();
	CPPUNIT_ASSERT( userParams.size() == 6 );

	ParameterNameValueDefaultRequiredTypeList::iterator userIt = userParams.begin();
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Species" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"Homo Sapien" );
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == true );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::TextParameterType);
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Individual" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"me" );
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == false );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::TextParameterType);
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"When" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"now" );
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == false );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::TimeParameterType);
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"How much" );
	const double eps = 1e-5;
	CPPUNIT_ASSERT_DOUBLES_EQUAL( (*userIt).GetValueFloat(), 1.0, eps);
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == false );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::FloatParameterType);
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Sex" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"" );
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == false );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::ControlMultiParameterType);

	std::list<std::wstring>::iterator cit = (*userIt).ControlMultiValues().begin();
	CPPUNIT_ASSERT( (*cit) == L"Male" );
	++cit;
	CPPUNIT_ASSERT( (*cit) == L"Other" );

	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 4 );
	cit = (*userIt).ControlledVocabulary().begin();
	CPPUNIT_ASSERT( (*cit) == L"Male" );
	++cit;
	CPPUNIT_ASSERT( (*cit) == L"Female" );
	++cit;
	CPPUNIT_ASSERT( (*cit) == L"Unknown" );
	++cit;
	CPPUNIT_ASSERT( (*cit) == L"Other" );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Age" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"30 something" );
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == false );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::TextParameterType);
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
	CPPUNIT_ASSERT (userIt == userParams.end() );

}

void ArrayFileReaderTest::testmethod_Read()
{
	CheckArrayData(false);
}

void ArrayFileReaderTest::testmethod_Read_header_data_only()
{
	CheckArrayData(true);
}

void ArrayFileReaderTest::testmethod_Read_when_file_does_not_exist()
{
	ArrayFileReader reader;
	ArrayData array;
	std::string name = TEST_NON_EXISTANT_FILE;
	CPPUNIT_ASSERT(reader.Read(name, array) == false);
}

void ArrayFileReaderTest::testmethod_Read_when_file_is_not_valid()
{
	ArrayFileReader reader;
	ArrayData array;
	std::string name = TEST_DATA_FILE_INVALID;
	CPPUNIT_ASSERT(reader.Read(name, array) == false);
}

void ArrayFileReaderTest::testmethod_IsFileType()
{
	std::string file;

	file = TEST_NON_EXISTANT_FILE;
	CPPUNIT_ASSERT(ArrayFileReader::IsFileType(file, ARRAY_SET_FILE_TYPE_IDENTIFIER) == false);

	file = TEST_DATA_FILE_INVALID;
	CPPUNIT_ASSERT(ArrayFileReader::IsFileType(file, ARRAY_SET_FILE_TYPE_IDENTIFIER) == false);

	file = TEST_DATA_FILE;
	CPPUNIT_ASSERT( ArrayFileReader::IsFileType(file, ARRAY_SET_FILE_TYPE_IDENTIFIER) == true );
}

void ArrayFileReaderTest::testmethod_DataTypeIdentifierStatic()
{
	std::string name = TEST_DATA_FILE;
	affymetrix_calvin_utilities::AffymetrixGuidType guid;
	guid = ArrayFileReader::DataTypeIdentifier(name);

	CPPUNIT_ASSERT(guid == ARRAY_SET_FILE_TYPE_IDENTIFIER);
}
