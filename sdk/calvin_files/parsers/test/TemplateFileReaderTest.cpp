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
#include "calvin_files/parsers/test/TemplateFileReaderTest.h"
//
#include "calvin_files/parsers/src/TemplateFileReader.h"
#include "calvin_files/template/src/TemplateId.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_template;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( TemplateFileReaderTest );

#define TEST_NON_EXISTANT_FILE "../data/file_does_not_exist";
#define TEST_DATA_FILE_INVALID "../data/test.file.data_header_only"
#define TEST_DATA_FILE "../data/MIAME.xml"

void TemplateFileReaderTest::setUp()
{
}

void TemplateFileReaderTest::tearDown()
{
}

void TemplateFileReaderTest::testCreation()
{
	TemplateFileReader reader;
	CPPUNIT_ASSERT(1);
}

void TemplateFileReaderTest::CheckTemplateData(bool headerOnly)
{
	TemplateFileReader reader;
	TemplateData temp;
	std::string name = TEST_DATA_FILE;
	CPPUNIT_ASSERT(reader.Read(name, temp, headerOnly) == true);

	CPPUNIT_ASSERT( temp.DataTypeIdentifier() == TEMPLATE_FILE_TYPE_IDENTIFIER );
	CPPUNIT_ASSERT(	temp.TemplateFileIdentifier() == std::string("987-987-987"));
	CPPUNIT_ASSERT( temp.CreationDateTime() == L"8/12/2005 9:00AM");
	CPPUNIT_ASSERT( temp.CreatedBy() == L"ljevon");

	if (headerOnly)
		return;

	ParameterNameValueDefaultRequiredTypeList &userParams = temp.UserAttributes();
	CPPUNIT_ASSERT( userParams.size() == 23 );

	ParameterNameValueDefaultRequiredTypeList::iterator userIt = userParams.begin();
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Genus" );
	CPPUNIT_ASSERT( (*userIt).GetDefaultValueText() == L"Unknown" );
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == true );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::TextParameterType);
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Species" );
	CPPUNIT_ASSERT( (*userIt).GetDefaultValueText() == L"" );
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == false );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::TextParameterType);
	++userIt;
	++userIt;
	++userIt;
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Sex" );
	CPPUNIT_ASSERT( (*userIt).GetDefaultValueText() == L"" );
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == false );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 4 );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::ControlSingleParameterType);
	std::list<std::wstring>::iterator cit = (*userIt).ControlledVocabulary().begin();
	CPPUNIT_ASSERT( (*cit) == L"Male" );
	++cit;
	CPPUNIT_ASSERT( (*cit) == L"Female" );
	++cit;
	CPPUNIT_ASSERT( (*cit) == L"Unknown" );
	++cit;
	CPPUNIT_ASSERT( (*cit) == L"Other" );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Age" );
	CPPUNIT_ASSERT( (*userIt).GetDefaultValueText() == L"" );
	CPPUNIT_ASSERT( (*userIt).RequiredFlag() == false );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*userIt).ValueType() == ParameterNameValueDefaultRequiredType::TextParameterType);
}

void TemplateFileReaderTest::testmethod_Read()
{
	CheckTemplateData(false);
}

void TemplateFileReaderTest::testmethod_Read_header_data_only()
{
	CheckTemplateData(true);
}

void TemplateFileReaderTest::testmethod_Read_when_file_does_not_exist()
{
	TemplateFileReader reader;
	TemplateData temp;
	std::string name = TEST_NON_EXISTANT_FILE;
	CPPUNIT_ASSERT(reader.Read(name, temp) == false);
}

void TemplateFileReaderTest::testmethod_Read_when_file_is_not_valid()
{
	TemplateFileReader reader;
	TemplateData temp;
	std::string name = TEST_DATA_FILE_INVALID;
	CPPUNIT_ASSERT(reader.Read(name, temp) == false);
}

void TemplateFileReaderTest::testmethod_IsFileType()
{
	std::string file;

	file = TEST_NON_EXISTANT_FILE;
	CPPUNIT_ASSERT(TemplateFileReader::IsFileType(file, TEMPLATE_FILE_TYPE_IDENTIFIER) == false);

	file = TEST_DATA_FILE_INVALID;
	CPPUNIT_ASSERT(TemplateFileReader::IsFileType(file, TEMPLATE_FILE_TYPE_IDENTIFIER) == false);

	file = TEST_DATA_FILE;
	CPPUNIT_ASSERT( TemplateFileReader::IsFileType(file, TEMPLATE_FILE_TYPE_IDENTIFIER) == true );
}

void TemplateFileReaderTest::testmethod_DataTypeIdentifierStatic()
{
	std::string name = TEST_DATA_FILE;
	affymetrix_calvin_utilities::AffymetrixGuidType guid;
	guid = TemplateFileReader::DataTypeIdentifier(name);

	CPPUNIT_ASSERT(guid == TEMPLATE_FILE_TYPE_IDENTIFIER);
}
