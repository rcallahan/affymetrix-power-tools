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
#include "calvin_files/writers/test/TemplateFileWriterTest.h"
//
#include "calvin_files/parsers/src/TemplateFileReader.h"
#include "calvin_files/template/src/TemplateData.h"
#include "calvin_files/template/src/TemplateId.h"
#include "calvin_files/writers/src/TemplateFileWriter.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_template;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( TemplateFileWriterTest );

void TemplateFileWriterTest::setUp()
{
}

void TemplateFileWriterTest::tearDown()
{
}

void TemplateFileWriterTest::testCreation()
{
	TemplateFileWriter writer;
	CPPUNIT_ASSERT(1);
}

void TemplateFileWriterTest::testmethod_WriteFile()
{
	TemplateData templ;

	templ.CreatedBy() = L"ljevon";
	templ.CreationDateTime() = L"datetime";
	templ.TemplateFileIdentifier() = "id";

	ParameterNameValueDefaultRequiredTypeList &userParams = templ.UserAttributes();
	ParameterNameValueDefaultRequiredType vparam;

	vparam.SetName(L"user-att-name-1") ;
	vparam.RequiredFlag() = false;
	vparam.HasDefault() = false;
	vparam.ValueType() = ParameterNameValueDefaultRequiredType::TextParameterType;
	vparam.ControlledVocabulary().clear() ;
	userParams.push_back(vparam);

	vparam.SetName(L"user-att-name-2") ;
	vparam.SetDefaultValueFloat(2.0);
	vparam.HasDefault() = true;
	vparam.RequiredFlag() = true;
	vparam.ValueType() = ParameterNameValueDefaultRequiredType::FloatParameterType;
	vparam.ControlledVocabulary().clear()  ;
	userParams.push_back(vparam);

	vparam.SetName(L"user-att-name-3") ;
	vparam.RequiredFlag() = false;
	vparam.HasDefault() = false;
	vparam.ValueType() = ParameterNameValueDefaultRequiredType::ControlMultiParameterType;
	vparam.ControlledVocabulary().push_back(L"one");
	vparam.ControlledVocabulary().push_back(L"two");
	vparam.ControlledVocabulary().push_back(L"three");
	userParams.push_back(vparam);


	// Write the file.
	const std::string TEST_FILE = "./test.file.templ";
	TemplateFileWriter writer;
	writer.Write(TEST_FILE, templ);

	// Read the file.
	TemplateFileReader reader;
	templ.Clear();
	reader.Read(TEST_FILE, templ);


	CPPUNIT_ASSERT(templ.CreatedBy() == L"ljevon");
	CPPUNIT_ASSERT(templ.CreationDateTime() == L"datetime");
	CPPUNIT_ASSERT(templ.TemplateFileIdentifier() == "id");

	CPPUNIT_ASSERT( templ.DataTypeIdentifier() == TEMPLATE_FILE_TYPE_IDENTIFIER );

	// User attributes
	CPPUNIT_ASSERT( templ.UserAttributes().size() == 3 );

	ParameterNameValueDefaultRequiredTypeList::iterator userAttIt = templ.UserAttributes().begin();

	CPPUNIT_ASSERT( (*userAttIt).GetName() == L"user-att-name-1" );
	CPPUNIT_ASSERT( (*userAttIt).GetValueText() == L"" );
	CPPUNIT_ASSERT( (*userAttIt).RequiredFlag() == false);
	CPPUNIT_ASSERT( (*userAttIt).ValueType() == ParameterNameValueDefaultRequiredType::TextParameterType);
	CPPUNIT_ASSERT( (*userAttIt).ControlledVocabulary().size() == 0 );
	++userAttIt;


	CPPUNIT_ASSERT( (*userAttIt).GetName() == L"user-att-name-2" );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( (*userAttIt).GetDefaultValueFloat(), 2, 0.0001);
	CPPUNIT_ASSERT( (*userAttIt).RequiredFlag() == true);
	CPPUNIT_ASSERT( (*userAttIt).ValueType() == ParameterNameValueDefaultRequiredType::FloatParameterType);
	CPPUNIT_ASSERT( (*userAttIt).ControlledVocabulary().size() == 0 );
	++userAttIt;


	CPPUNIT_ASSERT( (*userAttIt).GetName() == L"user-att-name-3" );
	CPPUNIT_ASSERT( (*userAttIt).RequiredFlag() == false);
	CPPUNIT_ASSERT( (*userAttIt).ValueType() == ParameterNameValueDefaultRequiredType::ControlMultiParameterType);
	CPPUNIT_ASSERT( (*userAttIt).ControlMultiValues().size() == 0 );
	std::list<std::wstring>::iterator cIt;
	CPPUNIT_ASSERT( (*userAttIt).ControlledVocabulary().size() == 3 );
	cIt = (*userAttIt).ControlledVocabulary().begin();
	CPPUNIT_ASSERT( (*cIt) == L"one" );
	++cIt;
	CPPUNIT_ASSERT( (*cIt) == L"two" );
	++cIt;
	CPPUNIT_ASSERT( (*cIt) == L"three" );
}
