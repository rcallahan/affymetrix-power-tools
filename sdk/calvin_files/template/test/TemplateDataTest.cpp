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
#include "calvin_files/template/test/TemplateDataTest.h"
//
#include "calvin_files/template/src/TemplateData.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_template;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( TemplateDataTest );

void TemplateDataTest::setUp()
{
}

void TemplateDataTest::tearDown()
{
}

void TemplateDataTest::testCreation()
{
	TemplateData temp;
	CPPUNIT_ASSERT(1);
}

void TemplateDataTest::testproperty_Properies()
{
	string id = "file_id";
	TemplateData temp;
	temp.TemplateFileIdentifier() = id;
	CPPUNIT_ASSERT(temp.TemplateFileIdentifier() == id);

	temp.DataTypeIdentifier() = "affymetrix-calvin-template";
	CPPUNIT_ASSERT(temp.DataTypeIdentifier() == "affymetrix-calvin-template");
	temp.CreationDateTime() = L"now";
	CPPUNIT_ASSERT(temp.CreationDateTime() == L"now");
	temp.CreatedBy() = L"me";
	CPPUNIT_ASSERT(temp.CreatedBy() == L"me");
}

void TemplateDataTest::testmethod_Clear()
{
	string id = "file_id";
	TemplateData temp;
	temp.UserAttributes().resize(4);
	temp.DataTypeIdentifier() = "affymetrix-calvin-template";
	temp.CreationDateTime() = L"now";
	temp.CreatedBy() = L"me";
	temp.Clear();

	CPPUNIT_ASSERT(temp.TemplateFileIdentifier() == "");
	CPPUNIT_ASSERT(temp.UserAttributes().size() == 0);
	CPPUNIT_ASSERT(temp.DataTypeIdentifier() == "");
	CPPUNIT_ASSERT(temp.CreationDateTime() == L"");
	CPPUNIT_ASSERT(temp.CreatedBy() == L"");
}

void TemplateDataTest::testproperty_UserAttributes()
{
	TemplateData temp;
	ParameterNameValueDefaultRequiredTypeList &params = temp.UserAttributes();
	ParameterNameValueDefaultRequiredType vparam;

	vparam.SetName(L"user-att-name-1");
	vparam.SetValueAscii("user-att-value-1");
	vparam.SetDefaultValueAscii("user-att-defaultvalue-1");
	vparam.RequiredFlag() = false;
	vparam.ControlledVocabulary().clear() ;
	params.push_back(vparam);

	vparam.SetName(L"user-att-name-2");
	vparam.SetValueText(L"user-att-value-2");
	vparam.SetDefaultValueText(L"user-att-defaultvalue-2");
	vparam.RequiredFlag() = true;
	vparam.ControlledVocabulary().clear()  ;
	params.push_back(vparam);

	vparam.SetName(L"user-att-name-3");
	vparam.SetValueText(L"user-att-value-3");
	vparam.SetDefaultValueText(L"user-att-defaultvalue-3");
	vparam.RequiredFlag() = false;
	vparam.ControlledVocabulary().push_back(L"control1");
	vparam.ControlledVocabulary().push_back(L"control2");
	params.push_back(vparam);

	CPPUNIT_ASSERT( temp.UserAttributes().size() == 3 );

	ParameterNameValueDefaultRequiredTypeList::iterator it = temp.UserAttributes().begin();

	CPPUNIT_ASSERT( (*it).GetName() == L"user-att-name-1" );
	CPPUNIT_ASSERT( (*it).GetValueAscii() == "user-att-value-1" );
	CPPUNIT_ASSERT( (*it).GetDefaultValueAscii() == "user-att-defaultvalue-1" );
	CPPUNIT_ASSERT( (*it).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*it).RequiredFlag() == false);
	++it;
	CPPUNIT_ASSERT( (*it).GetName() == L"user-att-name-2" );
	CPPUNIT_ASSERT( (*it).GetValueText() == L"user-att-value-2" );
	CPPUNIT_ASSERT( (*it).GetDefaultValueText() == L"user-att-defaultvalue-2" );
	CPPUNIT_ASSERT( (*it).ControlledVocabulary().size() == 0 );
	CPPUNIT_ASSERT( (*it).RequiredFlag() == true);
	++it;
	CPPUNIT_ASSERT( (*it).GetName() == L"user-att-name-3" );
	CPPUNIT_ASSERT( (*it).GetValueText() == L"user-att-value-3" );
	CPPUNIT_ASSERT( (*it).GetDefaultValueText() == L"user-att-defaultvalue-3" );
	CPPUNIT_ASSERT( (*it).ControlledVocabulary().size() == 2 );
	CPPUNIT_ASSERT( (*it).RequiredFlag() == false);

	std::list<std::wstring>::iterator controllIt = (*it).ControlledVocabulary().begin();
	CPPUNIT_ASSERT( (*controllIt) == L"control1" );
	++controllIt;
	CPPUNIT_ASSERT( (*controllIt) == L"control2" );
	++it;
	CPPUNIT_ASSERT (it == temp.UserAttributes().end() );
}


