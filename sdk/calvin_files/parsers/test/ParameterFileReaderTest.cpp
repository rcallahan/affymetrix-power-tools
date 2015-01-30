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
#include "calvin_files/parsers/src/ParameterFileReader.h"
//
#include <cppunit/extensions/HelperMacros.h>
//
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace std;

class ParameterFileReaderTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( ParameterFileReaderTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testmethod_Read );
	CPPUNIT_TEST ( testmethod_Read_when_file_does_not_exist );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testmethod_Read();
	void testmethod_Read_when_file_does_not_exist();
};


CPPUNIT_TEST_SUITE_REGISTRATION( ParameterFileReaderTest );

#define TEST_NON_EXISTANT_FILE "../data/file_does_not_exist";
#define TEST_DATA_FILE "../data/parameters.xml"

void ParameterFileReaderTest::setUp()
{
}

void ParameterFileReaderTest::tearDown()
{
}

void ParameterFileReaderTest::testCreation()
{
	ParameterFileReader reader;
	CPPUNIT_ASSERT(1);
}

void ParameterFileReaderTest::testmethod_Read()
{
	ParameterFileReader reader;
	ParameterFileData paramData;
	std::string name = TEST_DATA_FILE;
	CPPUNIT_ASSERT(reader.Read(name, paramData) == true);

    CPPUNIT_ASSERT(paramData.ParameterFileAttributes().company == L"Affymetrix");
    CPPUNIT_ASSERT(paramData.ParameterFileAttributes().userName == L"Srinivas");
    CPPUNIT_ASSERT(paramData.ParameterFileAttributes().contentVersion == L"1");

    CPPUNIT_ASSERT(paramData.ImplementationAttributes().name == L"Percentile Cell Average");
    CPPUNIT_ASSERT(paramData.ImplementationAttributes().version == L"0.1");
    CPPUNIT_ASSERT(paramData.ImplementationAttributes().executableFileName == L"percentile.exe");
    CPPUNIT_ASSERT(paramData.ImplementationAttributes().description == L"Percentile Average feature extraction algorithm for generating CEL files");

    CPPUNIT_ASSERT(paramData.Parameters().size() == 5);

    ParameterTypeList::iterator it;
    
    it = paramData.Parameters().begin();

    CPPUNIT_ASSERT(it->name==L"Percentile");
    CPPUNIT_ASSERT(it->displayName==L"Percentile");
    CPPUNIT_ASSERT(it->category==L"Percentile Parameters");
    CPPUNIT_ASSERT(it->isEditable==L"false");
    CPPUNIT_ASSERT(it->type==L"Int32");
    CPPUNIT_ASSERT(it->currentValue==L"75");
    CPPUNIT_ASSERT(it->minValue==L"0");
    CPPUNIT_ASSERT(it->maxValue==L"100");
    CPPUNIT_ASSERT(it->defaultValue==L"75");
    CPPUNIT_ASSERT(it->precision==L"0");
    CPPUNIT_ASSERT(it->maxLength==L"0");
    CPPUNIT_ASSERT(it->description==L"Percentile to use for cell intensity");
    CPPUNIT_ASSERT(it->index==L"-1");

    ++it;

    CPPUNIT_ASSERT(it->name==L"CellMargin");
    CPPUNIT_ASSERT(it->displayName==L"CellMargin");
    CPPUNIT_ASSERT(it->category==L"Percentile Parameters");
    CPPUNIT_ASSERT(it->isEditable==L"false");
    CPPUNIT_ASSERT(it->type==L"Int32");
    CPPUNIT_ASSERT(it->currentValue==L"2");
    CPPUNIT_ASSERT(it->minValue==L"0");
    CPPUNIT_ASSERT(it->maxValue==L"10");
    CPPUNIT_ASSERT(it->defaultValue==L"2");
    CPPUNIT_ASSERT(it->precision==L"0");
    CPPUNIT_ASSERT(it->maxLength==L"0");
    CPPUNIT_ASSERT(it->description==L"Cell Margin to use for cell intensity");
    CPPUNIT_ASSERT(it->index==L"-1");

    ++it;

    CPPUNIT_ASSERT(it->name==L"CellFileVersion");
    CPPUNIT_ASSERT(it->displayName==L"CellFileVersion");
    CPPUNIT_ASSERT(it->category==L"Percentile Parameters");
    CPPUNIT_ASSERT(it->isEditable==L"false");
    CPPUNIT_ASSERT(it->type==L"Int32");
    CPPUNIT_ASSERT(it->currentValue==L"3");
    CPPUNIT_ASSERT(it->minValue==L"3");
    CPPUNIT_ASSERT(it->maxValue==L"3");
    CPPUNIT_ASSERT(it->defaultValue==L"3");
    CPPUNIT_ASSERT(it->precision==L"0");
    CPPUNIT_ASSERT(it->maxLength==L"0");
    CPPUNIT_ASSERT(it->description==L"Version of Cell file to write");
    CPPUNIT_ASSERT(it->index==L"-1");

    ++it;

    CPPUNIT_ASSERT(it->name==L"OutlierLow");
    CPPUNIT_ASSERT(it->displayName==L"OutlierLow");
    CPPUNIT_ASSERT(it->category==L"Outlier Parameters");
    CPPUNIT_ASSERT(it->isEditable==L"false");
    CPPUNIT_ASSERT(it->type==L"Float");
    CPPUNIT_ASSERT(it->currentValue==L"1.004");
    CPPUNIT_ASSERT(it->minValue==L"1.004");
    CPPUNIT_ASSERT(it->maxValue==L"1.004");
    CPPUNIT_ASSERT(it->defaultValue==L"1.004");
    CPPUNIT_ASSERT(it->precision==L"3");
    CPPUNIT_ASSERT(it->maxLength==L"0");
    CPPUNIT_ASSERT(it->description==L"Outlier threshold low");
    CPPUNIT_ASSERT(it->index==L"-1");

    ++it;

    CPPUNIT_ASSERT(it->name==L"OutlierHigh");
    CPPUNIT_ASSERT(it->displayName==L"OutlierHigh");
    CPPUNIT_ASSERT(it->category==L"Outlier Parameters");
    CPPUNIT_ASSERT(it->isEditable==L"false");
    CPPUNIT_ASSERT(it->type==L"Float");
    CPPUNIT_ASSERT(it->currentValue==L"1.5");
    CPPUNIT_ASSERT(it->minValue==L"1.5");
    CPPUNIT_ASSERT(it->maxValue==L"1.5");
    CPPUNIT_ASSERT(it->defaultValue==L"1.5");
    CPPUNIT_ASSERT(it->precision==L"1");
    CPPUNIT_ASSERT(it->maxLength==L"0");
    CPPUNIT_ASSERT(it->description==L"Outlier threshold high");
    CPPUNIT_ASSERT(it->index==L"-1");

}

void ParameterFileReaderTest::testmethod_Read_when_file_does_not_exist()
{
	ParameterFileReader reader;
	ParameterFileData paramData;
	std::string name = TEST_NON_EXISTANT_FILE;
	CPPUNIT_ASSERT(reader.Read(name, paramData) == false);
}
