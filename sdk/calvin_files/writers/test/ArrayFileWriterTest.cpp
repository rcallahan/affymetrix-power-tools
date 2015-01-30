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
#include "calvin_files/writers/test/ArrayFileWriterTest.h"
//
#include "calvin_files/parsers/src/ArrayFileReader.h"
#include "calvin_files/writers/src/ArrayFileWriter.h"
//
#include <cmath>
#include <cstring>
#include <ostream>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_array;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( ArrayFileWriterTest );

void ArrayFileWriterTest::setUp()
{
}

void ArrayFileWriterTest::tearDown()
{
}

void ArrayFileWriterTest::testCreation()
{
	ArrayFileWriter writer;
	CPPUNIT_ASSERT(1);
}

void ArrayFileWriterTest::testmethod_WriteFile()
{
	ArrayData array1;

	array1.CreatedBy() = L"ljevon";
	array1.CreatedStep() = affymetrix_calvin_array::ScanningStep;
	array1.CreationDateTime() = L"datetime";
	array1.ArraySetFileIdentifier() = "id";
	array1.InitialProject() = L"proj1";

	array1.PhysicalArraysAttributes().resize(1);
	array1.PhysicalArraysAttributes()[0].Identifier() = "array1-id";
	array1.PhysicalArraysAttributes()[0].ArrayBarcode() = "123";
	array1.PhysicalArraysAttributes()[0].ArrayName() = "name";
	array1.PhysicalArraysAttributes()[0].Comment() = L"comment goes here";
	array1.PhysicalArraysAttributes()[0].CreatedBy() = L"me";
	array1.PhysicalArraysAttributes()[0].CreatedStep() = JobOrderServerStep;
	array1.PhysicalArraysAttributes()[0].CreationDateTime() = L"now";
	array1.PhysicalArraysAttributes()[0].MasterFile() = "master";
	array1.PhysicalArraysAttributes()[0].MasterFileId() = "456";
	array1.PhysicalArraysAttributes()[0].Media() = PlateOrStripMedia;
	array1.PhysicalArraysAttributes()[0].MediaRow() = 1;
	array1.PhysicalArraysAttributes()[0].MediaCol() = 2;
	array1.PhysicalArraysAttributes()[0].PatAssignment() = AffyBarcodeAssignment;
	array1.PhysicalArraysAttributes()[0].MediaFileGUID() = "mediaguid";
	array1.PhysicalArraysAttributes()[0].MediaFileName() = "mediafile";
	array1.PhysicalArraysAttributes()[0].LibraryPackageName() = "libpackage";

	ParameterNameValuePairVector &arrayParams = array1.PhysicalArraysAttributes()[0].Attributes();
	ParameterNameValuePair param;

	param.Name = L"array-att-name-1";
	param.Value = L"array-att-value-1";
	arrayParams.push_back(param);

	param.Name = L"array-att-name-2";
	param.Value = L"array-att-value-2";
	arrayParams.push_back(param);


	ParameterNameValueDefaultRequiredTypeList &userParams = array1.UserAttributes();
	ParameterNameValueDefaultRequiredType vparam;

	vparam.SetName(L"user-att-name-1") ;
	vparam.SetValueText(L"user-att-value-1");
	vparam.RequiredFlag() = false;
	vparam.HasDefault() = false;
	vparam.ValueType() = ParameterNameValueDefaultRequiredType::TextParameterType;
	vparam.ControlledVocabulary().clear() ;
	userParams.push_back(vparam);

	vparam.SetName(L"user-att-name-2") ;
	vparam.SetValueFloat(1.0);
	vparam.SetDefaultValueFloat(2.0);
	vparam.HasDefault() = true;
	vparam.RequiredFlag() = true;
	vparam.ValueType() = ParameterNameValueDefaultRequiredType::FloatParameterType;
	vparam.ControlledVocabulary().clear()  ;
	userParams.push_back(vparam);

	vparam.SetName(L"user-att-name-3") ;
	vparam.ControlMultiValues().push_back(L"one");
	vparam.ControlMultiValues().push_back(L"two");
	vparam.RequiredFlag() = false;
	vparam.HasDefault() = false;
	vparam.ValueType() = ParameterNameValueDefaultRequiredType::ControlMultiParameterType;
	vparam.ControlledVocabulary().push_back(L"one");
	vparam.ControlledVocabulary().push_back(L"two");
	vparam.ControlledVocabulary().push_back(L"three");
	userParams.push_back(vparam);


	// Write the file.
	const std::string TEST_FILE = "./test.file.array.xml";
	ArrayFileWriter writer;
	writer.Write(TEST_FILE, array1);

	// Read the file.
	ArrayFileReader reader;
	ArrayData array2;
	array2.Clear();
	reader.Read(TEST_FILE, array2);


	CPPUNIT_ASSERT(array2.CreatedBy() == L"ljevon");
	CPPUNIT_ASSERT(array2.CreatedStep() == affymetrix_calvin_array::ScanningStep);
	CPPUNIT_ASSERT(array2.CreationDateTime() == L"datetime");
	CPPUNIT_ASSERT(array2.ArraySetFileIdentifier() == "id");
	CPPUNIT_ASSERT(array2.InitialProject() == L"proj1");

	CPPUNIT_ASSERT( array2.DataTypeIdentifier() == ARRAY_SET_FILE_TYPE_IDENTIFIER );

	// Array attributes.
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes().size() == 1 );
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].Identifier() == "array1-id");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].ArrayBarcode() == "123");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].ArrayName() == "name");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].Comment() == L"comment goes here");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].CreatedBy() == L"me");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].CreatedStep() == JobOrderServerStep);
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].CreationDateTime() == L"now");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].MasterFile() == "master");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].MasterFileId() == "456");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].Media() == PlateOrStripMedia);
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].MediaRow() == 1);
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].MediaCol() == 2);
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].PatAssignment() == AffyBarcodeAssignment);
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].MediaFileGUID() == "mediaguid");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].MediaFileName() == "mediafile");
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].LibraryPackageName() == "libpackage");


	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].Attributes().size() == 2 );
	ParameterNameValuePairVector::iterator arrayAttIt = array2.PhysicalArraysAttributes()[0].Attributes().begin();
	CPPUNIT_ASSERT( (*arrayAttIt).Name == L"array-att-name-1" );
	CPPUNIT_ASSERT( (*arrayAttIt).Value == L"array-att-value-1" );
	++arrayAttIt;
	CPPUNIT_ASSERT( (*arrayAttIt).Name == L"array-att-name-2" );
	CPPUNIT_ASSERT( (*arrayAttIt).Value == L"array-att-value-2" );

	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].Attributes()[0].Name == L"array-att-name-1" );
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].Attributes()[0].Value == L"array-att-value-1" );
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].Attributes()[1].Name == L"array-att-name-2" );
	CPPUNIT_ASSERT( array2.PhysicalArraysAttributes()[0].Attributes()[1].Value == L"array-att-value-2" );



	// User attributes
	CPPUNIT_ASSERT( array2.UserAttributes().size() == 3 );

	ParameterNameValueDefaultRequiredTypeList::iterator userAttIt = array2.UserAttributes().begin();

	CPPUNIT_ASSERT( (*userAttIt).GetName() == L"user-att-name-1" );
	CPPUNIT_ASSERT( (*userAttIt).GetValueText() == L"user-att-value-1" );
	CPPUNIT_ASSERT( (*userAttIt).RequiredFlag() == false);
	CPPUNIT_ASSERT( (*userAttIt).ValueType() == ParameterNameValueDefaultRequiredType::TextParameterType);
	CPPUNIT_ASSERT( (*userAttIt).ControlledVocabulary().size() == 0 );
	++userAttIt;


	CPPUNIT_ASSERT( (*userAttIt).GetName() == L"user-att-name-2" );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( (*userAttIt).GetValueFloat(), 1, 0.0001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( (*userAttIt).GetDefaultValueFloat(), 2, 0.0001);
	CPPUNIT_ASSERT( (*userAttIt).RequiredFlag() == true);
	CPPUNIT_ASSERT( (*userAttIt).ValueType() == ParameterNameValueDefaultRequiredType::FloatParameterType);
	CPPUNIT_ASSERT( (*userAttIt).ControlledVocabulary().size() == 0 );
	++userAttIt;


	CPPUNIT_ASSERT( (*userAttIt).GetName() == L"user-att-name-3" );
	CPPUNIT_ASSERT( (*userAttIt).RequiredFlag() == false);
	CPPUNIT_ASSERT( (*userAttIt).ValueType() == ParameterNameValueDefaultRequiredType::ControlMultiParameterType);
	CPPUNIT_ASSERT( (*userAttIt).ControlMultiValues().size() == 2 );
	std::list<std::wstring>::iterator cIt;
	cIt = (*userAttIt).ControlMultiValues().begin();
	CPPUNIT_ASSERT( (*cIt) == L"one" );
	++cIt;
	CPPUNIT_ASSERT( (*cIt) == L"two" );
	CPPUNIT_ASSERT( (*userAttIt).ControlledVocabulary().size() == 3 );
	cIt = (*userAttIt).ControlledVocabulary().begin();
	CPPUNIT_ASSERT( (*cIt) == L"one" );
	++cIt;
	CPPUNIT_ASSERT( (*cIt) == L"two" );
	++cIt;
	CPPUNIT_ASSERT( (*cIt) == L"three" );

}
