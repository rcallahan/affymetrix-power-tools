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


#include "calvin_files/fusion/test/FusionArrayFileTest.h"
//
#include "calvin_files/fusion/src/FusionArrayFileReader.h"
//

using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_array;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_parameter;
using namespace std;


#define INVALID_FILE "../../parsers/data/test.file.text"
#define CALVIN_FILE "../../parsers/data/ArraySetFile.xml"
#define MAGE_FILE "../../../file/CPPTest/data/sample_and_exp_template.xml"
#define EXP_FILE "../../../file/CPPTest/data/full.EXP"

CPPUNIT_TEST_SUITE_REGISTRATION( FusionArrayFileReaderTest );

void FusionArrayFileReaderTest::setUp()
{
}

void FusionArrayFileReaderTest::tearDown()
{
}

/*! Test the ability to create the ArrayFileReader class. */
void FusionArrayFileReaderTest::testCreation()
{
	FusionArrayFileReader *reader = new FusionArrayFileReader;
	CPPUNIT_ASSERT(reader != NULL);
	delete reader;
}

/*! Test the handling of attempting to read an invalid array file. */
void FusionArrayFileReaderTest::testmethod_Read_invalid_file()
{
	FusionArrayFileReader reader;
	ArrayData data;
	CPPUNIT_ASSERT( reader.Read(INVALID_FILE, data) == false);
	CPPUNIT_ASSERT( data.ArraySetFileIdentifier() == "");
	CPPUNIT_ASSERT( data.DataTypeIdentifier() == "");
	CPPUNIT_ASSERT( data.PhysicalArraysAttributes().size() == 0);
	CPPUNIT_ASSERT( data.UserAttributes().size() == 0);
}

/*! Test the handling of attempting to read an array file that does not exist. */
void FusionArrayFileReaderTest::testmethod_Read_missing_file()
{
	FusionArrayFileReader reader;
	ArrayData data;
	CPPUNIT_ASSERT( reader.Read("no_file", data) == false);
	CPPUNIT_ASSERT( data.ArraySetFileIdentifier() == "");
	CPPUNIT_ASSERT( data.DataTypeIdentifier() == "");
	CPPUNIT_ASSERT( data.PhysicalArraysAttributes().size() == 0);
	CPPUNIT_ASSERT( data.UserAttributes().size() == 0);
}

/*! Test the ability to read a valid Calvin array file. */
void FusionArrayFileReaderTest::testmethod_Read_calvin_array_file()
{
	FusionArrayFileReader reader;
	ArrayData array;
	CPPUNIT_ASSERT( reader.Read(CALVIN_FILE, array) == true);

	CPPUNIT_ASSERT( array.DataTypeIdentifier() == ARRAY_SET_FILE_TYPE_IDENTIFIER );
	CPPUNIT_ASSERT(	array.ArraySetFileIdentifier() == std::string("432-432-432-432"));
	CPPUNIT_ASSERT( array.CreatedStep() == ArrayRegistrationStep);
	CPPUNIT_ASSERT( array.InitialProject() == L"my_project");
	CPPUNIT_ASSERT( array.CreationDateTime() == L"8/12/2005 9:00AM");
	CPPUNIT_ASSERT( array.CreatedBy() == L"ljevon");

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

/*! Test the ability to read an Data Transfer Tool formated array file. */
void FusionArrayFileReaderTest::testmethod_Read_dtt_array_file()
{
	FusionArrayFileReader reader;
	ArrayData data;
	CPPUNIT_ASSERT( reader.Read(MAGE_FILE, data) == true);

	CPPUNIT_ASSERT(	data.ArraySetFileIdentifier() == std::string(""));
	CPPUNIT_ASSERT( data.PhysicalArraysAttributes().size() == 1);
	CPPUNIT_ASSERT( data.PhysicalArraysAttributes()[0].Identifier() == "");
	const ParameterNameValuePairVector &arrayParams = data.PhysicalArraysAttributes()[0].Attributes();
	CPPUNIT_ASSERT( arrayParams.size() == 1 );
	ParameterNameValuePairVector::const_iterator arrayIt = arrayParams.begin();
	CPPUNIT_ASSERT( (*arrayIt).Name == L"Probe Array Type" );
	CPPUNIT_ASSERT( (*arrayIt).Value == L"Test3" );


	ParameterNameValueDefaultRequiredTypeList &userParams = data.UserAttributes();
	CPPUNIT_ASSERT( userParams.size() == 33 );
	ParameterNameValueDefaultRequiredTypeList::iterator userIt = userParams.begin();
	CPPUNIT_ASSERT( (*userIt).GetName() == L"GCOS Sample Name" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"Sample_with_Sample_Template" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Sample Template Name" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"MIAME Sample Information" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	userIt = userParams.end();
	--userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Grid Alignment Passed" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
}

/*! Test the ability to read an Micro Array Suit formated array file (EXP).  */
void FusionArrayFileReaderTest::testmethod_Read_mas_array_file()
{
	FusionArrayFileReader reader;
	ArrayData data;
	CPPUNIT_ASSERT( reader.Read(EXP_FILE, data) == true);

	CPPUNIT_ASSERT(	data.ArraySetFileIdentifier() == std::string(""));
	CPPUNIT_ASSERT( data.PhysicalArraysAttributes().size() == 1);
	CPPUNIT_ASSERT( data.PhysicalArraysAttributes()[0].Identifier() == "");
	const ParameterNameValuePairVector &arrayParams = data.PhysicalArraysAttributes()[0].Attributes();
	CPPUNIT_ASSERT( arrayParams.size() == 1 );
	ParameterNameValuePairVector::const_iterator arrayIt = arrayParams.begin();
	CPPUNIT_ASSERT( (*arrayIt).Name == L"Probe Array Type" );
	CPPUNIT_ASSERT( (*arrayIt).Value == L"Test3" );

	ParameterNameValueDefaultRequiredTypeList &userParams = data.UserAttributes();
	CPPUNIT_ASSERT( userParams.size() == 8 );
	ParameterNameValueDefaultRequiredTypeList::iterator userIt = userParams.begin();
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Chip Lot" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"123" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Operator" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"ljevon" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Sample Type" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"s_type" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Description" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"Demo data" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Project" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"s_proj" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Comments" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"None" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Solution Type" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"sol type" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
	CPPUNIT_ASSERT( (*userIt).GetName() == L"Solution Lot" );
	CPPUNIT_ASSERT( (*userIt).GetValueText() == L"123123" );
	CPPUNIT_ASSERT( (*userIt).ControlledVocabulary().size() == 0 );
	++userIt;
}
