////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

/**
 * @file   SnpModelConverterTest.cpp
 * @author Mybrid Spalding
 * @date   Fri Jan 28 16:06:52 PST 2011
 * 
 * @brief  
 * 
 * 
 */

#ifndef _SNP_MODEL_CONVERTER_TEST_H_
#define _SNP_MODEL_CONVERTER_TEST_H_

#include <iostream>
#include <string>
#include <vector>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "chipstream/SnpModelConverter.h"
#include "file/TsvFile/TsvFile.h"
#include "util/CPPTest/Setup.h"

class SnpModelConverterTest : public CppUnit::TestFixture {

public:
  CPPUNIT_TEST_SUITE( SnpModelConverterTest );
  CPPUNIT_TEST( GetModelHeaders );
  CPPUNIT_TEST_SUITE_END();

  void GetModelHeaders();

};

CPPUNIT_TEST_SUITE_REGISTRATION( SnpModelConverterTest );


void SnpModelConverterTest::GetModelHeaders() {

  std::map< std::string, std::string * > headerMap;
  std::string key, value;
  affx::TsvFile tsv;
  tsv.m_optAutoTrim = true;

  // VALID HEADER TYPES only
  NEGATIVE_TEST(SnpModelConverter::GetModelHeaders( "input/snp_model_converter_values_unit_test.txt",
                                                    "invalid-type", headerMap),std::exception);

  // BAD FILE
  NEGATIVE_TEST(SnpModelConverter::GetModelHeaders( "input/snp_model_converter_values_unit_test_file_not_found.txt",
                                                    "none", headerMap),std::exception);

  // MISSING FILE TYPE
  NEGATIVE_TEST(CPPUNIT_ASSERT(SnpModelConverter::GetModelHeaders( "input/snp_model_converter_values_unit_test_missing_file_type.txt",
                                                                   "none", headerMap)),std::exception);
  
  // MISSING CHIP TYPE HEADER, chip type is UNKNOWN_TYPE if data is valid.
  POSITIVE_TEST(CPPUNIT_ASSERT(SnpModelConverter::GetModelHeaders( "input/snp_model_converter_values_unit_test_missing_chip_type.txt",
                                                                   "none", headerMap)));
  
  // MISSING ALGORITHM TYPE HEADER, algorithm derived from data column headers
  POSITIVE_TEST(CPPUNIT_ASSERT(SnpModelConverter::GetModelHeaders( "input/snp_model_converter_values_unit_test_missing_algorithm_type.txt",
                                                                   "none", headerMap)));
  
  // MULTI-VALUE for GTC chip-type test.
  POSITIVE_TEST(SnpModelConverter::GetModelHeaders( "input/snp_model_converter_values_unit_test.txt",
                                      "GTC4.1", headerMap, ","));
  CPPUNIT_ASSERT((tsv.openTable("expected/snp_model_converter_test_results1.txt") ==
                  affx::TSV_OK));
  
  tsv.headersBegin();
  while ( tsv.headersNext(key, value ) == affx::TSV_OK ){
    CPPUNIT_ASSERT(headerMap.count(key));
    if ( *(headerMap[key]) != value ) {
      Verbose::out(1, key + ":" + "generated value=" + *(headerMap[key]) + "!= gold value=" + value);
    }
    CPPUNIT_ASSERT(*(headerMap[key]) == value);
  }

  tsv.clear();

  // MULTI-VALUE different order test for GTC chip-type, headers are FIFO
  POSITIVE_TEST(SnpModelConverter::GetModelHeaders( "input/snp_model_converter_values_unit_test2.txt",
                                      "GTC4.1", headerMap, ","));
  CPPUNIT_ASSERT((tsv.openTable("expected/snp_model_converter_test_results2.txt") ==
                  affx::TSV_OK));
  
  tsv.headersBegin();
  while ( tsv.headersNext(key, value ) == affx::TSV_OK ){
    CPPUNIT_ASSERT(headerMap.count(key));
    if ( *(headerMap[key]) != value ) {
      Verbose::out(1, key + ":" + "generated value=" + *(headerMap[key]) + "!= gold value=" + value);
    }
    CPPUNIT_ASSERT(*(headerMap[key]) == value);
  }

  tsv.clear();
  
}


#endif 
