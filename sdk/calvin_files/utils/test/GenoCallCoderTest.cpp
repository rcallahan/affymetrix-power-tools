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


/**
 * @file   GenoCallCoderTest.cpp
 * @author Ray Wheeler
 * @date   Mon Jun 16 13:54:36 PDT 2008
 * 
 * @brief  Test genotype call codes encoder/decoder
 */

//
#include "calvin_files/utils/src/GenoCallCoder.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <string>
//

using namespace std;

/**
 * @class RMATest
 * @brief cppunit class for testing conversion functions.
 */
class GenoCallCoderTest : public CPPUNIT_NS::TestFixture {

  CPPUNIT_TEST_SUITE( GenoCallCoderTest );
  CPPUNIT_TEST( testDecoding_v0 );
//   CPPUNIT_TEST( testDecoding_delimit_v0 );
  CPPUNIT_TEST( testDecoding_v1 );
  CPPUNIT_TEST( testDecoding_delimit_v1 );
  CPPUNIT_TEST( testEncoding_v0 );
//   CPPUNIT_TEST( testEncoding_delimit_v0 );
  CPPUNIT_TEST( testEncoding_v1 );
  CPPUNIT_TEST( testEncoding_delimit_v1 );
  CPPUNIT_TEST( testAbstractAlelleIntegerEncoding_v1 );
  CPPUNIT_TEST( testDecodeVec_v0 );
  CPPUNIT_TEST( testDecodeVec_v1 );
  CPPUNIT_TEST( testGetDecodeVec_v0 );
  CPPUNIT_TEST( testGetDecodeVec_v1 );
  CPPUNIT_TEST( testIsValidCode_v0 );
  CPPUNIT_TEST( testIsValidCode_v1 );
  CPPUNIT_TEST( testIsValidAllele_v0 );
//   CPPUNIT_TEST( testIsValidAllele_delimit_v0 );
  CPPUNIT_TEST( testIsValidAllele_v1 );
  CPPUNIT_TEST( testIsValidAllele_delimit_v1 );
  CPPUNIT_TEST( testIsHom_Het_v1 );
  CPPUNIT_TEST( testAbstractAlleleToReferenceAllele_delimit_v1 );
  CPPUNIT_TEST( testAbstractAlleleToReportAllele_delimit_v1 );
  CPPUNIT_TEST( testReferenceAlleleToReportAllele_delimit_v1 );
  CPPUNIT_TEST( testGenotypeCallNumToReferenceAllele_delimit_v1 );
  CPPUNIT_TEST( testGenotypeCallNumToReportAllele_delimit_v1 );
  CPPUNIT_TEST( testGetProbesetIds_delimit_v1 );
  CPPUNIT_TEST( testGetValidGenotypeCallCodes_delimit_v1 );
//   CPPUNIT_TEST( testFullCsvFile );
  CPPUNIT_TEST_SUITE_END();

public:

  
  vector<unsigned char> decodeTestVec_v0;
  vector<string> encodeTestVec_v0;
  //  string decodeVec_v1[];

  void setUp();
  void testDecoding_v0();
//   void testDecoding_delimit_v0();
  void testEncoding_v0();
//   void testEncoding_delimit_v0();
  void testDecoding_v1();
  void testDecoding_delimit_v1();
  void testEncoding_v1();
  void testEncoding_delimit_v1();
  void testAbstractAlelleIntegerEncoding_v1();
  void testDecodeVec_v0();
  void testDecodeVec_v1();
  void testEncodeVec_v0();
  void testEncodeVec_v1();
  void testGetDecodeVec_v0();  
  void testGetDecodeVec_v1();
  void testIsValidCode_v0();
  void testIsValidCode_v1();
  void testIsValidAllele_v0();
//   void testIsValidAllele_delimit_v0();
  void testIsValidAllele_v1();
  void testIsValidAllele_delimit_v1();
  void testIsHom_Het_v1();
  void testAbstractAlleleToReferenceAllele_delimit_v1();
  void testAbstractAlleleToReportAllele_delimit_v1();
  void testReferenceAlleleToReportAllele_delimit_v1();
  void testGenotypeCallNumToReferenceAllele_delimit_v1();
  void testGenotypeCallNumToReportAllele_delimit_v1();
  void testGetValidGenotypeCallCodes_delimit_v1();
  void testGetProbesetIds_delimit_v1();
//   void testFullCsvFile();
  
  vector<string> v0_abs_allele_string_vec;
  vector<string> v0_abs_allele_string_delimit_vec;
  vector<string> v1_abs_allele_string_vec;
  vector<string> v1_abs_allele_string_delimit_vec;
  vector<string> v1_ref_allele_string_delimit_vec;
  vector<string> v1_report_allele_string_delimit_vec;
  vector<unsigned char> v0_code_vec;
  vector<unsigned char> v1_code_vec;
  vector<bool> v1_code_hom;
  vector<bool> v1_code_het;
  vector<string> probeset_ids;
  vector<string> abstract_alleles;
  vector<string> reference_alleles;
  vector<string> report_alleles;
  vector<vector<unsigned char> > valid_genotype_call_codes;
  vector<vector<string> > reference_geno_calls_test;
  vector<vector<string> > report_geno_calls_test;
  string probeset_annot_csv_file;
  
};

CPPUNIT_TEST_SUITE_REGISTRATION( GenoCallCoderTest );

void GenoCallCoderTest::setUp() {

  v0_code_vec.push_back(6);
  v0_code_vec.push_back(7);
  v0_code_vec.push_back(8);
  v0_code_vec.push_back(9);
  v0_code_vec.push_back(10);
  v0_code_vec.push_back(11);
  
  v0_abs_allele_string_vec.push_back(std::string("AA"));	   
  v0_abs_allele_string_vec.push_back(string("BB"));	   
  v0_abs_allele_string_vec.push_back(string("AB"));	   
  v0_abs_allele_string_vec.push_back(string("AB_A"));	   
  v0_abs_allele_string_vec.push_back(string("AB_B"));	   
  v0_abs_allele_string_vec.push_back(string("NN"));
  
//   v0_abs_allele_string_delimit_vec.push_back(string("A/A"));	   
//   v0_abs_allele_string_delimit_vec.push_back(string("B/B"));	   
//   v0_abs_allele_string_delimit_vec.push_back(string("A/B"));	   
//   v0_abs_allele_string_delimit_vec.push_back(string("NN"));
  
  v1_code_vec.push_back(16);
  v1_code_vec.push_back(17);
  v1_code_vec.push_back(23);
  v1_code_vec.push_back(24);
  v1_code_vec.push_back(25);
  v1_code_vec.push_back(26);
  v1_code_vec.push_back(31);
  v1_code_vec.push_back(32);
  v1_code_vec.push_back(41);
  v1_code_vec.push_back(52);
  v1_code_vec.push_back(61);
  v1_code_vec.push_back(110);
  v1_code_vec.push_back(135);
  v1_code_vec.push_back(136);
  v1_code_vec.push_back(254);
  v1_code_vec.push_back(255);

  v1_abs_allele_string_vec.push_back(string("NoCall"));
  v1_abs_allele_string_vec.push_back(string("ZeroCopyNumber")); 
  //  v1_abs_allele_string_vec.push_back(string("0")); 
  v1_abs_allele_string_vec.push_back(string("F")); 
  v1_abs_allele_string_vec.push_back(string("n")); 
  v1_abs_allele_string_vec.push_back(string("AA")); 
  v1_abs_allele_string_vec.push_back(string("AB")); 
  v1_abs_allele_string_vec.push_back(string("An")); 
  v1_abs_allele_string_vec.push_back(string("BB")); 
  v1_abs_allele_string_vec.push_back(string("CF")); 
  v1_abs_allele_string_vec.push_back(string("nn"));
  v1_abs_allele_string_vec.push_back(string("ABC"));
  v1_abs_allele_string_vec.push_back(string("CDn")); 
  v1_abs_allele_string_vec.push_back(string("Fnn")); 
  v1_abs_allele_string_vec.push_back(string("nnn")); 
  v1_abs_allele_string_vec.push_back(string("PossibleRareAllele")); 
  v1_abs_allele_string_vec.push_back(string("NotAvailable")); 

  v1_abs_allele_string_delimit_vec.push_back(string("NoCall"));
  v1_abs_allele_string_delimit_vec.push_back(string("ZeroCopyNumber")); 
  //v1_abs_allele_string_delimit_vec.push_back(string("0")); 
  v1_abs_allele_string_delimit_vec.push_back(string("F")); 
  v1_abs_allele_string_delimit_vec.push_back(string("n")); 
  v1_abs_allele_string_delimit_vec.push_back(string("A/A")); 
  v1_abs_allele_string_delimit_vec.push_back(string("A/B")); 
  v1_abs_allele_string_delimit_vec.push_back(string("A/n")); 
  v1_abs_allele_string_delimit_vec.push_back(string("B/B")); 
  v1_abs_allele_string_delimit_vec.push_back(string("C/F")); 
  v1_abs_allele_string_delimit_vec.push_back(string("n/n"));
  v1_abs_allele_string_delimit_vec.push_back(string("A/B/C"));
  v1_abs_allele_string_delimit_vec.push_back(string("C/D/n")); 
  v1_abs_allele_string_delimit_vec.push_back(string("F/n/n")); 
  v1_abs_allele_string_delimit_vec.push_back(string("n/n/n")); 
  v1_abs_allele_string_delimit_vec.push_back(string("PossibleRareAllele")); 
  v1_abs_allele_string_delimit_vec.push_back(string("NotAvailable")); 


  v1_code_hom.push_back(0);
  v1_code_hom.push_back(1);
  v1_code_hom.push_back(1);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(1);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(1);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(0);
  v1_code_hom.push_back(0);

  v1_code_het.push_back(0);
  v1_code_het.push_back(0);
  v1_code_het.push_back(0);
  v1_code_het.push_back(0);
  v1_code_het.push_back(0);
  v1_code_het.push_back(1);
  v1_code_het.push_back(0);
  v1_code_het.push_back(0);
  v1_code_het.push_back(1);
  v1_code_het.push_back(0);
  v1_code_het.push_back(1);
  v1_code_het.push_back(0);
  v1_code_het.push_back(0);
  v1_code_het.push_back(0);
  v1_code_het.push_back(0);
  v1_code_het.push_back(0);

  probeset_ids.push_back(string("DMET3B10001"));
  probeset_ids.push_back(string("DMET3B10126"));
  probeset_ids.push_back(string("DMET3B12248"));
  probeset_ids.push_back(string("DMET3B12456"));
  probeset_ids.push_back(string("DMET3B13024"));
  probeset_ids.push_back(string("DMET3B20494"));

  abstract_alleles.push_back(string("A/B"));
  abstract_alleles.push_back(string("A/B/C/D"));
  abstract_alleles.push_back(string("A/B/C"));
  abstract_alleles.push_back(string("A/B/C/D/E/F"));
  abstract_alleles.push_back(string("A/B"));
  abstract_alleles.push_back(string("B/A"));

  reference_alleles.push_back(string("C/G"));
  reference_alleles.push_back(string("A/C/G/T"));
  reference_alleles.push_back(string("A/AGTGGGCACA/G"));
  reference_alleles.push_back(string("-/A/C/G/GC/T"));
  reference_alleles.push_back(string("TATATATATATA/TATATATATATATA"));
  reference_alleles.push_back(string("C/T"));

  report_alleles.push_back(string("C/G"));
  report_alleles.push_back(string("T/G/C/A"));
  report_alleles.push_back(string("T/TGTGCCCACT/C"));
  report_alleles.push_back(string("-/A/C/G/GC/T"));
  report_alleles.push_back(string("(TA)5or6/(TA)7or8"));
  report_alleles.push_back(string("C/T"));

  reference_geno_calls_test.resize(probeset_ids.size());
  reference_geno_calls_test[0].push_back("ZeroCopyNumber");
  reference_geno_calls_test[0].push_back("C");
  reference_geno_calls_test[0].push_back("G");
  reference_geno_calls_test[0].push_back("C/C");
  reference_geno_calls_test[0].push_back("C/G");
  reference_geno_calls_test[0].push_back("G/G");
  reference_geno_calls_test[1].push_back("ZeroCopyNumber");
  reference_geno_calls_test[1].push_back("A");
  reference_geno_calls_test[1].push_back("C");
  reference_geno_calls_test[1].push_back("G");
  reference_geno_calls_test[1].push_back("T");
  reference_geno_calls_test[1].push_back("A/A");
  reference_geno_calls_test[1].push_back("A/C");
  reference_geno_calls_test[1].push_back("A/G");
  reference_geno_calls_test[1].push_back("A/T");
  reference_geno_calls_test[1].push_back("C/C");
  reference_geno_calls_test[1].push_back("C/G");
  reference_geno_calls_test[1].push_back("C/T");
  reference_geno_calls_test[1].push_back("G/G");
  reference_geno_calls_test[1].push_back("G/T");
  reference_geno_calls_test[1].push_back("T/T");

  report_geno_calls_test.resize(probeset_ids.size());
  report_geno_calls_test[0].push_back("ZeroCopyNumber");
  report_geno_calls_test[0].push_back("C");
  report_geno_calls_test[0].push_back("G");
  report_geno_calls_test[0].push_back("C/C");
  report_geno_calls_test[0].push_back("C/G");
  report_geno_calls_test[0].push_back("G/G");
  report_geno_calls_test[1].push_back("ZeroCopyNumber");
  report_geno_calls_test[1].push_back("T");
  report_geno_calls_test[1].push_back("G");
  report_geno_calls_test[1].push_back("C");
  report_geno_calls_test[1].push_back("A");
  report_geno_calls_test[1].push_back("T/T");
  report_geno_calls_test[1].push_back("T/G");
  report_geno_calls_test[1].push_back("T/C");
  report_geno_calls_test[1].push_back("T/A");
  report_geno_calls_test[1].push_back("G/G");
  report_geno_calls_test[1].push_back("G/C");
  report_geno_calls_test[1].push_back("G/A");
  report_geno_calls_test[1].push_back("C/C");
  report_geno_calls_test[1].push_back("C/A");
  report_geno_calls_test[1].push_back("A/A");

  valid_genotype_call_codes.resize(probeset_ids.size());
  valid_genotype_call_codes[0].push_back(17);
  valid_genotype_call_codes[0].push_back(18);
  valid_genotype_call_codes[0].push_back(19);
  valid_genotype_call_codes[0].push_back(25);
  valid_genotype_call_codes[0].push_back(26);
  valid_genotype_call_codes[0].push_back(32);
  valid_genotype_call_codes[1].push_back(17);
  valid_genotype_call_codes[1].push_back(18);
  valid_genotype_call_codes[1].push_back(19);
  valid_genotype_call_codes[1].push_back(20);
  valid_genotype_call_codes[1].push_back(21);
  valid_genotype_call_codes[1].push_back(25);
  valid_genotype_call_codes[1].push_back(26);
  valid_genotype_call_codes[1].push_back(27);
  valid_genotype_call_codes[1].push_back(28);
  valid_genotype_call_codes[1].push_back(32);
  valid_genotype_call_codes[1].push_back(33);
  valid_genotype_call_codes[1].push_back(34);
  valid_genotype_call_codes[1].push_back(38);
  valid_genotype_call_codes[1].push_back(39);
  valid_genotype_call_codes[1].push_back(43);
  valid_genotype_call_codes[2].push_back(17);
  valid_genotype_call_codes[2].push_back(18);
  valid_genotype_call_codes[2].push_back(19);
  valid_genotype_call_codes[2].push_back(20);
  valid_genotype_call_codes[2].push_back(25);
  valid_genotype_call_codes[2].push_back(26);
  valid_genotype_call_codes[2].push_back(27);
  valid_genotype_call_codes[2].push_back(32);
  valid_genotype_call_codes[2].push_back(33);
  valid_genotype_call_codes[2].push_back(38);
  valid_genotype_call_codes[3].push_back(17);
  valid_genotype_call_codes[3].push_back(18);
  valid_genotype_call_codes[3].push_back(19);
  valid_genotype_call_codes[3].push_back(20);
  valid_genotype_call_codes[3].push_back(21);
  valid_genotype_call_codes[3].push_back(22);
  valid_genotype_call_codes[3].push_back(23);
  valid_genotype_call_codes[3].push_back(25);
  valid_genotype_call_codes[3].push_back(26);
  valid_genotype_call_codes[3].push_back(27);
  valid_genotype_call_codes[3].push_back(28);
  valid_genotype_call_codes[3].push_back(29);
  valid_genotype_call_codes[3].push_back(30);
  valid_genotype_call_codes[3].push_back(32);
  valid_genotype_call_codes[3].push_back(33);
  valid_genotype_call_codes[3].push_back(34);
  valid_genotype_call_codes[3].push_back(35);
  valid_genotype_call_codes[3].push_back(36);
  valid_genotype_call_codes[3].push_back(38);
  valid_genotype_call_codes[3].push_back(39);
  valid_genotype_call_codes[3].push_back(40);
  valid_genotype_call_codes[3].push_back(41);
  valid_genotype_call_codes[3].push_back(43);
  valid_genotype_call_codes[3].push_back(44);
  valid_genotype_call_codes[3].push_back(45);
  valid_genotype_call_codes[3].push_back(47);
  valid_genotype_call_codes[3].push_back(48);
  valid_genotype_call_codes[3].push_back(50);
  valid_genotype_call_codes[4].push_back(17);
  valid_genotype_call_codes[4].push_back(18);
  valid_genotype_call_codes[4].push_back(19);
  valid_genotype_call_codes[4].push_back(25);
  valid_genotype_call_codes[4].push_back(26);
  valid_genotype_call_codes[4].push_back(32);

  probeset_annot_csv_file = string("../data/genoCallCoder.testfile.csv");
}

void GenoCallCoderTest::testDecoding_v0() {
  //cout << "testDecoding_v0\n";
  GenoCallCoder v0;
  for (int i = 0; i < v0_code_vec.size(); i++) {
    CPPUNIT_ASSERT ( v0.genotypeCallNumToAbstractAllele( v0_code_vec[i] ) == v0_abs_allele_string_vec[i] );
  }
}

// void GenoCallCoderTest::testDecoding_delimit_v0() {
// //  cout << "testDecoding_v0\n";
//   GenoCallCoder v0(2,"UCHAR","0",'/');
//   for (int i = 0; i < v0_code_vec.size(); i++) {
//     CPPUNIT_ASSERT ( v0.genotypeCallNumToAbstractAllele( v0_code_vec[i] ) == v0_abs_allele_string_delimit_vec[i] );
//   }
// }

void GenoCallCoderTest::testEncoding_v0() {
  //cout << "testEncoding_v0\n";
  GenoCallCoder v0;
  for (int i = 0; i < v0_abs_allele_string_vec.size(); i++) {
    CPPUNIT_ASSERT ( v0.abstractAlleleToGenotypeCallNum(v0_abs_allele_string_vec[i]) == v0_code_vec[i] );
  }
}

// void GenoCallCoderTest::testEncoding_delimit_v0() {
// //  cout << "testEncoding_v0\n";
//   GenoCallCoder v0(2,"UCHAR","0",'/');
//   for (int i = 0; i < v0_abs_allele_string_delimit_vec.size(); i++) {
//     CPPUNIT_ASSERT ( v0.abstractAlleleToGenotypeCallNum(v0_abs_allele_string_delimit_vec[i]) == v0_code_vec[i] );
//   }
// }

void GenoCallCoderTest::testDecoding_v1() {
  //cout << "testDecoding_v1\n";
  GenoCallCoder v1(6, "UCHAR", "1.0");
  for (int i = 0; i < v1_code_vec.size(); i++) {
    CPPUNIT_ASSERT ( v1.genotypeCallNumToAbstractAllele(v1_code_vec[i]) == v1_abs_allele_string_vec[i] );
  }
}

void GenoCallCoderTest::testDecoding_delimit_v1() {
//  cout << "testDecoding_v1\n";
  GenoCallCoder v1(6, "UCHAR", "1.0",'/');
  for (int i = 0; i < v1_code_vec.size(); i++) {
    CPPUNIT_ASSERT ( v1.genotypeCallNumToAbstractAllele(v1_code_vec[i]) == v1_abs_allele_string_delimit_vec[i] );
  }
}

void GenoCallCoderTest::testEncoding_v1() {
//  cout << "testEncoding_v1\n";
  GenoCallCoder v1(6, string("UCHAR"), string("1.0"));
  for (int i = 0; i < v1_abs_allele_string_vec.size(); i++) {
    CPPUNIT_ASSERT ( v1.abstractAlleleToGenotypeCallNum(v1_abs_allele_string_vec[i]) == v1_code_vec[i] );
  }
}

void GenoCallCoderTest::testEncoding_delimit_v1() {
//  cout << "testEncoding_v1\n";
  GenoCallCoder v1(6, "UCHAR", "1.0", '/');
  for (int i = 0; i < v1_abs_allele_string_delimit_vec.size(); i++) {
    CPPUNIT_ASSERT ( v1.abstractAlleleToGenotypeCallNum(v1_abs_allele_string_delimit_vec[i]) == v1_code_vec[i] );
  }
}

void GenoCallCoderTest::testAbstractAlelleIntegerEncoding_v1() {
  GenoCallCoder v1(6, string("UCHAR"), string("1.0"));
  
  vector<int> F(1,5); 
  //cout << v1.abstractAlleleIntegersToGenotypeCallNum(F) << ":" << static_cast<int>(v1_code_vec[2])  << '\n';
  CPPUNIT_ASSERT ( v1.abstractAlleleIntegersToGenotypeCallNum(F) == v1_code_vec[2] );
  vector<int> AA(2,0);  
  CPPUNIT_ASSERT ( v1.abstractAlleleIntegersToGenotypeCallNum(AA) == v1_code_vec[4] );
  vector<int> AB(2,0);
  AB[1] = 1;
  CPPUNIT_ASSERT ( v1.abstractAlleleIntegersToGenotypeCallNum(AB) == v1_code_vec[5] );
  vector<int> BB(2,1);  
  CPPUNIT_ASSERT ( v1.abstractAlleleIntegersToGenotypeCallNum(BB) == v1_code_vec[7] );
  vector<int> CF(2,2);
  CF[1] = 5;
  CPPUNIT_ASSERT ( v1.abstractAlleleIntegersToGenotypeCallNum(CF) == v1_code_vec[8] );
  vector<int> ABC(3,0);
  ABC[1] = 1;
  ABC[2] = 2;
  CPPUNIT_ASSERT ( v1.abstractAlleleIntegersToGenotypeCallNum(ABC) == v1_code_vec[10] );
}

void GenoCallCoderTest::testDecodeVec_v0() {
//  cout << "testDecodeVec_v0\n";
  GenoCallCoder v0;
  vector<string> results = v0.genotypeCallNumVecToAbstractAlleleVec(v0_code_vec);
  CPPUNIT_ASSERT ( results == v0_abs_allele_string_vec );	   
}

void GenoCallCoderTest::testEncodeVec_v0() {
//  cout << "testEncodeVec_v0\n";
  GenoCallCoder v0;
  vector<unsigned char> results = v0.abstractAlleleVecToGenotypeCallNumVec(v0_abs_allele_string_vec);
  CPPUNIT_ASSERT ( results == v0_code_vec );
}

void GenoCallCoderTest::testDecodeVec_v1() {
//  cout << "testDecodeVec_v1\n";
  GenoCallCoder v1(6, string("UCHAR"), string("1.0"));
  vector<string> results = v1.genotypeCallNumVecToAbstractAlleleVec(v1_code_vec);
  CPPUNIT_ASSERT ( results == v1_abs_allele_string_vec);	 
}

void GenoCallCoderTest::testEncodeVec_v1() {
//  cout << "testEncodeVec_v1\n";
  GenoCallCoder v1(6, string("UCHAR"), string("1.0"));
  vector<unsigned char> results = v1.abstractAlleleVecToGenotypeCallNumVec(v1_abs_allele_string_vec);
  CPPUNIT_ASSERT ( results == v1_code_vec );
}

void GenoCallCoderTest::testGetDecodeVec_v0() {
//  cout << "testGetDecodeVec_v0\n";
  GenoCallCoder v0;
  vector<string> my_decodeVec(v0.getAbstractAlleleVector());
  for (int i = 0; i < v0_code_vec.size(); i++) {
    CPPUNIT_ASSERT ( my_decodeVec[v0_code_vec[i]] == v0_abs_allele_string_vec[i] );
  }  
}

void GenoCallCoderTest::testGetDecodeVec_v1() {
//  cout << "testGetDecodeVec_v1\n";
  GenoCallCoder v1(6, string("UCHAR"), string("1.0"));
  vector<string> my_decodeVec(v1.getAbstractAlleleVector());
  for (int i = 0; i < v1_code_vec.size(); i++) {
    CPPUNIT_ASSERT ( my_decodeVec[v1_code_vec[i]] == v1_abs_allele_string_vec[i] );
  }  
}

void GenoCallCoderTest::testIsValidCode_v0() {
//  cout << "testIsValidCode_v0\n";
  GenoCallCoder v0;
  for (int i = 0; i < v0_code_vec.size(); i++) {
    CPPUNIT_ASSERT ( v0.isValidCallNum(v0_code_vec[i]) );
  }
  CPPUNIT_ASSERT ( !v0.isValidCallNum(0) );
  CPPUNIT_ASSERT ( !v0.isValidCallNum(15) );
}

void GenoCallCoderTest::testIsValidAllele_v0() {
  //cout << "testIsValidAllele_v0\n";
  GenoCallCoder v0;
  for (int i = 0; i < v0_abs_allele_string_vec.size(); i++) {
    CPPUNIT_ASSERT ( v0.isValidAbstractAllele(v0_abs_allele_string_vec[i]) );
  }
  CPPUNIT_ASSERT ( !v0.isValidAbstractAllele(string("AC")) );
  CPPUNIT_ASSERT ( !v0.isValidAbstractAllele(string("NoCall")) );
}

// void GenoCallCoderTest::testIsValidAllele_delimit_v0() {
// //  cout << "testIsValidAllele_v0\n";
//   GenoCallCoder v0(2,"UCHAR","0",'/');
//   for (int i = 0; i < v0_abs_allele_string_delimit_vec.size(); i++) {
//     CPPUNIT_ASSERT ( v0.isValidAbstractAllele(v0_abs_allele_string_delimit_vec[i]) );
//   }
//   CPPUNIT_ASSERT ( !v0.isValidAbstractAllele(string("A/C")) );
//   CPPUNIT_ASSERT ( !v0.isValidAbstractAllele(string("N/C")) );
// }

void GenoCallCoderTest::testIsValidCode_v1() {
  //cout << "testIsValidCode_v1\n";
  GenoCallCoder v1(6, "UCHAR", "1.0", '/');
  for (int i = 0; i < v1_code_vec.size(); i++) {
    CPPUNIT_ASSERT ( v1.isValidCallNum( v1_code_vec[i]) );
  }
  CPPUNIT_ASSERT ( !v1.isValidCallNum(6) );
  CPPUNIT_ASSERT ( !v1.isValidCallNum(15) );
  CPPUNIT_ASSERT ( !v1.isValidCallNum(137) );
  CPPUNIT_ASSERT ( !v1.isValidCallNum(253) );
}

void GenoCallCoderTest::testIsValidAllele_v1() {
  //cout << "testIsValidAllele_v1\n";
  GenoCallCoder v1(6, string("UCHAR"), string("1.0"));
  for (int i = 0; i < v1_abs_allele_string_vec.size(); i++) {
    CPPUNIT_ASSERT ( v1.isValidAbstractAllele(v1_abs_allele_string_vec[i]) );
  }
  CPPUNIT_ASSERT ( !v1.isValidAbstractAllele(string("BA")) );
  CPPUNIT_ASSERT ( !v1.isValidAbstractAllele(string("NN")) );
}

void GenoCallCoderTest::testIsValidAllele_delimit_v1() {
  //cout << "testIsValidAllele_v1\n";
  GenoCallCoder v1(6, "UCHAR", "1.0", '/');
  for (int i = 0; i < v1_abs_allele_string_delimit_vec.size(); i++) {
    CPPUNIT_ASSERT ( v1.isValidAbstractAllele(v1_abs_allele_string_delimit_vec[i]) );
  }
  CPPUNIT_ASSERT ( !v1.isValidAbstractAllele(string("B/A")) );
  CPPUNIT_ASSERT ( !v1.isValidAbstractAllele(string("NN")) );
}

void GenoCallCoderTest::testIsHom_Het_v1() {
  GenoCallCoder v1(6, "UCHAR", "1.0", '/');
  for (int i = 0; i < v1_code_vec.size(); i++) {
//     cout << int(v1_code_vec[i]) << ':' << v1.isHet( v1_code_vec[i]) << ":" << v1_code_het[i] << ':' << v1.isHom(v1_code_vec[i]) << ':' << v1_code_hom[i] << '\n';
    CPPUNIT_ASSERT ( v1.isHom( v1_code_vec[i]) == v1_code_hom[i]);
    CPPUNIT_ASSERT ( v1.isHet( v1_code_vec[i]) == v1_code_het[i]);
  }

}

void GenoCallCoderTest::testAbstractAlleleToReferenceAllele_delimit_v1() {
  //cout << "testAbstractAlleleToReferenceAllele_delimit_v1\n";
  GenoCallCoder v1(6, "UCHAR", "1.0", '/', probeset_annot_csv_file);
  for (int i = 0; i < probeset_ids.size(); i++) {
//     cout << probeset_ids[i] + ':' + abstract_alleles[i] + ':' + reference_alleles[i] + ':' + v1.abstractAlleleToReferenceAllele(probeset_ids[i],abstract_alleles[i]) << '\n';
    CPPUNIT_ASSERT ( v1.abstractAlleleToReferenceAllele(probeset_ids[i],abstract_alleles[i]) == reference_alleles[i]);
  }
}


void GenoCallCoderTest::testAbstractAlleleToReportAllele_delimit_v1() {
  //cout << "testAbstractAlleleToReportAllele_delimit_v1\n";
  GenoCallCoder v1(6, "UCHAR", "1.0", '/', probeset_annot_csv_file);
  for (int i = 0; i < probeset_ids.size(); i++) {
    //cout << probeset_ids[i] + ':' + abstract_alleles[i] + ':' + report_alleles[i] + ':' + v1.abstractAlleleToReportAllele(probeset_ids[i],abstract_alleles[i]) << '\n';
    CPPUNIT_ASSERT ( v1.abstractAlleleToReportAllele(probeset_ids[i],abstract_alleles[i]) == report_alleles[i]);
  }
}

void GenoCallCoderTest::testReferenceAlleleToReportAllele_delimit_v1() {
  //cout << "testReferenceAlleleToReportAllele_delimit_v1\n";
  GenoCallCoder v1(6, "UCHAR", "1.0", '/', probeset_annot_csv_file);
  for (int i = 0; i < probeset_ids.size(); i++) {
    // cout << probeset_ids[i] + ':' + reference_alleles[i] + ':' + report_alleles[i] + ':' + v1.referenceAlleleToReportAllele(probeset_ids[i],reference_alleles[i]) << '\n';
    CPPUNIT_ASSERT ( v1.referenceAlleleToReportAllele(probeset_ids[i],reference_alleles[i]) == report_alleles[i]);
  }
}


void GenoCallCoderTest::testGenotypeCallNumToReferenceAllele_delimit_v1() {
  GenoCallCoder v1(6, "UCHAR", "1.0", '/', probeset_annot_csv_file);
  for (int i = 0; i < reference_geno_calls_test.size(); i++) {
    for (int j = 0; j < reference_geno_calls_test[i].size(); j++) {
      
      CPPUNIT_ASSERT ( v1.genotypeCallNumToReferenceAllele(probeset_ids[i], valid_genotype_call_codes[i][j]) == reference_geno_calls_test[i][j] );
      CPPUNIT_ASSERT ( v1.genotypeCallNumToReferenceAllele(probeset_ids[i], 16) == "NoCall" );
      CPPUNIT_ASSERT ( v1.genotypeCallNumToReferenceAllele(probeset_ids[i], 254) == "PossibleRareAllele" );
      CPPUNIT_ASSERT ( v1.genotypeCallNumToReferenceAllele(probeset_ids[i], 255) == "NotAvailable" );

    }
  }
}

void GenoCallCoderTest::testGenotypeCallNumToReportAllele_delimit_v1() {
  GenoCallCoder v1(6, "UCHAR", "1.0", '/', probeset_annot_csv_file);
  for (int i = 0; i < report_geno_calls_test.size(); i++) {
    for (int j = 0; j < report_geno_calls_test[i].size(); j++) {
      
      CPPUNIT_ASSERT ( v1.genotypeCallNumToReportAllele(probeset_ids[i], valid_genotype_call_codes[i][j]) == report_geno_calls_test[i][j] );
      CPPUNIT_ASSERT ( v1.genotypeCallNumToReportAllele(probeset_ids[i], 16) == "NoCall" );
      CPPUNIT_ASSERT ( v1.genotypeCallNumToReportAllele(probeset_ids[i], 254) == "PossibleRareAllele" );
      CPPUNIT_ASSERT ( v1.genotypeCallNumToReportAllele(probeset_ids[i], 255) == "NotAvailable" );

    }
  }
}

void GenoCallCoderTest::testGetValidGenotypeCallCodes_delimit_v1() {
  //cout << "testGetValidGenotypeCallCodes_delimit_v1\n";
  GenoCallCoder v1(6, "UCHAR", "1.0", '/', probeset_annot_csv_file);
  for (int i = 0; i < probeset_ids.size(); i++) {
    if (valid_genotype_call_codes[i].empty()) {
      break;
    }
//     cout << i << '\n';
//     vector<unsigned char> temp = v1.getValidGenotypeCallCodes(probeset_ids[i]);
//     for (int j = 0; j < temp.size(); j++) {
//       cout << temp.size() << ' ' << valid_genotype_call_codes[i].size() << '\n';
//       cout << static_cast<int>(temp[j]) << ':' << static_cast<int>(valid_genotype_call_codes[i][j]) << '\n';
//     }
    CPPUNIT_ASSERT ( v1.getValidGenotypeCallCodes(probeset_ids[i]) == valid_genotype_call_codes[i] );
  }
}

void GenoCallCoderTest::testGetProbesetIds_delimit_v1() {
  GenoCallCoder v1(6, "UCHAR", "1.0", '/', probeset_annot_csv_file);
//   vector<string> temp = v1.getProbesetIds();
//   for (int i = 0; i < probeset_ids.size(); i++) {
//     cout << temp[i] + "==" + probeset_ids[i] + '\n';
//   }
  CPPUNIT_ASSERT ( v1.getProbesetIds() == probeset_ids );
}

// void GenoCallCoderTest::testFullCsvFile() {
//   GenoCallCoder v1(6, "UCHAR", "1.0", '/', "/home/rwheel/affy/sdk/calvin_files/utils/data/DMET_Plus.r3.20081008.altered_a-a.dc_annot.csv");

//   CPPUNIT_ASSERT ( &v1 != NULL );

//   cout << v1.genotypeCallNumToReportAllele("AM_14631", 26) << '\n';
//   cout << v1.genotypeCallNumToReportAllele("AM_14633", 26) << '\n';
// }



/** 
 * Everbybodys favorite function. 
 * @param argc 
 * @param arg 
 * 
 * @return 
 */
// int main(int argc, char *arg[])
// {
//   // Get the top level suite from the registry
//   CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();

//   // Adds the test to the list of test to run
//   CppUnit::TextUi::TestRunner runner;
//   runner.addTest( suite );

//   // Change the default outputter to a compiler error format outputter
//   runner.setOutputter( new CppUnit::CompilerOutputter( &runner.result(),
//                                                        std::cerr ) );
//   // Run the tests.
//   bool wasSucessful = runner.run();

//   // Return error code 1 if the one of test failed.
//   return wasSucessful ? 0 : 1;
// }



