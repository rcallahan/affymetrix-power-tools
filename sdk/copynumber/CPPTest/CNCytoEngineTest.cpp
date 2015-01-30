////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
#include "copynumber/CNCytoEngine.h"
#include "copynumber/CPPTest/Setup.h"
//
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string> 
//
using namespace std;
/**
 * @class CNCytoEngineTest
 * @brief cppunit class for testing CNCytoEngine functions.
 * last change by vliber on 10/23/09
 */

class CNCytoEngineTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CNCytoEngineTest); 
  CPPUNIT_TEST(defaultDefineOptionsTest);
  CPPUNIT_TEST(addCelTest); 
  CPPUNIT_TEST(checkParametersTest); 
  CPPUNIT_TEST(checkOptionsTest); 
  CPPUNIT_TEST(runTest);
 
  
  CPPUNIT_TEST_SUITE_END();

public:  
  void defaultDefineOptionsTest();
  void addCelTest();
  void checkParametersTest();
  void checkOptionsTest();
  void runTest();
  
};

CPPUNIT_TEST_SUITE_REGISTRATION(CNCytoEngineTest );



void CNCytoEngineTest::defaultDefineOptionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNCytoEngineTest::defaultDefineOptionsTest****");	
	CNCytoEngine cnCyto;
	CPPUNIT_ASSERT(cnCyto.getOptBool("h")==false); 
	CPPUNIT_ASSERT(cnCyto.getOpt("explain")=="");
	CPPUNIT_ASSERT(cnCyto.getOptInt("verbose")==1); 
	CPPUNIT_ASSERT(cnCyto.getOptBool("version")==false);
	CPPUNIT_ASSERT(cnCyto.getOpt("config-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("xml-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOptBool("global-parameter-override")==false);
    CPPUNIT_ASSERT(cnCyto.getOptBool("keep-temp-reference-data")==false);
    // Not defined until checkOptions() is called
	//CPPUNIT_ASSERT(cnCyto.getOpt("snp-reference-input-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("snp-reference-output-file")=="");
    CPPUNIT_ASSERT(cnCyto.getOpt("gender-override-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("reference-input")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("reference-output")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("probe-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("cdf-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOptBool("force")==false);
    CPPUNIT_ASSERT(cnCyto.getOpt("qcc-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("qca-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("cel-files")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("chrX-probes")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("chrY-probes")=="");
    CPPUNIT_ASSERT(cnCyto.getOpt("genotype-call-override-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("reference-cels")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("cels")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("arrs")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("result-files")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("probeset-ids")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("annotation-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("gc-content-override-file")=="");
	CPPUNIT_ASSERT(cnCyto.getOptInt("gc-correction-bin-count")==25);
	CPPUNIT_ASSERT(cnCyto.getOpt("temp-dir")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("out-dir")==".");
	CPPUNIT_ASSERT(cnCyto.getOptBool("run-geno-qc")==false);
	CPPUNIT_ASSERT(cnCyto.getOptBool("adapter-type-normalization")==false);
    CPPUNIT_ASSERT(cnCyto.getOpt("analysis")=="");
    CPPUNIT_ASSERT(cnCyto.getOpt("wave-correction-reference-method")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("local-gc-background-correction-reference-method")=="none");
	CPPUNIT_ASSERT(cnCyto.getOpt("local-gc-background-intensity-adjustment-method")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("image-correction-intensity-adjustment-method")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("wave-correction-log2ratio-adjustment-method")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("allele-peaks-reporter-method")=="allele-peaks-reporter-method");
	CPPUNIT_ASSERT(cnCyto.getOpt("command-line")==""); 
	CPPUNIT_ASSERT(cnCyto.getOpt("exec-guid")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("program-name")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("program-company")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("program-version")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("program-cvs-id")=="");
	CPPUNIT_ASSERT(cnCyto.getOpt("version-to-report")=="");
	CPPUNIT_ASSERT(cnCyto.getOptInt("xChromosome")==24);
	CPPUNIT_ASSERT(cnCyto.getOptInt("yChromosome")==25);
    CPPUNIT_ASSERT(cnCyto.getOptBool("snp-qc-use-contrast")==false);
    CPPUNIT_ASSERT(cnCyto.getOpt("snp-qc-snp-list")=="");
    CPPUNIT_ASSERT(cnCyto.getOptDouble("snp-qc-k")==2.0);
    CPPUNIT_ASSERT(cnCyto.getOptDouble("snp-qc-em-threshold")==0.05);
    CPPUNIT_ASSERT(cnCyto.getOptDouble("snp-qc-bin-size")==0.04);
	CPPUNIT_ASSERT(cnCyto.getOptDouble("male-gender-ratio-cutoff")==1.3);
	CPPUNIT_ASSERT(cnCyto.getOptDouble("female-gender-ratio-cutoff")==1.0);
	CPPUNIT_ASSERT(cnCyto.getOptInt("reference-chromosome")==2);
	CPPUNIT_ASSERT(cnCyto.getOptDouble("xx-cutoff")==0.8); 
	CPPUNIT_ASSERT(cnCyto.getOptDouble("xx-cutoff-high")==1.07);
	CPPUNIT_ASSERT(cnCyto.getOptDouble("y-cutoff")==0.65);
}


void CNCytoEngineTest::addCelTest()
{
	Verbose::out(1, "****CNCytoEngineTest::addCelTest****");	
	CNCytoEngine cnCyto;
	cnCyto.addCel("./test_strCelFileName","./test_strSampleFileName","./test_strResultFileName");
	CPPUNIT_ASSERT(cnCyto.getOpt("cels")=="./test_strCelFileName");
	CPPUNIT_ASSERT(cnCyto.getOpt("arrs")=="./test_strSampleFileName");
	CPPUNIT_ASSERT(cnCyto.getOpt("result-files")=="./test_strResultFileName");
}	
void CNCytoEngineTest::checkParametersTest()
{
	Verbose::out(1, "****CNCytoEngineTest::checkParametersTest****");	

	//adapter-type-normalization cannot be run against the cyto2 arrays.
	CNCytoEngine cnCyto1;    
	cnCyto1.setOpt("cyto2","true");
	cnCyto1.setOpt("adapter-type-normalization","true");
	NEGATIVE_TEST(cnCyto1.checkOptions(),Except);
 
   	CNCytoEngine cnCyto2; 
	std::vector< std::string > arr_values;
	cnCyto2.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cnCyto2.setOpt("reference-output",OUTPUT + "/ref.a5");
   	arr_values.push_back(INPUT + "/Cyto/Beta1_M_01_Cyto_VH_20090202.ARR");
	POSITIVE_TEST(cnCyto2.setOpt("cels",arr_values));
	POSITIVE_TEST(cnCyto2.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto2.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto2.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
    POSITIVE_TEST(cnCyto2.setOpt("cychp-output","false"));
	POSITIVE_TEST(cnCyto2.setOpt("cnchp-output","false"));
	//Must specify either cychp-ouptut or cnchp-output.
	NEGATIVE_TEST(cnCyto2.checkOptions(),Except);
}
void CNCytoEngineTest::checkOptionsTest()
{
	Verbose::out(1, "****CNCytoEngineTest::checkOptionsTest****");	
	
	CNCytoEngine cnCyto1,cnCyto2;
	//Message: No cel files specified
	NEGATIVE_TEST(cnCyto1.checkOptions(), Except);
    cnCyto2.setOpt("out-dir","./Cyto");
    cnCyto2.setOpt("reference-input",INPUT + "/Cyto/ref.a5");
	cnCyto2.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	//FATAL ERROR: Please specify either a reference-input file or a reference-output file, but not both.
	NEGATIVE_TEST(cnCyto2.checkOptions(), Except);
	
	CNCytoEngine cnCyto3;
	cnCyto3.setOpt("reference-input",INPUT + "/Cyto/ref_test.a5");
	//FATAL ERROR: The reference-input specified either does not exist, or is not the correct type for the CNCytoEngine.
   	NEGATIVE_TEST(cnCyto3.checkOptions(), Except); 
    
	CNCytoEngine cnCyto3a;
	cnCyto3a.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	cnCyto3a.setOpt("cyto2","true");
	//FATAL ERROR: Please specify a valid probe-file when creating a reference.
   	NEGATIVE_TEST(cnCyto3a.checkOptions(), Except); 
	
    std::vector< std::string > values;
	CNCytoEngine cnCyto4;
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel");
	cnCyto4.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	POSITIVE_TEST(cnCyto4.setOpt("reference-cels",values));
	cnCyto4.setOpt("probe-file",INPUT+"/Cyto/CytogeneticsFocused_Array_probeSequences.txt");
	//FATAL ERROR: Specify reference-cels in combination with the reference-input option, not the reference-output option.
    NEGATIVE_TEST(cnCyto4.checkOptions(), Except);

	CNCytoEngine cnCyto4a;
	cnCyto4a.setOpt("snp-reference-input-file",OUTPUT + "/Cyto/CytogeneticsFocused_renamed.snpref_test.a5");
	//FATAL ERROR: The snp-reference-input-file specified either does not exist, or is not the correct type for the CNCytoEngine.
    NEGATIVE_TEST(cnCyto4a.checkOptions(), Except);
	
	
	CNCytoEngine cnCyto5;
	cnCyto5.setOpt("snp-reference-output-file",OUTPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	cnCyto5.setOpt("snp-reference-input-file","");
	//FATAL ERROR: When using the snp-reference-output-file option you must also use 
    //the snp-reference-input-file option.  This file must contain valid data for the
    //information component of each probeset.
	NEGATIVE_TEST(cnCyto5.checkOptions(), Except);

	CNCytoEngine cnCyto6;
	cnCyto6.setOpt("snp-reference-output-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	cnCyto6.setOpt("snp-reference-input-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	//FATAL ERROR: snp-reference-input-file cannot have the same name as the snp-reference-output-file
	NEGATIVE_TEST(cnCyto6.checkOptions(), Except);

	CNCytoEngine cnCyto7;
	cnCyto7.setOpt("snp-reference-output-file",OUTPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	cnCyto7.setOpt("snp-reference-input-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	cnCyto7.setOpt("reference-output","");
	//FATAL ERROR: To use the snp-reference-output-file option you must also use the reference-output option
	NEGATIVE_TEST(cnCyto7.checkOptions(), Except);

	CNCytoEngine cnCyto7a;
	cnCyto7a.setOpt("snp-reference-output-file",OUTPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	cnCyto7a.setOpt("snp-reference-input-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	cnCyto7a.setOpt("reference-output",INPUT + "/Cyto/ref.a5");
	//FATAL ERROR: To use the snp-reference-file-output option you must specify a gender-override-file containing the gender of the samples being analyzed
	NEGATIVE_TEST(cnCyto7a.checkOptions(), Except);

	CNCytoEngine cnCyto7b;
	cnCyto7b.setOpt("snp-reference-output-file",OUTPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	cnCyto7b.setOpt("snp-reference-input-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.snpref.a5");
	cnCyto7b.setOpt("reference-output",INPUT + "/Cyto/ref.a5");
	cnCyto7b.setOpt("gender-override-file",INPUT + "/Cyto/ref.a5");
	//FATAL ERROR: To use the snp-reference-file-output option you must specify a genotype-call-override-file containing the genotypes of all the SNP marker of all samples being analyzed.
	NEGATIVE_TEST(cnCyto7b.checkOptions(), Except);

	CNCytoEngine cnCyto8;
	cnCyto8.setOpt("cyto2","true");
	//FATAL ERROR: No cel files specified.
	NEGATIVE_TEST(cnCyto8.checkOptions(), Except);

	CNCytoEngine cnCyto8a;
	values.clear();
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18540_A01_01_CytoF_NN_20090121.cel");
	POSITIVE_TEST(cnCyto8a.setOpt("cels",values));
	cnCyto8a.setOpt("cyto2","true");
	cnCyto8a.setOpt("snp-reference-input-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.snpref_test.a5");
	//FATAL ERROR: The snp-reference-input-file specified either does not exist, or is not the correct type for the CNCytoEngine.
	NEGATIVE_TEST(cnCyto8a.checkOptions(), Except);
	values.clear();

	CNCytoEngine cnCyto9;
    CPPUNIT_ASSERT(cnCyto9.isOptDefined("cel-files"));
	CPPUNIT_ASSERT(cnCyto9.isOptDefined("cels"));
    cnCyto9.setOpt("cel-files",INPUT + "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel");
	cnCyto9.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	//FATAL ERROR: .\TsvFile\TsvFile.cpp:2952: Required column: 'cel_files' not found in file: '../../regression-data/data/copynumber-cyto/cppunit/input/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel'.
	NEGATIVE_TEST(cnCyto9.checkOptions(), Except);
	
	CNCytoEngine cnCyto10;
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel");
    POSITIVE_TEST(cnCyto10.isOptDefined("set-analysis-name"));
	POSITIVE_TEST(cnCyto10.getOpt("set-analysis-name")=="CN5");
	POSITIVE_TEST(cnCyto10.setOpt("cels",values));
	//Must specify a valid annotation-file.
	NEGATIVE_TEST(cnCyto10.checkOptions(), Except);

	CNCytoEngine cnCyto11;
	POSITIVE_TEST(cnCyto11.isOptDefined("annotation-file"));
	cnCyto11.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cnCyto11.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	POSITIVE_TEST(cnCyto11.setOpt("cels",values));
	//Must specify a Must specify a cdf-file or spf-file
	NEGATIVE_TEST(cnCyto11.checkOptions(), Except);
    
	
    CNCytoEngine cnCyto12;
	cnCyto12.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cnCyto12.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	POSITIVE_TEST(cnCyto12.setOpt("cels",values));
	POSITIVE_TEST(cnCyto12.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	//Must specify chrX-probes.
	NEGATIVE_TEST(cnCyto12.checkOptions(), Except);
	
	
    CNCytoEngine cnCyto13;
	cnCyto13.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cnCyto13.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	POSITIVE_TEST(cnCyto13.setOpt("cels",values));
	POSITIVE_TEST(cnCyto13.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto13.setOpt("chrX-probes",INPUT + "/Cyto/Cytogenetics.chrx.probes.txt"));
	//Must specify chrY-probes.
	NEGATIVE_TEST(cnCyto13.checkOptions(), Except);

	
    CNCytoEngine cnCyto13a;
	cnCyto13a.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cnCyto13a.setOpt("reference-output",OUTPUT + "/Cyto/bad_hdf_ref.a5");
	POSITIVE_TEST(cnCyto13a.setOpt("cels",values));
	POSITIVE_TEST(cnCyto13a.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto13a.setOpt("chrX-probes",INPUT + "/Cyto/Cytogenetics.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto13a.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	//FATAL ERROR: Filename: 'HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel' seen multiple times. Filenames must be unique
	NEGATIVE_TEST(cnCyto13a.checkOptions(), Except);
	
	CNCytoEngine cnCyto14; 
	cnCyto14.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cnCyto14.setOpt("reference-output",OUTPUT + "/ref.a5");
	POSITIVE_TEST(cnCyto14.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto14.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto14.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	POSITIVE_TEST(cnCyto14.setOpt("cels",values));
	CPPUNIT_ASSERT(cnCyto14.isOptDefined("run-geno-qc"));
	CPPUNIT_ASSERT(!cnCyto14.getOptBool("run-geno-qc"));
	CPPUNIT_ASSERT(cnCyto14.getOpt("qca-file")=="");
	CPPUNIT_ASSERT(cnCyto14.getOpt("qcc-file")=="");
    POSITIVE_TEST(cnCyto14.setOpt("run-geno-qc","true"));
	//FATAL ERROR: Please specify a qca-file.
	NEGATIVE_TEST(cnCyto14.checkOptions(), Except);

	CNCytoEngine cnCyto15; 
	cnCyto15.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cnCyto15.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	POSITIVE_TEST(cnCyto15.setOpt("cels",values));
	POSITIVE_TEST(cnCyto15.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto15.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto15.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	POSITIVE_TEST(cnCyto15.setOpt("run-geno-qc","true"));
	POSITIVE_TEST(cnCyto15.setOpt("qca-file",INPUT + "/Cyto/qca_test.txt"));
	//FATAL ERROR: Please specify a qcc-file.
	NEGATIVE_TEST(cnCyto15.checkOptions(),Except);
	
	CNCytoEngine cnCyto16; 
	cnCyto16.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cnCyto16.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18540_A01_01_CytoF_NN_20090121.cel");
	POSITIVE_TEST(cnCyto16.setOpt("cels",values));
	POSITIVE_TEST(cnCyto16.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto16.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto16.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
    std::vector< std::string > arr_values;
	arr_values.push_back(INPUT + "/Cyto/Beta1_M_01_Cyto_VH_20090202.ARR");
	POSITIVE_TEST(cnCyto16.setOpt("arrs",arr_values));
	//FATAL ERROR: When specifying ARR files, the number of ARR files specfified must match the number of CEL files specified .
	NEGATIVE_TEST(cnCyto16.checkOptions(),Except);

	arr_values.clear();
	values.clear();
	arr_values.push_back(INPUT + "/Cyto/Beta1_M_01_Cyto_VH_20090202.ARR");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18540_A01_01_CytoF_NN_20090121.cel");
	CNCytoEngine cnCyto17; 
	cnCyto17.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	cnCyto17.setOpt("reference-output",OUTPUT + "/Cyto/ref.a5");
	POSITIVE_TEST(cnCyto17.setOpt("cels",values));
	POSITIVE_TEST(cnCyto17.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto17.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto17.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
   	POSITIVE_TEST(cnCyto17.setOpt("arrs",arr_values));
	//Message: Different chiptypes or probe counts found. Only first 10 errors printed, see log for rest.
	NEGATIVE_TEST(cnCyto17.checkOptions(),Except);
		
}
void CNCytoEngineTest::runTest()
{
	Verbose::out(1, "****CNCytoEngineTest::runTest****");

	CNCytoEngine cnCyto1; 
	std::vector< std::string > values;
	cnCyto1.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	///@todo AW: cnCyto1.setState("create-reference","true");
	cnCyto1.setOpt("reference-output",OUTPUT + "/ref.a5");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel");
	POSITIVE_TEST(cnCyto1.setOpt("cels",values));
	POSITIVE_TEST(cnCyto1.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto1.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto1.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	CPPUNIT_ASSERT(cnCyto1.isOptDefined("run-geno-qc"));
	CPPUNIT_ASSERT(!cnCyto1.getOptBool("run-geno-qc"));
	CPPUNIT_ASSERT(!cnCyto1.getOptBool("force"));
	//FATAL ERROR: Different chiptypes or probe counts found. Only first 10 errors printed, see log for rest  
	NEGATIVE_TEST(cnCyto1.run(),Except);
	

	CNCytoEngine cnCyto1a; 
	cnCyto1a.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	///@todo AW: cnCyto1a.setState("create-reference","true");
	cnCyto1a.setOpt("reference-output",OUTPUT + "/ref.a5");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel");
	POSITIVE_TEST(cnCyto1a.setOpt("cels",values));
	POSITIVE_TEST(cnCyto1a.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto1a.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto1a.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	cnCyto1a.setOpt("force","true");
	//FATAL ERROR: Exception caught. Message is: At least six cel files are needed to get a valid CN Reference file.
	NEGATIVE_TEST(cnCyto1a.run(),Except);
	
    
	CNCytoEngine cnCyto2; 
	values.clear();
	cnCyto2.setOpt("reference-output",OUTPUT + "/ref.a5");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121_test.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18540_A01_01_CytoF_NN_20090121.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18547_A08_01_CytoF_NN_20090121.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18552_A12_01_CytoF_NN_20090121.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18555_B02_01_CytoF_NN_20090121.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18561_C10_01_CytoF_NN_20090121.cel");
	cnCyto2.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	POSITIVE_TEST(cnCyto2.setOpt("cels",values));
	POSITIVE_TEST(cnCyto2.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto2.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto2.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	cnCyto2.setOpt("adapter-type-normalization","false");
	cnCyto2.setOpt("out-dir",OUTPUT);
	//FATAL ERROR: Cannot read file ./input/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121_test.cel
	NEGATIVE_TEST(cnCyto2.run(),Except);

	
	CNCytoEngine cnCyto4;
	values.clear();
	cnCyto4.setOpt("reference-output",OUTPUT + "/ref.a5");
	cnCyto4.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18526_C09_01_CytoF_NN_20090121.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18540_A01_01_CytoF_NN_20090121.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18547_A08_01_CytoF_NN_20090121.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18552_A12_01_CytoF_NN_20090121.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18555_B02_01_CytoF_NN_20090121.cel");
	values.push_back(INPUT + "/Cyto/HapMap-As_NA18561_C10_01_CytoF_NN_20090121.cel");
	POSITIVE_TEST(cnCyto4.setOpt("cels",values));
	POSITIVE_TEST(cnCyto4.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array_test.CDF"));
	POSITIVE_TEST(cnCyto4.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto4.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	cnCyto4.setOpt("adapter-type-normalization","false");
	cnCyto4.setOpt("out-dir",OUTPUT);
	cnCyto4.setOpt("force","true");
	//FATAL ERROR: Error opening CDF file to read header ./input/Cyto/CytogeneticsFocused_Array_test.CDF.
	NEGATIVE_TEST(cnCyto4.run(),Except);

	
	CNCytoEngine cnCyto5;
	cnCyto5.setOpt("reference-output",OUTPUT + "/ref.a5");
	cnCyto5.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	POSITIVE_TEST(cnCyto5.setOpt("cels",values));
	POSITIVE_TEST(cnCyto5.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto5.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes_test.txt"));
	POSITIVE_TEST(cnCyto5.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	cnCyto5.setOpt("adapter-type-normalization","false");
	cnCyto5.setOpt("out-dir",OUTPUT);
	cnCyto5.setOpt("force","true");
	//FATAL ERROR: .\TsvFile\TsvFile.cpp:2352: open: Could not open file: './input/Cyto/CytogeneticsFocused.chrx.probes_test.txt' to read.
	NEGATIVE_TEST(cnCyto5.run(),Except);

	
	CNCytoEngine cnCyto6;
	cnCyto6.setOpt("reference-output",OUTPUT + "/ref.a5");
	cnCyto6.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	POSITIVE_TEST(cnCyto6.setOpt("cels",values));
	POSITIVE_TEST(cnCyto6.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto6.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto6.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes_test.txt"));
	cnCyto6.setOpt("adapter-type-normalization","false");
	cnCyto6.setOpt("out-dir",OUTPUT);
	cnCyto6.setOpt("force","true");
	//FATAL ERROR: .\TsvFile\TsvFile.cpp:2352: open: Could not open file: './input/Cyto/CytogeneticsFocused.chry.probes_test.txt' to read.
	NEGATIVE_TEST(cnCyto6.run(),Except);

    
	CNCytoEngine cnCyto7;
	cnCyto7.setOpt("reference-output",OUTPUT + "/ref.a5");
	cnCyto7.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot_empty.a5");
	POSITIVE_TEST(cnCyto7.setOpt("cels",values));
	POSITIVE_TEST(cnCyto7.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto7.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto7.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	cnCyto7.setOpt("adapter-type-normalization","false");
	cnCyto7.setOpt("out-dir",OUTPUT);
	cnCyto7.setOpt("force","true");
	//FATAL ERROR: Exception caught. Message is: ERROR: SQLiteCode: 1, Message: Failed to prepare SQL statement.
	NEGATIVE_TEST(cnCyto7.run(),Except);

    
	CNCytoEngine cnCyto8;
	cnCyto8.setOpt("reference-output",OUTPUT + "/ref.a5");
	cnCyto8.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	POSITIVE_TEST(cnCyto8.setOpt("cels",values));
	POSITIVE_TEST(cnCyto8.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array_empty.CDF"));
	POSITIVE_TEST(cnCyto8.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto8.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	cnCyto8.setOpt("adapter-type-normalization","false");
	cnCyto8.setOpt("out-dir",OUTPUT);
	cnCyto8.setOpt("force","true");
	//FATAL ERROR: Error opening CDF file to read header ./input/Cyto/CytogeneticsFocused_Array_empty.CDF: 
	NEGATIVE_TEST(cnCyto8.run(),Except);

   
	CNCytoEngine cnCyto9;
	cnCyto9.setOpt("reference-output",OUTPUT + "/ref.a5");
	cnCyto9.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	POSITIVE_TEST(cnCyto9.setOpt("cels",values));
	POSITIVE_TEST(cnCyto9.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto9.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes_empty.txt"));
	POSITIVE_TEST(cnCyto9.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
	cnCyto9.setOpt("adapter-type-normalization","false");
	cnCyto9.setOpt("out-dir",OUTPUT);
	cnCyto9.setOpt("force","true");
	 //FATAL ERROR: Unable to find probe_id column in file './input/Cyto/CytogeneticsFocused.chrx.probes_empty.txt'
	NEGATIVE_TEST(cnCyto9.run(),Except);


	
	CNCytoEngine cnCyto10;
	cnCyto10.setOpt("reference-output",OUTPUT + "/ref.a5");
	cnCyto10.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	POSITIVE_TEST(cnCyto10.setOpt("cels",values));
	POSITIVE_TEST(cnCyto10.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto10.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto10.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes_empty.txt"));
	cnCyto10.setOpt("adapter-type-normalization","false");
	cnCyto10.setOpt("out-dir",OUTPUT);
	cnCyto10.setOpt("force","true");
	//FATAL ERROR: Unable to find probe_id column in file './input/Cyto/CytogeneticsFocused.chry.probes_empty.txt'
	NEGATIVE_TEST(cnCyto10.run(),Except);

	CNCytoEngine cnCyto11;
	cnCyto11.setOpt("reference-output",OUTPUT + "/ref.a5");
	cnCyto11.setOpt("annotation-file",INPUT + "/Cyto/CytogeneticsFocused_renamed.annot.a5");
	POSITIVE_TEST(cnCyto11.setOpt("cels",values));
	POSITIVE_TEST(cnCyto11.setOpt("cdf-file",INPUT + "/Cyto/CytogeneticsFocused_Array.CDF"));
	POSITIVE_TEST(cnCyto11.setOpt("chrX-probes",INPUT + "/Cyto/CytogeneticsFocused.chrx.probes.txt"));
	POSITIVE_TEST(cnCyto11.setOpt("chrY-probes",INPUT + "/Cyto/CytogeneticsFocused.chry.probes.txt"));
    POSITIVE_TEST(cnCyto11.setOpt("probe-file",INPUT + "/Cyto/CytogeneticsFocused_Array_probeSequences.txt"));
    
	cnCyto11.setOpt("adapter-type-normalization","false");
	cnCyto11.setOpt("out-dir",OUTPUT);
	cnCyto11.setOpt("force","true");

        //APT-663: Test case is taking 13 minutes and needs moved to regression.
        Verbose::out(1,"APT-663: Test case cnCyto11 commented out due to 13 minute run time.");
#if 0        
	POSITIVE_TEST(cnCyto11.run());
    
	CPPUNIT_ASSERT(cnCyto11.getOptInt("sample-count")==6);
	CPPUNIT_ASSERT(cnCyto11.getOpt("snp-reference-input-file")=="");
	CPPUNIT_ASSERT(cnCyto11.getOpt("snp-reference-output-file")=="");
	CPPUNIT_ASSERT(cnCyto11.getOptBool("log2ratio-hdf5-output")==false);
	CPPUNIT_ASSERT(cnCyto11.getOptBool("log2ratio-text-output")==false);
	CPPUNIT_ASSERT(cnCyto11.getOptBool("reference-text-output")==false);
 	CPPUNIT_ASSERT(cnCyto11.getOptBool("create-snp-reference")==false);
 	CPPUNIT_ASSERT(cnCyto11.getOptBool("create-reference")==true);
 	CPPUNIT_ASSERT(cnCyto11.getOpt("copynumber-reference-file")=="");
 	CPPUNIT_ASSERT(cnCyto11.getOptInt("recommended-reference-sample-count")==44);
 	CPPUNIT_ASSERT(cnCyto11.getOptInt("recommended-reference-female-sample-count")==20);
 	CPPUNIT_ASSERT(cnCyto11.getOptInt("reference-sample-count")==6);
 	CPPUNIT_ASSERT(cnCyto11.getOptInt("reference-female-sample-count")==5);
 	CPPUNIT_ASSERT(cnCyto11.getOptInt("reference-male-sample-count")==1);
 	CPPUNIT_ASSERT(cnCyto11.getOptInt("reference-unknown-sample-count")==0);
 	CPPUNIT_ASSERT(cnCyto11.getOptInt("probe-count")==432964);
#endif

}
	

	

