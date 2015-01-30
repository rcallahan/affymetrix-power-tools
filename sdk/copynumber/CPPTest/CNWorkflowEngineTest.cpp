////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

//
#include "chipstream/SelfCreate.cpp"
#include "chipstream/apt-geno-qc/GenoQC.cpp"
#include "chipstream/apt-probeset-genotype/ProbesetGenotypeEngine.cpp"
#include "copynumber/CNAnalysisEngine.h"
#include "copynumber/CNLog2RatioData.h"
#include "copynumber/CNWorkflowEngine.h"
#include "copynumber/CPPTest/Setup.h"

#include "util/Fs.h"
#include "util/PgOptions.h"
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
/**
 * @class CNWorkflowEngineTest
 * @brief cppunit class for testing CNWorkflowEngine functions.
 * last change by vliber on 03/23/09
 */

class CNWorkflowEngineTest : public CNWorkflowEngine, public CppUnit::TestFixture
{
	
	CPPUNIT_TEST_SUITE(CNWorkflowEngineTest);
    CPPUNIT_TEST(defineOptionsTest);
    CPPUNIT_TEST(addCelTest);
    CPPUNIT_TEST(defineStateTest);
    CPPUNIT_TEST(checkOptionsTest);
	
	
	CPPUNIT_TEST_SUITE_END();

public:
   void defineOptionsTest();
   void addCelTest();
   void defineStateTest();
   void checkOptionsTest();  
   void setUp();
  
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNWorkflowEngineTest );

void CNWorkflowEngineTest::setUp()
{
  Fs::ensureWriteableDirPath(OUTPUT);
  Fs::fileCopy(INPUT + "/CNReferenceOriginal.a5",OUTPUT + "/CNReference.a5");
}

void CNWorkflowEngineTest::defineOptionsTest()
{
	cout<<endl;
	Verbose::out(1, "****CNWorkflowEngineTest::defineOptionsTest****");
	CNWorkflowEngine m_objCNWorkflowEngine;
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("help")==false);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("explain")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("verbose")==1);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("version")==false);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("explain")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("xml-file")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("reference-input")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("reference-output")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("cdf-file")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("force")==false);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("qcc-file")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("qca-file")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("cel-files")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("special-snps")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("chrX-probes")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("chrY-probes")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("target-sketch")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("use-feat-eff")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("read-models-brlmmp")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("out-dir")==".");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("analysis")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("adapter-type-normalization")==true);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("normalization-type")==1);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("adapter-parameters")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("brlmmp-parameters")=="");
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("allele-peaks-reporter-method")=="allele-peaks-reporter-method");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("command-line")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("exec-guid")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("program-name")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("program-company")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("program-version")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("program-cvs-id")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("version-to-report")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("mem-usage")==0);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("run-geno-qc")==false);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("run-probeset-genotype")==true);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("prior-size")==10000);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("cels")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("arrs")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("result-files")=="");
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("gc-correction-bin-count")==25);
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("probeset-ids")=="");
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("annotation-file")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("xChromosome")==24);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("yChromosome")==25);
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("delete-files")==false);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("log2-input")==false);
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("gc-content-override-file")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("geno-qc-file")=="");
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("cyto2")==false);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("array-name")=="");
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("set-analysis-name")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("text-output")==false);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("cnchp-output")==true);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptBool("cychp-output")==false);
}


void CNWorkflowEngineTest::addCelTest()
{
	Verbose::out(1, "****CNWorkflowEngineTest::addCelTest****");
	CNWorkflowEngine m_objCNWorkflowEngine;
	m_objCNWorkflowEngine.addCel("./test_strCelFileName","./test_strSampleFileName","./test_strResultFileName");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("cels")=="./test_strCelFileName");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("arrs")=="./test_strSampleFileName");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("result-files")=="./test_strResultFileName");
}
void CNWorkflowEngineTest::defineStateTest()
{
    /* State not defined till check options called */
    /*
	Verbose::out(1, "****CNWorkflowEngineTest::defineStateTest****");
	CNWorkflowEngine m_objCNWorkflowEngine;
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("time-start")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("time-end")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("free-mem-at-start")=="0");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("gender-method")=="cn-probe-chrXY-ratio");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("genotype-analysis")=="");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOpt("qmethod-spec")=="");
    CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("recommended-reference-sample-count")==44);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("recommended-reference-female-sample-count")==15);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("recommended-reference-male-sample-count")==15);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("reference-sample-count")==0);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("reference-female-sample-count")==0);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("reference-male-sample-count")==0);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptInt("reference-unknown-sample-count")==0);
	CPPUNIT_ASSERT(m_objCNWorkflowEngine.getOptDouble("yTarget")==-0.56747);
    */
}
void CNWorkflowEngineTest::checkOptionsTest()
{
	//FATAL ERROR: Must specify an output directory.
	Verbose::out(1, "****CNWorkflowEngineTest::checkOptionsTest****");
	CNWorkflowEngine m_objCNWorkflowEngine1;
    m_objCNWorkflowEngine1.setOpt("out-dir","");
    NEGATIVE_TEST(m_objCNWorkflowEngine1.checkOptions(), Except);

     //FATAL ERROR: Please specify either a reference-input file or a reference-output file, but not both.
	CNWorkflowEngine m_objCNWorkflowEngine2;
	m_objCNWorkflowEngine2.setOpt("out-dir","./Cyto");
    m_objCNWorkflowEngine2.setOpt("reference-input",INPUT + "/CNReference.a5");
	m_objCNWorkflowEngine2.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	NEGATIVE_TEST(m_objCNWorkflowEngine2.checkOptions(), Except);

    //FATAL ERROR: Please specify either a reference-input file or a reference-output
	CNWorkflowEngine m_objCNWorkflowEngine3;
	m_objCNWorkflowEngine3.setOpt("out-dir","./Cyto");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine3.isOptDefined("reference-input"));
	CPPUNIT_ASSERT(m_objCNWorkflowEngine3.isOptDefined("reference-output"));
	CPPUNIT_ASSERT(m_objCNWorkflowEngine3.isOptDefined("out-dir"));
	NEGATIVE_TEST(m_objCNWorkflowEngine3.checkOptions(), Except);     

    //FATAL ERROR:  The reference-input specified either does not exist, or is not the correct type for the CNWorkflowEngine.
	CNWorkflowEngine m_objCNWorkflowEngine4;
	m_objCNWorkflowEngine4.setOpt("reference-input",INPUT + "/CNReference_1.a5");
	m_objCNWorkflowEngine4.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine4.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	m_objCNWorkflowEngine4.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs");
	m_objCNWorkflowEngine4.setOpt("chrX-probes",INPUT + "/GenomeWideSNP_6.chrXprobes");
	m_objCNWorkflowEngine4.setOpt("chrY-probes",INPUT + "/GenomeWideSNP_6.chrYprobes");
	std::vector<std::string> values;
	values.push_back(INPUT + "/NA06985_GW6_C.CEL");
	m_objCNWorkflowEngine4.setOpt("cels",values);
	NEGATIVE_TEST(m_objCNWorkflowEngine4.checkOptions(),Except);

    //FATAL ERROR: Option create-reference cannot be found in the options for this engine.
	CNWorkflowEngine m_objCNWorkflowEngine5;
	m_objCNWorkflowEngine5.setOpt("reference-input",INPUT + "CNReference.a5");
    // cannot check this till after checkOptions is called
    //CPPUNIT_ASSERT(m_objCNWorkflowEngine5.getOptBool("create-reference")==false);
	//NEGATIVE_TEST(m_objCNWorkflowEngine5.getOpt("create-reference"),Except);
    
	
	//FATAL ERROR: The cels option and the cel-files option cannot both be used together. Choose one or the other.
	CNWorkflowEngine m_objCNWorkflowEngine6;
	CPPUNIT_ASSERT(m_objCNWorkflowEngine6.isOptDefined("cel-files"));
	CPPUNIT_ASSERT(m_objCNWorkflowEngine6.isOptDefined("cels"));
    m_objCNWorkflowEngine6.setOpt("cel-files",INPUT + "NA06985_GW6_C.CEL");
	m_objCNWorkflowEngine6.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	POSITIVE_TEST(m_objCNWorkflowEngine6.setOpt("cels",values));
	NEGATIVE_TEST(m_objCNWorkflowEngine6.checkOptions(), Except);
	
	CNWorkflowEngine m_objCNWorkflowEngine7;
    CPPUNIT_ASSERT(m_objCNWorkflowEngine7.isOptDefined("set-analysis-name"));
	CPPUNIT_ASSERT(m_objCNWorkflowEngine7.getOpt("set-analysis-name")=="");
	m_objCNWorkflowEngine7.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	POSITIVE_TEST(m_objCNWorkflowEngine7.setOpt("cels",values));
	POSITIVE_TEST(m_objCNWorkflowEngine7.getOpt("set-analysis-name")=="CN5");
	//FATAL ERROR: Must specify an annotation-file.
	NEGATIVE_TEST(m_objCNWorkflowEngine7.checkOptions(), Except);
    

	CNWorkflowEngine m_objCNWorkflowEngine8;
    m_objCNWorkflowEngine8.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	POSITIVE_TEST(m_objCNWorkflowEngine8.setOpt("cels",values));
	CPPUNIT_ASSERT(m_objCNWorkflowEngine8.isOptDefined("annotation-file"));
	m_objCNWorkflowEngine8.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	POSITIVE_TEST(m_objCNWorkflowEngine8.setOpt("cels",values));
	//FATAL ERROR: Must specify a CDF file
	NEGATIVE_TEST(m_objCNWorkflowEngine8.checkOptions(), Except); 

	//Must specify special-snps
    CNWorkflowEngine m_objCNWorkflowEngine8a;
	m_objCNWorkflowEngine8a.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	m_objCNWorkflowEngine8a.setOpt("cels",values);
	m_objCNWorkflowEngine8a.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	POSITIVE_TEST(m_objCNWorkflowEngine8a.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf"));
	NEGATIVE_TEST(m_objCNWorkflowEngine8a.checkOptions(), Except);

	//Must specify chrX-probes
    CNWorkflowEngine m_objCNWorkflowEngine9;
	m_objCNWorkflowEngine9.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	m_objCNWorkflowEngine9.setOpt("cels",values);
	m_objCNWorkflowEngine9.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine9.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	POSITIVE_TEST(m_objCNWorkflowEngine9.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs"));
	NEGATIVE_TEST(m_objCNWorkflowEngine9.checkOptions(), Except);
	
	//Must specify chrY-probes
    CNWorkflowEngine m_objCNWorkflowEngine10;
	m_objCNWorkflowEngine10.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	m_objCNWorkflowEngine10.setOpt("cels",values);
	m_objCNWorkflowEngine10.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine10.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	m_objCNWorkflowEngine10.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs");
	POSITIVE_TEST(m_objCNWorkflowEngine10.setOpt("chrX-probes",INPUT + "/GenomeWideSNP_6.chrXprobes"));
	NEGATIVE_TEST(m_objCNWorkflowEngine10.checkOptions(), Except);
	
	//FATAL ERROR: Must not specify target-sketch. It will be read from the Reference.
    CNWorkflowEngine m_objCNWorkflowEngine11;
    // Cannot check this till after checkOptions is called
	//POSITIVE_TEST(m_objCNWorkflowEngine11.getOpt("create-reference"));
	POSITIVE_TEST(m_objCNWorkflowEngine11.setOpt("target-sketch","1.0"));
	m_objCNWorkflowEngine11.setOpt("reference-input",INPUT + "/CNReference.a5");
	m_objCNWorkflowEngine11.setOpt("cels",values);
	m_objCNWorkflowEngine11.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine11.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	m_objCNWorkflowEngine11.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs");
	m_objCNWorkflowEngine11.setOpt("chrX-probes",INPUT + "/GenomeWideSNP_6.chrXprobes");
	m_objCNWorkflowEngine11.setOpt("chrY-probes",INPUT + "/GenomeWideSNP_6.chrYprobes");
	NEGATIVE_TEST(m_objCNWorkflowEngine11.checkOptions(), Except);
	
    //FATAL ERROR: Must not specify use-feat-eff. It will be read from the Reference.
    CNWorkflowEngine m_objCNWorkflowEngine11a;
	POSITIVE_TEST(m_objCNWorkflowEngine11a.setOpt("use-feat-eff","false"));
	m_objCNWorkflowEngine11a.setOpt("reference-input",INPUT + "/CNReference.a5");
	m_objCNWorkflowEngine11a.setOpt("cels",values);
	m_objCNWorkflowEngine11a.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine11a.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	m_objCNWorkflowEngine11a.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs");
	m_objCNWorkflowEngine11a.setOpt("chrX-probes",INPUT + "/GenomeWideSNP_6.chrXprobes");
	m_objCNWorkflowEngine11a.setOpt("chrY-probes",INPUT + "/GenomeWideSNP_6.chrYprobes");
    NEGATIVE_TEST(m_objCNWorkflowEngine11a.checkOptions(), Except);

	//FATAL ERROR: Must not specify read-model-brlmmp. It will be read from the Reference.
    CNWorkflowEngine m_objCNWorkflowEngine11b; 
	POSITIVE_TEST(m_objCNWorkflowEngine11b.setOpt("read-models-brlmmp","false"));
	m_objCNWorkflowEngine11b.setOpt("reference-input",INPUT + "/CNReference.a5");
	m_objCNWorkflowEngine11b.setOpt("cels",values);
	m_objCNWorkflowEngine11b.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine11b.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	m_objCNWorkflowEngine11b.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs");
	m_objCNWorkflowEngine11b.setOpt("chrX-probes",INPUT + "/GenomeWideSNP_6.chrXprobes");
	m_objCNWorkflowEngine11b.setOpt("chrY-probes",INPUT + "/GenomeWideSNP_6.chrYprobes");
    NEGATIVE_TEST(m_objCNWorkflowEngine11b.checkOptions(), Except);
    
    
	//FATAL ERROR: Please specify a qca-file
    CNWorkflowEngine m_objCNWorkflowEngine12;
	m_objCNWorkflowEngine12.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	m_objCNWorkflowEngine12.setOpt("cels",values);
	m_objCNWorkflowEngine12.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine12.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	m_objCNWorkflowEngine12.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs");
	m_objCNWorkflowEngine12.setOpt("chrX-probes",INPUT + "/GenomeWideSNP_6.chrXprobes");
	m_objCNWorkflowEngine12.setOpt("chrY-probes",INPUT + "/GenomeWideSNP_6.chrYprobes");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine12.isOptDefined("run-geno-qc"));
   	m_objCNWorkflowEngine12.setOpt("run-geno-qc","true");
    NEGATIVE_TEST(m_objCNWorkflowEngine12.checkOptions(), Except);
	
	//FATAL ERROR: Please specify a qcc-file.
    CNWorkflowEngine m_objCNWorkflowEngine13;
	m_objCNWorkflowEngine13.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	m_objCNWorkflowEngine13.setOpt("cels",values);
	m_objCNWorkflowEngine13.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine13.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	m_objCNWorkflowEngine13.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs");
	m_objCNWorkflowEngine13.setOpt("chrX-probes",INPUT + "/GenomeWideSNP_6.chrXprobes");
	m_objCNWorkflowEngine13.setOpt("chrY-probes",INPUT + "/GenomeWideSNP_6.chrYprobes");
	m_objCNWorkflowEngine13.setOpt("run-geno-qc","true");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine13.isOptDefined("qca-file"));
   	m_objCNWorkflowEngine13.setOpt("qca-file",INPUT + "/testing_qca.txt");
	NEGATIVE_TEST(m_objCNWorkflowEngine13.checkOptions(), Except);
    
	//FATAL ERROR: When specifying ARR files, the number of ARR files specfified must match the number of CEL files specified
    CNWorkflowEngine m_objCNWorkflowEngine14;
	m_objCNWorkflowEngine14.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	m_objCNWorkflowEngine14.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine14.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	m_objCNWorkflowEngine14.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs");
	m_objCNWorkflowEngine14.setOpt("chrX-probes",INPUT + "/GenomeWideSNP_6.chrXprobes");
	m_objCNWorkflowEngine14.setOpt("chrY-probes",INPUT + "/GenomeWideSNP_6.chrYprobes");
	m_objCNWorkflowEngine14.setOpt("cels",values);
	std::vector<std::string> arr;
	arr.push_back(INPUT + "/one.arr");
	arr.push_back(INPUT + "/one.arr");
	CPPUNIT_ASSERT(m_objCNWorkflowEngine14.isOptDefined("arrs"));
	m_objCNWorkflowEngine14.setOpt("arrs",arr);
	NEGATIVE_TEST(m_objCNWorkflowEngine14.checkOptions(), Except);
	
	//FATAL ERROR: When specifying CHP files, the number of chp file names specified with result-files must match the number of CEL files specified
    CNWorkflowEngine m_objCNWorkflowEngine15;
	m_objCNWorkflowEngine15.setOpt("reference-output",OUTPUT + "/CNReference.a5");
	m_objCNWorkflowEngine15.setOpt("annotation-file",INPUT + "/GenomeWideSNP_test_6.na23.annot_small.csv");
	m_objCNWorkflowEngine15.setOpt("cdf-file",INPUT + "/GenomeWideSNP_6.cdf");
	m_objCNWorkflowEngine15.setOpt("special-snps",INPUT + "/GenomeWideSNP_6.specialSNPs");
	m_objCNWorkflowEngine15.setOpt("chrX-probes",INPUT + "/GenomeWideSNP_6.chrXprobes");
	m_objCNWorkflowEngine15.setOpt("chrY-probes",INPUT + "/GenomeWideSNP_6.chrYprobes");
	m_objCNWorkflowEngine15.setOpt("cels",values);
	std::vector<std::string> chip;
	chip.push_back(INPUT + "/one.chip");
	chip.push_back(INPUT + "/two.chip");
	chip.push_back(INPUT + "/three.chip");
	POSITIVE_TEST(m_objCNWorkflowEngine15.isOptDefined("result-files"));
	m_objCNWorkflowEngine15.setOpt("result-files",arr);
	NEGATIVE_TEST(m_objCNWorkflowEngine15.checkOptions(), Except);
	values.clear();
		
}
