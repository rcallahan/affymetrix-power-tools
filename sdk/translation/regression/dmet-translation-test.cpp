////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
 * @file   dmet-translation-test.cpp
 * @author Mybrid Spalding
 * @date   Wed Jun  4 08:25:55 PDT 2008
 * @brief  Regression tests for apt-dmet-translation.cpp.
 */

//
#include "translation/RegressionExperimentReport.h"
#include "translation/RunTimeEnvironment.h"
#include "translation/regression/ATDRegression.h"
#include "translation/regression/DMET2Check.h"
#include "translation/regression/DMET3TestCaseCheck.h"
#include "translation/regression/DMET3TestCaseTableModel.h"
//
#include "calvin_files/utils/src/FileUtils.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/FsTestDir.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/RegressionTest.h"
#include "util/Util.h"
//
#include "pcrecpp.h"
//
#include <cassert>
#include <sstream>

//
#define NDEBUG

using namespace std;

/*****************************************************************************/
/**
 * displayHelpAndExit:
 * Synopsis: Calls the PGOptions usage method.
 *
 * @param opts - PgOptions for printing the usage.
 * return - exits program.
 *
 */
/*****************************************************************************/
void displayHelpAndExit(PgOptions * opts)
{

  assert(opts);

  opts->usage();

  cout << endl << endl;

  exit(0);

}
// end displayHelpAndExit
/*****************************************************************************/
/**
 * getOptions:
 * Synopsis: Fill the PGOptions data structure with our options. Note that
 * this program is intended to run with no options via "make regression".
 * If there are options they should be passed in via the Makefile.
 *
 * @param name1 - description
 * @param name1 - description
 * @return - description
 */
/*****************************************************************************/
void getOptions(int argc, char * argv[],   RunTimeEnvironment & rte,
                const FsTestDir & testDir)
{

  assert(argv);

  PgOptions  * pgOpts        = new PgOptions();

  pgOpts->defineOption("i", "input-dir", PgOpt::STRING_OPT,
                       "Input regression data directory", "");

  // Parse options and store into the opts data structure.
  pgOpts->parseArgv(argv);

  rte.m_adtOpts.m_progName        = pgOpts->getProgName();
  
  int i                      = rte.m_adtOpts.m_progName.rfind("/") ;
  rte.m_adtOpts.m_progName        = rte.m_adtOpts.m_progName.substr(i > 0 ? i + 1 : 0);

  rte.m_adtOpts.m_inputDir
  = pgOpts->get("input-dir").empty() ? INPUT_DIR : pgOpts->get("input-dir");

  rte.m_adtOpts.m_outputDir       = testDir.asString();

  rte.m_adtOpts.m_verbosity       = ADT_VERBOSE_EXCEPTION;

  Verbose::setLevel(rte.m_adtOpts.m_verbosity);

  return;

}
// end getOptions
/*****************************************************************************/
/*****************************************************************************/
/**
 * getDMET2AlgorithmChecks:
 * Synopsis: Get the DME2 Genotype files and their corresponding DMET2 files
 *           for regression comparision
 *
 *
 *
 *
 * @param rte - the run time environment
 * @param cmd -  for return
 * @return - The vector of RegressionCheck objects as well as a command to run them.
 */
/*****************************************************************************/
map< string, vector< RegressionCheck *> > getDMET2AlgorithmChecks(RunTimeEnvironment & rte)
{

  map< string, vector< RegressionCheck *> > checks;

  list <string> searchFilesList
  = affymetrix_calvin_utilities::
    FileUtils::ListFiles(rte.m_adtOpts.m_inputDir.c_str(), "txt");

  list<string>::iterator itS;

  pcrecpp::RE dmet2RE("DMET2_Genotypes_Short");

  for (itS = searchFilesList.begin(); itS != searchFilesList.end(); itS++) {

    string genoFile = *itS;

    if (dmet2RE.PartialMatch(genoFile)) {

      // Replace the input directory with the output directory.
      string genoFileRoot = Fs::basename(genoFile);

      Verbose::out(ADT_VERBOSE_TMI, genoFile);

      string dmet2MarkerFile = genoFile + "." + DMET2_REGRESSION_MARKER_FILE_EXT;
      string dmet3MarkerFile = Fs::join(rte.m_adtOpts.m_outputDir,
                                        genoFileRoot+"."+RegressionExperimentReport::getDmet3MarkerFileNameExt());

      string dmet2HaplotypeFile = genoFile + "." + DMET2_REGRESSION_HAPLOTYPE_FILE_EXT;
      string dmet3HaplotypeFile = Fs::join(rte.m_adtOpts.m_outputDir,
                                           genoFileRoot+"."+RegressionExperimentReport::getDmet3HaplotypeFileNameExt());

      if (! Fs::isReadable(dmet2MarkerFile)) {
        APT_ERR_ABORT(genoFile + ": missing DMET2 marker regression file.");
      }
      if (! Fs::isReadable(dmet2HaplotypeFile)) {
        APT_ERR_ABORT(genoFile + ": missing DMET2 haplotype regression file.");
      }
      // HACK: in reality the Marker and Haplotype files are created
      // in the same call. The RegressionTest object doesn't support this,
      // so we'll invoke the program twice and pretend the invocations were
      // different.
      
      string cmd;
      cmd = Fs::join(".","apt-dmet-translation");
      cmd = cmd + " -v 2 ";
      cmd = cmd + " --marker-report --dmet2-calling -r 1 -o " + rte.m_adtOpts.m_outputDir  + " -i " ;
      cmd = cmd + rte.m_adtOpts.m_inputDir.c_str() + " -g " + genoFile;

      checks[cmd].push_back(new DMET2MarkerCheck(rte, dmet2MarkerFile, dmet3MarkerFile, genoFile));
      checks[cmd].push_back(new DMET2HaplotypeCheck(rte, dmet2HaplotypeFile, dmet3HaplotypeFile , genoFile));

    }

  }

  return checks;
}

// end getDMET2AlgorithmChecks
/*****************************************************************************/
/*****************************************************************************/
/**
 * runDmet2AlgorithmRegression:
 * Synopsis: Run algorithem regression 
 *
 *
 * @param name1 - description
 * @param name1 - description
 * @return - description
 */
/*****************************************************************************/
bool runDmet2AlgorithmRegression(RunTimeEnvironment & rte)
{

  bool allPass = true;
  string cmd;

  map< string, vector< RegressionCheck *> > checks = getDMET2AlgorithmChecks(rte);
  map< string, vector< RegressionCheck *> >::iterator itSVRC;

  Verbose::out(ADT_VERBOSE_NORMAL, "Running DMET2 tests...");

  int i = 0;
  for (itSVRC = checks.begin(); itSVRC != checks.end(); itSVRC++) {
    i++;

    Verbose::out(ADT_VERBOSE_NORMAL, itSVRC->first);
    std::string name = "DMET2-" + ToStr(i);
    RegressionTest dmet2Test(name.c_str(), itSVRC->first.c_str(), itSVRC->second);

    int pass = dmet2Test.pass();

    if (!pass) {
      Verbose::out(ADT_VERBOSE_NORMAL, dmet2Test.getErrorMsg(), false);
    }

    allPass = allPass && pass;

  }

  if (!allPass) {
    Verbose::out(ADT_VERBOSE_NORMAL, "Failed DMET2 algorithm tests!");
  } else {
    Verbose::out(ADT_VERBOSE_NORMAL, "Passed DMET2 algorithm tests!");
  }

  return allPass;

}
// end runDmet2AlgorithmRegression
/*****************************************************************************/
/*****************************************************************************/
/**
 * getDMET3ReportChecksFilePath:
 * Synopsis:
 * Helper function for getDMET3ReportChecks.
 *
 * The file path can either be
 * 1.) absolute or relative to the current path (AS IS).
 *
 * 2.) In the corresponding test case directory. Take the test case
 * file, strip the extension and look in the corresponding directory.
 * If the file exists in this test case directory then return the path
 * relative to the regression directory. 
 * 
 * *
 * @param tesCase - the test case CSV file
 * @param pathToCheck - the file from the test case CSV to get the path for.
 * @param pathToUse   - the return path
 * 
 * @return - true of the file exists, false otherwise, pathToUse is filled in.
 */
/*****************************************************************************/
static bool getDMET3ReportChecksFilePath( string testCase, const string & pathToCheck, string & pathToUse, bool returnInvalidPathToUse = false) {


  if ( Fs::isReadable(pathToCheck) ) {
    pathToUse = pathToCheck;
    return true;
  }
  
  string testCaseDir = testCase.substr(0,testCase.rfind("."));

  if ( Fs::isReadable(Fs::join(TEST_DATA_REGRESSION_DIR,pathToCheck))) {
    pathToUse = Fs::join(TEST_DATA_REGRESSION_DIR,pathToCheck);
    return true;
  }

  if ( Fs::isReadable(Fs::join(testCaseDir,pathToCheck))) {
    pathToUse = Fs::join(testCaseDir,pathToCheck);
    return true;
  }

  if ( returnInvalidPathToUse ) {
    pathToUse = Fs::join(TEST_DATA_REGRESSION_DIR  + testCaseDir,pathToCheck);
  }
    
  return false;

}
// end getDMET3ReportChecksFilePath
/*****************************************************************************/
/*****************************************************************************/
/**
 * getDMET3ReportChecksInitializeFiles
 * Synopsis:
 * Helper function for getDMET3ReportChecks.
 *
 * Initialize the files from the CSV file by possible adding
 * the test case directory. Make sure all files exist or abort. 
 *
 * 
 * *
 * @param testCase - the test case CSV file
 * @param pathToCheck - the file from the test case CSV to get the path for.
 * @param pathToUse   - the return path
 * 
 * @return - true of the file exists, false otherwise, pathToUse is filled in.
 */
/*****************************************************************************/
static void getDMET3ReportChecksInitializeFiles( string testCase, const vector< string > testCaseRow, DMET3ReportTestCheckFiles & returnFiles ) {

  
  if ( ! getDMET3ReportChecksFilePath(testCase, testCaseRow[DMET3TestCaseTableModel::CHP_FILES], returnFiles.m_chpFiles) ) {
    APT_ERR_ABORT(testCaseRow[DMET3TestCaseTableModel::CHP_FILES] +  string(" file not found.") );
  }
  if ( ! getDMET3ReportChecksFilePath(testCase, testCaseRow[DMET3TestCaseTableModel::TRANSLATION_FILE], returnFiles.m_translationFile) ) {
    APT_ERR_ABORT(testCaseRow[DMET3TestCaseTableModel::TRANSLATION_FILE] +  string(" file not found.") );
  }
  if ( ! getDMET3ReportChecksFilePath(testCase, testCaseRow[DMET3TestCaseTableModel::ANNOTATION_FILE], returnFiles.m_annotationFile) ) {
    APT_ERR_ABORT(testCaseRow[DMET3TestCaseTableModel::ANNOTATION_FILE] +  string(" file not found.") );
  }
  if ( !testCaseRow[DMET3TestCaseTableModel::MARKER_FILE].empty() &&
       ! getDMET3ReportChecksFilePath(testCase, testCaseRow[DMET3TestCaseTableModel::MARKER_FILE], returnFiles.m_markerFile)) {
    APT_ERR_ABORT(testCaseRow[DMET3TestCaseTableModel::MARKER_FILE] +  string(" file not found.") );
  }
  if ( !testCaseRow[DMET3TestCaseTableModel::OVERRIDE_FILE].empty() &&
       !getDMET3ReportChecksFilePath(testCase, testCaseRow[DMET3TestCaseTableModel::OVERRIDE_FILE], returnFiles.m_overrideFile)) {
    APT_ERR_ABORT(testCaseRow[DMET3TestCaseTableModel::OVERRIDE_FILE] +  string(" file not found.") );
  }

  vector< std::string > reports;
  reports.push_back("comprehensive");
  reports.push_back("uncalled");
  reports.push_back("summary");

  for ( int i =0; i < reports.size(); i++ ) {
    
    // BUILD EXPECTED RESULT FILES
    if ( pcrecpp::RE(reports[i]).PartialMatch(testCaseRow[DMET3TestCaseTableModel::EXPECTED_RESULT_REPORTS]) ) {
      string test = testCaseRow[DMET3TestCaseTableModel::BASE_NAME] +  string("_expected_results_") + reports[i] + string(".rpt");
      string result;
      if ( ! getDMET3ReportChecksFilePath( testCase, test, result, true ) ) {
        APT_ERR_ABORT( result + string(" file not found.") );
      }
      returnFiles.m_resultFiles[reports[i]] = result;
    }
  }
  
}
// end getDMET3ReportChecksInitializeFiles
/*****************************************************************************/
/**
 * getDMET3ReportChecks:
 * Synopsis: Get the DMET3 test case CSV files. 
 *
 *
 *
 * @param rte - the run time environment
 * 
 * @return - The vector of RegressionCheck objects 
 */
/*****************************************************************************/
map< string, vector< RegressionCheck *> > getDMET3ReportChecks(RunTimeEnvironment & rte) {
  string cmd;
  
  map<  string, vector< RegressionCheck *> >checks;
  
  string comprehensive_dir = Fs::join(rte.m_adtOpts.m_inputDir,"DMET3");

  /// @todo add this to Fs.
  list <string> searchFilesList
  = affymetrix_calvin_utilities::
    FileUtils::ListFiles(comprehensive_dir.c_str(), "csv");

  list<string>::iterator sflIt1 = searchFilesList.begin();
  
  for (; sflIt1 != searchFilesList.end(); ++sflIt1) {
    Verbose::out(ADT_VERBOSE_NORMAL,string("DMET3 Test Case File: ") + *sflIt1);

    DMET3TestCaseTableModel testCases(rte, *sflIt1 );

    for (int i =0; i < testCases.size(); i++ ) {

      DMET3ReportTestCheckFiles tcf;
      
      getDMET3ReportChecksInitializeFiles( *sflIt1, testCases.m_rows[i], tcf);
      
      cmd = Fs::join(".",TRANSLATION_PROGRAM);
      //cmd = cmd + string(" -v -1 ");
      if ( pcrecpp::RE("[xX][mM][lL]$").PartialMatch(tcf.m_chpFiles) ) {
        cmd = cmd + string(" --xml-file ") + tcf.m_chpFiles;
      }
      else {
        cmd = cmd + string(" -b ") + testCases.m_rows[i][DMET3TestCaseTableModel::BASE_NAME];
        cmd = cmd + string(" -o ") + rte.m_adtOpts.m_outputDir + string(" ") ;
        cmd = cmd + string(" -a ") + tcf.m_annotationFile;
        cmd = cmd + string(" -t ") + tcf.m_translationFile;
        if ( pcrecpp::RE("[tT][Xx][tT]$").PartialMatch(tcf.m_chpFiles) ) {
          cmd = cmd + string(" -E ") + tcf.m_chpFiles;
        }
        else {
          cmd = cmd + string(" -e ") + tcf.m_chpFiles;
        }
        if ( !tcf.m_markerFile.empty() ) {
          cmd = cmd + string(" -m ") + tcf.m_markerFile;
        }
        if ( !tcf.m_overrideFile.empty() ) {
          cmd = cmd + string(" -n ") + tcf.m_overrideFile;
        }
      }
      
      cerr << "CMD: " << endl;
      checks[cmd].push_back(new DMET3TestCaseCheck(rte, *sflIt1, testCases.m_rows[i], tcf ));
    }
  }


  return checks;
  
}
// end getDMET3ReportChecks
/*****************************************************************************/
/*****************************************************************************/
/**
 * runDmet3ReportRegression:
 * Synopsis: 
 * 
 * 
 * 
 *
 * @param rte - the run time envrionment 
 * 
 * @return - true if all tests passed, false otherwise. 
 */
/*****************************************************************************/
bool runDmet3ReportRegression(RunTimeEnvironment & rte) {

  bool allPass = true;
  string cmd;

  map< string, vector< RegressionCheck *> > checks = getDMET3ReportChecks(rte);
  map< string, vector< RegressionCheck *> >::iterator itSVRC;

  Verbose::out(ADT_VERBOSE_NORMAL, string("Running DMET3 Report Test Cases (")  + ToStr(checks.size()) + string(")..."));

  int i=0;
  for (itSVRC = checks.begin(); itSVRC != checks.end(); itSVRC++) {
    i++;

    Verbose::out(ADT_VERBOSE_NORMAL, itSVRC->first);

    std::string name = "DMET2-" + ToStr(i);
    RegressionTest dmet3TestCase(name.c_str(), itSVRC->first.c_str(), itSVRC->second);

    int pass = dmet3TestCase.pass();

    allPass = allPass && pass;

  }

  if (!allPass) {
    Verbose::out(ADT_VERBOSE_NORMAL, "Failed DMET3 report tests!");
  } else {
    Verbose::out(ADT_VERBOSE_NORMAL, "Passed DMET3 report tests!");
  }

  return allPass;
  

  return true;
  
}
// end runDmet3ReportRegression
/*****************************************************************************/
/*****************************************************************************/
/**
 * main
 * Synopsis: main control loop.
 *
 * @param argc - command line
 * @param argv - command line
 * @return - 0 ok, 1 on error
 */
/*****************************************************************************/
int main(int argc, char* argv[])
{
  try {

    bool okToContinue = true;

    string timeStamp(Util::getTimeStamp());

    FsTestDir testDir;
    testDir.setTestDir("translation", true);

    RunTimeEnvironment  *rte = new RunTimeEnvironment();

    getOptions(argc, argv, *rte, testDir);
    
    rte->initializeRunTimeEnvironment();
    
    Err::setThrowStatus(true);
    
    string verboseLogName = Fs::join(rte->m_adtOpts.m_outputDir,"dmet-translation-test.log");
    
    std::ofstream *verboseOFStream = new std::ofstream;
    
    if (!Fs::isWriteableDir(rte->m_adtOpts.m_outputDir)) {
      if (Fs::mkdirPath(rte->m_adtOpts.m_outputDir, false) != APT_OK) {
        APT_ERR_ABORT("Can't make or write to directory: " + rte->m_adtOpts.m_outputDir);
      }
    }
    
    Fs::mustOpenToWrite(*(verboseOFStream), verboseLogName.c_str());
    
    ADT_VERBOSE_ENUM logVerbosity = ADT_VERBOSE_EXCEPTION;
    if (rte->m_adtOpts.m_verbosity > ADT_VERBOSE_EXCEPTION) {
        logVerbosity = (ADT_VERBOSE_ENUM) rte->m_adtOpts.m_verbosity;
    }
    
    LogStream *verboseLogStream = new
    LogStream(logVerbosity,  verboseOFStream);
    
    Verbose::pushMsgHandler(verboseLogStream);
    Verbose::pushProgressHandler(verboseLogStream);
    Verbose::pushWarnHandler(verboseLogStream);
    Verbose::setLevel(rte->m_adtOpts.m_verbosity);
    
    timeStamp = Util::getTimeStamp();
    
    rte->m_programTime.begin();
    Verbose::out(ADT_VERBOSE_NORMAL, "[" + rte->m_programName  + " " + timeStamp + "]  BEGIN ");
    
    
    okToContinue = runDmet2AlgorithmRegression(*rte) && okToContinue;
    okToContinue = runDmet3ReportRegression(*rte) && okToContinue;
    
    
    timeStamp = Util::getTimeStamp();
    
    if (!okToContinue) {
        Verbose::out(ADT_VERBOSE_NORMAL, "Failed regression tests!");
    } else {
        Verbose::out(ADT_VERBOSE_NORMAL, "Passed ALL regression tests!");
    }
    
    Verbose::out(ADT_VERBOSE_NORMAL, "[" + rte->m_programName  + " " + timeStamp + "]  END ");
    
    delete rte;
    rte = NULL;
    
    if (!okToContinue) {
        exit(1);
    }
    
    exit(0);
    
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
// end main
/*****************************************************************************/

