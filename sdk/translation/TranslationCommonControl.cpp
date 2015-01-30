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
 * @file   TranslationCommonControl.cpp
 * @author Mybrid Spalding
 * @date   Tue Jul  8 09:26:28 PDT 2008
 * @brief  Controller code common to both command line and console controllers. 
 */

//
#include "calvin_files/utils/src/FileUtils.h"
#include "file/TsvFile/TsvFile.h"
#include "translation/ADTOptions.h"
#include "translation/CopyNumberTableModel.h"
#include "translation/ExperimentFileTableModel.h"
#include "translation/GenotypeOverrideTableModel.h"
#include "translation/GenotypeTableModel.h"
#include "translation/MarkerListModel.h"
#include "translation/SampleInfoTableModel.h"
#include "translation/TranslationCommonControl.h"
#include "translation/TranslationEngine.h"
#include "translation/TranslationTableModel.h"
#include "util/Err.h"    
#include "util/Fs.h"
#include "util/Verbose.h"
//
#include "pcrecpp.h"
//
#include <cstring>
#include <string>
//
using namespace affx;
using namespace std;

const std::string TRANSLATION_OPT_GENO_DMET2_FILE = \
    "(DMET2: 1 of 2 experiment files) => YYYYMMDD_*_DMET[23]_Genotypes*.txt";
const std::string EXCEPTION_MISSING_TRANSLATION_TABLE
= "Invalid option: translation file is required (i.e. -t, --translate-file option).";

const std::string EXCEPTION_MISSING_EXPERIMENT_CHP_FILE
= "Invalid option: experiment CHP files is required (i.e. -e, --experiment-list option).";

const std::string HAPLOTYPE_REPORT_OPTION_FIRST = "Reporting only first named haplotype among non-differentiable haplotype names in translation file.";

const std::string HAPLOTYPE_REPORT_OPTION_DUPLICATE =  "Reporting a combined haplotype name based on all non-differentiable haplotype names in translation file.";


/*****************************************************************************/
/**
 * TranslationCommonControl::defineCommonOptions
 * Synopsis:
 *
 *  Define the space of possible PgOptions for command line or console. 
 *  
 *
 * @param opts - newly allocated PgOptions from either console or command line.
 *
 */
/*****************************************************************************/
void TranslationCommonControl::defineCommonOptions(TranslationEngine *engine)
{


  APT_ERR_ASSERT(engine, "");

  // Program Control
  engine->defineOptionSection("Program Control");
  engine->defineOption("", "dmet2-calling", PgOpt::BOOL_OPT,
                     "=> use DMET2 logic for making allele calls." , "false");
  engine->defineOption("", "enforce-complete-haplotypes", PgOpt::BOOL_OPT,
                     "=> Enforce that experiment data has a complete set of markers for any haplotype set and abort when incomplete haplotype sets are detected." , "false");
  engine->defineOption("", "ignore-unknown-alleles", PgOpt::BOOL_OPT,
                     "=> Ignore unknown alleles for a multi-allelic marker where the allele is not specified in the translation file. " , "false");
  engine->defineOption("", "ignore-report-allele", PgOpt::BOOL_OPT,
                     "=> Use the alleles from the translation file for reporting and not the report allele in the annotation file. ", "false");

  engine->defineOption("", "audit", PgOpt::STRING_OPT,
                     "=> Audit various components of translation, don't translate.\nAvailable audit types are:\n * annotation - audit annotation file probe sets against the translation file.", "");

  engine->defineOption("", "use-first-dup-allele-def", PgOpt::BOOL_OPT,
                     "=> If marker list filtering removes markers needed to differntiate among multiple possible haplotypes then report only the first named haplotype in the Translation file.", "false");


  // Input Section
  engine->defineOptionSection("Input Options");

  engine->defineOption("a", "annotation-file", PgOpt::STRING_OPT,
                     "=> file.csv (probe set annotation file, CSV)", "");

  engine->defineOption("c", "copy-dmet2-number-file", PgOpt::STRING_OPT,
                     "(DMET2: 2 of 2 experiment files) => YYYYMMDD_*_DMET[23]_cn_*.txt", "");
  engine->defineOption("e", "experiment-list", PgOpt::STRING_OPT,
                     "=> file1.chp;file2.chp;... OR file1.chp:file2.chp:..", "");
  engine->defineOption("E", "experiment-list-file", PgOpt::STRING_OPT,
                     "=> file with listing of CHP files. (Commented out lines starting with '#' are allowed).", "");
  engine->defineOption("g", "geno-dmet2-type-file", PgOpt::STRING_OPT,
                     TRANSLATION_OPT_GENO_DMET2_FILE, "");
  engine->defineOption("i", "input-dir", PgOpt::STRING_OPT,
                     "=> Directory to read input files from when file is not specified (-c, -g, -t_). The current directory is assumed when not specified." , "");
  engine->defineOption("m", "marker-list-file", PgOpt::STRING_OPT,
                     "=> marker-list-file.*.txt file is a simple list of probe sets, one per line, to be used for translation.", "");
  engine->defineOption("n", "genotype-override-file", PgOpt::STRING_OPT,
                     "=> file with genotype override base calls in TSV format." , "");

  engine->defineOption("s", "sample-file", PgOpt::STRING_OPT,
                     "=> file.txt (sample info file, TSV)" , "");


  engine->defineOption("t", "translate-file", PgOpt::STRING_OPT,
                     "=> DMET3_TTTable_vYYYYMMDD_*.txt (DMET2_TTTable_vYYYYMMDD_*.txt for DMET2 testing)", "");

  engine->defineOption("", "metabolizer-file", PgOpt::STRING_OPT,
					 "The metabolizer bin file for translating allele calls to metabolizer genotypes", "");
  
  // Output Section
  engine->defineOptionSection("Output Options");
  engine->defineOption("b", "base-name-prefix", PgOpt::STRING_OPT,
                     "=> Base name prefix of all output reports.", "");

  // Reporting Section
  engine->defineOptionSection("Reporting Options");
  engine->defineOption("p", "profile", PgOpt::BOOL_OPT,
                     " => Profile reporting.", "false");

  engine->defineOption("r", "regression", PgOpt::INT_OPT,
                     " => Regression log level.", "0");

  engine->defineOption("", "marker-report", PgOpt::BOOL_OPT,
                     " => Output DMET2 Marker Report.", "false");

  engine->defineOption("", "summary-report-sort", PgOpt::BOOL_OPT,
                     "=> Sort the summary report rather than use the default translation table order.", "true");

  engine->defineOption("u", "uncalled-report-all-markers", PgOpt::BOOL_OPT,
                     "=> Report all markers in the uncalled report, not just the uncalled markers.","false");
  

  return;

}
// end TranslationCommonControl::defineCommonOptions
/*****************************************************************************/
/*****************************************************************************/
/**
 * Synopsis: Convert the Copy Number Report, Genotype Short Report
 * and Translation Tablem input files into in memory models.
 * Note that all the input paramaters are for return. 
 *
 * 
 * @param rte         - single instance run-time environment data model
 * @param returnMlm   - marker list single instance data model 
 * @param returnCnrm  - DMET2 only copy number list single instance data model
 * @param returnGtm   - Genotype data (CHP) single instance data model. 
 * @param returnTtm   - Translation table data single instance data mode.
 * @param returnSitm  - Sample info single instance data model.
 * @param returnGotm  - Genotype override single instance data model. 
 * 
 * @return true - all models where successfully allocated on the heap. 
 * 
 */
/*****************************************************************************/
bool TranslationCommonControl::initializeInputFileModels(RunTimeEnvironment &rte,
    MarkerListModel             **returnMlm,
    CopyNumberTableModel        **returnCnrm,
    GenotypeTableModel          **returnGtm,
    TranslationTableModel       **returnTtm,
    SampleInfoTableModel        **returnSitm,
    GenotypeOverrideTableModel  **returnGotm)
{

  std::stringstream msgSStr;
  APT_ERR_ASSERT(returnCnrm && returnGtm && returnTtm, "");

  bool okToContinue = true;

  APT_ERR_ASSERT(!*returnCnrm && !*returnGtm && !*returnTtm && !*returnMlm, "");


  // TRANSLATION TABLE
  *returnTtm = new TranslationTableModel(rte, rte.m_adtOpts.m_inputTTableFile, rte.m_adtOpts.m_inputTTableType);
  msgSStr <<  "Input translation file records read:      " << (*returnTtm)->size() ;

  Verbose::out(rte.m_currentVerbosity, msgSStr.str());
  msgSStr.str("");


  // COPY NUMBER, DMET2 ONLY
  if (rte.m_adtOpts.m_dmet2Calling) {
    *returnCnrm = new CopyNumberTableModel(rte, rte.m_adtOpts.m_inputCopyFile);
    msgSStr <<  "Input copy number records read:     " << (*returnCnrm)->size();
    Verbose::out(rte.m_currentVerbosity, msgSStr.str());
  } else {
    *returnCnrm = new CopyNumberTableModel();
  }

  msgSStr.str("");

  // MARKER LIST
  // Translation Table Filter.
  // Modify the Translation Table and filter the appropriate probe sets
  // and to merge the allele names for shared haplotype calls.
  if (*returnTtm && (! rte.m_adtOpts.m_inputMarkerListFile.empty())) {
    *returnMlm = new MarkerListModel(rte, rte.m_adtOpts.m_inputMarkerListFile);
  } else if (*returnTtm && (rte.m_adtOpts.m_inputProbeSetVector.size() > 0)) {
    *returnMlm = new MarkerListModel(rte);
  }
  
  if (*returnMlm) {
    msgSStr << "Input marker list, probe set records read: "  << (*returnMlm)->size();
    Verbose::out(rte.m_currentVerbosity, msgSStr.str());

    okToContinue = okToContinue && (*returnTtm)->probeSetFilter(rte, **returnMlm);
    (*returnMlm)->describeVerbose(rte);
  }

  // TRANSLATION TABLE GENE MARKER INDEX
  // Generate the complete list of gene markers for those operations
  // that require gene marker sets.
  // Must be done after the MARKER LIST has filtered the translation table.
  if ( *returnTtm ) {
    (*returnTtm)->generateGeneMarkerList();
  }

  msgSStr.str("");

  // GENOTYPE OVERRIDE FILE
  if (returnGotm && (!rte.m_adtOpts.m_inputGenotypeOverrideFile.empty())) {

    *returnGotm = new GenotypeOverrideTableModel(rte, rte.m_adtOpts.m_inputGenotypeOverrideFile, **returnTtm);
    msgSStr << "Input genotype override records read:      " << (*returnGotm)->size();
    Verbose::out(rte.m_currentVerbosity, msgSStr.str());
  }

  msgSStr.str("");

  // GENOTYPE
  if (rte.m_adtOpts.m_streamType == ADT_EXPERIMENT_STREAM_TYPE_TSV) {
    *returnGtm = new GenotypeTableModel(rte, rte.m_adtOpts.m_inputGenoFile, returnGotm ? *returnGotm : NULL);
  } else {
    *returnGtm = new GenotypeTableModel(rte, returnGotm ? *returnGotm : NULL);
  }

  // SAMPLE INFO FILE
  if (returnSitm && (! rte.m_adtOpts.m_inputSampleFile.empty())) {

    *returnSitm = new SampleInfoTableModel(rte, rte.m_adtOpts.m_inputSampleFile);
    msgSStr << "Input sample info records read:            " << (*returnSitm)->size();
    Verbose::out(rte.m_currentVerbosity, msgSStr.str());
  } else if (returnSitm && (rte.m_adtOpts.m_sampleTable.size() != 0)) {
    *returnSitm = new SampleInfoTableModel(rte);
    msgSStr << "Input sample info records read:            " << (*returnSitm)->size();
    Verbose::out(rte.m_currentVerbosity, msgSStr.str());
  }

  msgSStr.str("");

  return okToContinue;

}
// end TranslationCommonControl::initializeInputFileModels
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationCommonControl::setCommonADTOptions:
 * Synopsis:
 * A quasie ADTOptions constructor. Perhaps it should be?
 * Copies all the passed in options to the ADTOptions single instance
 * object that itself is stashed in the rte. 
 *
 *
 * @param adtOpts - the ADTOptions to set up.
 * @param ptOpts  - from command line or console.
 *
 */
/*****************************************************************************/
void TranslationCommonControl::setCommonADTOptions(ADTOptions & adtOpts, class TranslationEngine *engine)
{
  // Copy into our friendly ADTOptions class.
  adtOpts.m_progName        = engine->getProgName();
  int i                      = adtOpts.m_progName.rfind("/") ;
  adtOpts.m_progName        = adtOpts.m_progName.substr(i ? i + 1 : 0);

  adtOpts.m_inputDir        = engine->getOpt("input-dir");

  adtOpts.m_dmet2Calling    = engine->getOptBool("dmet2-calling");
  adtOpts.m_enforceCompleteHaplotypeGroup
  = engine->getOptBool("enforce-complete-haplotypes");

  adtOpts.m_ignoreReportAllele = engine->getOptBool("ignore-report-allele");

  adtOpts.m_inputAnnotationFile
  = engine->getOpt("annotation-file");

  adtOpts.m_ignoreUnknownAlleles
  = engine->getOptBool("ignore-unknown-alleles");
  adtOpts.m_useFirstDupAlleleDef
  = engine->getOptBool("use-first-dup-allele-def");

  adtOpts.m_audit            = engine->getOpt("audit");

  
  adtOpts.m_inputExperimentList
  = engine->getOpt("experiment-list");
  adtOpts.m_inputExperimentListFile
  = engine->getOpt("experiment-list-file");


  adtOpts.m_inputTTableFile = engine->getOpt("translate-file");

  adtOpts.m_inputTTableType = TranslationTableModel::getTranslationTableFileType(adtOpts.m_inputTTableFile);

  adtOpts.m_inputGenoFile   = engine->getOpt("geno-dmet2-type-file");

  adtOpts.m_inputCopyFile   = engine->getOpt("copy-dmet2-number-file");

  adtOpts.m_inputMarkerListFile  = engine->getOpt("marker-list-file");

  if (! adtOpts.m_inputMarkerListFile.empty()) {
    if (adtOpts.m_useFirstDupAlleleDef) {
      adtOpts.m_useFirstDupAlleleDefHeaderText =  HAPLOTYPE_REPORT_OPTION_FIRST;
    } else {
      adtOpts.m_useFirstDupAlleleDefHeaderText =  HAPLOTYPE_REPORT_OPTION_DUPLICATE;
    }
  }

  adtOpts.m_inputGenotypeOverrideFile
  = engine->getOpt("genotype-override-file");

  adtOpts.m_inputSampleFile    = engine->getOpt("sample-file");

  adtOpts.m_outputDir = engine->getOpt("out-dir").empty()?".":engine->getOpt("out-dir");

  adtOpts.m_logFile       = engine->getOpt("log-file");

  adtOpts.m_outputReportPrefix = engine->getOpt("base-name-prefix");

  adtOpts.m_verbosity = engine->getOptInt("verbose") ? engine->getOptInt("verbose") : 1;

  adtOpts.m_profile = engine->getOptBool("profile");

  adtOpts.m_regression = engine->getOptInt("regression") ? engine->getOptInt("regression") : 0;

  adtOpts.m_markerReport = engine->getOptBool("marker-report");

  adtOpts.m_summaryReportSort = engine->getOptBool("summary-report-sort");

  adtOpts.m_uncalledReportAllMarkers = engine->getOptBool("uncalled-report-all-markers");
  

  return;

}
// end TranslationCommonControl::setCommonADTOptions
/*****************************************************************************/
/*****************************************************************************/
/**
 * _validateInputExperimentFilesFromListDMET3:
 * Synopsis:
 *
 * Helper function to _validateInputFilesDMET3, broken out for clarity.
 *
 * This routine reads a listing of files from a list passed in at the
 * command line, i.e. ( -e file1,file2,...,filen
 *
 *
 * @param rte     - the single instance run time environment that houses the ADTOptions object as m_adtOpts. 
 * @params pgOopts - command line options object
 *
 * @return - true if all input data was valid
 */
/*****************************************************************************/
static bool _validateInputExperimentFilesFromListDMET3(RunTimeEnvironment & rte, TranslationEngine* engine)
{

  bool validOptions = true;

  // Create the experiment file list from a command line list.
  pcrecpp::RE extensionRE("(\\.chp(?:[:,]|$))",
                          pcrecpp::RE_Options().set_caseless(true));

  pcrecpp::StringPiece inputTest(rte.m_adtOpts.m_inputExperimentList);
  int numFiles = 0;
  std::string experimentFile;

  while (extensionRE.FindAndConsume(&inputTest, &experimentFile)) {
    numFiles++;
  }

  if ((numFiles == 0) || !(inputTest.empty())) {
    cerr << endl << "Invalid option: " << rte.m_adtOpts.m_inputExperimentList << " invalid CHP experiment file name detected." << endl;
    validOptions  = false;
  } else if (numFiles == 1) {
    // Just one file.
    if (! Fs::fileExists(rte.m_adtOpts.m_inputExperimentList.c_str())) {
      cerr << endl << "Invalid option: " << rte.m_adtOpts.m_inputExperimentList << " file not found.\n";
      validOptions = false;
    } else {
      rte.m_adtOpts.m_inputExperimentFiles.push_back(rte.m_adtOpts.m_inputExperimentList);
    }
  } else {

    pcrecpp::RE colonSeparatorRE("\\.chp:", pcrecpp::RE_Options().set_caseless(true));
    pcrecpp::RE commaSeparatorRE("\\.chp,",
                                 pcrecpp::RE_Options().set_caseless(true));

    bool isColonSeparated     = false;
    bool isCommaSeparated = false;
    std::string separator;

    if (colonSeparatorRE.PartialMatch(rte.m_adtOpts.m_inputExperimentList)) {
      isColonSeparated = true;
      separator = ':';
    }
    if (commaSeparatorRE.PartialMatch(rte.m_adtOpts.m_inputExperimentList)) {
      isCommaSeparated = true;
      separator = ",";
    }

    if (isColonSeparated && isCommaSeparated) {
      cerr << endl << "Invalid input file options: colon and comma separators cannot be mixed in the 'experiment-list'." << endl;
      validOptions = false;
    } else if (!isColonSeparated && !isCommaSeparated) {
      cerr << endl <<  "Invalid input file options: invalid 'experiment-list'." << endl;
      validOptions = false;

    } else if (validOptions) {
      pcrecpp::RE endFileRE("([^" + separator + "]+)");

      pcrecpp::StringPiece input(rte.m_adtOpts.m_inputExperimentList);

      while (endFileRE.FindAndConsume(&input, &experimentFile)) {

        if (! Fs::fileExists(experimentFile.c_str())) {
          validOptions = false;
          cerr << endl << "Invalid option: " << experimentFile << " 'experiment-list' file not found."  << endl;
        }
        rte.m_adtOpts.m_inputExperimentFiles.push_back(experimentFile);

      }
      if (!input.empty()) {
        cerr << endl << "Invalid option: " << rte.m_adtOpts.m_inputExperimentList << " invalid CHP experiment file name detected. " << endl;
        validOptions = false;
      }
    }

  }

  return validOptions;

}
// end _validateInputExperimentFilesFromListDMET3
/*****************************************************************************/
/*****************************************************************************/
/**
 * _validateInputExperimentFilesFromFileDMET3:
 * Synopsis:
 *
 * Helper function to _validateInputFilesDMET3, broken out for clarity.
 *
 * This routine reads a listing of files from a file, i.e. -E chp_files.txt.
 *
 * @param rte     - the single instance run time environment that houses the ADTOptions object as m_adtOpts. 
 * @param pgOopts - command line options object
 *
 * @return - true if all input files were valid.
 */
/*****************************************************************************/
bool _validateInputExperimentFilesFromFileDMET3(RunTimeEnvironment & rte, TranslationEngine * engine)
{


  bool validOptions = true;
  TsvFile tsv;
  set<std::string> duplicateCheck;

  if (!  Fs::fileExists(rte.m_adtOpts.m_inputExperimentListFile.c_str())) {
    cerr << "Invalid option: " << rte.m_adtOpts.m_inputExperimentListFile << " file not found."  << endl;
    validOptions = false;
  } else {
    ExperimentFileTableModel eftm(rte, rte.m_adtOpts.m_inputExperimentListFile);

    rte.m_adtOpts.m_inputExperimentFiles = eftm.m_experimentFiles;
    validOptions = eftm.validateFiles();
  }

  return validOptions;

}
// end TranslationCommonControl::_validateInputExperimentFilesFromFileDMET3
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationCommonControl::_validateInputExperimentFilesFromVectorDMET3:
 * Synopsis:
 *
 * See validateInputFiles, broken out for clarity.
 *
 * This routine reads a listing of files from a std::vector passed in by the
 * console.
 *
 *
 * @param rte     - the single instance run time environment that houses the ADTOptions object as m_adtOpts. 
 * @param pgOopts - command line options object
 *
 * @return - true if all input std::vector data was valid.
 */
/*****************************************************************************/
bool _validateInputExperimentFilesFromVectorDMET3(RunTimeEnvironment & rte, TranslationEngine* engine)
{


  bool validOptions = true;
  set<std::string> duplicateCheck;

  pcrecpp::RE extensionRE("\\.chp$",
                          pcrecpp::RE_Options().set_caseless(true));
  std::stringstream msgSStr;

  std::vector< std::string > experimentFiles = engine->getOptVector("experiment-list-vector");

  for (int i = 0; i < experimentFiles.size(); i++) {

    if (!extensionRE.PartialMatch(experimentFiles[i])) {
      msgSStr << experimentFiles[i] << ": invalid CHP file name." << endl;
      validOptions = false;
    } else if (!  Fs::fileExists(experimentFiles[i].c_str())) {
      msgSStr << experimentFiles[i] << ": file not found." << endl;
      validOptions = false;
    } else {
      if (duplicateCheck.count(experimentFiles[i]) > 0) {
        cerr << experimentFiles[i] << ": duplicate file found, ignoring.\n";
      } else {
        duplicateCheck.insert(experimentFiles[i]);
        rte.m_adtOpts.m_inputExperimentFiles.push_back(experimentFiles[i]);
      }
    }
  }


  if (! validOptions) {
    Verbose::out(ADT_VERBOSE_NORMAL, msgSStr.str());
  } else if (rte.m_adtOpts.m_inputExperimentFiles.size() == 0) {
    cerr << "Invalid option: " << rte.m_adtOpts.m_inputExperimentListFile << " no valid CHP files found." << endl;
    validOptions = false;
  }

  return validOptions;

}
// end _validateInputExperimentFilesFromVectorDMET3
/*****************************************************************************/
/*****************************************************************************/
/**
 * _validateInputFilesDMET3:
 * Synopsis:
 *
 * Helper function to  _validateInputFiles, broken out for clarity.
 *
 * ADDITIONAL NOTES: The list of files can be either a colon or
 * comma separated list. We check by looking for ".chp," or
 * ".chp:" in the list (when more than file exists).
 *
 *
 * @param rte     - the single instance run time environment that houses the ADTOptions object as m_adtOpts. 
 * @params pgOopts - command line options object
 *
 * @return - true if all input files were valid.
 */
/*****************************************************************************/
static bool _validateInputFilesDMET3(RunTimeEnvironment & rte, TranslationEngine *engine)
{

  bool validOptions = true;


  // HACK: change this once the DMET3 table is finalized
  if ((rte.m_adtOpts.m_inputTTableType == ADT_TRANSLATION_TABLE_TYPE_DMET2) &&
      (!rte.m_adtOpts.m_inputAnnotationFile.empty())) {
    cerr << "Invalid option: DMET2 translation file passed for DMET3 calling.\n";
    validOptions = false;
  }

  if (rte.m_adtOpts.m_inputAnnotationFile.empty()) {
    cerr << "Invalid option: annotation file is required (i.e. -a, --annotation-file option)" << endl;
    validOptions = false;
  } else if (!  Fs::fileExists(rte.m_adtOpts.m_inputAnnotationFile.c_str())) {
    cerr << "Invalid option: " << rte.m_adtOpts.m_inputAnnotationFile << " probe set annotation file not found." << endl;
    validOptions = false;
  }

  if (rte.m_adtOpts.m_inputTTableFile.empty()) {
    cerr << EXCEPTION_MISSING_TRANSLATION_TABLE << endl;
    rte.m_invalidOptionsMask |= ADT_INVALID_OPTION_MASK_MISSING_TRANSLATION_TABLE;
    validOptions = false;
  } else if (rte.m_adtOpts.m_inputTTableType == ADT_TRANSLATION_TABLE_TYPE_INVALID) {
    validOptions = false;
    cerr << "Invalid option: " << rte.m_adtOpts.m_inputTTableFile << ": not a DMET3 tranlsation file."  << endl;
  } else {
    if (!  Fs::fileExists(rte.m_adtOpts.m_inputTTableFile.c_str())) {
      validOptions = false;
      cerr << rte.m_adtOpts.m_inputTTableFile << ": input 'translate-file' file not found."  << endl;
    }


  }

  // Validate the CHP Genotype experiment files.

  // Add to the CHP experimet file list from a file.
  if (rte.m_adtOpts.m_audit.empty()) {
    if (engine->getOpt("experiment-list-file") != "") {
      validOptions = validOptions && _validateInputExperimentFilesFromFileDMET3(rte, engine);
    } else  if (engine->getOpt("experiment-list") != "") {
      validOptions = validOptions && _validateInputExperimentFilesFromListDMET3(rte, engine);
    } else if (engine->getOpt("experiment-list-vector") != "") {
      validOptions = validOptions && _validateInputExperimentFilesFromVectorDMET3(rte, engine);
    } else if (engine->getOpt("geno-dmet2-type-file") != "") {
      // Using Genotype short report style input.
      if (!  Fs::fileExists(rte.m_adtOpts.m_inputGenoFile.c_str())) {
        validOptions = false;
        cerr << rte.m_adtOpts.m_inputGenoFile << ": input genotype short report experiment file not found."  << endl;
      }
    } else  {
      cerr << EXCEPTION_MISSING_EXPERIMENT_CHP_FILE << endl;
      rte.m_invalidOptionsMask |= ADT_INVALID_OPTION_MASK_MISSING_EXPERIMENT_CHP;
      validOptions = false;
    }
  }

  return validOptions;

}
// end _validateInputFilesDMET3
/*****************************************************************************/
/*****************************************************************************/
/**
 * _validateInputFilesDMET2:
 * Synopsis:
 *
 * Helper function to _validateInputFilesDMET2, broken out for clarity.
 *
 *
 * @param rte     - the single instance run time environment that houses the ADTOptions object as m_adtOpts. 
 * @params pgOopts - command line options object
 *
 * @return - true if all input files were valid.
 */
/*****************************************************************************/
static bool _validateInputFilesDMET2(RunTimeEnvironment & rte, TranslationEngine * engine)
{

  bool validOptions = true;

  // INPUT DIRECTORY
  if (! Fs::isReadableDir(rte.m_adtOpts.m_inputDir.c_str())) {
    validOptions = false;
    cerr << "Invalid option: directory \"" << rte.m_adtOpts.m_inputDir + "\" input directory not found.\n" ;
  } else if (! Fs::isReadableDir(rte.m_adtOpts.m_inputDir.c_str())) {
    validOptions = false;
    cerr << "Invalid option: directory " << rte.m_adtOpts.m_inputDir + "\" input directory read permission denied.\n" ;
  }

  pcrecpp::RE reCN("(?:(?:(?:^)|(?:[\\./]))\\d{8,8}_[\\d\\w]+_DMET[23]_cn_*)");
  pcrecpp::RE reGS("(?:(?:(?:^)|(?:[\\./]))\\d{8,8}_[\\d\\w]+_DMET[23]_Genotypes_Short_qc_*)");
  pcrecpp::RE reTT("((^)|([/]))DMET2_TTable_v\\d{8,8}_[^\\.]+\\.txt$");
  pcrecpp::RE reCurrentDirectory("\\.");


  std::list <std::string> searchFilesList
  = affymetrix_calvin_utilities::
    FileUtils::ListFiles(rte.m_adtOpts.m_inputDir.c_str(), "txt");

  unsigned int numFoundFiles = 0;
  std::string *testFiles[] = { &(rte.m_adtOpts.m_inputCopyFile),
                          &(rte.m_adtOpts.m_inputGenoFile),
                          &(rte.m_adtOpts.m_inputTTableFile)
                        };

  const char * testOptionNames[] = { "copy-dmet2-number-file",
                                     "geno-dmet2-type-file",
                                     "translate-file",
                                   };

  pcrecpp::RE *testRegexes[] = { &reCN, &reGS, &reTT };

  std::list<std::string>::iterator sflIt1 = searchFilesList.begin();

  for (; sflIt1 != searchFilesList.end(); ++sflIt1) {
    for (int i = 0; i < 3; i++) {
      if (testFiles[i]->empty()   &&
          (testRegexes[i]->PartialMatch(*sflIt1))) {
        *(testFiles[i]) = *sflIt1 ;
        numFoundFiles++;
      }
    }
  }

  if ((numFoundFiles == 0) &&
      (engine->getOpt(testOptionNames[0]) == "") &&
      (engine->getOpt(testOptionNames[1]) == "") &&
      (engine->getOpt(testOptionNames[2]) == "") ) {

    if (engine->getOpt("input-dir") == "") {
      APT_ERR_ABORT("\nCurrent directory contains no input files.");
    }
    cerr << "Invalid option: directory \"" << rte.m_adtOpts.m_inputDir << "\" does not contain any experiment files or the required DMET3 input files." << endl;
    return false;
  }

  const char *testFileTypes[] = { "input copy number",
                                  "geno type short report",
                                  "translation file",
                                };
  const std::string TRANSLATION_OPT_COPY_DMET2_FILE =                   \
      "(DMET2: 2 of 2 experiment files) => YYYYMMDD_*_DMET[23]_cn_*.txt";
  const std::string TRANSLATION_OPT_TRANSLATE_DMET2_FILE =      \
      "=> DMET2_TTTable_vYYYYMMDD_*.txt";

  const std::string * testFileNames[] = {
    &TRANSLATION_OPT_COPY_DMET2_FILE,
    &TRANSLATION_OPT_GENO_DMET2_FILE,
    &TRANSLATION_OPT_TRANSLATE_DMET2_FILE,
  };

  for (int i = 0; i < 3; i++) {
    if (testFiles[i]->empty()) {
      validOptions = false;
      if (reCurrentDirectory.FullMatch(rte.m_adtOpts.m_inputDir)) {
        cerr << "Invalid option: current directory";
      } else {
        cerr << "Invalid option: directory \"" << rte.m_adtOpts.m_inputDir << "\"";
      }
      cerr << " is missing " << testFileTypes[i] << " file ";
      cerr << *testFileNames[i] << "." << endl;
    } else if (!  Fs::fileExists(testFiles[i]->c_str())) {
      validOptions = false;
      cerr << *testFiles[i] << ": input " << testFileTypes[i] << " file not found."  << endl;
    } else if ((engine->getOpt(testOptionNames[i]) == "") && !testRegexes[i]->PartialMatch(*testFiles[i])) {
      validOptions = false;
      cerr << "Invalid option: ";
      cerr << *testFiles[i] << ": invalid input " << testFileTypes[i];
      cerr << " file name, " ;
      cerr << "expecting " << *testFileNames[i] ;
      cerr << "." << endl;
    }
  }

  rte.m_adtOpts.m_inputTTableType = TranslationTableModel::getTranslationTableFileType(rte.m_adtOpts.m_inputTTableFile);

  return validOptions;

}
// end _validateInputFilesDMET2
/*****************************************************************************/
/*****************************************************************************/
/**
 * _validateInputFiles:
 *
 * Synopsis: Helper to "validateOptions", broken out for clarity.
 * Validate user specified input files.
 *
 * This routine accomodates:
 * 1.) DMET2 files for DMET2 testing.
 * 2.) DMET3 files for DMET3 analysis
 * 3.) DMET2 files for DMET3 testing.
 * 4.) DMET2 files for DMET3 analysis*
 *
 * *The 4th case supports end users creating text files of data they
 * want analyzed. The text file format will be that of DMET2. Both
 * the Genotype Short Report and Copy Number Report files are required
 * for DMET2.
 *
 * The "-g" and "-c" options are exclusive to DMET2 input and it
 * is an error to use them in combination with "-e".
 *
 * The "-t" option for the translation table is required for both DMET2 and
 * DMET3, but the table formats are different and hence a DMET2 translation
 * table will need modified to perform DMET3 analysis. 
 *
 *
 * @param rte     - the single instance run time environment that houses the ADTOptions object as m_adtOpts. 
 * @params pgOopts - command line options object
 *
 * @return - true if all input files were valid.
 */
/*****************************************************************************/
static bool _validateInputFiles(RunTimeEnvironment & rte, TranslationEngine *engine)
{

  APT_ERR_ASSERT(engine, "");

  bool validOptions = true;


  // Assume the default call behavior is for DMET3. If none of the
  // DMET2 or DMET3 input file options are passed the
  // give missing DMET3 input file message.

  bool isDMET3Call = (!rte.m_adtOpts.m_dmet2Calling ||
                      (engine->getOpt("experiment-list") != "") ||
                      (engine->getOpt("experiment-list-file") != "") ||
                      (engine->getOpt("experiment-list-vector") != ""));

  bool isDMET2Call = (rte.m_adtOpts.m_dmet2Calling ||
                      (engine->getOpt("input-dir") != "") ||
                      (engine->getOpt("copy-dmet2-number-file") != ""));


  if (!isDMET3Call && !isDMET2Call) {
    cerr << "Invalid option: an experiment list of CHP files is required." << endl;
    validOptions = false;
  } else if (isDMET2Call && !isDMET3Call) {
    validOptions = _validateInputFilesDMET2(rte, engine);
    rte.m_adtOpts.m_streamType = ADT_EXPERIMENT_STREAM_TYPE_TSV;
  } else if (isDMET3Call && !isDMET2Call) {
    validOptions = _validateInputFilesDMET3(rte, engine);
  } else {
    cerr  << "Invalid option: use \"--dmet2-calling\" with DMET2 specific input files or the DMET2 input directory." << endl;
    validOptions = false;
  }

  if (engine->getOpt("geno-dmet2-type-file") != "") {
    rte.m_adtOpts.m_streamType = ADT_EXPERIMENT_STREAM_TYPE_TSV;
  } else {
    rte.m_adtOpts.m_streamType = ADT_EXPERIMENT_STREAM_TYPE_CHP;
  }

  // PROBE SET FILTER

  if (!rte.m_adtOpts.m_inputMarkerListFile.empty()) {

    if (!  Fs::fileExists(rte.m_adtOpts.m_inputMarkerListFile.c_str())) {
      validOptions = false;
      cerr << "Invalid option: " <<  rte.m_adtOpts.m_inputMarkerListFile << ": input marker list filter file not found."  << endl;
    }
  }


  // SAMPLE INFO FILE
  if (! rte.m_adtOpts.m_inputSampleFile.empty()) {
    if (!  Fs::fileExists(rte.m_adtOpts.m_inputSampleFile.c_str())) {
      validOptions = false;
      cerr << "Invalid option: " <<  rte.m_adtOpts.m_inputSampleFile << ": input sample information file not found."  << endl;
    }

  }

  // GENOTYPE OVERRIDE FILE
  if (!rte.m_adtOpts.m_inputGenotypeOverrideFile.empty()) {
    if (!  Fs::fileExists(rte.m_adtOpts.m_inputGenotypeOverrideFile.c_str())) {
      validOptions = false;
      cerr << "Invalid option: " <<  rte.m_adtOpts.m_inputGenotypeOverrideFile << ": input genotype override file not found."  << endl;
    }
  }

  return validOptions;

}
// end _validateInputFiles
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationCommonControl::validateCommonOptions
 * Synopsis:
 *
 * Validate options after they have been filled in by the call the
 * setCommonOptions. If any errors are detected during validation
 * then an exception is thrown via Err::Abort. 
 *
 * @param rte     - the single instance run time environment that houses the ADTOptions object as m_adtOpts previously set in setCommonOptions. . 
 * @param engine  - to be validated.
 *
 */
/*****************************************************************************/
void TranslationCommonControl::validateCommonOptions(RunTimeEnvironment & rte , TranslationEngine *engine)
{

  APT_ERR_ASSERT(engine, "");

  bool validOptions = true;

  // INPUT FILES
  validOptions = _validateInputFiles(rte, engine);

  // OUTPUT DIRECTORY
  if (Fs::isReadableDir(rte.m_adtOpts.m_outputDir.c_str())) {

    if (!Fs::isWriteableDir(rte.m_adtOpts.m_outputDir)) {
      validOptions = false;
      cerr  << "Invalid option: directory \"" << rte.m_adtOpts.m_outputDir << "\" output directory write permission denied." << endl;
    } else if (!Fs::isReadableDir(rte.m_adtOpts.m_outputDir)) {
      validOptions = false;
      cerr << "Invalid option: directoy \"" <<  rte.m_adtOpts.m_outputDir << "\" output directory read permission denied." << endl;
    }
  } else {

    // Fs::mkdir Throws fatal exception.
    Fs::mkdirPath(rte.m_adtOpts.m_outputDir, false);
    Verbose::out(1, rte.m_adtOpts.m_outputDir + ": output directory created.\n");
  }

  if (! validOptions) {

    APT_ERR_ABORT(rte.m_adtOpts.m_progName + ": invalid options detected.");

  }

}
// end TranslationCommonControl::validateCommonOptions
/*****************************************************************************/
