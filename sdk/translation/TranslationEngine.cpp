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
 * @file   TranslationEngine.cpp
 * @author Alan Williams
 * @date   Mon Jun 23 14:57:34 PDT 2008
 *
 * @brief Analysis engine for DMET 3.0 Translation Tool
 */
#include "translation/TranslationEngine.h"
//
#include "translation/ADTOptions.h"
#include "translation/CopyNumberTableModel.h"
#include "translation/ExperimentReport.h"
#include "translation/ExperimentResults.h"
#include "translation/GenotypeOverrideTableModel.h"
#include "translation/GenotypeTableModel.h"
#include "translation/MarkerExperimentReport.h"
#include "translation/MarkerListModel.h"
#include "translation/RegressionExperimentReport.h"
#include "translation/SampleInfoTableModel.h"
#include "translation/TranslationAudit.h"
#include "translation/TranslationCommonControl.h"
#include "translation/TranslationExperimentReport.h"
#include "translation/TranslationTable.h"
#include "translation/TranslationTableModel.h"
#include "translation/MetabolizerPhenotypingEngine.h"
//
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "util/Err.h"    
#include "util/Fs.h"    
#include "util/Verbose.h"
//
#include "pcrecpp.h"
//
#include <cstring>
#include <ctime>
#include <string>
#include <vector>
//

char *CONSOLE_ARGV[] = {(char*) "apt-dmet-translation", NULL };

using namespace affymetrix_calvin_exceptions;
using namespace std;
using namespace affymetrix;

TranslationEngine::Reg TranslationEngine::reg;

TranslationEngine * TranslationEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == TranslationEngine::EngineName())
		return (TranslationEngine *)engine;
	return NULL;
}

void TranslationEngine::checkOptionsImp() {
  setLibFileOpt("annotation-file");
  setLibFileOpt("marker-list-file");
  setLibFileOpt("translate-file");
  setLibFileOpt("metabolizer-file");
}

/*****************************************************************************/
/**
 * _defineOptions:
 * Synopsis:
 *
 * Helper function for the TranslationEngine::TranslationEngine constructor. 
 *
 *   Call TranslationCommonControl::defineCommonOptions and then define
 * the translate engine specific engines. It is required to define them
 * out side of the common controller so that command line help does not
 * display help for these options.
 *
 * @param adtOpts - the ATDOption object to validate
 * @param engine  - required for complex validation only.
 *
 * @return - Throws APT_ERR_ABORT exception on error.
 */
/*****************************************************************************/
static void _defineOptions(unsigned char controllerMask, TranslationEngine *engine)
{


  engine->setUsage("Allele translatione engine for the generation of allele translation reports.");

  TranslationCommonControl::defineCommonOptions(engine);

  // WARNING: The m_controllerMask if conditions are deliberately separate and not
  // coded as "if/else" pairs in order to afford testing both
  // code paths simultaneously. DO NOT make these "if/else".

  engine->defOptMult("", "experiment-list-vector", PgOpt::STRING_OPT,
                       "=> Not a command line option.", "");

  engine->defOptMult("", "marker-list-vector", PgOpt::STRING_OPT,
                       "=> Not a command line option.", "");

  engine->defOptMult("", "sample-table", PgOpt::INT_OPT,
                       "=> Not a command line option.", "0");

  engine->defineOptionSection("Internal Use Only");
  engine->defineOption("", "prototype-chp-files", PgOpt::BOOL_OPT,
                         "=> For internal use only.", "false" );
  engine->defineOption("", "salt", PgOpt::STRING_OPT, "Salt to add to the MD5 values.", "salt");
  
  return;

}
// end _defineOptions
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationEngine::TranslationEngine:
 * Synopsis:
 *
 *   Default constructor that primarily defines the options to use.
 *
 */
/*****************************************************************************/
TranslationEngine::TranslationEngine(char *argv[],
                                     unsigned char controllerMask,
                                     void (*messageHandlerCallback)(ADTOptions & adtOpts))
    : m_controllerMask(controllerMask)
{

  APT_ERR_ASSERT(m_controllerMask & (C_CONSOLE | C_CMDLINE), "");

  _defineOptions(m_controllerMask, this);

  parseArgv(argv == NULL ? CONSOLE_ARGV : argv);

  m_messageHandlerCallback = messageHandlerCallback;

  return;

}
// end TranslationEngine::TranslationEngine()
/*****************************************************************************/
/*****************************************************************************/
/**
 * _validateOptions:
 * Synopsis:
 *
 * Helper function for TranslationEngine::run. 
 *
 *   Call TranslateCommonControl::validateCommonOptions. The
 *   validation for translate engine specific options are handled in the
 * common controller. If this ever changes put the specific changes here.
 *
 * @param rte     - the single instance run time environment that houses the ADTOptions object as m_adtOpts.
 * @param engine  - required for complex validation only.
 *
 * @return - Throws APT_ERR_ABORT exception on error.
 */
/*****************************************************************************/
static void _validateOptions(RunTimeEnvironment & rte, TranslationEngine *engine)
{

  TranslationCommonControl::validateCommonOptions(rte, engine);

}
// end _validateOptions
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationEngine::_setADTOptions:
 * Synopsis:
 * Helper function for TranslationEngine::run.
 * 
 * Call TranslateCommonControl::setADTOptions and then set the translate
 * specific ATDOption values.
 *
 *  We need to set the specific options here so that when command help
 * is called the options are not visible.
 *
 *
 * @param adtOpts                 - the ATDOption object to set
 * @param messageHandlerCallback  - controller specific callback for setting controller independent options. 
 *
 */
/*****************************************************************************/
void TranslationEngine::_setADTOptions(ADTOptions & adtOpts,
                          void (*messageHandlerCallback)(ADTOptions *adtOpts))
{

  TranslationCommonControl::setCommonADTOptions(adtOpts, this);

    // WARNING: The if conditions are deliberately separate and not
  // coded as "if/else" pairs in order to afford testing both
  // code paths simultaneously. DO NOT make m_controllerMask "if/else".

  // CONSOLE SPECIFIC

  if (m_controllerMask & C_CONSOLE) {
    adtOpts.m_inputExperimentListVector = getOptVector("experiment-list-vector");

    if (getOptVector("marker-list-vector").size() > 1) {
      adtOpts.m_inputProbeSetVector = getOptVector("marker-list-vector");
      if (adtOpts.m_useFirstDupAlleleDef) {
        adtOpts.m_useFirstDupAlleleDefHeaderText =  HAPLOTYPE_REPORT_OPTION_FIRST;
      } else {
        adtOpts.m_useFirstDupAlleleDefHeaderText =  HAPLOTYPE_REPORT_OPTION_DUPLICATE;
      }
    }

    adtOpts.m_inputSampleTable   = getOptInt("sample-table");

    if (adtOpts.m_inputSampleTable > 0) {
      int numSampleRows = adtOpts.m_inputSampleTable;

      for (int row = 0; row < numSampleRows; row++) {
        adtOpts.m_sampleTable.push_back(getOptVector("sample-table::" + ToStr(row)));
      }
    }


  }

  if (m_controllerMask & C_CMDLINE) {
      adtOpts.m_prototypeCHPFiles = getOptBool("prototype-chp-files");
  }

  return;

}
// TranslationEngine::_setADTOptions
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationEngine::run
 * Synopsis:
 *
 *   Common option-less factory interface for running an engine.
 *
 *  Arguably this method should be broken into the following sections.
 *  Then again doing so would require passing 20-30 state variables so
 *  maybe not.
 *  1.) runInitialize
 *  2.) runDoAnalysis
 *  3.) runCleanup
 *
 */
/*****************************************************************************/
void TranslationEngine::runImp()
{

  //runInitialize
  
  // INPUT DATA: all the single instance classes are allocated here.
  // A single instance class may be instantiated from either a file
  // or a PgOptions std::vector. In either case the data is stashed in
  // a common TSV file model.
  MarkerListModel             *mlm  = NULL;
  CopyNumberTableModel        *cntm = NULL;
  GenotypeTableModel          *gtm  = NULL;
  GenotypeOverrideTableModel  *gotm = NULL;
  TranslationTableModel       *ttm  = NULL;
  SampleInfoTableModel        *sitm = NULL;
  TranslationTable            *tt   = NULL;

  std::stringstream                exceptionSStr;
  bool                        okToExit = true;
  bool						  okClose = true;
  bool						  reportsClosed = false;
  float                       previousCallSecond = 0.0;
  std::vector<ExperimentReport *>  report;
  std::string                      timeStamp;

  try {

    _setADTOptions(m_rte.m_adtOpts);

    if (m_messageHandlerCallback != NULL) {
      (*m_messageHandlerCallback)(m_rte.m_adtOpts);
    }

    // COMMAND LINE SPECIFIC
    if (m_controllerMask & C_CMDLINE) {

      timeStamp =  Util::getTimeStamp();
      m_rte.m_programTime.begin();

      Verbose::out(ADT_VERBOSE_NORMAL, "[" + m_rte.m_programName  + " " + timeStamp + "]  BEGIN ");
    }

    _validateOptions(m_rte, this);

    // AUDIT, don't translate. 
    if (! m_rte.m_adtOpts.m_audit.empty()) {
      TranslationAudit ta(m_rte);
      exit(ta.audit(m_rte));
    }

    m_rte.initializeRunTimeEnvironment();

    // Aggregate as many validation problems similar to how a compiler
    // attempts to report all problems. Once all validation error messages
    // have been written to standard error, thrown the exception. 
    if (! TranslationCommonControl::initializeInputFileModels(m_rte, &mlm, &cntm, &gtm, &ttm, &sitm, &gotm))  {
      APT_ERR_ABORT("Errors detected initializing input files.");
    }

    ADT_VERBOSE_ENUM headerLevel = ADT_VERBOSE_EXCEPTION;
    if ( m_controllerMask & C_CONSOLE ) {
      headerLevel = ADT_VERBOSE_NORMAL;
    }
    
    Verbose::out(headerLevel, "For research use only. Not for diagnostic purposes.");
    Verbose::out(headerLevel, "Program=" +  m_rte.m_adtOpts.m_progName);
    Verbose::out(headerLevel, "Version=" +  ADT_VERSION);
    Verbose::out(headerLevel, "Date=" + Util::getTimeStamp());
    Verbose::out(headerLevel, "TranslationFile=" + Fs::basename(m_rte.m_adtOpts.m_inputTTableFile));
    Verbose::out(headerLevel, "MarkerList=" +  Fs::basename(m_rte.m_adtOpts.m_inputMarkerListFile));
    Verbose::out(headerLevel, "HaplotypeReportOption=" +  m_rte.m_adtOpts.m_useFirstDupAlleleDefHeaderText);
    Verbose::out(headerLevel, "GenotypeOverrideFile=" + Fs::basename(m_rte.m_adtOpts.m_inputGenotypeOverrideFile));

    if (m_controllerMask & C_CMDLINE) {
      if (m_rte.m_adtOpts.m_inputExperimentFiles.size() > 0) {
        std::stringstream msgSStr;
        msgSStr << "Input CHP experiment files     count:      " << m_rte.m_adtOpts.m_inputExperimentFiles.size();
        Verbose::out(ADT_VERBOSE_NORMAL, msgSStr.str());
      }
    }

    // TRANSLATION TABLE: single instance business object that only stores data
    // derived from the the single instance data models, including
    // the translation table model. Specifically the translation table
    // business object stores business objects CallSets of CallElements
    // for the translation table in the form of GeneCall business objects. 
    tt = new TranslationTable(m_rte, *ttm , *cntm);

    // REPORT: factory, initialize reports in our report factory. 
    if (m_rte.m_adtOpts.m_regression > 0) {
      report.push_back(new RegressionExperimentReport());
    }
    if (m_rte.m_adtOpts.m_markerReport) {
      report.push_back(new MarkerExperimentReport());
    }

    // TRANSLATION REPORT: the three reports for DMET3 of comprehensive,
    // summary and uncalled are so similar they are implemented
    // as a generice report object simply called TranslationExperimentReport
	TranslationExperimentReport *ter = new TranslationExperimentReport(sitm, gotm);
    report.push_back(ter);

    // runDoAnalysis
    int totalMarkersRead = 0;
    int markersRead;
    int totalExperiments = m_rte.m_adtOpts.m_inputExperimentFiles.size();

    // How to report progress
    int progressStatus = 0; // >0 means we have a progress meter, <= 0 no progress meter

    if ((m_controllerMask & C_CMDLINE) && !(m_controllerMask & C_CONSOLE)) {
      if (Verbose::getParam().m_Verbosity == ADT_VERBOSE_NORMAL) {
          progressStatus = 1; //fine resolution progress meter
          Verbose::progressBegin(1, "Translating allele calls", 1 , 1, 1);
          Verbose::out(ADT_VERBOSE_NORMAL, "\n");
      } else {
          progressStatus = -1; //no progress meter but we still need to handle some internal state info
      }
    } else if (!(m_controllerMask & C_CMDLINE) && (m_controllerMask & C_CONSOLE)) {
      Verbose::progressBegin(ADT_VERBOSE_NORMAL, "Translating allele calls", totalExperiments , 1, totalExperiments);
      progressStatus = 3; // course progress meter
    }

    if ((m_controllerMask & C_CMDLINE) && !(m_controllerMask & C_CONSOLE)) {
      if (m_rte.m_adtOpts.m_profile) {
        m_rte.m_profiles["readNextExperiment"] = new Profile();
        m_rte.m_profiles["call_loop1"] = new Profile();
        m_rte.m_profiles["call_loop2"] = new Profile();
        m_rte.m_profiles["call_loop3"] = new Profile();
      }
      m_rte.m_profiles["call_loop"]    = new Profile();
      previousCallSecond = m_rte.m_profiles["call_loop"]->getElapsedSeconds();
    }


    // MAIN LOOP: stream each experiment into memory and output the
    // report results. 
    while ((markersRead = gtm->getNextExperiment(m_rte, *ttm, tt->m_geneExperimentCopyNumberCall)) > 0) {

      totalMarkersRead += markersRead;

      std::map< std::string, ExperimentResults*> er;

      if (m_controllerMask & C_CMDLINE) {
        Verbose::out(ADT_VERBOSE_NORMAL, gtm->m_experimentName, false);
      }

      // Do the dirty deed and get the allele translations.
      set< std::string >::iterator itS;

      for (itS = gtm->m_experimentGenes.begin();
           itS != gtm->m_experimentGenes.end(); itS++) {

        gtm->m_geneName = *itS;

        if (er.find(gtm->m_experimentName) == er.end()) {
          er[gtm->m_experimentName] = new ExperimentResults();
          er[gtm->m_experimentName]->m_experiment = gtm->m_experimentName;
        }

        if(progressStatus == 1 || progressStatus == -1) {
          m_rte.m_profiles["call_loop"]->begin();
          if((progressStatus == 1) && 
             ((m_rte.m_profiles["call_loop"]->getElapsedSeconds() - previousCallSecond) > 1)) {
            previousCallSecond = m_rte.m_profiles["call_loop"]->getElapsedSeconds();
            Verbose::progressStep(1);
          }
        }

        er[gtm->m_experimentName]->m_geneResults[gtm->m_geneName] = new ExperimentGeneResults(gtm->m_geneName, gtm->m_experimentName, tt->m_geneCopyNumberZeroCallSet[gtm->m_geneName].m_name);

        tt->translateExperimentGene(m_rte, *gtm, *ttm, *cntm, er[gtm->m_experimentName]->m_geneResults[gtm->m_geneName]);

        if (gtm->m_rows.size() > 0) {
          er[gtm->m_experimentName]->m_sample
          = gtm->m_rows[0][GT_SAMPLE_INDEX];
          er[gtm->m_experimentName]->m_experimentGuid = gtm->getCHPGuid();
        }

        er[gtm->m_experimentName]->m_markerCallCount += er[gtm->m_experimentName]->m_geneResults[gtm->m_geneName]->m_markerCallCount;

        if ((m_controllerMask & C_CMDLINE) && !(m_controllerMask & C_CONSOLE)) {
          m_rte.m_profiles["call_loop"]->end();
        }

      } // foreach experiment gene, translate

      for (size_t j = 0; j < report.size(); j++) {
        if (!(*report[j]).generate(m_rte, *ttm, er)) {
          Verbose::out(ADT_VERBOSE_NORMAL, "Report failed: " + (*report[j]).name());
        }
      }

      //cerr << gtm->m_experimentName <<  endl;
      if (progressStatus == 3) {
        Verbose::progressStep(1);
      }

      if (m_controllerMask & C_CMDLINE) {
        Verbose::out(ADT_VERBOSE_NORMAL, "...input markers read:         " + ToStr(totalMarkersRead) + "(" + ToStr(markersRead) + ")");
      }

      std::map< std::string, ExperimentResults*>::iterator it;

      it = er.begin();
      for (; it != er.end(); it++) {
        it->second->Clear();
      }

    } // foreach experiment, translate each gene.

	// runCleanup
    if (progressStatus > 0) {
      Verbose::progressEnd(1, "Done");
    }
	
	std::string alleleFile = ter->getComprehensiveFileName();
	for (int j = report.size() - 1; reportsClosed == false && j >= 0 ; j--) {
		okClose = (*report[j]).close(!okToExit) && okClose;
		delete report[j];
		report.erase(report.begin() + j);
	}
	reportsClosed = true;
    
    if (m_controllerMask & C_CMDLINE) {
      Verbose::out(ADT_VERBOSE_NORMAL, "Input experiments read:                    " + ToStr(totalExperiments));
      Verbose::out(ADT_VERBOSE_NORMAL, "Input markers read:                        " + ToStr(totalMarkersRead));
    }

	// Run the metabolizer genotyping
	std::string metabolizerFile = getOpt("metabolizer-file");
	if (metabolizerFile.empty() == false)
	{
		std::string outputMetabolizerFile = Fs::join(m_rte.m_adtOpts.m_outputDir, m_rte.m_adtOpts.m_outputReportPrefix + "_phenotype.rpt");
		MetabolizerPhenotypingEngine mg;
		mg.setOpt("metabolizer-file", metabolizerFile);
		mg.setOpt("allele-file", alleleFile);
		mg.setOpt("output-file", outputMetabolizerFile);
		mg.setOpt("salt", getOpt("salt"));
		mg.setOpt("md5-file", Fs::join(m_rte.m_adtOpts.m_outputDir, m_rte.m_adtOpts.m_outputReportPrefix + "_md5.rpt"));
		mg.setOpt("program-name", getOpt("program-name"));
		mg.setOpt("program-version", getOpt("program-version"));
		mg.checkOptions();
		mg.run();
		Verbose::out(ADT_VERBOSE_NORMAL, "Metabolizer file: " + outputMetabolizerFile);
	}

  } catch (const Except &  e) {
    okToExit = false;
    exceptionSStr << "Message is: " << e.what();

  } catch (CalvinException &ce) {
    okToExit = false;
    exceptionSStr << "Affymetrix Calvin library has thrown an exception. ";
    exceptionSStr << "Description: '" + StringUtils::ConvertWCSToMBS(ce.Description()) + "'";
  } catch (const std::exception &e) {
    okToExit = false;
    exceptionSStr << "Exception caught. " "Most problems with this program are memory related. ";
    exceptionSStr << "Message is: " << e.what();
  } catch (...) {
    okToExit = false;
    exceptionSStr << "Unknown exception caught (most problems with this program are memory related).";
  }

  // !!! WARNING !!!
  // Beyond here no calls to Verbose::out or Verbose::warn can be
  // made if okToExit is false and C_CMDLINE is false.
  // This is because Vebose::out will
  // throw an exception in Windows. All messaging beyond this
  // point must be conditioned upon okToExit being true
  // if C_CMDLINE is false;

  for (int j = report.size() - 1; reportsClosed == false && j >= 0 ; j--) {
    okClose = (*report[j]).close(!okToExit) && okClose;
    delete report[j];
    report.erase(report.begin() + j);
  }

  if (!okClose) {
    okToExit = false;
    exceptionSStr << "Errors occurred in creating final report files.";
  }

  if (m_controllerMask & C_CMDLINE) {
    if (okToExit) {
      m_rte.m_programTime.end();

      if (!(m_controllerMask & C_CONSOLE)) {
        if (m_rte.m_adtOpts.m_profile) {
          m_rte.profilesReport(m_rte.m_programTime.getElapsedSeconds());
        }
        Verbose::out(ADT_VERBOSE_NORMAL, "Elasped time: " +
                     m_rte.m_programTime.getElapsedFormatedString() + ".");
      }
    }
    timeStamp = Util::getTimeStamp();
    Verbose::out(ADT_VERBOSE_NORMAL, "[" + m_rte.m_programName  + " " + timeStamp + "]  END ");

  }


  if (mlm     != NULL)  delete mlm;
  if (cntm    != NULL)  delete cntm;
  if (gtm     != NULL)  delete gtm;
  if (gotm    != NULL)  delete gotm;
  if (ttm     != NULL)  delete ttm;
  if (sitm    != NULL)  delete sitm;
  if (tt      != NULL)  delete tt;

  // Rethrow the exception
  if (!okToExit) {
    if ((m_controllerMask & C_CMDLINE)) {
      if (m_controllerMask & C_CONSOLE) {
        APT_ERR_ABORT("");
      }
      exit(1);
    }
    APT_ERR_ABORT(exceptionSStr.str());
  }


  return;
}
// end TranslationEngine::run
/*****************************************************************************/
