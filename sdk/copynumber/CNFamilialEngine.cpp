////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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
   @file   CNFamilialEngine.cpp

   @brief Analyzes trios and duos of cychp files.
*/

/*
  Brief design notes:

*/


#include "copynumber/CNFamilialEngine.h"
//
#include "copynumber/CNAnalysisMethodLOH.h"
#include "copynumber/CNFamilialAnalysisMethodFactory.h"
//
#include "broadutil/BroadException.h"
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/EngineUtil.h"
#include "file/TsvFile/TsvFile.h"

#include "util/Err.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/PgOptions.h"
#include "util/Verbose.h"
//
#include "../external/newmat/myexcept.h"
//
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//
#ifndef WIN32
#include <unistd.h>
#endif /* WIN32 */

CNFamilialEngine::Reg CNFamilialEngine::reg;

CNFamilialEngine * CNFamilialEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == CNFamilialEngine::EngineName())
        return (CNFamilialEngine *)engine;
    return NULL;
}

/**
 * @brief Constructor
 */
CNFamilialEngine::CNFamilialEngine() {
    m_pstdMethods = NULL;
    clear();
    defineStdMethods();
    defineOptions();
}

/**
 * @brief Destructor
 */
CNFamilialEngine::~CNFamilialEngine() {
    clear();
}

/**
 * @brief Initialize the data and free any memory allocated.
 */
void CNFamilialEngine::clear()
{
    if (m_pstdMethods != NULL) {delete m_pstdMethods; m_pstdMethods = NULL;}
    m_cychpIndex.clear();
    m_cychpMother.clear();
    m_cychpFather.clear();
    m_vCNFamilialAnalysisMethods.deleteAll();
    m_vCNFamilialReporters.deleteAll();
}

/**
 * @brief Define any standard methods for the engine
 */
void CNFamilialEngine::defineStdMethods() {
//    m_stdMethods["brlmm"] = "quant-norm.sketch=50000,pm-only,brlmm.transform=ccs.K=4";
}

/**
 * @brief Dump the standard methods to the output stream specified
 * @param std::ostream& - The specified output stream
 */
void CNFamilialEngine::printStandardMethods(std::ostream &out) {
    if (m_pstdMethods == NULL) {return;}
    std::map<std::string,std::string>::iterator iter;
    out << std::endl << "Standard Methods:" << std::endl;
      unsigned int maxLength = 0;
      for(iter = m_pstdMethods->begin(); iter != m_pstdMethods->end(); iter++) {
        if(iter->first.size() > maxLength)
          maxLength = (unsigned int)iter->first.size();
      }
      for(iter = m_pstdMethods->begin(); iter != m_pstdMethods->end(); iter++) {
        unsigned int currentLength = 0;
        out << " '" << iter->first << "' ";
        currentLength = (unsigned int)iter->first.size();
        while(currentLength < maxLength + 1) {
            out.put(' ');
          currentLength++;
        }
        out << iter->second << std::endl;
      }
}

/**
 * @brief Define the options used by this engine
 */
void CNFamilialEngine::defineOptions() {

  defineOptionSection("Input Options");
  defineOption("", "index-cychp-file", PgOpt::STRING_OPT,
                     "Cychp file specifying the sample to process as the index.",
                     "");
  defineOption("", "mother-cychp-file", PgOpt::STRING_OPT,
                     "Cychp file specifying the sample to process as the mother.",
                     "");
  defineOption("", "father-cychp-file", PgOpt::STRING_OPT,
                     "Cychp file specifying the sample to process as the father.",
                     "");
  defineOption("", "allele-frequency-file", PgOpt::STRING_OPT,
                     "Text file specifying allele frequencies.",
                     "");

  defineOptionSection("Output Options");
  defineOption("","familial-file", PgOpt::STRING_OPT,
                     "The file name to output the familial data (in calvin format) to.",
                     "");
  defineOption("","text-output", PgOpt::BOOL_OPT,
                     "Output data in ASCII text format in addition to calvin format.",
                     "true");

  defineOptionSection("Analysis Options");

  defOptMult("a", "analysis", PgOpt::STRING_OPT,
                   "String representing analysis pathway desired.",
                   "");

  defineOptionSection("Advanced Options"); 
  defineOption("", "xChromosome", PgOpt::INT_OPT, "X Chromosome.", "24");
  defineOption("", "yChromosome", PgOpt::INT_OPT, "Y Chromosome.", "25");

  defineOptionSection("Misc Options");
  defineOption("", "explain", PgOpt::STRING_OPT,
                    "Explain a particular operation (i.e. --explain brlmm or --explain brlmm-p).",
                    "");
}

/**
 * @brief Define the states used by this engine
 */
void CNFamilialEngine::defineStates() { }

/**
 * @brief Make sure that our options are sane. Call Err::errAbort if not.
 */
void CNFamilialEngine::checkOptionsImp() {

  defineStates();

  if(getOpt("explain") != "") { explain(); exit(0); }

  std::string outDir = getOpt("out-dir");
  std::vector<std::string> analysis = getOptVector("analysis");

  if (!Fs::fileExists(getOpt("index-cychp-file"))) {
    error("Must provide a valid index-cychp-file.");
  }

  if ((!Fs::fileExists(getOpt("father-cychp-file"))) && (!Fs::fileExists(getOpt("mother-cychp-file")))) {
    error("Must provide a valid father-cychp-file and/or a valid mother-cychp-file.");
  }

  if (!Fs::fileExists(getOpt("allele-frequency-file"))) {
    error("Must provide a valid allele-frequency-file.");
  }

  if (getOpt("familial-file") == "")
  {
    error("Must provide a familial-file.");
  }

  // setup output folder
  Fs::ensureWriteableDirPath(outDir, false);

  if (getOpt("familial-file") != "") {
    if ( !Fs::dirExists(Fs::dirname(getOpt("familial-file"))) ) {
      Fs::mkdirPath(Fs::dirname(getOpt("familial-file")), false);
    }
  }

  createAnalysis();
}

/**
 * @brief Loop over all familial analyses and store calvin group and dataset names needed to run them
 */
void CNFamilialEngine::prepareGroupsDatasets()
{
    using std::string;
    using std::set;

    // Loop over all analyses and store the (set-theoretic) union
    // of all cychp groups/datasets required by them. Only those
    // groups/datasets will be read in later.
    for (int iIndex = 0; (iIndex < (int)m_vCNFamilialAnalysisMethods.size()); iIndex++)
    {
        const groupDatasets_t groupDatasets = m_vCNFamilialAnalysisMethods[iIndex]->getGroupDatasets();
        for (groupDatasets_t::const_iterator gIter = groupDatasets.begin(); gIter != groupDatasets.end(); ++gIter)
        {
            for (set<string>::const_iterator sIter = gIter->second.begin(); sIter != gIter->second.end(); ++sIter)
            {
                m_groupDatasetsAllFamilial[gIter->first].insert(*sIter);
            }
        }
    }
}

/**
 * @brief Run the specifed analysis
 */
void CNFamilialEngine::runImp() {
  // make this all one big try catch to make sure destructors are
  // called if exceptions thrown.
  try { // outer try for handling the exception.
    try { // inner try for memory clean up.
        createAnalysis();
        prepareGroupsDatasets();
        loadCychpFiles();
        for (int iIndex = 0; (iIndex < (int)m_vCNFamilialAnalysisMethods.size()); iIndex++)
        {
            m_vCNFamilialAnalysisMethods[iIndex]->run();
            Verbose::out(1, "*");
        }
        for (int iIndex = 0; (iIndex < (int)m_vCNFamilialReporters.size()); iIndex++)
        {
            m_vCNFamilialReporters[iIndex]->setup(*this, m_cychpIndex, m_cychpMother, m_cychpFather, m_vCNFamilialAnalysisMethods);
            m_vCNFamilialReporters[iIndex]->run();
            Verbose::out(1, "*");
        }
    } // inner try end
    catch (...) {
      // Free up any allocated memory here...
        clear();
      // re-throw the unchanged exception.
      throw;
    }
  } // outer try end
  /* When things go wrong see if we can die gracefully here. */
  catch(Except &e) {
    Verbose::out(0,"");
    error("Exception caught. "
                  "Message is: " + ToStr(e.what()));
  }
  catch(const std::bad_alloc &e) {
    Verbose::out(0,"");
    error("std::bad_alloc caught. "
                  "The application has run out of memory. "
                  "Message is: " + ToStr(e.what()));
  }
  catch(affymetrix_calvin_exceptions::CalvinException &ce) {
    Verbose::out(0,"");
    error("affymetrix_calvin_exceptions::CalvinException caught. "
                   "Affymetrix GeneChip Command Console library has thrown an exception. "
                  "Message is: " + affymetrix_calvin_utilities::StringUtils::ConvertWCSToMBS(ce.Description()));
  }
  catch(BroadException &e) {
    Verbose::out(0,"");
    error("BroadException caught. "
                  "Message is: '" + ToStr(e.m_msg) + "' source file: '" + ToStr(e.m_sourcefile) +
                  ":" + ToStr(e.m_sourceline) + "'");
  }
  catch(const std::exception &e) {
    Verbose::out(0,"");
    error("std::exception caught. "
                  "Message is: " + ToStr(e.what()));
  }
  catch(const BaseException &e) {
    Verbose::out(0,"");
    error("newmat BaseException caught. "
                  "Message is: " + ToStr(e.what()));
  }
  catch(...) {
    Verbose::out(0,"");
    error("Unknown exception caught. "
                  "No message is available.");
  }
  // free any allocated memory here...
  clear();
}

/**
 * @brief Create the specified analysis
 */
void CNFamilialEngine::createAnalysis()
{
    Verbose::out(1, "CNFamilialEngine::createAnalysis()");
    affymetrix_calvin_parameter::ParameterNameValueType param;
    std::vector<std::string>& analysis = getOptVector("analysis");
    if (analysis.empty())
    {
        analysis.push_back("paternity");
        //analysis.push_back("segment-overlap");
        //analysis.push_back("genotype-concordance");
        //analysis.push_back("genotype-discordance");
        //analysis.push_back("cn-neutral-loh-concordance");
        //analysis.push_back("cn-loss-loh-concordance");
        //analysis.push_back("hetero-upd");
        //analysis.push_back("iso-upd");
        //analysis.push_back("denovo-cn");
        //analysis.push_back("hemizygous-parent-of-origin");
    }
    CNFamilialAnalysisMethod::getParams().clear();
    CNFamilialAnalysisMethodFactory factory;
    m_vCNFamilialAnalysisMethods.deleteAll();
    m_vCNFamilialReporters.deleteAll();
    for (unsigned int i = 0; i < analysis.size(); i++)
    {
        CNFamilialAnalysisMethod* pMethod = NULL;
        if (analysis[i] == "") {throw(Except("Invalid (blank) analysis specification"));}
        Verbose::out(1, analysis[i]);
        pMethod = factory.CNFamilialAnalysisMethodForString(analysis[i]);
        pMethod->setup(m_cychpIndex, m_cychpMother, m_cychpFather, getOpt("allele-frequency-file"));
        m_vCNFamilialAnalysisMethods.push_back(pMethod);
    }
    if (getOpt("familial-file") != "")
    {
        CNFamilialReporter* pReporter = NULL;
        pReporter = factory.CNFamilialReporterForString("familial-output");
        m_vCNFamilialReporters.push_back(pReporter);
    }
}

/**
 * @brief Load up the data from the cychp files into memory.
 */
void CNFamilialEngine::loadCychpFiles()
{
    Verbose::out(1, "*");
    // load index
    AffxString strIndexCychpFileName = getOpt("index-cychp-file");
    if (!m_cychpIndex.readFile(strIndexCychpFileName, m_groupDatasetsAllFamilial, "index", "index-cychp-file")) {error("Cannot read cychp file: " + strIndexCychpFileName);}
    int iMarkerCount = m_cychpIndex.getCychpProbeSetsCopyNumbers().getCount();
    if (iMarkerCount == 0) {error("No entries found in the index-cychp-file for the ProbeSets\\CopyNumber data set.");}

    // load mother
    AffxString strMotherCychpFileName = getOpt("mother-cychp-file");
    if (strMotherCychpFileName != "")
    {
        if (!m_cychpMother.readFile(strMotherCychpFileName, m_groupDatasetsAllFamilial, "mother", "mother-cychp-file")) {error("Cannot read cychp file: " + strMotherCychpFileName);}
        if (m_cychpMother.getCychpHeader().getGender() != affx::Female) {
            error("The mother-cychp-file does not contain female gender data");
        }
        if (m_cychpMother.getCychpProbeSetsCopyNumbers().getCount() != iMarkerCount)
        {
            error("The number of entries in the mother-cychp-file ProbeSets\\CopyNumber data set does not equal the number of entries in the index-cychp-file ProbeSets\\CopyNumber data set.");
        }
    }

    // load father
    AffxString strFatherCychpFileName = getOpt("father-cychp-file");
    if (strFatherCychpFileName != "")
    {
        if (!m_cychpFather.readFile(strFatherCychpFileName, m_groupDatasetsAllFamilial, "father", "father-cychp-file")) {error("Cannot read cychp file: " + strFatherCychpFileName);}
        if (m_cychpFather.getCychpHeader().getGender() != affx::Male) {
            error("The father-cychp-file does not contain male gender data");
        }
        if (m_cychpFather.getCychpProbeSetsCopyNumbers().getCount() != iMarkerCount)
        {
            error("The number of entries in the father-cychp-file ProbeSets\\CopyNumber data set does not equal the number of entries in the index-cychp-file ProbeSets\\CopyNumber data set.");
        }
    }
    checkFamilialConsistency();
    Verbose::out(1, "*");
}

// Make sure the cychp trio use the same CN reference and annotation file
void CNFamilialEngine::checkFamilialConsistency()
{
    if (getOpt("mother-cychp-file") != "") {
        checkFamilialConsistencyPair(m_cychpIndex, m_cychpMother);
    }
    if (getOpt("father-cychp-file") != "") {
        checkFamilialConsistencyPair(m_cychpIndex, m_cychpFather);
    }
}

void CNFamilialEngine::checkFamilialConsistencyPair(CNCychp& cychp1, CNCychp& cychp2)
{
    if (cychp1.getCychpHeader().getCNReferenceFileName() != cychp2.getCychpHeader().getCNReferenceFileName())
    {
        error("CN reference files mismatch in CYCHP files. "
              "File '" + ToStr(cychp1.getFileName()) + "' lists CN reference '" + ToStr(cychp1.getCychpHeader().getCNReferenceFileName()) + "', "
              "File '" + ToStr(cychp2.getFileName()) + "' lists CN reference '" + ToStr(cychp2.getCychpHeader().getCNReferenceFileName()) + "'");
    }
    if (cychp1.getCychpHeader().getAnnotationFileName() != cychp2.getCychpHeader().getAnnotationFileName())
    {
        error("Annotation files mismatch in CYCHP files. "
              "File '" + ToStr(cychp1.getFileName()) + "' lists annotation file '" + ToStr(cychp1.getCychpHeader().getAnnotationFileName()) + "', "
              "File '" + ToStr(cychp2.getFileName()) + "' lists annotation file '" + ToStr(cychp2.getCychpHeader().getAnnotationFileName()) + "'");
    }
    if (cychp1.getCychpHeader().getArrayType() != cychp2.getCychpHeader().getArrayType())
    {
        error("Array type mismatch in CYCHP files. "
              "File '" + ToStr(cychp1.getFileName()) + "' lists array type '" + ToStr(cychp1.getCychpHeader().getArrayType()) + "', "
              "File '" + ToStr(cychp2.getFileName()) + "' lists array type '" + ToStr(cychp2.getCychpHeader().getArrayType()) + "'");
    }
    if (cychp1.getCychpHeader().getDbsnpDate() != cychp2.getCychpHeader().getDbsnpDate())
    {
        error("dbsnp date mismatch in CYCHP files. "
              "File '" + ToStr(cychp1.getFileName()) + "' lists dbsnp date '" + ToStr(cychp1.getCychpHeader().getDbsnpDate()) + "', "
              "File '" + ToStr(cychp2.getFileName()) + "' lists dbsnp date '" + ToStr(cychp2.getCychpHeader().getDbsnpDate()) + "'");
    }
    if (cychp1.getCychpHeader().getDbsnpVersion() != cychp2.getCychpHeader().getDbsnpVersion())
    {
        error("dbsnp version mismatch in CYCHP files. "
              "File '" + ToStr(cychp1.getFileName()) + "' lists dbsnp date '" + ToStr(cychp1.getCychpHeader().getDbsnpVersion()) + "', "
              "File '" + ToStr(cychp2.getFileName()) + "' lists dbsnp date '" + ToStr(cychp2.getCychpHeader().getDbsnpVersion()) + "'");
    }
    if (cychp1.getCychpHeader().getNetaffxAnnotDate() != cychp2.getCychpHeader().getNetaffxAnnotDate())
    {
        error("Netaffx annotation date mismatch in CYCHP files. "
              "File '" + ToStr(cychp1.getFileName()) + "' lists dbsnp date '" + ToStr(cychp1.getCychpHeader().getNetaffxAnnotDate()) + "', "
              "File '" + ToStr(cychp2.getFileName()) + "' lists dbsnp date '" + ToStr(cychp2.getCychpHeader().getNetaffxAnnotDate()) + "'");
    }
    if (cychp1.getCychpHeader().getNetaffxBuild() != cychp2.getCychpHeader().getNetaffxBuild())
    {
        error("Netaffx build number mismatch in CYCHP files. "
              "File '" + ToStr(cychp1.getFileName()) + "' lists dbsnp date '" + ToStr(cychp1.getCychpHeader().getNetaffxBuild()) + "', "
              "File '" + ToStr(cychp2.getFileName()) + "' lists dbsnp date '" + ToStr(cychp2.getCychpHeader().getNetaffxBuild()) + "'");
    }
}


void CNFamilialEngine::explain() {
    CNFamilialAnalysisMethodFactory factory;
    std::vector<SelfDoc>docs = factory.getDocs();

    for(unsigned int i = 0; i < docs.size(); i++)
    {
        if(docs[i].getDocName() == getOpt("explain"))
        {
            SelfDoc::printExplanation(docs[i],cout);
        }
    }
    docs.clear();
}

void CNFamilialEngine::extraHelp() {
    CNFamilialAnalysisMethodFactory factory;
    printStandardMethods(cout);
    EngineUtil::printSelfDocs("Data transformations:", factory.getDocs());
}
