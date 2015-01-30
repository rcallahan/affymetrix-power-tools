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
 * @file CNWaveEngine.cpp
 *
 * @brief This file contains the CNWaveEngine class members.
 */

//
#include "copynumber/CNWaveEngine.h"
//
#include "copynumber/CNAnalysisMethodFactory.h"
//
#include "broadutil/BroadException.h"
#include "calvin_files/exception/src/ExceptionBase.h"
#include "chipstream/EngineUtil.h"
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "file5/File5_Tsv.h"
#include "util/Fs.h"
#include "util/TmpFileFactory.h"
//
#include "../external/newmat/myexcept.h"
//
CNWaveEngine::Reg CNWaveEngine::reg;

CNWaveEngine * CNWaveEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == CNWaveEngine::EngineName())
        return (CNWaveEngine *)engine;
    return NULL;
}

/**
 * @brief Constructor
 */
CNWaveEngine::CNWaveEngine()
{
    m_pGenderCalls = NULL;
    m_pstdMethods = NULL;
    defineStdMethods();
    defineOptions();
}

/**
 * @brief Destructor
 */
CNWaveEngine::~CNWaveEngine()
{
    clear();
}

/**
 * @brief Initialize the data and free any memory allocated.
 */
void CNWaveEngine::clear()
{
    if (m_pstdMethods != NULL) {delete m_pstdMethods; m_pstdMethods = NULL;}
    // Clean up static data.
    CNAnalysisMethod::getCelFileParams()->clear();
    CNAnalysisMethod::getParams()->clear();
}

/**
 * @brief Define any standard methods used by the engine
 */
void CNWaveEngine::defineStdMethods() {
//    m_stdMethods["brlmm"] = "quant-norm.sketch=50000,pm-only,brlmm.transform=ccs.K=4";
}

/**
 * @brief Print any standard methods used by the engine
 * @param std::ostream& - The specified output stream
 */
void CNWaveEngine::printStandardMethods(std::ostream &out) {
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
void CNWaveEngine::defineOptions() {

  defineOptionSection("Input Options");
  defineOption("", "cn-reference-input", PgOpt::STRING_OPT,
                    "The CN Reference input to the process.",
                    "");
  defineOption("", "config-file", PgOpt::STRING_OPT,
                    "The configuration file name as passed from GTC or the Cyto Browser.",
                    "");
  defineOption("", "cychp-files", PgOpt::STRING_OPT,
                     "Text file specifying cychp files to process, one per line with the "
                     "first line being 'cychp_files'.",
                     "");

  defineOptionSection("Output Options");
  defineOption("", "cn-reference-output", PgOpt::STRING_OPT,
                    "Output CN reference file name.", "");

  defineOptionSection("Analysis Options");
  defineOption("a", "analysis", PgOpt::STRING_OPT,
                    "String representing analysis pathway desired.",
                    "additional-waves-reference-method");

  defineOptionSection("Misc Options");
  defineOption("", "explain", PgOpt::STRING_OPT,
                    "Explain a particular operation (i.e. --explain additional-waves-reference-method).",
                    "");

  defineOptionSection("Advanced Options");
  defineOption("","xChromosome", PgOpt::INT_OPT,
                    "X Chromosome",
                    "24");
  defineOption("","yChromosome", PgOpt::INT_OPT,
                    "Y Chromosome",
                    "25");

  defineOptionSection("Engine Options (Not used on command line)");
  defOptMult("", "cychps", PgOpt::STRING_OPT,
                     "CYCHP files to process (cannot be used with the cychp-files option).",
                     "");
}

/**
 * @brief Add a cel file specification to the options
 * @param const std::string& - The cel file name
 * @param const std::string& - The sample (arr) file name
 * @param const std::string& - The result (cychp) file name
 */
void CNWaveEngine::addCychp(const std::string& strCychpFileName)
{
    pushOpt("cychps", strCychpFileName);
}

/**
 * @brief Define the states used by this engine
 */
void CNWaveEngine::defineStates() {
  defineOption("", "temp-reference-file", PgOpt::STRING_OPT,
                    "A temporary file use to store intermediate data.",
                    "");
}


void CNWaveEngine::checkParameters()
{
}

void getCychpFiles(std::vector<std::string> &cychpFiles, BaseEngine *engine) {
    if(engine->getOpt("cychp-files")!="") {
      affx::TsvFile tsv;
#ifdef WIN32
      tsv.m_optEscapeOk = false;
#endif
      std::string cychpFilesFile = engine->getOpt("cychp-files");
      string file;
      tsv.bind(0, "cychp_files", &file, affx::TSV_BIND_REQUIRED);
      if(tsv.open(cychpFilesFile) != affx::TSV_OK) {
        Err::errAbort("Couldn't open cychp-files file: " + cychpFilesFile);
      }
      tsv.rewind();
      while(tsv.nextLevel(0) == affx::TSV_OK) {
        cychpFiles.push_back(file);
      }
      tsv.close();
      Verbose::out(1, "Read " + ToStr(cychpFiles.size()) + " cychp files from: " + Fs::basename(cychpFilesFile));
    } else {
      cychpFiles = engine->getOptVector("cychps");
      if(cychpFiles.size() == 0) {
        for(vector<const char *>::size_type i = 0; i < engine->getArgCount(); i++)
            cychpFiles.push_back(engine->getArg(i));
      }
    }

#ifdef WIN32
  // Windows doesn't do wildcard expansion in windows, let the user know.
  if(cychpFiles.size() > 0 && cychpFiles[0].find("\\*") != string::npos) {
    Err::errAbort("Wildcard ('*') expansion not available in windows, please use --cychp-files option.");
  }
#endif

}

/**
 * @brief Make sure that our options are sane. Call Err::errAbort if not.
 */
void CNWaveEngine::checkOptionsImp() {

    defineStates();

    // file is not used by the engine
    //setLibFileOpt("config-file");

    if(getOpt("explain") != "") { explain(); exit(0); }

    checkParameters();
    Verbose::out(4, "Parameters are ok.");

    if (getOpt("cn-reference-output") == "") {error("Must specify an cn-reference-output.");}
    if (getOpt("temp-dir") == "") { setOpt("temp-dir", Fs::join(getOpt("out-dir"),"temp")); }

    if (getOpt("cn-reference-input") != "")
    {
        if ((!isCyto2Reference(getOpt("cn-reference-input"))) && (!isCopyNumberReference(getOpt("cn-reference-input")))) {error("The cn-reference-input specified either does not exist, or is not the correct type for the CNWaveEngine.");}
    }

    if (getOpt("cn-reference-input") == getOpt("cn-reference-output")) {error("The cn-reference-output must not be the same as the cn-reference-input");}

    /* Read in cychp file list from other file if specified. */
    vector<string> cychpFiles;
    getCychpFiles(cychpFiles, this);
    if(cychpFiles.size() == 0)
        Err::errAbort("No CYCHP files specified.");
    setOpt("cychps",cychpFiles);

    if (Fs::isReadable(getOpt("cn-reference-input")))  {
      if (affx::File5_File::isHdf5file(getOpt("cn-reference-input"))) {Verbose::out(1, "Input reference is in HDF5 format.");}
      else {error("cn-reference-input is not in HDF5 file format: " + getOpt("cn-reference-input"));}
    }
    else {error("Cannot open cn-reference-input: " + getOpt("cn-reference-input"));}
    for (unsigned int uiIndex = 0; (uiIndex < getOptVector("cychps").size()); uiIndex++)
    {
        if (!Fs::fileExists(getOptVector("cychps")[uiIndex]))
        {
            Err::errAbort("A file specified as a CYCHP input does not exist:\t" + getOptVector("cychps")[uiIndex]);
        }
    }
}

/**
 * @brief Run the specified analysis
 */
void CNWaveEngine::runImp() {
  // make this all one big try catch to make sure destructors are
  // called if exceptions thrown.
  try { // outer try for handling the exception.
    try { // inner try for memory clean up.
      Fs::ensureWriteableDirPath(getOpt("out-dir"), false);
      Fs::ensureWriteableDirPath(getOpt("temp-dir"), false);

        modifyReference();

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

void CNWaveEngine::modifyReference()
{
    Verbose::out(1, "*");
    CNAnalysisMethodFactory amFactory(true);
    CNAnalysisMethod* am = NULL;
    am = amFactory.CNAnalysisMethodForString(getOpt("analysis"));
    am->setEngine(this);
    am->run();
    delete am;
    Verbose::out(1, "*");
}


// This method is called by ChAS without having called checkOptions.
// In short, best not to refer to Engine options/state in this method
AffxString CNWaveEngine::getAnnotationParameter(const AffxString& strFileName, const AffxString& strParameterName)
{
    AffxString strParameterValue;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        AffxString str;
        try
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            if (group5 != NULL)
            {
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if (tsv5 != NULL)
            {
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &str);
                int iIndex = str.indexOf("=");
                if (str.startsWith("#%" + strParameterName + "=")) {strParameterValue = str.substring(iIndex + 1);}
            }
            tsv5->close();
            delete tsv5;
            }
            group5->close();
            delete group5;
            }
            file5.close();
        } catch(...) {return "";}
    }
    return strParameterValue;
}

// This method is called by ChAS without having called checkOptions.
// In short, best not to refer to Engine options/state in this method
bool CNWaveEngine::isCyto2Reference(const AffxString& strReferenceFileName)
{
    int i = 0;
    try {i = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "Cyto2", "MedianSignals");}
    catch (...) {i = 0;}
    return (i > 0);
}

bool CNWaveEngine::isCopyNumberReference(const AffxString& strReferenceFileName)
{
    int i = 0;
    try {i = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "CopyNumber", "Reference");}
    catch (...) {i = 0;}
    return (i > 0);
}

void CNWaveEngine::explain() {
    CNAnalysisMethodFactory factory (true);
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

void CNWaveEngine::extraHelp() {
    CNAnalysisMethodFactory factory(true);
    printStandardMethods(cout);
    EngineUtil::printSelfDocs("Data transformations:", factory.getDocs());
}
