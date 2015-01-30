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
 * @file CNCorrelationEngine.cpp
 *
 * @brief This file contains the CNCorrelationEngine class members.
 */

//
#include "copynumber/CNCorrelationEngine.h"
//
#include "copynumber/CNCychp.h"
#include "copynumber/CNProbeSet.h"
//
#include "broadutil/BroadException.h"
#include "calvin_files/exception/src/ExceptionBase.h"
#include "chipstream/EngineUtil.h"
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "util/AffxMultiDimensionalArray.h"
#include "util/Fs.h"
#include "util/TmpFileFactory.h"
//
#include "../external/newmat/myexcept.h"

CNCorrelationEngine::Reg CNCorrelationEngine::reg;

CNCorrelationEngine * CNCorrelationEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == CNCorrelationEngine::EngineName())
        return (CNCorrelationEngine *)engine;
    return NULL;
}

/**
 * @brief Constructor
 */
CNCorrelationEngine::CNCorrelationEngine()
{
    m_pGenderCalls = NULL;
    m_pstdMethods = NULL;
    defineStdMethods();
    defineOptions();
}

/**
 * @brief Destructor
 */
CNCorrelationEngine::~CNCorrelationEngine()
{
    clear();
}

/**
 * @brief Initialize the data and free any memory allocated.
 */
void CNCorrelationEngine::clear()
{
    if (m_pstdMethods != NULL) {delete m_pstdMethods; m_pstdMethods = NULL;}
}

/**
 * @brief Define any standard methods used by the engine
 */
void CNCorrelationEngine::defineStdMethods() {
//    m_stdMethods["brlmm"] = "quant-norm.sketch=50000,pm-only,brlmm.transform=ccs.K=4";
}

/**
 * @brief Print any standard methods used by the engine
 * @param std::ostream& - The specified output stream
 */
void CNCorrelationEngine::printStandardMethods(std::ostream &out) {
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
void CNCorrelationEngine::defineOptions() {

  defineOptionSection("Input Options");
  defineOption("", "cychp-files", PgOpt::STRING_OPT,
                     "Text file specifying cychp files to process, one per line with the "
                     "first line being 'cychp_files'.",
                     "");
  defineOption("","start-file-number", PgOpt::INT_OPT,
                    "Start File Number",
                    "-1");
  defineOption("","end-file-number", PgOpt::INT_OPT,
                    "End File Number",
                    "-1");

  defineOptionSection("Output Options");
  defineOption("", "a5-output", PgOpt::STRING_OPT,
                    "A5 output correlation matrix file.", "");
  defineOption("", "text-output", PgOpt::STRING_OPT,
                    "Text output correlation matrix file.", "");

  defineOptionSection("Analysis Options");

  defineOptionSection("Misc Options");
  defineOption("", "explain", PgOpt::STRING_OPT,
                    "Explain a particular operation (i.e. --explain additional-Correlations-reference-method).",
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
                     "CYCHP files to process.",
                     "");
}

/**
 * @brief Define the states used by this engine
 */
void CNCorrelationEngine::defineStates() {
}


void CNCorrelationEngine::checkParameters()
{
}

void getCychpInputFiles(const std::string& str, std::vector<std::string> &cychpFiles, BaseEngine *engine) {
    if(engine->getOpt(str)!="") {
      affx::TsvFile tsv;
#ifdef WIN32
      tsv.m_optEscapeOk = false;
#endif
      std::string cychpFilesFile = engine->getOpt(str);
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
        Err::errAbort("Please specify the " + str + " parameter.");
    }
}

/**
 * @brief Make sure that our options are sane. Call Err::errAbort if not.
 */
void CNCorrelationEngine::checkOptionsImp() {

    defineStates();

    // file is not used by the engine
    //setLibFileOpt("config-file");

    if(getOpt("explain") != "") { explain(); exit(0); }

    checkParameters();
    Verbose::out(4, "Parameters are ok.");

    /* Read in cychp file list from other file if specified. */
    vector<string> cychpFiles;
    getCychpInputFiles("cychp-files", cychpFiles, this);
    if(cychpFiles.size() == 0)
        Err::errAbort("No CYCHP files specified.");
    setOpt("cychps",cychpFiles);
    cychpFiles.clear();

    for (unsigned int uiIndex = 0; (uiIndex < getOptVector("cychps").size()); uiIndex++)
    {
        if (Fs::fileExists(getOptVector("cychps")[uiIndex])==false) {
            Err::errAbort("A file specified as a CYCHP input does not exist:\t" + getOptVector("cychps")[uiIndex]);
        }
    }
}

/**
 * @brief Run the specified analysis
 */
void CNCorrelationEngine::runImp() {
  // make this all one big try catch to make sure destructors are
  // called if exceptions thrown.
  try { // outer try for handling the exception.
    try { // inner try for memory clean up.
      calculateCorrelation();
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

void CNCorrelationEngine::explain() {
}

void CNCorrelationEngine::extraHelp() {
    printStandardMethods(cout);
}

void CNCorrelationEngine::calculateCorrelation()
{
    if (getOpt("a5-output") == "") {Err::errAbort("a5-output must be specified, text-output is optional.");}
    unsigned int uiCychpCount1 = getOptVector("cychps").size();
    unsigned int uiCychpCount2 = (getOptInt("end-file-number") - getOptInt("start-file-number") + 1);

    if (getOptInt("start-file-number") < 1) {Err::errAbort("start-file-number must be greater than or equal to 1.");}
    if (getOptInt("end-file-number") < 1) {Err::errAbort("end-file-number must be greater than or equal to 1.");}
    if (getOptInt("start-file-number") > uiCychpCount1) {Err::errAbort("start-file-number must be less than or equal to " + ::getInt(uiCychpCount1) + ".");}
    if (getOptInt("end-file-number") > uiCychpCount1) {Err::errAbort("end-file-number must be less than or equal to " + ::getInt(uiCychpCount1) + ".");}
    if (getOptInt("end-file-number") < getOptInt("start-file-number")) {Err::errAbort("end-file-number must be greater than or equal to start-file-number.");}

    AffxMultiDimensionalArray<double> mx(uiCychpCount2, uiCychpCount1);
    Verbose::out(1, "*");
    CNCychp cychp;
    AffxArray<CNCychpProbeSetsCopyNumber> v;

    CNCychpProbeSetsCopyNumber* pv = new CNCychpProbeSetsCopyNumber[3000000];
    v.resize(3000000);
    for (int iIndex = 0; (iIndex < 3000000); iIndex++)
    {
        v[iIndex] = &pv[iIndex];
    }
    int iProbeSetCount = -1;
    int iX = getOptInt("xChromosome");
    int iY = getOptInt("yChromosome");

    // Load up the Log2Ratios in a5 format.
    Verbose::out(1, "*");
    Verbose::out(1, "Loading Log2Ratios...");
    affx::File5_File file5;
    affx::File5_Group* group5;
    affx::File5_Tsv* tsv5 = NULL;

    AffxString strFileName = getOpt("a5-output");
    try
    {
        file5.open(strFileName, affx::FILE5_RW);
        group5 = file5.openGroup("Cyto2", affx::FILE5_REPLACE);
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    try
    {
        for (unsigned int uiIndex = (getOptInt("start-file-number") - 1); (uiIndex < getOptInt("end-file-number")); uiIndex++)
        {
            AffxString strFileName2 = getOptVector("cychps").at(uiIndex);
            cychp.readFile(strFileName2, v, false);
            tsv5 = group5->openTsv("Log2Ratios-2-" + ::getInt(uiIndex), affx::FILE5_REPLACE);
            tsv5->defineColumn(0, 0, "Log2Ratios", affx::FILE5_DTYPE_DOUBLE);
            if (uiIndex == (getOptInt("start-file-number") - 1)) {iProbeSetCount = 0;}
            int iCount = 0;
            for (int iProbeSetIndex = 0; (iProbeSetIndex < v.size()); iProbeSetIndex++)
            {
                if ((v[iProbeSetIndex]->getChromosome() != 0) && (v[iProbeSetIndex]->getChromosome() != iX) && (v[iProbeSetIndex]->getChromosome() != iY))
                {
                    if (uiIndex == (getOptInt("start-file-number") - 1)) {iProbeSetCount++;}
                    iCount++;
                    tsv5->set_d(0, 0, v[iProbeSetIndex]->getLog2Ratio());
                    tsv5->writeLevel(0);
                }
            }
            tsv5->close();
            delete tsv5;
            if (iCount != iProbeSetCount) {Err::errAbort("Probe Set counts do not match.");}
        }
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}

    // Process the log2 Ratios.
    std::vector<double> vMeans1(uiCychpCount1);
    std::vector<double> vMeans2(uiCychpCount1);

    std::vector<double> vLog2Ratios1(iProbeSetCount);
    std::vector<double> vLog2Ratios2(iProbeSetCount);
    Verbose::out(1, "*");
    Verbose::out(1, "Calculating Correlations...");
    for (unsigned int uiIndex = 0; (uiIndex < uiCychpCount1); uiIndex++)
    {
        AffxString strFileName1 = getOptVector("cychps").at(uiIndex);
        cychp.readFile(strFileName1, v, false);
        int i = 0;
        for (int iProbeSetIndex = 0; (iProbeSetIndex < v.size()); iProbeSetIndex++)
        {
            if ((v[iProbeSetIndex]->getChromosome() != 0) && (v[iProbeSetIndex]->getChromosome() != iX) && (v[iProbeSetIndex]->getChromosome() != iY))
            {
                if (i >= iProbeSetCount) {Err::errAbort("Probe Set counts do not match.");}
                vLog2Ratios1[i] = v[iProbeSetIndex]->getLog2Ratio();
                i++;
            }
        }
        if (i != iProbeSetCount) {Err::errAbort("Probe Set counts do not match.");}
        for (unsigned int uiIndex2 = (getOptInt("start-file-number") - 1); (uiIndex2 <= (getOptInt("end-file-number") - 1)); uiIndex2++)
        {
            if (uiIndex == uiIndex2) {mx.set((uiIndex2 - (getOptInt("start-file-number") - 1)), uiIndex, 1.0); continue;}
            if (uiIndex > uiIndex2) {continue;}
            tsv5 = group5->openTsv("Log2Ratios-2-" + ::getInt(uiIndex2), affx::FILE5_OPEN);
            int i = 0;
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &vLog2Ratios2[i]);
                i++;
            }
            tsv5->close();
            delete tsv5;
            if (vMeans1[uiIndex] == 0)
            {
                double dSum = 0;
                for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
                {
                    dSum += vLog2Ratios1[iProbeSetIndex];
                }
                vMeans1[uiIndex] = (dSum / (double)iProbeSetCount);
            }
            if (vMeans2[uiIndex2] == 0)
            {
                double dSum = 0;
                for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
                {
                    dSum += vLog2Ratios2[iProbeSetIndex];
                }
                vMeans2[uiIndex2] = (dSum / (double)iProbeSetCount);
            }
            double dNumerator = 0;
            double d1 = 0;
            double d2 = 0;
            for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
            {
                dNumerator += (vLog2Ratios1[iProbeSetIndex] - vMeans1[uiIndex]) * (vLog2Ratios2[iProbeSetIndex] - vMeans2[uiIndex2]);
                d1 += (vLog2Ratios1[iProbeSetIndex] - vMeans1[uiIndex]) * (vLog2Ratios1[iProbeSetIndex] - vMeans1[uiIndex]);
                d2 += (vLog2Ratios2[iProbeSetIndex] - vMeans2[uiIndex2]) * (vLog2Ratios2[iProbeSetIndex] - vMeans2[uiIndex2]);
            }
            mx.set((uiIndex2 - (getOptInt("start-file-number") - 1)), uiIndex, (dNumerator / sqrt(d1) / sqrt(d2)));
        }
    }
    group5->close();
    delete group5;
    file5.close();

    delete[] pv;

  Fs::rm(strFileName,false);

    Verbose::out(1, "*");
    Verbose::out(1, "Writing out Correlations...");

    try
    {
        file5.open(strFileName, affx::FILE5_RW);
        group5 = file5.openGroup("Cyto2", affx::FILE5_REPLACE);
        tsv5 = group5->openTsv("Correlations", affx::FILE5_REPLACE);
        for (int iIndex = 0; (iIndex < uiCychpCount2); iIndex++)
        {
            tsv5->defineColumn(0, iIndex, "Correlation_" + ::getInt(iIndex + 1), affx::FILE5_DTYPE_DOUBLE);
        }
        for (int iIndex = 0; (iIndex < uiCychpCount1); iIndex++)
        {
            for (int iIndex2 = 0; (iIndex2 < uiCychpCount2); iIndex2++)
            {
                tsv5->set_d(0, iIndex2, mx.get(iIndex2, iIndex));
            }
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;
        group5->close();
        delete group5;
        file5.close();
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}

    strFileName = getOpt("text-output");
    if (!strFileName.empty())  {
      affx::TsvFile tsv;
      tsv.writeOpen(strFileName);
      AffxString str;
      for (int iIndex = 0; (iIndex < uiCychpCount1); iIndex++) {
        for (int iIndex2 = 0; (iIndex2 < uiCychpCount2); iIndex2++) {
          tsv.set(0,iIndex2, ::getDouble(mx.get(iIndex2, iIndex), 6));
        }
        tsv.writeLevel(0);
      }
      tsv.clear();

    }
}

