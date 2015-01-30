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
 * @file CNReferenceEngine.cpp
 *
 * @brief This file contains the CNReferenceEngine class members.
 */

//
#include "chipstream/EngineUtil.h"
#include "copynumber/CNAnalysisEngine.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNReferenceEngine.h"
#include "util/AffxBinaryFile.h"
#include "util/Fs.h"
//

CNReferenceEngine::Reg CNReferenceEngine::reg;

CNReferenceEngine * CNReferenceEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == CNReferenceEngine::EngineName())
        return (CNReferenceEngine *)engine;
    return NULL;
}

/**
 * @brief Constructor
 */
CNReferenceEngine::CNReferenceEngine()
{
    m_pEngine = this;
    m_iProbeSetCount = 0;
    m_iExperimentCount = 0;
    m_iProbeSetsToProcess = 0;
    m_iExperimentsToProcess = 0;
    defineOptions();
}

/**
 * @brief Destructor
 */
CNReferenceEngine::~CNReferenceEngine()
{
    clear();
}

void CNReferenceEngine::clear()
{
    m_data.clear();
}

/**
 * @brief Define the options used by the log2 ratio engine.
 */
void CNReferenceEngine::defineOptions()
{
  defineOptionSection("Input Options");
  defineOption("","expr-summary-file", PgOpt::STRING_OPT,
                     "Expression Summary Table file.",
                     "");
  defineOption("","genotype-calls-file", PgOpt::STRING_OPT,
                     "Genotype Calls Table file.",
                     "");
  defineOption("","genotype-confidences-file", PgOpt::STRING_OPT,
                     "Genotype Confidences Table file.",
                     "");
  defineOption("","genotype-report-file", PgOpt::STRING_OPT,
                     "Genotype Report Table file.",
                     "");
//  defOptMult("", "genotype-chp-file", PgOpt::STRING_OPT, "Genotype CHP file.", "");
//  defineOption("","chrRefId", PgOpt::INT_OPT,
//                     "Reference Chromosome",
//                     "2");
  defineOption("","reference-file", PgOpt::STRING_OPT, "Copy Number Reference file.", "");

  defineOptionSection("Output Options");
  defineOption("","reference-text-output", PgOpt::BOOL_OPT,
                     "Output the reference-file in acii text format.",
                     "false");

  defineOptionSection("Analysis Options");
  defineOption("", "log2-input", PgOpt::BOOL_OPT, "Input Allele Summaries are in log2.", "false");
  defineOption("", "adapter-type-normalization", PgOpt::BOOL_OPT,
                    "Adapter Type Normalization option. true = perform adapter type normalization.",
                    "true");

  defineOptionSection("Advanced Options");
  defineOption("", "set-analysis-name", PgOpt::STRING_OPT, "Analysis name to use as prefix for output files.", "");

  defineOptionSection("Execution Control Options");
  defineOption("", "mem-usage", PgOpt::INT_OPT,
                     "How many MB of memory to use for this run.",
                     "0");

  defineSharedOptions();
}

void CNReferenceEngine::defineStates() {
  defineSharedStates();
}

/**
* @brief Define options used by both this engine and the CNWorkflow engine
*/
void CNReferenceEngine::defineSharedOptions()
{
    m_pEngine->defineOptionSection("Additional CNReferenceEngine Options");
    m_pEngine->defineOption("", "probeset-ids", PgOpt::STRING_OPT,
                     "Tab delimited file with column 'probeset_id' specifying probesets to summarize.",
                     "");
    m_pEngine->defineOption("","annotation-file", PgOpt::STRING_OPT,
                    "NetAffx Annotation file.",
                    "");
    m_pEngine->defineOption("","xChromosome", PgOpt::INT_OPT,
                    "X Chromosome",
                    "24");
    m_pEngine->defineOption("","yChromosome", PgOpt::INT_OPT,
                    "Y Chromosome",
                    "25");
}

/**
* @brief Define state used by both this engine and the CNWorkflow engine
*/
void CNReferenceEngine::defineSharedStates()
{
    if(!m_pEngine->isOptDefined("sample-count")) {
        m_pEngine->defineOption("","sample-count", PgOpt::INT_OPT,
                        "sample-count",
                        "0");
        m_pEngine->defineOption("","create-reference", PgOpt::BOOL_OPT,
                        "create-reference",
                        "false");
    }
}

/**
 * @brief Check to see if the options used by the copy number engine are valid.
 */
void CNReferenceEngine::checkOptionsImp()
{
    defineStates();

    setLibFileOpt("reference-file");
    m_pEngine->setLibFileOpt("probeset-ids");
    m_pEngine->setLibFileOpt("annotation-file");

    if (m_pEngine == this)
    {
        if (getOpt("set-analysis-name") == "") {setOpt("set-analysis-name", "CN5");}
        if (m_pEngine->getOpt("expr-summary-file") == "") {throw(Except("Must specify a expr-summary-file."));}
        if (m_pEngine->getOpt("genotype-calls-file") == "") {throw(Except("Must specify a genotype-calls-file."));}
        if (m_pEngine->getOpt("genotype-confidences-file") == "") {throw(Except("Must specify a genotype-confidences-file."));}
        if (m_pEngine->getOpt("genotype-report-file") == "") {throw(Except("Must specify a genoptype-report-file."));}
        if (m_pEngine->getOpt("reference-file") == "") {throw(Except("Must specify a reference-file."));}
    }
    if (m_pEngine->getOpt("annotation-file") == "") {throw(Except("Must specify a netaffx annotation-file."));}

    AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
    if (m_pEngine == this)
    {
        m_pEngine->setOpt("create-reference", "true");
        m_pEngine->setOpt("reference-file", strReferenceFileName);
    }
}

/**
 * @brief Run the copy number ratio engine. The copy number ratio engine will inturn call the copy number engine.
 * @param CopyNumberOptions* - The options to run with.
 * @return int - 0 if successful else 1
 */
void CNReferenceEngine::runImp()
{
    try
    {
        m_data.setEngine(m_pEngine);
        uint64_t memAvail = ::getDouble(m_pEngine->getOpt("free-mem-at-start"));
        if (memAvail >= (2000 * MEGABYTE))
        {
            m_pEngine->setOpt("mem-usage", "2000");
        }
        else
        {
            int iMemUsage = (memAvail / MEGABYTE);
            m_pEngine->setOpt("mem-usage", ::getInt(iMemUsage));
        }

        m_data.loadAnnotation(true);
        m_iProbeSetCount = m_data.getProbeSets()->getCount();
        m_iExperimentCount = m_data.loadExperiments(m_pEngine->getOpt("expr-summary-file"));
        if (m_iExperimentCount == 0) {throw(Except("No experiments to process."));}
        if (m_iProbeSetCount == 0) {throw(Except("No probesets to process."));}
        if (!determineMemoryUsage()) {return;}
        if (!processData()) {return;}

        m_pEngine->setOpt("sample-count", ::getInt(m_iExperimentCount));
        if (!m_data.writeModelFile(m_pEngine->getOpt("reference-file"))) {return;}

        if (m_pEngine->isOptDefined("wave-correction-reference-method"))
        {
            if ((m_pEngine->getOpt("wave-correction-reference-method") != "none") && (m_pEngine->getOpt("wave-correction-reference-method") != ""))
            {
                CNAnalysisMethodFactory amFactory;
                CNAnalysisMethod* am = NULL;
                am = amFactory.CNAnalysisMethodForString(m_pEngine->getOpt("wave-correction-reference-method"));
                am->setEngine(m_pEngine);
                am->run();
                delete am;
            }
        }
        if ((m_pEngine->isOptDefined("temp-reference-file")) && (m_pEngine->getOpt("temp-reference-file") != ""))
        {
            Fs::rm(m_pEngine->getOpt("temp-reference-file"));
        }
        clear();

    }
    catch(...) {m_data.clear(); throw;}
}

/**
 * Guesstimate the maximum number of probe sets we can load for all experiments.
 * Guesstimate the maximum number of experiments we can load for all probe sets.
 * @return bool - true if successful
 */
bool CNReferenceEngine::determineMemoryUsage()
{
    int iMemUsage = m_pEngine->getOptInt("mem-usage");
    int iMBUsedForProbeSetList = (int)((77 * m_iProbeSetCount) / 1000000.0);
    int iMBUsedForExperimentList = (int)((58 * m_iExperimentCount) / 1000000.0);

    m_iProbeSetsToProcess = (int)((((iMemUsage - iMBUsedForProbeSetList - iMBUsedForExperimentList) * 1000000.0) / 9.0) / m_iExperimentCount);
    if (m_iProbeSetsToProcess <= 0) {throw(Except("Not enough memory to run application. ProbeSetsToProcess = " + ::getInt(m_iProbeSetsToProcess)));}
    if (m_iProbeSetsToProcess > m_iProbeSetCount) {m_iProbeSetsToProcess = m_iProbeSetCount;}

    m_iExperimentsToProcess = (int)((((iMemUsage - iMBUsedForProbeSetList - iMBUsedForExperimentList) * 1000000.0) / 9.0) / m_iProbeSetCount);
    if (m_iExperimentsToProcess <= 0) {throw(Except("Not enough memory to run application. ExperimentsToProcess = " + ::getInt(m_iExperimentsToProcess)));}
    if (m_iExperimentsToProcess > m_iExperimentCount) {m_iExperimentsToProcess = m_iExperimentCount;}

    return true;
}

/**
 * @brief Process the data.
 * @return bool - true if successful
 */
bool CNReferenceEngine::processData()
{
    if ((m_pEngine->isOptDefined("temp-reference-file")) && (m_pEngine->getOpt("temp-reference-file") != ""))
    {
        bool bFailed = true;
        int iMB = 80;
        while ((bFailed) && (iMB > 0))
        {
            try
            {
            Verbose::out(3, "AAlleleEstimates Allocation = " + ::getInt(m_iProbeSetCount * m_iExperimentsToProcess * 4));
            m_data.getAAlleleEstimates()->initialize(m_iExperimentsToProcess, m_iProbeSetCount);
            Verbose::out(3, "BAlleleEstimates Allocation = " + ::getInt(m_iProbeSetCount * m_iExperimentsToProcess * 4));
            m_data.getBAlleleEstimates()->initialize(m_iExperimentsToProcess, m_iProbeSetCount);
            Verbose::out(3, "Genotype Calls Allocation = " + ::getInt(m_iProbeSetCount * m_iExperimentsToProcess * 1));
            m_data.getGenotypeCalls()->initialize(m_iExperimentsToProcess, m_iProbeSetCount);
            Verbose::out(3, "Done with Allocations.");
            bFailed = false;
            } catch(...) {bFailed = true;}
            if (bFailed)
            {
                Verbose::out(3, "Memory allocations failed. Trying again...");
                iMB -= 10;
                m_iProbeSetsToProcess  =int(Min(double(m_iProbeSetsToProcess)  , ((iMB * 1000000.0) / m_iExperimentCount / 4.0)));
                m_iExperimentsToProcess=int(Min(double(m_iExperimentsToProcess), ((iMB * 1000000.0) / m_iProbeSetCount   / 4.0)));
            }
        }
        if (iMB <= 0) {throw(Except("Memory allocations failed."));}
        if ( !Fs::dirExists(Fs::dirname(m_pEngine->getOpt("temp-reference-file"))) ) {
          Fs::mkdirPath(Fs::dirname(m_pEngine->getOpt("temp-reference-file")), false);
        }
        Fs::rm(m_pEngine->getOpt("temp-reference-file"), false);
        affx::File5_File file5;
        file5.open(m_pEngine->getOpt("temp-reference-file"), affx::FILE5_CREATE);
        affx::File5_Group* pGroup5 = file5.openGroup("Intensities", affx::FILE5_RW);
        int iExperimentIndex = 0;
        while (iExperimentIndex < m_iExperimentCount)
        {
            if (!m_data.loadAlleleEstimatesByExperiment(m_pEngine->getOpt("expr-summary-file"), iExperimentIndex, m_iExperimentsToProcess)) {return false;}

            for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + m_iExperimentsToProcess), m_iExperimentCount)); iColIndex++)
            {
                affx::File5_Tsv* pTsv5 = pGroup5->openTsv("Signals-" + ToStr(iColIndex), affx::FILE5_REPLACE);
                pTsv5->defineColumn(0, 0, "A", affx::FILE5_DTYPE_DOUBLE);
                pTsv5->defineColumn(0, 1, "B", affx::FILE5_DTYPE_DOUBLE);
                for (int iIndex = 0; (iIndex < m_data.getProbeSets()->getCount()); iIndex++)
                {
                    pTsv5->set_d(0, 0, m_data.getAAlleleEstimates()->get(iColIndex - iExperimentIndex, iIndex));
                    pTsv5->set_d(0, 1, m_data.getBAlleleEstimates()->get(iColIndex - iExperimentIndex, iIndex));
                    pTsv5->writeLevel(0);
                }
                pTsv5->close();
                delete pTsv5;
            }

            iExperimentIndex += m_iExperimentsToProcess;
        }
        pGroup5->close();
        delete pGroup5;
        file5.close();

        // Free up some memory.
        m_data.getAAlleleEstimates()->initialize(1);
        m_data.getBAlleleEstimates()->initialize(1);
        m_data.getGenotypeCalls()->initialize(1);
    }

    int iProbeSetIndex = 0;

    // Read the gender calls from the report file.
    m_data.loadGenotypeReport(m_pEngine->getOpt("genotype-report-file"));

    // Calculate the median signals for each probe sets.
    bool bFailed = true;
    int iMB = 80;
    while ((bFailed) && (iMB > 0))
    {
        try
        {
        Verbose::out(3, "AAlleleEstimates Allocation = " + ::getInt(m_iProbeSetsToProcess * m_iExperimentCount * 4));
        m_data.getAAlleleEstimates()->initialize(m_iExperimentCount, m_iProbeSetsToProcess);
        Verbose::out(3, "BAlleleEstimates Allocation = " + ::getInt(m_iProbeSetsToProcess * m_iExperimentCount * 4));
        m_data.getBAlleleEstimates()->initialize(m_iExperimentCount, m_iProbeSetsToProcess);
        Verbose::out(3, "Genotype Calls Allocation = " + ::getInt(m_iProbeSetsToProcess * m_iExperimentCount * 1));
        m_data.getGenotypeCalls()->initialize(m_iExperimentCount, m_iProbeSetsToProcess);
        Verbose::out(3, "Done with Allocations.");
        bFailed = false;
        } catch(...) {bFailed = true;}
        if (bFailed)
        {
            Verbose::out(3, "Memory allocations failed. Trying again...");
            iMB -= 10;
            m_iProbeSetsToProcess  =int(Min(double(m_iProbeSetsToProcess)  , ((iMB * 1000000.0) / m_iExperimentCount / 4.0)));
            m_iExperimentsToProcess=int(Min(double(m_iExperimentsToProcess), ((iMB * 1000000.0) / m_iProbeSetCount   / 4.0)));
        }
    }
    if (iMB <= 0) {throw(Except("Memory allocations failed."));}
    iProbeSetIndex = 0;
    while (iProbeSetIndex < m_iProbeSetCount)
    {
        if (!m_data.loadAlleleEstimatesByProbeSet(m_pEngine->getOpt("expr-summary-file"), iProbeSetIndex, m_iProbeSetsToProcess)) {return false;}
        if (!m_data.loadGenotypeCallsByProbeSet(m_pEngine->getOpt("genotype-calls-file"), iProbeSetIndex, m_iProbeSetsToProcess)) {return false;}
        m_data.calculateMedianSignal(iProbeSetIndex, m_iProbeSetsToProcess);
        m_data.calculateXYMedianSignal(iProbeSetIndex, m_iProbeSetsToProcess);
        m_data.calculateMedianSignalsPerGenotype(iProbeSetIndex, m_iProbeSetsToProcess);
        iProbeSetIndex += m_iProbeSetsToProcess;
    }

    // Free up some memory.
    m_data.getAAlleleEstimates()->initialize(1);
    m_data.getBAlleleEstimates()->initialize(1);
    m_data.getGenotypeCalls()->initialize(1);
    return true;
}
