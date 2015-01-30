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

/**
 * @file CNLog2RatioEngine.cpp
 *
 * @brief This file contains the CNLog2RatioEngine class members.
 */

#include "copynumber/CNLog2RatioEngine.h"
//
#include "chipstream/EngineUtil.h"
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "util/AffxBinaryFile.h"
#include "util/Fs.h"
#include "util/Guid.h"
//
CNLog2RatioEngine::Reg CNLog2RatioEngine::reg;

CNLog2RatioEngine * CNLog2RatioEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == CNLog2RatioEngine::EngineName())
        return (CNLog2RatioEngine *)engine;
    return NULL;
}

/**
 * @brief Constructor
 */
CNLog2RatioEngine::CNLog2RatioEngine()
{
    m_pEngine = this;
    m_SharedStateIsDefined = false;
    m_iProbeSetCount = 0;
    m_iExperimentCount = 0;
    m_iProbeSetsToProcess = 0;
    m_iExperimentsToProcess = 0;
    defineOptions();
}

/**
 * @brief Destructor
 */
CNLog2RatioEngine::~CNLog2RatioEngine()
{
    clear();
}

/**
 * @brief Free memory allocated for this object
 */
void CNLog2RatioEngine::clear()
{
    m_data.clear();
    m_objCNAnalysisEngine.clear();
}

/**
 * @brief Define the options used by the log2 ratio engine.
 */
void CNLog2RatioEngine::defineOptions()
{
  defineOptionSection("Input Options");
  defineOption("", "probeset-ids", PgOpt::STRING_OPT,
                     "Tab delimited file with column 'probeset_id' specifying probesets to summarize.",
                     "");
  // while this isn't defined as a "shared" option, we assume it could be
  defineOption("","annotation-file", PgOpt::STRING_OPT,
                     "NetAffx Annotation file.",
                     "");
  // while this isn't defined as a "shared" option, we assume it could be
  defineOption("","expr-summary-file", PgOpt::STRING_OPT,
                     "Expression Summary Table file.",
                     "");
  // while this isn't defined as a "shared" option, we assume it could be
  defineOption("","genotype-calls-file", PgOpt::STRING_OPT,
                     "Genotype Calls Table file.",
                     "");
  // while this isn't defined as a "shared" option, we assume it could be
  defineOption("","genotype-confidences-file", PgOpt::STRING_OPT,
                     "Genotype Confidences Table file.",
                     "");
  // while this isn't defined as a "shared" option, we assume it could be
  defineOption("","genotype-report-file", PgOpt::STRING_OPT,
                     "Genotype Report Table file.",
                     "");
  // while this isn't defined as a "shared" option, we assume it could be
  defineOption("", "reference-file", PgOpt::STRING_OPT, "Copy Number Reference file.", "");

  defineOptionSection("Output Options");
  defineOption("", "call-copynumber-engine", PgOpt::BOOL_OPT, "Call the copy number engine to output CNCHP files.", "true");
  defineOption("", "log2ratio-hdf5-output", PgOpt::BOOL_OPT, "Write intermediate output in HDF5 format.", "false");
  defineOption("", "log2ratio-text-output", PgOpt::BOOL_OPT, "Write intermediate output in ASCII Text format.", "false");

  defineOptionSection("Analysis Options");

  defineOptionSection("Advanced Options");
  defineOption("","xChromosome", PgOpt::INT_OPT, "X Chromosome", "24");
  defineOption("","yChromosome", PgOpt::INT_OPT, "Y Chromosome", "25");
  defineOption("", "gc-correction-bin-count", PgOpt::INT_OPT, "The number of bins to use for GC content.", "25");


  defineOptionSection("Execution Control Options");
  defineOption("", "mem-usage", PgOpt::INT_OPT,
                    "How many MB of memory to use for this run.",
                    "0");

  defineSharedOptions();
}

void CNLog2RatioEngine::defineStates() {
    defineOption("", "create-reference", PgOpt::BOOL_OPT,
                    "create-reference",
                    "false");
  defineSharedStates();
}

/**
 * @brief Define options used by both the CNLog2RatioEngine and the CNWorkflowEngine
 */
void CNLog2RatioEngine::defineSharedOptions()
{
    m_pEngine->defineOptionSection("Additional CNLog2RatioEngine Options");
    m_pEngine->defOptMult("a", "analysis", PgOpt::STRING_OPT,
                   "String representing analysis pathway desired.",
                   "");
    m_pEngine->defineOption("", "delete-files", PgOpt::BOOL_OPT,
                    "Delete extra output files after the run has completed.",
                    "false");
    m_pEngine->defineOption("", "log2-input", PgOpt::BOOL_OPT, "Input Allele Summaries are in log2.", "false");

    m_pEngine->defineOption("", "gc-content-override-file", PgOpt::STRING_OPT, "Input file used to override the GC content read from the annotation file (Two columns with header line, ProbeSetName/GCContent).", "");

    m_objCNAnalysisEngine.setEngine(m_pEngine);
    m_objCNAnalysisEngine.defineSharedOptions();
}

void CNLog2RatioEngine::defineSharedStates()
{
    if(!m_SharedStateIsDefined) {
        m_objCNAnalysisEngine.setEngine(m_pEngine);
        m_objCNAnalysisEngine.defineSharedStates();
        m_SharedStateIsDefined = true;
    }
}


/**
 * @brief Check to see if the options used by the copy number engine are valid.
 */
void CNLog2RatioEngine::checkOptionsImp()
{
    defineStates();

    setLibFileOpt("probeset-ids");
    setLibFileOpt("annotation-file");
    setLibFileOpt("reference-file");
    m_pEngine->setLibFileOpt("gc-content-override-file");

    std::vector<std::string>& analysis = m_pEngine->getOptVector("analysis");

    if(analysis.size() == 0) {
        analysis.push_back("log2-ratio");
        analysis.push_back("allelic-difference");
    }
    if (m_pEngine == this)
    {
        if (getOpt("set-analysis-name") == "") {setOpt("set-analysis-name", "CN5");}
        if (getOpt("expr-summary-file") == "") {throw(Except("Must specify a expr-summary-file."));}
        if (getOpt("genotype-calls-file") == "") {throw(Except("Must specify a genotype-calls-file."));}
        if (getOpt("genotype-confidences-file") == "") {throw(Except("Must specify a genotype-confidences-file."));}
        if (getOpt("genotype-report-file") == "") {throw(Except("Must specify a genoptype-report-file."));}
    }
    if (m_pEngine->getOpt("annotation-file") == "") {throw(Except("Must specify a netaffx annotation-file."));}
    if (m_pEngine->getOpt("reference-file") == "") {throw(Except("Must specify a reference-file."));}

    AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");

    if (Fs::fileExists(strReferenceFileName))
    {
        if (affx::File5_File::isHdf5file(strReferenceFileName)) {Verbose::out(1, "Input reference is in HDF5 format.");}
        else {Verbose::out(1, "Input reference is in text format.");}
        if ((affx::File5_File::isHdf5file(strReferenceFileName)) && (getOptBool("log2ratio-text-output")))
        {
            Verbose::out(1, "WARNING: Input reference-file is in HDF5 format, overriding log2ratio-text-output to log2ratio-hdf5-output.");
            setOpt("log2ratio-hdf5-output", "true");
            setOpt("log2ratio-text-output", "false");
        }
        if ((!affx::File5_File::isHdf5file(strReferenceFileName)) && (getOptBool("log2ratio-hdf5-output")))
        {
            Verbose::out(1, "WARNING: Input reference-file is in text format, overriding log2ratio-hdf5-output to log2ratio-text-output.");
            setOpt("log2ratio-hdf5-output", "false");
            setOpt("log2ratio-text-output", "true");
        }
    }

    m_objCNAnalysisEngine.checkOptions();
}

/**
 * @brief Run the copy number ratio engine. The copy number ratio engine will inturn call the copy number engine.
 */
void CNLog2RatioEngine::runImp()
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

        m_objCNAnalysisEngine.setEngine(m_pEngine);
        m_objCNAnalysisEngine.createAnalysis();

        m_data.loadAnnotation(true);
        if (!m_data.loadModelFile(m_pEngine->getOpt("reference-file"))) {throw(Except("No reference file to process."));}

        m_iProbeSetCount = m_data.getProbeSets()->getCount();
        m_iExperimentCount = m_data.loadExperiments(m_pEngine->getOpt("expr-summary-file"));
        if (m_pEngine->isOptDefined("geno-qc-file"))
        {
            m_data.loadQCMetrics(m_pEngine->getOpt("geno-qc-file"));
        }
        if (m_iExperimentCount == 0) {throw(Except("No experiments to process."));}

        if (!determineMemoryUsage()) {return;}
        if (!processData()) {return;}

        m_data.clear();

    }
    catch(...)
    {
        m_data.clear();
        throw;
    }
}

/**
 * Guesstimate the maximum number of probe sets we can load for all experiments.
 * Guesstimate the maximum number of experiments we can load for all probe sets.
 * @return bool - true if successful
 */
bool CNLog2RatioEngine::determineMemoryUsage()
{
    Verbose::out(3, "CNLog2RatioEngine::determineMemoryUsage");
    bool bCyto2 = ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")));
    int iMemUsage = m_pEngine->getOptInt("mem-usage");
    int iMBUsedForPreAllocation = (int)((20 * m_iProbeSetCount) / 1000000.0);
    if (bCyto2) {iMBUsedForPreAllocation = (int)((8 * m_iProbeSetCount) / 1000000.0);}
    int iMBUsedForProbeSetList = (int)((77 * m_iProbeSetCount) / 1000000.0);
    int iMBUsedForExperimentList = (int)((58 * m_iExperimentCount) / 1000000.0);
    int iMBUsedForSnpArray = 140;
    int iMBFudgeFactor = 25;

    m_iProbeSetsToProcess = (int)((((iMemUsage - iMBUsedForPreAllocation - iMBUsedForProbeSetList - iMBUsedForExperimentList - iMBUsedForSnpArray - iMBFudgeFactor) * 1000000.0) / 21.0) / m_iExperimentCount);
    if (m_iProbeSetsToProcess <= 0) {throw(Except("Not enough memory to run application. ProbeSetsToProcess = " + ::getInt(m_iProbeSetsToProcess)));}
    if (m_iProbeSetsToProcess > m_iProbeSetCount) {m_iProbeSetsToProcess = m_iProbeSetCount;}

    m_iExperimentsToProcess = (int)((((iMemUsage - iMBUsedForPreAllocation - iMBUsedForProbeSetList - iMBUsedForExperimentList - iMBUsedForSnpArray - iMBFudgeFactor) * 1000000.0) / 21.0) / m_iProbeSetCount);
    if (m_iExperimentsToProcess <= 0) {throw(Except("Not enough memory to run application. ExperimentsToProcess = " + ::getInt(m_iExperimentsToProcess)));}
    if (m_iExperimentsToProcess > m_iExperimentCount) {m_iExperimentsToProcess = m_iExperimentCount;}
    Verbose::out(1, ::getInt(m_iExperimentsToProcess) + " ExperimentsToProcess out of " + ::getInt(m_iExperimentCount));
    return true;
}

/**
 * @brief Process the data.
 * @return bool - true if successful
 */
bool CNLog2RatioEngine::processData()
{
    Verbose::out(3, "CNLog2RatioEngine::processData");
    m_data.setLog2RatioEngine(true);
    int iExperimentIndex = 0;

    // Read the gender calls from the report file.
    m_data.loadGenotypeReport(m_pEngine->getOpt("genotype-report-file"));

    // Calculate the allelic differences and write them in transposed binary format.
//    Verbose::out(3, "ExperimentsToProcess = " + ::getInt(m_iExperimentsToProcess) + ", ProbeSetCount = " + ::getInt(m_iProbeSetCount));
    bool bFailed = true;
    int iMB = 80;
    while ((bFailed) && (iMB > 0))
    {
        try
        {
        Verbose::out(3, "AAlleleEstimates Allocation = " + ::getInt(m_iExperimentsToProcess * m_iProbeSetCount * 4));
        m_data.getAAlleleEstimates()->initialize(m_iExperimentsToProcess, m_iProbeSetCount);
        Verbose::out(3, "BAlleleEstimates Allocation = " + ::getInt(m_iExperimentsToProcess * m_iProbeSetCount * 4));
        m_data.getBAlleleEstimates()->initialize(m_iExperimentsToProcess, m_iProbeSetCount);
        Verbose::out(3, "GenotypeCalls Allocation = " + ::getInt(m_iExperimentsToProcess * m_iProbeSetCount * 1));
        m_data.getGenotypeCalls()->initialize(m_iExperimentsToProcess, m_iProbeSetCount);
        Verbose::out(3, "GenotypeConfidences Allocation = " + ::getInt(m_iExperimentsToProcess * m_iProbeSetCount * 4));
        m_data.getGenotypeConfidences()->initialize(m_iExperimentsToProcess, m_iProbeSetCount);
        Verbose::out(3, "Done with Allocations.");
        bFailed = false;
        } catch(...) {bFailed = true;}
        if (bFailed)
        {
            Verbose::out(3, "Memory allocations failed. Trying again...");
            iMB -= 10;
            m_iExperimentsToProcess = (int)Min((double)m_iExperimentsToProcess, ((iMB * 1000000.0) / m_iProbeSetCount / 4.0));
        }
    }
    if (iMB <= 0) {throw(Except("Memory allocations failed."));}
    iExperimentIndex = 0;
    Verbose::progressBegin(1, "CNLog2RatioEngine::processData() ", m_iExperimentCount, 1, m_iExperimentCount);
    while (iExperimentIndex < m_iExperimentCount)
    {
        if (!m_data.loadAlleleEstimatesByExperiment(m_pEngine->getOpt("expr-summary-file"), iExperimentIndex, m_iExperimentsToProcess)) {return false;}
        if (!m_data.loadGenotypeCallsByExperiment(m_pEngine->getOpt("genotype-calls-file"), iExperimentIndex, m_iExperimentsToProcess)) {return false;}
        if (!m_data.loadGenotypeConfidencesByExperiment(m_pEngine->getOpt("genotype-confidences-file"), iExperimentIndex, m_iExperimentsToProcess)) {return false;}

        for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + m_iExperimentsToProcess), m_data.getExperiments()->getCount())); iColIndex++)
        {
            Verbose::out(1, "Processing sample: " + Fs::basename(m_data.getExperiments()->getAt(iColIndex)->getExperimentName()));
            transferProbeSetData(iColIndex - iExperimentIndex);
            if (getOptBool("call-copynumber-engine")) {callCNAnalysisEngine(iColIndex);}
            if (getOptBool("log2ratio-text-output")) {writeOutputText(iColIndex);}
            if (getOptBool("log2ratio-hdf5-output")) {writeOutputHdf5(iColIndex);}
            Verbose::progressStep(1);
        }

        iExperimentIndex += m_iExperimentsToProcess;
    }
    Verbose::progressEnd(1, "Done");
    m_data.dumpGenotypeReport(m_pEngine->getOpt("genotype-report-file"));

    // Free up some memory.
    m_data.getAAlleleEstimates()->initialize(1);
    m_data.getBAlleleEstimates()->initialize(1);
    m_data.getGenotypeCalls()->initialize(1);
    m_data.getGenotypeConfidences()->initialize(1);

    if (getOptBool("delete-files")) {deleteFiles();}

    return true;
}

/**
 * @brief Initialize the ProbeSet Data for the next run of analysis.
 * @param int - The experiment index
 */
void CNLog2RatioEngine::transferProbeSetData(int iExperimentIndex)
{
    Verbose::out(3, "CNLog2RatioEngine::transferProbeSetDataToSnps");
    CNProbeSetArray* parProbeSets = m_data.getProbeSets();

    for (int iIndex = 0; (iIndex < parProbeSets->getCount()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = parProbeSets->getAt(iIndex);
        pobjProbeSet->setAAlleleSignal(m_data.getAAlleleEstimates()->get(iExperimentIndex, iIndex));
        pobjProbeSet->setBAlleleSignal(m_data.getBAlleleEstimates()->get(iExperimentIndex, iIndex));
        pobjProbeSet->setGenotypeCall(m_data.getGenotypeCalls()->get(iExperimentIndex, iIndex));
        pobjProbeSet->setGenotypeConfidence(m_data.getGenotypeConfidences()->get(iExperimentIndex, iIndex));
        pobjProbeSet->setCNState(0);
        pobjProbeSet->setSmoothedLog2Ratio(0);
        pobjProbeSet->setLoh(0);
    }
}

/**
 * @brief Run the CNAnalysisEngine
 * @param int - The experiment index
 */
void CNLog2RatioEngine::callCNAnalysisEngine(int iExperimentIndex)
{
    CNExperiment* pobjExperiment = m_data.getExperiments()->getAt(iExperimentIndex);
    m_objCNAnalysisEngine.process(*pobjExperiment, *m_data.getProbeSetsAlternateSort());
}

/**
 * @brief Ouptut the SNP data in ASCII text format.
 * @param int - The experiment index to process.
 * @return bool - true if successful
 */
bool CNLog2RatioEngine::writeOutputText(int iExperimentIndex)
{
    Verbose::out(3, "CNLog2RatioEngine::writeOutputText");
    CNExperiment* pobjExperiment = m_data.getExperiments()->getAt(iExperimentIndex);
    AffxString strFileName;
    if (m_pEngine->getOpt("set-analysis-name") != "")
    {
        strFileName = Fs::join(m_pEngine->getOpt("out-dir"),
                           m_pEngine->getOpt("set-analysis-name") + "." + pobjExperiment->getExperimentName() + ".txt");
    }
    else
    {
        strFileName = Fs::join(m_pEngine->getOpt("out-dir"),
                           pobjExperiment->getExperimentName() + ".txt");
    }
    affx::TsvFile tsv;

    // write out the options used.
    tsv.addHeader("guid", affxutil::Guid::GenerateNewGuid());

    // dump options
    vector<string> optionNames;
    m_pEngine->getOptionNames(optionNames,1);
    for(int i=0; i< optionNames.size(); i++) {
      std::string name = optionNames[i];
      std::vector<std::string> vals = m_pEngine->getOptVector(name,1);
      if (vals.size() > 1) {
        for(int j=0; j< vals.size(); j++) {
          tsv.addHeader("affymetrix-algorithm-param-apt-" + name, vals[j]);
        }
      }
      else {
        tsv.addHeader("affymetrix-algorithm-param-apt-" + name,  m_pEngine->getOpt(name,1));
      }
    }
    tsv.writeTsv_v1(strFileName);

    // Normally we don't have tables concatenated in this fashion.
    // There is a header row, record row, header row, records.
    // TSV does not support this so just write out the ad hoc lines.
    std::stringstream lineSStr;
    lineSStr << "experiment\tMadDiffCN\tIqr\tMeanAbsRle\tMedianAutosomeMedian\tGender\tGCCorrectionSize" << std::endl;
    tsv.write_str(lineSStr.str());
    lineSStr.str("");
    lineSStr << pobjExperiment->getExperimentName() + "\t" + ::getDouble(pobjExperiment->getMadDiffCN(), 6) + "\t" + ::getDouble(pobjExperiment->getIqr(), 6) + "\t" + ::getDouble(pobjExperiment->getMeanAbsRle(), 6) + "\t" + ::getDouble(pobjExperiment->getMedianAutosomeMedian(), 6) + "\t" + pobjExperiment->getCNCallGender() + "\t" + ::getDouble(pobjExperiment->getGCCorrectionMetric(), 6) << std::endl;
    tsv.write_str(lineSStr.str());
    lineSStr.str("");
    lineSStr << "probeset_id\tChromosome\tPosition\tLog2Ratio\tAllelicDifference\tGenotypeCall\tGenotypeConfidence\tChrXPar1\tChrXPar2" << std::endl;
    tsv.write_str(lineSStr.str());
    
    for (int iProbeSetIndex = 0, iCidx =0 ; (iProbeSetIndex < m_data.getProbeSets()->getCount()); iProbeSetIndex++, iCidx = 0)
    {
        CNProbeSet* pobjProbeSet = m_data.getProbeSets()->getAt(iProbeSetIndex);
        tsv.set(0, iCidx++, pobjProbeSet->getProbeSetName()); 
        tsv.set(0, iCidx++, ::getInt(pobjProbeSet->getChromosome())); 
        tsv.set(0, iCidx++, ::getInt((int)pobjProbeSet->getPosition())); 
        tsv.set(0, iCidx++, ::getDouble(pobjProbeSet->getLog2Ratio(), 10));
        tsv.set(0, iCidx++, ::getDouble(pobjProbeSet->getAllelicDifference(), 10));
        tsv.set(0, iCidx++, ::getInt((int)pobjProbeSet->getGenotypeCall())); 
        tsv.set(0, iCidx++, ::getDouble(pobjProbeSet->getGenotypeConfidence(), 10)); 
        tsv.set(0, iCidx++, ::getInt((int)pobjProbeSet->isPseudoAutosomalRegion()));
        tsv.set(0, iCidx++, "0"); 
        tsv.writeLevel(0);
    }

    tsv.clear();
    return true;
}

/**
 * @brief Ouptut the SNP data in HDF5 fformat.
 * @param int - The experiment index to process.
 * @return bool - true if successful
 */
bool CNLog2RatioEngine::writeOutputHdf5(int iExperimentIndex)
{
    Verbose::out(3, "CNLog2RatioEngine::writeOutputHdf5");
    CNExperiment* pobjExperiment = m_data.getExperiments()->getAt(iExperimentIndex);
    AffxString strFileName;
    if (m_pEngine->getOpt("set-analysis-name") != "")
    {
        strFileName = Fs::join(m_pEngine->getOpt("out-dir"),
                           m_pEngine->getOpt("set-analysis-name") + "." + pobjExperiment->getExperimentName() + ".a5");
    }
    else
    {
        strFileName = Fs::join(m_pEngine->getOpt("out-dir"),
                           pobjExperiment->getExperimentName() + ".a5");
    }
    try
    {
        affx::File5_File file5;
        file5.open(strFileName, affx::FILE5_REPLACE);

        affx::File5_Tsv* tsv5 = file5.openTsv("CNLog2Ratio_Parameters", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "Parameter", affx::FILE5_DTYPE_STRING, 1024);
        tsv5->set_string(0, 0, "#%guid=" + affxutil::Guid::GenerateNewGuid()); tsv5->writeLevel(0);

        // dump options
        vector<string> optionNames;
        m_pEngine->getOptionNames(optionNames,1);
        for(int i=0; i< optionNames.size(); i++) {
            std::string name = optionNames[i];
            std::vector<std::string> vals = m_pEngine->getOptVector(name,1);
            if (vals.size() > 1)
            {
                for(int j=0; j< vals.size(); j++)
                {
                    std::string val = vals[j];
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-" + name + "=" + val); tsv5->writeLevel(0);
                }
            }
            else
            {
                std::string val = m_pEngine->getOpt(name,1);
                tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-" + name + "=" + val); tsv5->writeLevel(0);
            }
        }

        tsv5->close();
        delete tsv5;

        tsv5 = file5.openTsv("CNLog2Ratio_Experiments", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "experiment", affx::FILE5_DTYPE_STRING, pobjExperiment->getExperimentName().length());
        tsv5->defineColumn(0, 1, "MadDiffCN", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 2, "Iqr", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 3, "MeanAbsRle", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 4, "MedianAutosomeMedian", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 5, "Gender", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 6, "GCCorrectionSize", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->set_string(0, 0, pobjExperiment->getExperimentName());
        tsv5->set_d(0, 1, pobjExperiment->getMadDiffCN());
        tsv5->set_d(0, 2, pobjExperiment->getIqr());
        tsv5->set_d(0, 3, pobjExperiment->getMeanAbsRle());
        tsv5->set_d(0, 4, pobjExperiment->getMedianAutosomeMedian());
        tsv5->set_string(0, 5, pobjExperiment->getCNCallGender());
        tsv5->set_d(0, 6, pobjExperiment->getGCCorrectionMetric());
        tsv5->writeLevel(0);
        tsv5->close();
        delete tsv5;

        int iMaxProbeSetNameLength = 0;
        for (int iProbeSetIndex = 0; (iProbeSetIndex < m_data.getProbeSets()->getCount()); iProbeSetIndex++)
        {
            CNProbeSet* pobjProbeSet = m_data.getProbeSets()->getAt(iProbeSetIndex);
            iMaxProbeSetNameLength = Max(iMaxProbeSetNameLength, (int)pobjProbeSet->getProbeSetName().length());
        }
        tsv5 = file5.openTsv("CNLog2Ratio_" + pobjExperiment->getExperimentName(), affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "probeset_id", affx::FILE5_DTYPE_STRING, iMaxProbeSetNameLength);
        tsv5->defineColumn(0, 1, "Chromosome", affx::FILE5_DTYPE_INT, 0);
        tsv5->defineColumn(0, 2, "Position", affx::FILE5_DTYPE_INT, 0);
        tsv5->defineColumn(0, 3, "Log2Ratio", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 4, "AllelicDifference", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 5, "GenotypeCall", affx::FILE5_DTYPE_INT, 0);
        tsv5->defineColumn(0, 6, "GenotypeConfidence", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 7, "ChrXPar1", affx::FILE5_DTYPE_INT, 0);
        tsv5->defineColumn(0, 8, "ChrXPar2", affx::FILE5_DTYPE_INT, 0);
        for (int iProbeSetIndex = 0; (iProbeSetIndex < m_data.getProbeSets()->getCount()); iProbeSetIndex++)
        {
            CNProbeSet* pobjProbeSet = m_data.getProbeSets()->getAt(iProbeSetIndex);
            tsv5->set_string(0, 0, pobjProbeSet->getProbeSetName());
            tsv5->set_i(0, 1, (int)pobjProbeSet->getChromosome());
            tsv5->set_i(0, 2, pobjProbeSet->getPosition());
            tsv5->set_d(0, 3, pobjProbeSet->getLog2Ratio());
            tsv5->set_d(0, 4, pobjProbeSet->getAllelicDifference());
            tsv5->set_i(0, 5, (int)pobjProbeSet->getGenotypeCall());
            tsv5->set_d(0, 6, pobjProbeSet->getGenotypeConfidence());
            tsv5->set_i(0, 7, (pobjProbeSet->isPseudoAutosomalRegion() ? 1 : 0));
            tsv5->set_i(0, 8, 0);
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;
        file5.close();
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    return true;
}

/**
 * @brief Delete temporary files
 */
void CNLog2RatioEngine::deleteFiles()
{
    try{::H5close();} catch(...) {} // Make sure HDF5 releases the files so that we can delete them.
    Fs::rm(Fs::dirname(Fs::join(m_pEngine->getOpt("expr-summary-file"),
                                                     m_pEngine->getOpt("set-analysis-name") + ".normalized-summary.a5")), false);
    Fs::rm(m_pEngine->getOpt("expr-summary-file"), false);
    Fs::rm(m_pEngine->getOpt("genotype-calls-file"), false);
    Fs::rm(m_pEngine->getOpt("genotype-confidences-file"), false);
    Fs::rm(m_pEngine->getOpt("genotype-report-file"), false);
    Fs::rm(Fs::join(m_pEngine->getOpt("out-dir"),m_pEngine->getOpt("expr-summary-file") + ".CopyNumber.Report.txt"), false);
}
