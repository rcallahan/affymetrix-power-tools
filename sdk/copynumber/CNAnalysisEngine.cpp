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
 * @file CNAnalysisEngine.cpp
 *
 * @brief This file contains the CNAnalysisEngine class members.
 */

#include "copynumber/CNAnalysisEngine.h"
//
#include "copynumber/CNAnalysisMethodCN.h"
#include "copynumber/CNAnalysisMethodCNCyto2.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNLog2RatioData.h"
#include "copynumber/CNReporterCychp.h"
//
#include "chipstream/HomHiLoCelListener.h"
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "file5/File5_Tsv.h"
//

using namespace std;

CNAnalysisEngine::Reg CNAnalysisEngine::reg;

CNAnalysisEngine * CNAnalysisEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == CNAnalysisEngine::EngineName())
        return (CNAnalysisEngine *)engine;
    return NULL;
}

/**
 * Constructor. Set up the CN Quant methods, and CN reporters we want to use.
 */
CNAnalysisEngine::CNAnalysisEngine()
{
    m_pEngine = this;

    // Reporters
    m_vCNReporters.push_back(new CNReporterCnchp());
    m_vCNReporters.push_back(new CNReporterCychp());

    defineOptions();
    m_bComponentRun = false;
}

/**
 * Destructor. Destroy the CN Quant methods, and CN reporters.
 */
CNAnalysisEngine::~CNAnalysisEngine()
{
    clear();
}

/**
 * Release the memory used by the engine.
 */
void CNAnalysisEngine::clear()
{
    m_vCNAnalysisMethods.deleteAll();
    m_vCNReporters.deleteAll();
}

/**
 * Define the options used by the engine.
 */
void CNAnalysisEngine::defineOptions()
{
  defineOptionSection("Input Options");
  defOptMult("", "log2ratio-file", PgOpt::STRING_OPT, "Log2Ratio file.", "");

  defineOptionSection("Analysis Options");
  defOptMult("a", "analysis", PgOpt::STRING_OPT,
                     "String representing analysis pathway desired.",
                     "");
  defineOption("", "gc-correction-bin-count", PgOpt::INT_OPT,
                     "The number of bins to use for GC content.",
                     "25");

  defineOptionSection("Advanced Options");
  defineOption("", "xChromosome", PgOpt::INT_OPT, "X Chromosome.", "24");
  defineOption("", "yChromosome", PgOpt::INT_OPT, "Y Chromosome.", "25");

  defineOptionSection("Execution Control Options");
  defineOption("", "mem-usage", PgOpt::INT_OPT,
                    "How many MB of memory to use for this run.",
                    "0");

  defineSharedOptions();
}

/**
 * Define options used by this engine, and also the CNWorkflowEngine.
 */
void CNAnalysisEngine::defineSharedOptions()
{
    m_pEngine->defineOptionSection("Additional CNAnalysisEngine Options");
    m_pEngine->defineOption("", "geno-qc-file", PgOpt::STRING_OPT, "The file output from GenoQC.", "");
    m_pEngine->defineOption("", "cyto2", PgOpt::BOOL_OPT, "Processing CYTO2 chip.", "false");
    m_pEngine->defineOption("", "array-name", PgOpt::STRING_OPT, "Array name or type to use.", "");
    m_pEngine->defineOption("", "set-analysis-name", PgOpt::STRING_OPT, "Analysis name to use as prefix for output files.", "");
    m_pEngine->defineOption("", "text-output", PgOpt::BOOL_OPT, "Output data in ASCII text format in addition to calvin format.", "false");
    for (int iIndex = 0; (iIndex < (int)m_vCNReporters.size()); iIndex++)
    {
        m_vCNReporters[iIndex]->defineOptions(*m_pEngine);
    }
}



/**
 * Define options used by this engine, and also the CNWorkflowEngine.
 */
void CNAnalysisEngine::defineSharedStates()
{
    if(!m_pEngine->isOptDefined("beta-cn-calibrate")) {
        m_pEngine->defineOption("", "beta-cn-calibrate", PgOpt::DOUBLE_OPT, "beta-cn-calibrate", "0");
        m_pEngine->defineOption("", "beta-X-cn-calibrate", PgOpt::DOUBLE_OPT, "beta-X-cn-calibrate", "0");
        m_pEngine->defineOption("", "beta-Y-cn-calibrate", PgOpt::DOUBLE_OPT, "beta-Y-cn-calibrate", "0");
        m_pEngine->defineOption("", "alpha-cn-calibrate", PgOpt::DOUBLE_OPT, "alpha-cn-calibrate", "1");
        m_pEngine->defineOption("", "alpha-X-cn-calibrate", PgOpt::DOUBLE_OPT, "alpha-X-cn-calibrate", "1");
        m_pEngine->defineOption("", "alpha-Y-cn-calibrate", PgOpt::DOUBLE_OPT, "alpha-Y-cn-calibrate", "1");
        m_pEngine->defineOption("", "gaussian-smooth-exp", PgOpt::INT_OPT, "gaussian-smooth-exp", "1");
        m_pEngine->defineOption("", "cancer", PgOpt::INT_OPT, "The state of cancer processing.", "0");
        m_pEngine->defineOption("", "hmmCN_state", PgOpt::STRING_OPT, "Used for calculating CN calibration parameters.", "");
        m_pEngine->defineOption("", "hmmCN_mu", PgOpt::STRING_OPT, "Used for calculating CN calibration parameters.", "");
    }
}

// Create all the specified analysis.
void CNAnalysisEngine::createAnalysis()
{
    Verbose::out(1, "CNAnalysisEngine::createAnalysis()");
    affymetrix_calvin_parameter::ParameterNameValueType param;
    std::vector<std::string>& analysis = m_pEngine->getOptVector("analysis");
    // Check for trailing "." (usually due to an extra space in the analysis string)
    for (int i = 0; i < analysis.size(); i++) {
        if (analysis[i][analysis[i].length() - 1] == '.') {
            Err::errAbort(std::string("Trailing '.' found in analysis string '") + analysis[i] + "'");
        }
    }
    m_vCNAnalysisMethods.deleteAll();
    CNAnalysisMethodFactory amFactory;
    for (unsigned int i = 0; i < analysis.size(); i++)
    {
        CNAnalysisMethod* am = NULL;
        if (analysis[i] == "") {throw(Except("Invalid (blank) analysis specification"));}

        am = amFactory.CNAnalysisMethodForString(analysis[i]);
        am->setEngine(m_pEngine);
        if (am->getName() == "cancer") {m_pEngine->setOpt("cancer", "1");}
        m_vCNAnalysisMethods.push_back(am);
    }
    if (needToCalculateCalibrationParameters()) {
        calculateCalibrationParameters();
    }
    else {
        CalibrationParams dAlphaCNCalibrate = retrieveCalibrationParams();
        m_pEngine->setOpt("alpha-cn-calibrate",::getDouble(dAlphaCNCalibrate.alphaCNCalibration, 6));
        m_pEngine->setOpt("alpha-X-cn-calibrate",::getDouble(dAlphaCNCalibrate.alphaCNCalibration_X, 6));
        m_pEngine->setOpt("alpha-Y-cn-calibrate",::getDouble(dAlphaCNCalibrate.alphaCNCalibration_Y, 6));
        m_pEngine->setOpt("beta-cn-calibrate",::getDouble(dAlphaCNCalibrate.betaCNCalibration, 6));
        m_pEngine->setOpt("beta-X-cn-calibrate",::getDouble(dAlphaCNCalibrate.betaCNCalibration_X, 6));
        m_pEngine->setOpt("beta-Y-cn-calibrate",::getDouble(dAlphaCNCalibrate.betaCNCalibration_Y, 6));
    }
}

bool CNAnalysisEngine::needToCalculateCalibrationParameters()
{
    std::string calibrateStr = m_pEngine->getOpt("cn-calibrate-parameters");
    return calibrateStr == "";
}

CalibrationParams CNAnalysisEngine::retrieveCalibrationParams()
{
    CalibrationParams params;

    // Default values for calibration parameters
    params.alphaCNCalibration = 0.564278;
    params.alphaCNCalibration_X = 0.619453;
    params.alphaCNCalibration_Y = 0.494620;
    params.betaCNCalibration = 1.0;
    params.betaCNCalibration_X = 1.0;
    params.betaCNCalibration_Y = 1.0;

    std::string calibrateStr = m_pEngine->getOpt("cn-calibrate-parameters");
    string name;
    map<string, string> param;
    SelfCreate::fillInNameParam(calibrateStr, name, param);
    map<string, string>::iterator iter;

    iter = param.find("alpha-cn-calibrate");
    if (iter != param.end())
    {
        params.alphaCNCalibration = Convert::toDouble(iter->second);
    }
    iter = param.find("alpha-X-cn-calibrate");
    if (iter != param.end())
    {
        params.alphaCNCalibration_X = Convert::toDouble(iter->second);
    }
    iter = param.find("alpha-Y-cn-calibrate");
    if (iter != param.end())
    {
        params.alphaCNCalibration_Y = Convert::toDouble(iter->second);
    }
    iter = param.find("beta-cn-calibrate");
    if (iter != param.end())
    {
        params.betaCNCalibration = Convert::toDouble(iter->second);
    }
    iter = param.find("beta-X-cn-calibrate");
    if (iter != param.end())
    {
        params.betaCNCalibration_X = Convert::toDouble(iter->second);
    }
    iter = param.find("beta-Y-cn-calibrate");
    if (iter != param.end())
    {
        params.betaCNCalibration_Y = Convert::toDouble(iter->second);
    }
    return params;
}

/**
 * @brief Calculate the parameters used to calibrate the log2 ratio to
 * the CN, and store the parameters in the engine's state.
 *
 * This information seems to be recovered by a getOpt call in
 * CNReporterCychp::loadChromosomes
 *
 * This should be moved out of this class.
 */

void CNAnalysisEngine::calculateCalibrationParameters()
{
    std::vector<int> vCN;
    std::vector<double> vMu;
    AffxByteArray ba;
    ba.assign(m_pEngine->getOpt("hmmCN_state"));
    ba.replace(",", " ");
    vCN.resize(ba.parameterCount());
    for (int k=0; k<ba.parameterCount(); k++)
    {
        vCN[k] = ba.getParameter(k+1).parseInt();
    }
    ba.assign(m_pEngine->getOpt("hmmCN_mu"));
    ba.replace(",", " ");
    vMu.resize(ba.parameterCount());
    for (int k=0; k<ba.parameterCount(); k++)
    {
        vMu[k] = ba.getParameter(k+1).parseDouble();
    }
    if (vCN.size() != vMu.size()) {throw(Except("Cannot calculate CN calibration parameters."));}
    double dBetaCNCalibrate = log2(2.0);
    m_pEngine->setOpt("beta-cn-calibrate",::getDouble(dBetaCNCalibrate, 6));
    m_pEngine->setOpt("beta-X-cn-calibrate",::getDouble(dBetaCNCalibrate, 6));
    m_pEngine->setOpt("beta-Y-cn-calibrate",::getDouble(dBetaCNCalibrate, 6));
    int iCount = 0;

    for (int k=0; (k<(int)vCN.size()); k++)
    {
        if (vCN[k] > 0) {iCount++;}
    }

    double dAlphaCNCalibrate = 1;
    if (iCount > 0)
    {
        AffxMultiDimensionalArray<double> mx(iCount, 2);
        int iIndex = 0;
        for (int k=0; k<(int)vCN.size(); k++)
        {
            if (vCN[k] > 0)
            {
                mx.set(iIndex, 0, log2((double)vCN[k]) - dBetaCNCalibrate);
                mx.set(iIndex, 1, vMu[k]);
                iIndex++;
            }
        }
        dAlphaCNCalibrate = mx.regress();
    }
    m_pEngine->setOpt("alpha-cn-calibrate",::getDouble(dAlphaCNCalibrate, 6));
    m_pEngine->setOpt("alpha-X-cn-calibrate",::getDouble(dAlphaCNCalibrate, 6));
    m_pEngine->setOpt("alpha-Y-cn-calibrate",::getDouble(dAlphaCNCalibrate, 6));
}

/**
 * Check to see if the options sepcified are valid.
 */
void CNAnalysisEngine::checkOptionsImp()
{
    defineSharedStates();

    if (m_pEngine == this)
    {
        if (getOpt("set-analysis-name") == "") {setOpt("set-analysis-name", "CN5");}
        setOpt("geno-qc-file", getOpt("geno-qc-file"));
    }
    if (m_pEngine->isOptDefined("log2ratio-file"))
    {
        std::vector<std::string> vLog2RatioFileNames = m_pEngine->getOptVector("log2ratio-file");
        if (vLog2RatioFileNames.size() == 0) {throw(Except("Invalid (missing) log2ratio-file specification"));}
        for (unsigned int i = 0; i < vLog2RatioFileNames.size(); i++)
        {
            if (vLog2RatioFileNames[i] == "") {throw(Except("Invalid (blank) log2ratio-file specification"));}
        }
    }
    std::vector<std::string>& analysis = m_pEngine->getOptVector("analysis");

    if(analysis.size() == 0) {
        m_bComponentRun = true;
        analysis.push_back("gaussian-smooth");
        analysis.push_back("cn-state");
        analysis.push_back("cn-gender");
        analysis.push_back("cn-segment");
        analysis.push_back("loh");
        analysis.push_back("loh-segment");
        analysis.push_back("cn-neutral-loh");
        analysis.push_back("normal-diploid");
        analysis.push_back("mosaicism");
    //        analysis.push_back("no-call");
    }
    createAnalysis();
    for (int iIndex = 0; (iIndex < (int)m_vCNReporters.size()); iIndex++)
    {
        m_vCNReporters[iIndex]->checkOptions(*m_pEngine);
    }
}

/**
 * Run the engine.
 */
void CNAnalysisEngine::runImp()
{
    AffxString str;
    int i = 0;
    int i1 = 0;
    int i2 = 0;
    double d = 0;
    affx::TsvFile tsv;
    CNLog2RatioData data;

    data.setEngine(m_pEngine);
    std::vector<std::string> vLog2RatioFileNames = m_pEngine->getOptVector("log2ratio-file");
    for (int iIndex = 0; (iIndex < (int)vLog2RatioFileNames.size()); iIndex++)
    {
        AffxString strFileName = vLog2RatioFileNames[iIndex];
        if (affx::File5_File::isHdf5file(strFileName))
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);

            affx::File5_Tsv* tsv5 = file5.openTsv("CNLog2Ratio_Parameters", affx::FILE5_OPEN);
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                AffxString str;
                tsv5->get(0, 0, &str);
            }
            tsv5->close();
            delete tsv5;

            tsv5 = file5.openTsv("CNLog2Ratio_Experiments", affx::FILE5_OPEN);
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                CNExperiment* pobjExperiment = new CNExperiment();
                tsv5->get(0, 0, &str); pobjExperiment->setExperimentName(str);
                tsv5->get(0, 1, &d); pobjExperiment->setMadDiffCN((float)d);
                tsv5->get(0, 2, &d); pobjExperiment->setIqr((float)d);
                tsv5->get(0, 3, &d); pobjExperiment->setMeanAbsRle((float)d);
                tsv5->get(0, 4, &d); pobjExperiment->setMedianAutosomeMedian((float)d);
                tsv5->get(0, 5, &str); pobjExperiment->setRawIntensityRatioGenderFromString(str);
                tsv5->get(0, 6, &d); pobjExperiment->setGCCorrectionMetric((float)d);
                int iSearchIndex = data.getExperiments()->binarySearch(*pobjExperiment, 0);
                if (iSearchIndex == -1)
                {
                    data.getExperiments()->add(pobjExperiment);
                    data.getExperiments()->quickSort(0); // ExperimentName
                }
            }
            tsv5->close();
            delete tsv5;
            file5.close();
        }
        else {
          tsv.m_optAutoTrim = true;
          if (tsv.open(strFileName) != affx::TSV_OK) {
            throw(Except("Cannot open file: " + strFileName));
          }
          int ridx = 0;
          std::vector< AffxString > cols;
          cols.resize( tsv.getColumnCount(0) );
          for ( int i = 0; i < tsv.getColumnCount(0) ; i++ ) {
            tsv.bind(0, i, &(cols[i]), affx::TSV_BIND_REQUIRED );
          }
          while ( tsv.nextLevel(0) == affx::TSV_OK ) {
            int cidx = 0;
            if ((ridx++ == 0) && (cols[0].startsWith("experiment")) ) {
              throw(Except("The input file does not appear to be in the correct format: " + strFileName));
            }
            CNExperiment* pobjExperiment = new CNExperiment;
            pobjExperiment->setExperimentName(cols[cidx++]);
            pobjExperiment->setMadDiffCN((float)AffxByteArray(cols[cidx++]).parseDouble());
            pobjExperiment->setIqr((float)AffxByteArray(cols[cidx++]).parseDouble());
            pobjExperiment->setMeanAbsRle((float)AffxByteArray(cols[cidx++]).parseDouble());
            pobjExperiment->setMedianAutosomeMedian((float)AffxByteArray(cols[cidx++]).parseDouble());
            AffxString strGender = cols[cidx++];
            pobjExperiment->setRawIntensityRatioGenderFromString(strGender);
            pobjExperiment->setGCCorrectionMetric((float)AffxByteArray(cols[cidx++]).parseDouble());
            int iSearchIndex = data.getExperiments()->binarySearch(*pobjExperiment, 0);
            if (iSearchIndex == -1)
            {
                data.getExperiments()->add(pobjExperiment);
                data.getExperiments()->quickSort(0); // ExperimentName
            }
            else {
                delete pobjExperiment;
            }
            tsv.close();
          }
        } // file5 vs tab delimited
    }
    if (m_pEngine->isOptDefined("geno-qc-file"))
    {
        data.loadQCMetrics(m_pEngine->getOpt("geno-qc-file"));
    }
    CNProbeSetArray vProbeSets;
    vProbeSets.reserve(3000000);
    for (int iExperimentIndex = 0; (iExperimentIndex < data.getExperiments()->getCount()); iExperimentIndex++)
    {
        // For each experiment, load the data into the ProbeSet array and write it out.
        CNExperiment* pobjExperiment = data.getExperiments()->getAt(iExperimentIndex);
        vProbeSets.deleteAll();
        for (int iIndex = 0; (iIndex < (int)vLog2RatioFileNames.size()); iIndex++)
        {
            AffxString strFileName = vLog2RatioFileNames[iIndex];
            if (affx::File5_File::isHdf5file(strFileName))
            {
                affx::File5_File file5;
                file5.open(strFileName, affx::FILE5_OPEN_RO);
                affx::File5_Tsv* tsv5 = file5.openTsv("CNLog2Ratio_" + pobjExperiment->getExperimentName(), affx::FILE5_OPEN);
                while (tsv5->nextLine() == affx::FILE5_OK)
                {
                    CNProbeSet* pobjProbeSet = new CNProbeSet();
                    tsv5->get(0, 0, &str); pobjProbeSet->setProbeSetName(str);
                    tsv5->get(0, 1, &i); pobjProbeSet->setChromosome((unsigned char)i);
                    tsv5->get(0, 2, &i); pobjProbeSet->setPosition(i);
                    tsv5->get(0, 3, &d); pobjProbeSet->setLog2Ratio((float)d);
                    tsv5->get(0, 4, &d); pobjProbeSet->setAllelicDifference((float)d);
                    tsv5->get(0, 5, &i); pobjProbeSet->setGenotypeCall((char)i);
                    tsv5->get(0, 6, &d); pobjProbeSet->setGenotypeConfidence((float)d);
                    tsv5->get(0, 7, &i1);
                    tsv5->get(0, 8, &i2);
                    if ((i1 == 1) || (i2 == 1)) {pobjProbeSet->setPseudoAutosomalRegion(true);}
                    vProbeSets.add(pobjProbeSet);
                }
                tsv5->close();
                delete tsv5;
                file5.close();
            }
            else {
              affx::TsvFile tsv;
              tsv.m_optAutoTrim = true;
              tsv.open(strFileName);
              std::string strExperimentName;
              tsv.bind(0,"experiment", &strExperimentName, affx::TSV_BIND_REQUIRED);


              tsv.nextLevel(0);
              
              if (strExperimentName == pobjExperiment->getExperimentName()) {
                AffxString strCol;
                int iCidx = 0;
                tsv.nextLevel(0);
                tsv.get(0,0,strCol);
                if (!strCol.startsWith("probeset_id")) {throw(Except("The input file does not appear to be in the correct format: " + strFileName));}
                while (tsv.nextLevel(0) == affx::TSV_OK ) {
                  iCidx = 0;
                  CNProbeSet* pobjProbeSet = new CNProbeSet();
                  tsv.get(0,iCidx++, strCol);
                  pobjProbeSet->setProbeSetName(strCol);
                  tsv.get(0,iCidx++, strCol);
                  pobjProbeSet->setChromosome((char)AffxByteArray(strCol).parseInt());
                  tsv.get(0,iCidx++, strCol);
                  pobjProbeSet->setPosition(AffxByteArray(strCol).parseInt());
                  tsv.get(0,iCidx++, strCol);
                  pobjProbeSet->setLog2Ratio((float)AffxByteArray(strCol).parseDouble());
                  tsv.get(0,iCidx++, strCol);
                  pobjProbeSet->setAllelicDifference((float)AffxByteArray(strCol).parseDouble());
                  tsv.get(0,iCidx++, strCol);
                  pobjProbeSet->setGenotypeCall((char)AffxByteArray(strCol).parseInt());
                  tsv.get(0,iCidx++, strCol);
                  pobjProbeSet->setGenotypeConfidence((float)AffxByteArray(strCol).parseDouble());
                  tsv.get(0,iCidx++, strCol);
                  bool bPAR1 = AffxByteArray(strCol).parsebool();
                  tsv.get(0,iCidx++, strCol);
                  bool bPAR2 = AffxByteArray(strCol).parsebool();
                  if ((bPAR1) || (bPAR2))   {
                    pobjProbeSet->setPseudoAutosomalRegion(true);
                  }
                  vProbeSets.add(pobjProbeSet);
                }
                // Look for duplicates.
                    vProbeSets.quickSort(0); // ProbeSetName
                    for (int iProbeSetIndex = 1; (iProbeSetIndex < vProbeSets.getCount()); iProbeSetIndex++)
                    {
                        CNProbeSet* pPrev = vProbeSets.getAt(iProbeSetIndex - 1);
                        CNProbeSet* pNext = vProbeSets.getAt(iProbeSetIndex);
                        if (pPrev->getProbeSetName() == pNext->getProbeSetName())
                        {
                            if ((pPrev->getChromosome() == pNext->getChromosome()) && (pPrev->getPosition() == pNext->getPosition()))
                            {
                                Verbose::out(1, "WARNING: Duplicate ProbeSet found. Removing the one found in file: " + strFileName);
                                vProbeSets.deleteAt(iProbeSetIndex);
                                iProbeSetIndex--;
                            }
                            else
                            {
                                throw(Except("Annotation mismatch for duplicate ProbeSet. Experiment: " + strExperimentName + ", ProbeSet: " + pPrev->getProbeSetName()));
                            }
                        }
                    }
                }
                tsv.clear();
            }
        }

        vProbeSets.quickSort(1); // Chromosome, Position, ProbeSetName
        process(*pobjExperiment, vProbeSets);
    }
    data.clear();
}

/**
 * Process the experiment. Call the plugged in CN Quant methods and then call the plugged in CN reporters.
 * @param CNExperiment& - The experiment being processed.
 * @param CNProbeSetArray& - The vector of ProbeSets associated with this experiment.
 */
void CNAnalysisEngine::process(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets, CNProbeArray* pvProbes)
{
	// Calculations
	if (vProbeSets.size() > 0)
	{
		for (int iIndex = 0; (iIndex < (int)m_vCNAnalysisMethods.size()); iIndex++)
		{
			if (m_pEngine->getOptInt("cancer") == 1)
			{
				if (m_vCNAnalysisMethods[iIndex]->getName() == "cancer")
				{
					m_vCNAnalysisMethods[iIndex]->setup(objExperiment, vProbeSets, pvProbes);
					m_vCNAnalysisMethods[iIndex]->run();
				}
			}
			else
			{
				if (m_vCNAnalysisMethods[iIndex]->getName() != "cancer")
				{
					m_vCNAnalysisMethods[iIndex]->setup(objExperiment, vProbeSets, pvProbes);
					m_vCNAnalysisMethods[iIndex]->run();
				}
			}
		}
	}
	if (m_pEngine->getOptInt("cancer") == 1)
	{
		m_pEngine->setOpt("cancer", "2");
		return;
	}
	postProcessing(objExperiment, vProbeSets);
	// Ouput data
	for (int iIndex = 0; (iIndex < (int)m_vCNReporters.size()); iIndex++)
	{
		m_vCNReporters[iIndex]->setup(*m_pEngine, objExperiment, vProbeSets, m_vCNAnalysisMethods);
		m_vCNReporters[iIndex]->run();
	}
}

/**
 * Calculate some sample metrics and store them in the experiment object.
 * @param CNExperiment& - The experiment being processed.
 * @param CNProbeSetArray& - The vector of ProbeSets associated with this experiment.
 */
void CNAnalysisEngine::postProcessing(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets)
{
    const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
	float fConfidenceThreshold = CNAnalysisMethod::getConfidenceThreshold(m_pEngine->getOpt("brlmmp-parameters"));
    int iSnpTotalCount = 0;
    int iHomTotalCount = 0;
    int iHetTotalCount = 0;
    for (int iIndex = 0; (iIndex < (int)vProbeSets.size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
        if (pobjProbeSet->processAsSNP())
        {
			if (isCytoScanHD && (pobjProbeSet->getGenotypeConfidence() >= fConfidenceThreshold)) {continue;}
            if ((pobjProbeSet->getGenotypeCall() == 0) || (pobjProbeSet->getGenotypeCall() == 2))
            {
                iHomTotalCount++;
            }
            if (pobjProbeSet->getGenotypeCall() == 1)
            {
                iHetTotalCount++;
            }
            iSnpTotalCount++;
        }
    }
    AffxMultiDimensionalArray<float> vCnState((int)vProbeSets.size());
    for (int iIndex = 0; (iIndex < (int)vProbeSets.size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
        int iXChromosome = m_pEngine->getOptInt("xChromosome");
        int iYChromosome = m_pEngine->getOptInt("yChromosome");
        float fCalibratedLog2Ratio = 0.0;
        if (pobjProbeSet->getChromosome() < iXChromosome) {
            fCalibratedLog2Ratio = pobjProbeSet->getCalibratedLog2Ratio(m_pEngine->getOptDouble("alpha-cn-calibrate"), m_pEngine->getOptDouble("beta-cn-calibrate"));
        }
        else if (pobjProbeSet->getChromosome() == iXChromosome) {
            fCalibratedLog2Ratio = pobjProbeSet->getCalibratedLog2Ratio(m_pEngine->getOptDouble("alpha-X-cn-calibrate"), m_pEngine->getOptDouble("beta-X-cn-calibrate"));
        }
        else if (pobjProbeSet->getChromosome() == iYChromosome) {
            fCalibratedLog2Ratio = pobjProbeSet->getCalibratedLog2Ratio(m_pEngine->getOptDouble("alpha-Y-cn-calibrate"), m_pEngine->getOptDouble("beta-Y-cn-calibrate"));
        }
        vCnState.set(iIndex, fCalibratedLog2Ratio);
    }
    objExperiment.setMedianCnState(vCnState.median());
    objExperiment.setHomFrequency((double)iHomTotalCount / (double)iSnpTotalCount);
    objExperiment.setHetFrequency((double)iHetTotalCount / (double)iSnpTotalCount);

    // Calculate SNPQC and RawSNPQC
    if (objExperiment.isSNPQCset()) {return;}
    if (!m_pEngine->isOptDefined("snp-qc-use-contrast")) {return;}
    if (!m_pEngine->isOptDefined("snp-qc-k")) {return;}
    if (!m_pEngine->isOptDefined("snp-qc-em-threshold")) {return;}
    if (!m_pEngine->isOptDefined("snp-qc-bin-size")) {return;}
    AffxArray<AffxString> arSNPIds;
    if ((m_pEngine->isOptDefined("snp-qc-snp-list")) && (m_pEngine->getOpt("snp-qc-snp-list") != ""))
    {
      affx::TsvFile tsv;
      tsv.m_optAutoTrim = true;
      tsv.openTable(m_pEngine->getOpt("snp-qc-snp-list"));
      std::string strSnp;
      while (tsv.nextLevel(0) == affx::TSV_OK) {
        tsv.get(0,0,strSnp);
        arSNPIds.add(new AffxString(strSnp));
      }
      tsv.close();

    }
    arSNPIds.quickSort(0);

    const double log2 = log(2.0);

    std::vector<double> vData;
    if (m_pEngine->getOptBool("cyto2")) {
        if (m_pEngine->getOptBool("snp-qc-use-contrast")) {
            if (arSNPIds.empty()) {
                for (int iIndex = 0; iIndex < vProbeSets.getCount(); iIndex++) {
                    CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
                    if (!pobjProbeSet->processAsVisualization()) {
                        continue;
                    }
                    vData.push_back(pobjProbeSet->getIntensityContrast(m_pEngine->getOptDouble("snp-qc-k")));
                }
            } else {
                for (int iIndex = 0; iIndex < vProbeSets.getCount(); iIndex++) {
                    CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
                    if (!pobjProbeSet->processAsVisualization()) {
                        continue;
                    }
                    AffxString str = pobjProbeSet->getProbeSetName();
                    if (arSNPIds.binarySearch(str, 0) != -1)
                    {
                        vData.push_back(pobjProbeSet->getIntensityContrast(m_pEngine->getOptDouble("snp-qc-k")));
                    }
                }
            }
        } else {
            if (arSNPIds.empty()) {
                for (int iIndex = 0; iIndex < vProbeSets.getCount(); iIndex++) {
                    CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
                    if (pobjProbeSet->getInformation()==0.0) {
                        continue;
                    }
                    if (!pobjProbeSet->processAsVisualization()) {
                        continue;
                    }
                    vData.push_back(pobjProbeSet->getSCAR());
                }
            } else {
                for (int iIndex = 0; iIndex < vProbeSets.getCount(); iIndex++) {
                    CNProbeSet* pobjProbeSet = vProbeSets.at(iIndex);
                    if (pobjProbeSet->getInformation()==0.0) {
                        continue;
                    }
                    if (!pobjProbeSet->processAsVisualization()) {
                        continue;
                    }
                    AffxString str = pobjProbeSet->getProbeSetName();
                    if (arSNPIds.binarySearch(str, 0) != -1)
                    {
                        vData.push_back(pobjProbeSet->getSCAR());
                    }
                }
            }
        }
    } else {
        // Make a copy of pointers to probe sets sorted by name
        CNProbeSetArray resortedProbeSets = vProbeSets;
        resortedProbeSets.quickSort(0);

        AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
        affx::File5_File file5;
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;

        file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
        group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
        tsv5 = group5->openTsv("CN5.snp-posteriors", affx::FILE5_OPEN);

        while (tsv5->nextLine() == affx::FILE5_OK) {
            string probeSetName;
            double clusterAAmean;
            double clusterABmean;
            double clusterBBmean;

            tsv5->get(0, 0, &probeSetName);
            tsv5->get(0, 2, &clusterAAmean);
            tsv5->get(0, 6, &clusterABmean);
            tsv5->get(0, 10, &clusterBBmean);

            if (!arSNPIds.empty()) {
                AffxString str(probeSetName);
                if (arSNPIds.binarySearch(str, 0) == -1) {
                    continue;
                }
            }
            CNProbeSet searchProbeSet;
            searchProbeSet.setProbeSetName(probeSetName);
            pair<CNProbeSetArray::iterator, CNProbeSetArray::iterator> range = equal_range(
                                                                                    resortedProbeSets.begin(),
                                                                                    resortedProbeSets.end(),
                                                                                    &searchProbeSet,
                                                                                    ProbeSetNameCompare()
                                                                                    );
            if (range.first == range.second) {
                continue;
            }

            CNProbeSet* foundProbeSet = *range.first;
            if (foundProbeSet->processAsVisualization())
            {
                double log2ratio = log(foundProbeSet->getAAlleleSignal()/foundProbeSet->getBAlleleSignal())/log2;
                double scaledLog2Ratio = log2ratio - clusterABmean;
                scaledLog2Ratio /= 0.5*(abs(clusterAAmean - clusterABmean) + abs(clusterBBmean - clusterABmean));
                vData.push_back(scaledLog2Ratio);
            }
        }
        tsv5->close();
        delete tsv5;

        group5->close();
        delete group5;
        file5.close();
    }
    arSNPIds.deleteAll();

    for (int i = 0; (i < vData.size()); i++)
    {
        if (vData[i] != vData[i])
        {
            vData[i] = 0;
        }
    }

    if (m_pEngine->getOptBool("cyto2")) {
        objExperiment.setSNPQC(HomHiLoCelListener::computeStatistic(vData, m_pEngine->getOptDouble("snp-qc-em-threshold"), m_pEngine->getOptDouble("snp-qc-bin-size")));
    } else {
        objExperiment.setSNPQC(HomHiLoCelListener::computeSNPQC(vData));
    }
    objExperiment.setIsSNPQCset(true);

    Verbose::out(1, "*");
    Verbose::out(1, "CNQC  (MAPD) = " + ToStr(objExperiment.getMadDiffCN()));
    if (m_pEngine->getOptBool("cyto2")) {
        if (m_pEngine->getOptBool("snp-qc-use-contrast"))
        {
            Verbose::out(1, "SNPQC (Contrast) = " + ToStr(objExperiment.getSNPQC()));
        }
        else 
        {
            Verbose::out(1, "SNPQC (PVQC) = " + ToStr(objExperiment.getSNPQC()));
        }
    } else {
        Verbose::out(1, "SNPQC (EMQC) = " + ToStr(objExperiment.getSNPQC()));
        Verbose::out(1, "RawSNPQC = " + ToStr(objExperiment.getRawSNPQC()));
    }
    Verbose::out(1, "*");
}
