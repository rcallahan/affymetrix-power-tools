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
 * @file CNCytoEngine.cpp
 *
 * @brief This file contains the CNCytoEngine class members.
 */

#include "copynumber/CNCytoEngine.h"
//
#include "copynumber/Annotation.h"
#include "copynumber/CNAnalysisEngine.h"
#include "copynumber/CNAnalysisMethodChipstream.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNAnalysisMethodReference.h"
#include "copynumber/CNAnalysisMethodCovariateParams.h"
//
#include "broadutil/BroadException.h"
#include "calvin_files/exception/src/ExceptionBase.h"
#include "chipstream/CnProbeGenderCelListener.h"
#include "chipstream/EngineUtil.h"
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "util/Fs.h"
#include "util/TmpFileFactory.h"
//
#include "../external/newmat/myexcept.h"
//
#include <sstream>
//


CNCytoEngine::Reg CNCytoEngine::reg;

CNCytoEngine * CNCytoEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == CNCytoEngine::EngineName())
        return (CNCytoEngine *)engine;
    return NULL;
}

/**
 * @brief Constructor
 */
CNCytoEngine::CNCytoEngine()
{
    m_pstdMethods = NULL;
    defineStdMethods();
    defineOptions();
}

/**
 * @brief Destructor
 */
CNCytoEngine::~CNCytoEngine()
{
    clear();
}

/**
 * @brief Initialize the data and free any memory allocated.
 */
void CNCytoEngine::clear()
{
    if (m_pstdMethods != NULL) {delete m_pstdMethods; m_pstdMethods = NULL;}
    // Clean up static data.
    CNExperiment::getQCMetricColumnNames()->deleteAll();
    CNAnalysisMethod::getCelFileParams()->clear();
    Annotation::getParams()->clear();
    CNAnalysisMethod::getParams()->clear();
}

/**
 * @brief Define any standard methods used by the engine
 */
void CNCytoEngine::defineStdMethods() {
//    m_stdMethods["brlmm"] = "quant-norm.sketch=50000,pm-only,brlmm.transform=ccs.K=4";
}

/**
 * @brief Print any standard methods used by the engine
 * @param std::ostream& - The specified output stream
 */
void CNCytoEngine::printStandardMethods(std::ostream &out) {
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
void CNCytoEngine::defineOptions() {

  defineOptionSection("Input Options");
  defineOption("", "chipstream", PgOpt::STRING_OPT,
                    "String representing chipstream parameters. This includes normal-diploid analysis parameters .",
                    "chipstream");
  defineOption("", "check-input-files", PgOpt::BOOL_OPT,
                     "Does an upfront check of the syntax and content of the input files: gender-override, genotype-call-override, snp-reference-input files.",
                     "true");
  defineOption("", "doDualNormalization", PgOpt::BOOL_OPT,
                     "The usual default action is to normalize SNP and CN probeset separately.  This option allows a single normalization set for CytoScanHD chips.  ",
                     "true");
  defineOption("", "config-file", PgOpt::STRING_OPT,
                    "The configuration file name as passed from GTC or the Cyto Browser.",
                    "");
  defineOption("", "antigenomic-probe-file", PgOpt::STRING_OPT,
                    "The probe indexes (ProbeID - 1) of the antigenomic probes on the array.",
                    "");
  defineOption("", "snp-reference-input-file", PgOpt::STRING_OPT,
                    "Input SNP reference file name.", "");
  defineOption("", "reference-input", PgOpt::STRING_OPT, "Input reference file name.", "");
  defineOption("", "probe-file", PgOpt::STRING_OPT,
                     "File defining probe sequences and locations.",
                     "");
  defineOption("", "cdf-file", PgOpt::STRING_OPT,
                     "File defining probe sets.",
                     "");
  defineOption("", "spf-file", PgOpt::STRING_OPT,
                     "spf format file defining probe sets.",
                     "");
  defineOption("", "qcc-file", PgOpt::STRING_OPT,
                     "File defining QC probesets.",
                     "");
  defineOption("", "qca-file", PgOpt::STRING_OPT,
                     "File defining QC analysis methods.",
                     "");
  defineOption("", "cel-files", PgOpt::STRING_OPT,
                     "Text file specifying cel files to process, one per line with the "
                     "first line being 'cel_files'.",
                     "");
  defineOption("", "chrX-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrX. "
                     "Used for copy number probe chrX/Y ratio gender calling.",
                     "");
  defineOption("", "chrY-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrY. "
                     "Used for copy number probe chrX/Y ratio gender calling.",
                     "");
  defineOption("", "normal-diploid-files-file", PgOpt::STRING_OPT,
                   "Text file specifying normal-diploid probeset files. First line must be 'normal_diploid_files'.",
                     "");
  defOptMult("", "reference-cels", PgOpt::STRING_OPT,
                     "'Normal' CEL file(s) to process when doing paired analysis of Cancer vs. Normal.",
                     "");
  defineOption("","annotation-file", PgOpt::STRING_OPT,
                     "NetAffx Annotation database file.",
                     "");
  defOptMult("", "normal-diploid-files", PgOpt::STRING_OPT,
                     "Normal Diploid probeset file names.",
                     "");
  defineOption("", "sketch-size", PgOpt::INT_OPT,
                    "The number number of point to be used for a sketch.",  "50000");
  defineOption("","covariates-file", PgOpt::STRING_OPT,
                     "External covariates file.",
                     "");
  defineOption("", "minSegSeparation", PgOpt::INT_OPT,
                    "Value used to skip over centromere in LOH.",  "1000000000");

  defineOptionSection("Output Options");
  defineOption("", "snp-reference-output-file", PgOpt::STRING_OPT,
                    "Output SNP reference file name.",
                    "");
  defineOption("", "reference-output", PgOpt::STRING_OPT,
                    "Output reference file name.",
                    "");
  defineOption("", "file-name-prefix", PgOpt::STRING_OPT,
                    "CYCHP file name prefix.",
                    "");
  defineOption("", "file-name-suffix", PgOpt::STRING_OPT,
                    "CYCHP file name suffix.",
                    "");
  defineOption("", "file-name-ext", PgOpt::STRING_OPT,
                    "CYCHP file name extension.",
                    "cychp");

  defineOptionSection("Analysis Options");
  defineOption("", "run-geno-qc", PgOpt::BOOL_OPT,
                     "Run the GenoQC engine.",
                     "false");
  defineOption("", "adapter-type-normalization", PgOpt::BOOL_OPT,
                    "Adapter Type Normalization option. true = perform adapter type normalization.",
                    "false");
  defOptMult("a", "analysis", PgOpt::STRING_OPT,
                    "String representing analysis pathway desired.",
                    "");
  defineOption("", "wave-correction-reference-method", PgOpt::STRING_OPT,
                    "String representing wave correction pathway desired.",
                    "");
  defineOption("", "local-gc-background-correction-reference-method", PgOpt::STRING_OPT,
                    "String representing local gc correction background pathway desired.",
                    "none");
  defineOption("", "local-gc-background-intensity-adjustment-method", PgOpt::STRING_OPT,
                    "String representing local gc correction background intensity adjustment "
                    "pathway desired.",
                    "");
  defineOption("", "image-correction-intensity-adjustment-method", PgOpt::STRING_OPT,
                    "String representing image correction intensity adjustment pathway desired.",
                    "");
  defineOption("", "log2ratio-adjustment-method", PgOpt::STRING_OPT,
                    "String representing high pass filter parameters for modification of log2ratios.",
                    "");
  defineOption("", "wave-correction-log2ratio-adjustment-method", PgOpt::STRING_OPT,
                    "String representing wave correction log2ratio adjustment pathway desired.",
                    "");
  defineOption("", "allele-peaks-reporter-method", PgOpt::STRING_OPT,
                    "String representing allele peaks reporter pathway desired.",
                    "allele-peaks-reporter-method");
  defineOption("", "hmm-priors", PgOpt::BOOL_OPT,
                    "Read HMM priors from the reference file.",
                    "false");
  defineOption("", "hmm-means-check", PgOpt::BOOL_OPT,
                    "Check HMM priors in the reference file for consistency.",
                    "false");
  defineOption("", "brlmmp-parameters", PgOpt::STRING_OPT,
                     "Parameters to use when running brlmmp.",
                     "");
  defineOption("", "signal-adjustment-covariates", PgOpt::STRING_OPT,
                    "String representing the covariate-based signal adjustment pathway desired.",
                    "");
  defineOption("", "lr-adjustment-covariates", PgOpt::STRING_OPT,
                    "String representing the covariate-based log2 ratio adjustment pathway desired.",
                    "");
  defineOption("", "allele-peaks-adjustment-covariates", PgOpt::STRING_OPT,
                   "String representing the covariate-based allele peaks adjustment pathway desired.",
                   "");

  defineOptionSection("QC Options");
  defineOption("", "snp-qc-use-contrast", PgOpt::BOOL_OPT,
                    "SNP QC use intensity contrast.", "false");
  defineOption("", "snp-qc-snp-list", PgOpt::STRING_OPT,
                    "Input file containing a list of SNP Ids to be used in calculating the SNPQC value.",
                    "");
  defineOption("", "snp-qc-k", PgOpt::DOUBLE_OPT,
                    "SNP QC K value.", "2.0");
  defineOption("", "snp-qc-em-threshold", PgOpt::DOUBLE_OPT,
                    "SNP QC EM Threshold.", "0.05");
  defineOption("", "snp-qc-bin-size", PgOpt::DOUBLE_OPT,
                    "SNP QC bin size.", "0.04");

  defineOptionSection("Misc Options");
  defineOption("", "explain", PgOpt::STRING_OPT,
                    "Explain a particular operation (i.e. --explain cn-state or --explain loh).",
                    "");

  defineOptionSection("Advanced Options");
  defineOption("", "probeset-ids", PgOpt::STRING_OPT,
                     "Tab delimited file with column 'probeset_id' specifying probesets to summarize.",
                     "");
  defineOption("", "global-parameter-override", PgOpt::BOOL_OPT,
                    "Global parameters in loh-cyto analysis input string will override those "
                    "given in the snp-reference-input-file.", "false" );
  defineOption("", "keep-temp-reference-data", PgOpt::BOOL_OPT,
                    "Set to true, this option will keep the final signal and intensity values computed "
                    "while determining the CN reference file.  Warning: It may be large.", "false" );
  defineOption("", "keep-intermediate-data", PgOpt::BOOL_OPT,
                    "Set to true, this option will keep all, intensity values computed "
                    "while invoking any intensity adjustment method.", "false" );
  defineOption("", "keep-intermediate-data-local", PgOpt::BOOL_OPT,
                    "Set to true, this option will keep all, intensity values computed "
                    "while invoking any intensity adjustment method. This is a duplicate option for local testing.", "false" );
  defineOption("", "genotype-call-override-file", PgOpt::STRING_OPT,
                    "Input file containing genotype calls to be used in place of brlmmp-p calls.", "");
  defineOption("", "gender-override-file", PgOpt::STRING_OPT,
                    "A file containing externally computed genders. It is used in reference "
                    "and snp-reference generation", "");
  defineOption("", "gc-content-override-file", PgOpt::STRING_OPT,
                    "Input file used to override the GC content read from the annotation "
                    "files (Two columns with header line, ProbeSetName/GCContent).", "");
  defineOption("", "gc-correction-bin-count", PgOpt::INT_OPT,
                    "The number of bins to use for GC content.", "25");
  defineOption("","xChromosome", PgOpt::INT_OPT,
                    "X Chromosome",
                    "24");
  defineOption("","yChromosome", PgOpt::INT_OPT,
                    "Y Chromosome",
                    "25");
  defineOption("", "male-gender-ratio-cutoff", PgOpt::DOUBLE_OPT,
                    "Male gender ratio cutoff", "1.3");
  defineOption("", "female-gender-ratio-cutoff", PgOpt::DOUBLE_OPT,
                    "Female gender ratio cutoff", "1.0");
  defineOption("", "reference-chromosome", PgOpt::INT_OPT,
                    "Reference chromosome", "2");
  defineOption("", "xx-cutoff", PgOpt::DOUBLE_OPT,
                    "XX cutoff", "0.8");
  defineOption("", "xx-cutoff-high", PgOpt::DOUBLE_OPT,
                    "XX cutoff high", "1.07");
  defineOption("", "y-cutoff", PgOpt::DOUBLE_OPT, "Y cutoff", "0.65");
  defineOption("", "warning-message-limit", PgOpt::INT_OPT, "Set the limit on the number of warning message produced.", "10");
//  defineOption("", "ad-limit", PgOpt::DOUBLE_OPT, "Allelic Difference Limit", "0.2");
//  defineOption("", "xx-cutoff-ad", PgOpt::DOUBLE_OPT, "XX cutoff Allelic Difference", "0.066");
//  defineOption("", "xx-cutoff-override", PgOpt::DOUBLE_OPT, "XX cutoff Override", "0.9");

  defineOption("", "waviness-block-size", PgOpt::INT_OPT,
                    "marker count", "50");
  defineOption("", "waviness-genomic-span", PgOpt::INT_OPT,
                    "genomic segment length", "0");
  defineOption("", "cn-calibrate-parameters", PgOpt::STRING_OPT,
                    "SmoothSignal calibration parameters", "");

  defineOption("", "use-old-kdensity-function", PgOpt::BOOL_OPT,
                   "For testing only, revert to the old slow kdensity method for density estimation.", "false");

  defineOptionSection("Engine Options (Not used on command line)");
  defOptMult("", "cels", PgOpt::STRING_OPT,
                     "CEL files to process.",
                     "");
  defOptMult("", "arrs", PgOpt::STRING_OPT,
                     "ARR files to process. Must be paired with cels.",
                     "");
  defOptMult("", "result-files", PgOpt::STRING_OPT,
                     "CYCHP files to output. Must be paired with cels.",
                     "");

  CNAnalysisEngine objCNAnalysisEngine;
  objCNAnalysisEngine.setEngine(this);
  objCNAnalysisEngine.defineSharedOptions();
}

/**
 * @brief Add a cel file specification to the options
 * @param const std::string& - The cel file name
 * @param const std::string& - The sample (arr) file name
 * @param const std::string& - The result (cychp) file name
 */
void CNCytoEngine::addCel(const std::string& strCelFileName, const std::string& strSampleFileName, const std::string& strResultFileName)
{
    pushOpt("cels", strCelFileName);
    pushOpt("arrs", strSampleFileName);
    pushOpt("result-files", strResultFileName);
}

/**
 * @brief Define the states used by this engine
 */
void CNCytoEngine::defineStates() {

    defineOption("", "cytoscan-hd", PgOpt::BOOL_OPT, "The array is the CytoScanHD product.", "false");
	setOpt("cytoscan-hd", ((Annotation::isCytoScanHD(getOpt("annotation-file"))) ? "true" : "false"));
    defineOption("", "num-rows", PgOpt::INT_OPT, "The number of rows on the chip.", "-1");
    defineOption("", "num-cols", PgOpt::INT_OPT, "The number of cols on the chip.", "-1");
    defineOption("", "reference-wave-count", PgOpt::INT_OPT, "wave-count-used", "-1");
    defineOption("", "sample-count", PgOpt::INT_OPT, "sample-count", "0");
    defineOption("", "reference-file", PgOpt::STRING_OPT, "reference-file", "");
    defineOption("", "create-snp-reference", PgOpt::BOOL_OPT, "create-snp-reference", "false");
    defineOption("", "create-reference", PgOpt::BOOL_OPT, "create-reference", "false");
    defineOption("", "copynumber-reference-file", PgOpt::STRING_OPT, "copynumber-reference-file", "");
    defineOption("", "reference-text-output", PgOpt::BOOL_OPT, "The reference file is in ascii text format.", "false");
    defineOption("", "recommended-reference-sample-count", PgOpt::INT_OPT, "The recommended number of samples for a Copy Number reference.", "44");
    defineOption("", "recommended-reference-female-sample-count", PgOpt::INT_OPT, "The recommended number of female samples for a Copy Number reference.", "20");
    defineOption("", "recommended-reference-male-sample-count", PgOpt::INT_OPT, "The recommended number of male samples for a Copy Number reference.", "20");
    defineOption("", "reference-sample-count", PgOpt::INT_OPT, "The actual number of samples for the Copy Number reference.", "0");
    defineOption("", "reference-female-sample-count", PgOpt::INT_OPT, "The actual number of female samples for the Copy Number reference.", "0");
    defineOption("", "reference-male-sample-count", PgOpt::INT_OPT, "The actual number of male samples for the Copy Number reference.", "0");
    defineOption("", "reference-unknown-sample-count", PgOpt::INT_OPT, "The actual number of unknown gender samples for the Copy Number reference.", "0");
    defineOption("", "probe-count", PgOpt::INT_OPT, "The actual number of probes on the chip.", "0");
    defineOption("", "probe-rows-count", PgOpt::INT_OPT, "The number of probe rows on the chip.", "0");
    defineOption("", "probe-columns-count", PgOpt::INT_OPT, "The number of probe columns on the chip.", "0");
    defineOption("", "probeset-count", PgOpt::INT_OPT, "The number of probesets on the chip.", "-1");
    defineOption("", "log2ratio-hdf5-output", PgOpt::BOOL_OPT, "Write intermediate output in HDF5 format.", "false");
    defineOption("", "log2ratio-text-output", PgOpt::BOOL_OPT, "Write intermediate output in ASCII Text format.", "false");
    defineOption("", "yTarget", PgOpt::DOUBLE_OPT, "Store yTarget for sharing between modules.", "NaN");
    defineOption("", "array-size", PgOpt::INT_OPT, "Array size for determining which array we are running full or focused.", "2015");
    defineOption("", "temp-reference-file", PgOpt::STRING_OPT, "temp-reference-file", "");
    defineOption("", "pdnn-reference-values-exist", PgOpt::BOOL_OPT, "A boolean value set to true if pdnn parameters exist in the CN reference-input file.", "false");
    defineOption("", "chrX-probes-filtered", PgOpt::STRING_OPT, "Temp file for filtered chr-X probes.", "");
    defineOption("", "chrY-probes-filtered", PgOpt::STRING_OPT, "Temp file for filtered chr-Y probes.", "");
    defineOption("", "disable-covariates-file-warning", PgOpt::BOOL_OPT, "Indicates that the covariates file warning should be disabled to avoid redundant warnings.", "false");
    defineOption("", "AD-modulation-no-single-peak", PgOpt::BOOL_OPT, "Indicates whether allelic differences have been changed by the AD modulation correction method", "false");

    CNAnalysisEngine objCNAnalysisEngine;
    objCNAnalysisEngine.setEngine(this);
    objCNAnalysisEngine.defineSharedStates();
}

class Analysis
{
public:
    int Order;
    std::string Name;

    AffxString getName()
    {
        AffxString str = Name;
        int iFindIndex = str.indexOf(".");
        if (iFindIndex != std::string::npos)
        {
            str = str.substring(0, iFindIndex);
        }
        return str;
    }

    int compareTo(Analysis& that, int iCompareCode)
    {
        int iCompareResult = 0;
        switch (iCompareCode)
        {
        case 0:
            iCompareResult = AffxArray<int>::compare(Order, that.Order);
            break;
        }
        return iCompareResult;
    }

    template<int k> struct ComparePred {
        bool operator()(const Analysis* lhs, const Analysis* rhs) const {
            Err::errAbort("CNCytoEngine/Analysis: ComparePred instantiated with an invalid compare code = " + ToStr(k));
            return false;
        }
    };
};

template<> struct Analysis::ComparePred<0> {
    bool operator()(const Analysis* lhs, const Analysis* rhs) const {
        return AffxArray<int>::compare(lhs->Order, rhs->Order) < 0;
    }
};

void CNCytoEngine::checkParameters()
{
    defineStates();
    if (getOpt("reference-input") != "")
    {
        if ( (getOpt("local-gc-background-intensity-adjustment-method") != "")  && (getOpt("local-gc-background-intensity-adjustment-method") != "none") )
        {
                if(!predictedIntensitiesExist(getOpt("reference-input")))
                {
                    Err::errAbort("PDNN predicted intensity values are missing from input reference file. To use the pdnn intensity adjustment while analyzing samples you must use a reference input file where these values were created.  To create these values use the the option --local-gc-background-correction-reference-method pdnn-reference-method  when generating the reference input file. You may be using an out-dated input reference file.");
                }
                else
                {
                    setOpt("pdnn-reference-values-exist", "true");
                }
        }
    } 

    if ((getOptBool("cyto2")) && (getOptBool("adapter-type-normalization")))
    {
        Err::errAbort("adapter-type-normalization cannot be run against the cyto2 arrays.");
    }

    std::vector<std::string>& analysis = getOptVector("analysis");

    if(analysis.size() == 0) {
        if (getOptBool("cyto2"))
        {
            analysis.push_back("log2-ratio-cyto2");
            analysis.push_back("kernel-smooth");
            analysis.push_back("cn-cyto2");
            analysis.push_back("cn-cyto2-gender");
            analysis.push_back("cn-segment");
            analysis.push_back("loh-cyto2");
            analysis.push_back("mosaicism");
            analysis.push_back("allele-peaks");
        }
        else
        {
			if (getOptBool("cytoscan-hd"))
			{
				analysis.push_back("log2-ratio");
				analysis.push_back("allelic-difference");
				analysis.push_back("kernel-smooth");
				analysis.push_back("cn-cyto2");
				analysis.push_back("cn-cyto2-gender");
				analysis.push_back("cn-segment");
				analysis.push_back("lohCytoScan");
				analysis.push_back("loh-segment");
				analysis.push_back("cn-neutral-loh");
				analysis.push_back("genotype");
			}
			else
			{
				analysis.push_back("log2-ratio");
				analysis.push_back("allelic-difference");
				analysis.push_back("gaussian-smooth");
				analysis.push_back("cn-state");
				analysis.push_back("cn-gender");
				analysis.push_back("cn-segment");
				analysis.push_back("loh");
				analysis.push_back("loh-segment");
				analysis.push_back("cn-neutral-loh");
				analysis.push_back("normal-diploid");
				analysis.push_back("mosaicism");
			}
        }
    }

    // Make consistency checks
    checkConsistency(analysis);

    // Insure proper analysis order.
    bool bCyto2LOH = false;
    AffxArray<Analysis> ar;
    for (int i = 0; (i < analysis.size()); i++)
    {
        Analysis* p = new Analysis;
        p->Name = analysis[i];
        p->Order = CNSegment::getSegmentType(p->getName());
        ar.add(p);
    }
    ar.quickSort(0);
    for (int i = 0; (i < analysis.size()); i++)
    {
        Analysis* p = ar.getAt(i);
        analysis[i] = p->Name;
        if (p->getName() == "loh-cyto2") {bCyto2LOH = true;}
    }
    ar.deleteAll();


    CNAnalysisEngine objCNAnalysisEngine;
    objCNAnalysisEngine.setEngine(this);
    objCNAnalysisEngine.checkOptions();
    objCNAnalysisEngine.clear();

    if ((!getOptBool("cychp-output")) && (!getOptBool("cnchp-output")))
    {
        Err::errAbort("Must specify either cychp-ouptut or cnchp-output.");
    }

    CovariateParams::checkParams(*this);
}

/**
 * @brief Make sure that our options are sane. Call Err::errAbort if not.
 */
void CNCytoEngine::checkOptionsImp() {

    //defineStates();

    // file is not used by the engine
    //setLibFileOpt("config-file");
    setLibFileOpt("snp-reference-input-file");
    setLibFileOpt("reference-input");
    setLibFileOpt("cdf-file");
    setLibFileOpt("spf-file");
    setLibFileOpt("qcc-file");
    setLibFileOpt("qca-file");
    setLibFileOpt("chrX-probes");
    setLibFileOpt("chrY-probes");
    setLibFileOpt("annotation-file");
    setLibFileOpt("probeset-ids");
    setLibFileOpt("gc-content-override-file");
    setLibFileOpt("snp-qc-snp-list");
    setLibFileOpt("probe-file");
    setLibFileOpt("covariates-file");

    if(getOpt("explain") != "") { explain(); exit(0); }

    if (getOpt("cel-files") != "" && getOpt("cels") != "") {
        if (getArgCount() > 0) {
            error("CEL files cannot be specified both as command line arguments and as lists with \"cel-files\"/\"cels\" options");
        }
    }

    checkParameters();
    Verbose::out(4, "Parameters are ok.");

    if(getOptBool("cyto2"))
    {
        if(!getOptBool("doDualNormalization"))
        {
            error("You cannot turn off dual normalization when the cyto2 input parameter is true.  This program will not allow single normalization on 2.7M or 310K chips.  It may only be turned off when analyzing snp6 and snp7 chips."  );
        }
    }

    if (getOpt("out-dir") == "") {error("Must specify an output directory.");}
    if (getOpt("temp-dir") == "") {
      setOpt("temp-dir", Fs::join(getOpt("out-dir"),"temp"));
    }

    std::vector<std::string>& analysis = getOptVector("analysis");
    if(analysis.size() == 0) {
        if (getOptBool("cyto2"))
        {
            analysis.push_back("log2-ratio-cyto2");
            analysis.push_back("kernel-smooth");
            analysis.push_back("cn-cyto2");
            analysis.push_back("cn-cyto2");
            analysis.push_back("cn-cyto2-gender");
            analysis.push_back("cn-segment");
            analysis.push_back("loh-cyto2");
            analysis.push_back("mosaicism");
            analysis.push_back("allele-peaks");
        }
        else
        {
            analysis.push_back("log2-ratio");
            analysis.push_back("allelic-difference");
            analysis.push_back("gaussian-smooth");
            analysis.push_back("cn-state");
            analysis.push_back("cn-gender");
            analysis.push_back("cn-segment");
            analysis.push_back("loh");
            analysis.push_back("loh-segment");
            analysis.push_back("cn-neutral-loh");
            analysis.push_back("normal-diploid");
            analysis.push_back("mosaicism");
            analysis.push_back("cn-snp7");
        }
    }

    // Insure proper analysis order.
    bool bCyto2LOH = false;
    AffxArray<Analysis> ar;
    for (int i = 0; (i < analysis.size()); i++)
    {
        Analysis* p = new Analysis;
        p->Name = analysis[i];
        p->Order = CNSegment::getSegmentType(p->getName());
        ar.add(p);
    }
    ar.quickSort(0);
    for (int i = 0; (i < analysis.size()); i++)
    {
        Analysis* p = ar.getAt(i);
        analysis[i] = p->Name;
        if (p->getName() == "loh-cyto2") {bCyto2LOH = true;}
    }
    ar.deleteAll();

    if (getOptBool("keep-temp-reference-data") )
    {
            setOpt("keep-temp-reference-data","true");
    }
    if (getOptBool("keep-intermediate-data") )
    {
            setOpt("keep-intermediate-data","true");
    }

    if ((getOpt("reference-input") != "") && (getOpt("reference-output") != ""))
    {
        error("Please specify either a reference-input file or a reference-output file, but not both.");
    }

    if (getOpt("reference-input") != "")
    {
        if (getOptBool("cyto2"))
        {
            if (!isCyto2Reference(getOpt("reference-input"))) {error("The reference-input specified either does not exist, or is not the correct type for the CNCytoEngine.");}
        }
        else
        {
            if (!isCN5Reference(getOpt("reference-input"))) {error("The reference-input specified either does not exist, or is not the correct type for the CNCytoEngine.");}
        }

        if(getOptBool("doDualNormalization"))
        {
            if (!isDualNormalizeReference(getOpt("reference-input"))) 
            {
                 error("The reference-input was generated with the doDualNormalization set to false. This is incompatible with the present run where the doDualNormalization option has been set to true. Either regenerate the reference or set the doDualNormalization option to false.");
            }
        }
        else
        {
            if (isDualNormalizeReference(getOpt("reference-input"))) 
            {
                error("The reference-input was generated with the doDualNormalization set to true. This is incompatible with the present run where the doDualNormalization option has been set to false. Either regenerate the reference or set the doDualNormalization option to true.");
            }
        }

        setOpt("create-reference", "false");
        setOpt("reference-file", getOpt("reference-input"));
    }

    if (getOpt("reference-output") != "")
    {
        // if (getOptBool("cyto2"))
        if (isOptDefined("local-gc-background-correction-reference-method")      &&
            getOpt("local-gc-background-correction-reference-method") != "none"  && 
            getOpt("local-gc-background-correction-reference-method") != "")
        {
            if ((getOpt("probe-file") == "") || (!Fs::fileExists(getOpt("probe-file"))))
            {
                error("Please specify a valid probe-file when creating a reference.");
            }
        }
        if (getOptVector("reference-cels").size() > 0)
        {
            error("Specify reference-cels in combination with the reference-input option, not the reference-output option.");
        }
        setOpt("create-reference", "true");
        setOpt("reference-file", getOpt("reference-output"));
    }

    if (getOpt("snp-reference-input-file") != "")
        {

        if (!isCyto2SnpReference(getOpt("snp-reference-input-file")))
                {
                        error("The snp-reference-input-file specified either does not exist, or is not the correct type for the CNCytoEngine.");
                }
        }

    if (getOpt("snp-reference-output-file") != "")
    {
            if (getOpt("snp-reference-input-file") == ""){
            error("When using the snp-reference-output-file option you must also use the snp-reference-input-file option.  This file must contain valid data for the information component of each probeset.");
        }
            if (getOpt("snp-reference-output-file") == getOpt("snp-reference-input-file")){
                        error("snp-reference-input-file cannot have the same name as the snp-reference-output-file");
                }
            if (getOpt("reference-output") == ""){
                        error("To use the snp-reference-output-file option you must also use the reference-output option");
                }
            if (getOpt("reference-input") != ""){
                        error("To use the snp-reference-file-output option you cannot set the reference-input option");
                }
        setOpt("create-snp-reference", "true");
        setOpt("snp-reference-output-file", getOpt("snp-reference-output-file"));
            if (getOpt("gender-override-file") == ""){
                        error("To use the snp-reference-file-output option you must specify a gender-override-file containing the gender of the samples being analyzed");
                }
            if (getOpt("genotype-call-override-file") == ""){
                        error("To use the snp-reference-file-output option you must specify a genotype-call-override-file containing the genotypes of all the SNP marker of all samples being analyzed.");
                }

    }
    if (bCyto2LOH)
    {
        if (getOpt("snp-reference-input-file") == "")
        {
            Verbose::warn(1, "A snp-reference-input-file must be set to perform an loh-cyto2 analysis");
        }
        else
        {
            if ((!Fs::fileExists(getOpt("snp-reference-input-file"))) || (!affx::File5_File::isHdf5file(getOpt("snp-reference-input-file"))) || (affx::File5_Tsv::getFileTsvLineCount(getOpt("snp-reference-input-file"), "SNPReference", "SNPReference") == -1))
            {
                error("A valid snp-reference-input-file must be set to perform an loh-cyto2 analysis");
            }
        }
    }

    /* Read in cel file list from other file if specified. */
    vector<string> celFiles;
    EngineUtil::getCelFiles(celFiles, this);
    if(celFiles.size() == 0)
        Err::errAbort("No cel files specified.");
    setOpt("cels",celFiles);

//MG temporary code for Bitao
        if(getOpt("normal-diploid-files-file")!="")
        {

                std::string normalDiploidFilesFile = getOpt("normal-diploid-files-file");
                std::vector<std::string> normalDiploidFiles;
                affx::TsvFile::extractColToVec(normalDiploidFilesFile,"normal_diploid_files",&normalDiploidFiles, affx::TSV_BIND_REQUIRED);
                setOpt("normal-diploid-files", normalDiploidFiles);
        }
//

    if (getOpt("set-analysis-name") == "") {setOpt("set-analysis-name", "CN5");}
    if (getOptBool("run-geno-qc"))
    {
        if (getOpt("reference-output") != "") {
            if (isOptDefined("geno-qc-file") && getOpt("geno-qc-file") != "") {
                setOpt("geno-qc-file", Fs::join(getOpt("out-dir"),getOpt("geno-qc-file")));
            } else {
                string refFileName = Fs::basename(getOpt("reference-output"));
                setOpt("geno-qc-file", Fs::join(getOpt("out-dir"),refFileName + "." + getOpt("set-analysis-name") + ".GenoQC.txt"));
            }
        }
        if (getOpt("reference-input") != "") {
            if (isOptDefined("geno-qc-file") && getOpt("geno-qc-file") != "") {
                setOpt("geno-qc-file", Fs::join(getOpt("out-dir"),getOpt("geno-qc-file")));
            } else {
                wstring wDateTime = DateTime::GetCurrentDateTime().ToString();
                AffxString dateTime = StringUtils::ConvertWCSToMBS(wDateTime);
                dateTime.replace(':', '-');
                setOpt("geno-qc-file", Fs::join(getOpt("out-dir"),dateTime + "." + getOpt("set-analysis-name") + ".GenoQC.txt"));
            }
        }
    }
//    if ((popts->m_iCnChpVersion < 1) || (popts->m_iCnChpVersion > 2)) {popts->errAbort("Invalid CNCHP version: " + ::getInt(popts->m_iCnChpVersion));}
    if ((getOpt("annotation-file") == "") || (!Fs::fileExists(getOpt("annotation-file"))))
    {
        error("Must specify a valid annotation-file.");
    }
    /*
    if ((getOpt("annotation-file") == "") || ((getOpt("netaffx-snp-annotation-file") == "") && (getOpt("netaffx-cn-annotation-file") == "")))
    {
        error("Must specify an annotation-file.");
    }
    */
    if (getOpt("cdf-file") == "" && getOpt("spf-file") == "") {error("Must specify a cdf-file or spf-file.");}
//    if (getOpt("special-snps") == "") {error("Must specify special-snps.");}
    if (getOpt("chrX-probes") == "") {error("Must specify chrX-probes.");}
    if (getOpt("chrY-probes") == "") {error("Must specify chrY-probes.");}

    bool bFileExists = Fs::fileExists(getOpt("reference-file"));
    
    if (bFileExists)
    {
        if (affx::File5_File::isHdf5file(getOpt("reference-file"))) {Verbose::out(1, "Input reference is in HDF5 format.");}
        if (!getOptBool("create-reference"))
        {
            if (!affx::File5_File::isHdf5file(getOpt("reference-file"))){error("Text format reference-file is not supported in apt-copynumber-Cyto.");}
            AffxString strAnnotationFileName = getAnnotationFileName(getOpt("reference-file"));
            Verbose::out(1, "CN Reference File\'s annotation-file " + strAnnotationFileName);
            if (getOptBool("cyto2"))
            {
                AffxString strAnnotationVersion = getAnnotationVersion(getOpt("reference-file"));
                Verbose::out(1, "MAJOR PROGRESS UPDATE: CN Reference annotation-version to be used is:" + strAnnotationVersion);
                AffxString strAnnotationArrayType = getAnnotationArrayType(getOpt("reference-file"));
                Verbose::out(1, "CN Reference annotation-array-type " + strAnnotationArrayType);
            }
            if (strAnnotationFileName != "")
            {
/* MG  put a warning here
                if (Fs::basename(strAnnotationFileName) != Fs::basename(getOpt("annotation-file")))
                {
                    error("Annotation file name does not match the reference file.");
                }
*/
                int iNormalizationType = -1;
                int iSampleCount = -1;
                bool bAdapterTypeNormalization = false;
                if (getNormalizationTypes(getOpt("reference-file"), iNormalizationType, bAdapterTypeNormalization, iSampleCount))
                {
                    if(iSampleCount != -1)
                    {
                        setOpt("sample-count",ToStr(iSampleCount));
                    }
                    if (bAdapterTypeNormalization != getOptBool("adapter-type-normalization"))
                    {
                        error("Selected Adapter Type Normalization does not match reference.");
                    }
                } else {error("Cannot open reference-file: " + getOpt("reference-file"));}
            }
            else
            {
                error("Cannot retreive annotation-file name from the reference file.");
            }
        }
    }
    else if (!getOptBool("create-reference")) {error("Cannot open reference-file: " + getOpt("reference-file"));}

    if (getOptBool("run-geno-qc"))
    {
        if (getOpt("qca-file") == "") {error("Please specify a qca-file.");}
        if (getOpt("qcc-file") == "") {error("Please specify a qcc-file.");}
    }

    if (getOptVector("arrs").size() > 0)
    {
        if (getOptVector("arrs").size() != getOptVector("cels").size())
        {
          error("When specifying ARR files, the number of ARR files specified must match the number of CEL files specified .");
        }
        for (unsigned int uiIndex = 0; (uiIndex < getOptVector("arrs").size()); uiIndex++)
        {
            AffxString strCelFileName = getOptVector("cels")[uiIndex];
            AffxString strArrFileName = getOptVector("arrs")[uiIndex];
            if (strArrFileName == "")
            {
                AffxString strTempFileName = strCelFileName;
                if (strTempFileName.toLowerCase().endsWith(".cel"))
                {
                    strArrFileName = strCelFileName.substring(0, (unsigned int)strCelFileName.length() - 4) + ".ARR";
                }
                else
                {
                    strArrFileName = strCelFileName + ".ARR";
                }
                if (Fs::fileExists(strArrFileName)) {getOptVector("arrs")[uiIndex] = strArrFileName;}
            }
            else
            {
                if (!Fs::fileExists(strArrFileName)) {error("Cannot find ARR file: " + strArrFileName);}
            }
            // TODO: Check to see that ARR guid matches CEL file header.
        }
    }

    if (getOptVector("result-files").size() > 0)
    {
        if (getOptVector("result-files").size() != getOptVector("cels").size())
        {
          error("When specifying CYCHP files, the number of ARR files specified must match the number of CEL files specified .");
        }
        for (unsigned int uiIndex = 0; (uiIndex < getOptVector("arrs").size()); uiIndex++)
        {
            AffxString strCychpFileName = getOptVector("result-files")[uiIndex];
            if (strCychpFileName != "")
            {
              if ( !Fs::dirExists(Fs::dirname(strCychpFileName)) ) {
                Fs::mkdirPath(Fs::dirname(strCychpFileName), false);
              }
            }
        }
    }

    checkChipType();

    if (getOptBool("create-reference"))
    {
        if (getOptVector("cels").size() < 6) {error("At least six cel files are needed to get a valid CN Reference file.");}
    }
    if ((getOptBool("cyto2")) && (!getOptBool("create-reference")))
    {
        if (getOpt("snp-qc-snp-list") == "")
        {
            error("snp-qc-snp-list must be specified when running cyto2 arrays.");
        }
        else if (!Fs::fileExists(getOpt("snp-qc-snp-list")))
        {
            error("snp-qc-snp-list must exist when specified.");
        }
    }

    if ((getOpt("local-gc-background-correction-reference-method") != "none") && (getOpt("local-gc-background-correction-reference-method") != ""))
    {
        if ((getOpt("local-gc-background-intensity-adjustment-method") == "none") || (getOpt("local-gc-background-intensity-adjustment-method") == ""))
        {
            Err::errAbort("When using the --local-gc-background-correction-reference-method, the --local-gc-background-intensity-adjustment-method must also be defined.");
        }
    }

    if ((getOpt("local-gc-background-intensity-adjustment-method") != "none") && (getOpt("local-gc-background-intensity-adjustment-method") != ""))
    {
        CNAnalysisMethodFactory amFactory;
        CNAnalysisMethod* am = NULL;
        am = amFactory.CNAnalysisMethodForString(getOpt("local-gc-background-intensity-adjustment-method"));
        delete am;
    }
    if ((getOpt("log2ratio-adjustment-method") != "none") && (getOpt("log2ratio-adjustment-method") != ""))
    {
        CNAnalysisMethodFactory amFactory;
        CNAnalysisMethod* am = NULL;
        am = amFactory.CNAnalysisMethodForString(getOpt("log2ratio-adjustment-method"));
        delete am;
    }
    if ((getOpt("image-correction-intensity-adjustment-method") != "none") && (getOpt("image-correction-intensity-adjustment-method") != ""))
    {
        CNAnalysisMethodFactory amFactory;
        CNAnalysisMethod* am = NULL;
        am = amFactory.CNAnalysisMethodForString(getOpt("image-correction-intensity-adjustment-method"));
        delete am;
    }
    if (getOptBool("create-reference"))
    {
        if ((getOpt("wave-correction-reference-method") != "none") && (getOpt("wave-correction-reference-method") != ""))
        {
            CNAnalysisMethodFactory amFactory;
            CNAnalysisMethod* am = NULL;
            am = amFactory.CNAnalysisMethodForString(getOpt("wave-correction-reference-method"));
            delete am;
        }
    }
    if (!getOptBool("create-reference"))
    {
        if ((getOpt("wave-correction-log2ratio-adjustment-method") != "none") && (getOpt("wave-correction-log2ratio-adjustment-method") != ""))
        {
            CNAnalysisMethodFactory amFactory;
            CNAnalysisMethod* am = NULL;
            am = amFactory.CNAnalysisMethodForString(getOpt("wave-correction-log2ratio-adjustment-method"));
            delete am;
        }
    }

   if(getOptBool("check-input-files"))
    {
            Verbose::out(1, "The --check-input-files default of true has been invoked.  The snp-reference-input, genotype-call-override and gender-override  files will be checked for syntax and content.  This process may take some time.  This checking may be turned of by setting the --check-input-file option to false.");
        if(getOpt("snp-reference-output-file") != "")
        {

            // By the existence of the snp-reference-output-file we know we are doing a "create snp-reference file" run.
            // Thus we check the genotype-call-override-file, gender-override, and snp-reference-input files.

            // Check the snp-reference-input file.
            if(getOpt("snp-reference-input-file") != "")
            {
                Verbose::out(1, "Checking the --snp-reference-input-file file for syntax and content.");
                if(! checkSnpReferenceInputFile(getOpt("snp-reference-input-file")))
                {
                    error("Invalid format/content found in the  --snp-reference-input-file");
                }
            } else
            {
                Verbose::warn(1, "No snp-reference-input-file found during --check-input-files.");
            }

            // Check the genotype-call-override-file.
            if (getOpt("genotype-call-override-file") == "")
            {
                Verbose::warn(1, "No genotype-call-override file found during --check-input-files.");
            }
            else
            {
                // We need to have the snp-reference-input-file to get the probesets with non-zero information content.
                // Thus we need to check for the existence of this file before beginning the check.
                if(getOpt("snp-reference-input-file") != "")
                {
                    Verbose::out(1, "Checking the --genotype-call-override-file file for syntax and content.");
                    if(! checkGenotypeCallOverrideFile( getOpt("genotype-call-override-file"),
                                                        getOptVector("cels"),
                                                        getOpt("snp-reference-input-file") )
                      )
                     {
                         error("Invalid format/content found in the  --genotype-call-override-file");
                     }
                }
                else
                {
                    Verbose::warn(1, "No snp-reference-input-file found during --check-input-files.");
                }
            }

            // Check the gender-override-file.
            if (getOpt("gender-override-file") == "")
            {
                Verbose::warn(1, "No gender-override-file found during --check-input-files.");
            }
            else
            {
                Verbose::out(1, "Checking the --gender-override-file file for syntax and content.");
                if(! checkGenderOverrideFile(getOpt("gender-override-file"), getOptVector("cels")))
                 {
                     error("Invalid format/content found in the  --gender-override-file");
                 }
            }
        }
    }
    else
    {
            Verbose::warn(1, "The --check-input-files has been disabled.   The snp-reference-input, genotype-call-override and gender-override  files will NOT be checked for syntax and content at the start of the program.  Any syntax and content errors in these files will not be detected until they are used in the analyses.  This may be several hours/days into the run.");
    }
}

bool CNCytoEngine::checkSnpReferenceInputFile(std::string strSnpReferenceFileName){

    bool bSuccessful=true;

    affx::File5_File file5;
    if(  file5.open(strSnpReferenceFileName, affx::FILE5_OPEN_RO) != affx::FILE5_OK)
    {
        bSuccessful=false;
        Verbose::warn(1, "Could not open snp-reference-input-file " + strSnpReferenceFileName);
    }
    else
    {
        affx::File5_Group* group5 = file5.openGroup("SNPReference", affx::FILE5_OPEN);
        if( group5->open("SNPReference", affx::FILE5_OPEN)!= affx::FILE5_OK)
        {
            bSuccessful=false;
            Verbose::warn(1, "Could not open group SNPReference in the snp-reference-input-file " + strSnpReferenceFileName);
        }
        else
        {
            affx::File5_Tsv* tsv5 = group5->openTsv("EMSNParameters", affx::FILE5_OPEN);
            if( tsv5->open("EMSNParameters", affx::FILE5_OPEN)!= affx::FILE5_OK)
            {
                bSuccessful=false;
                Verbose::warn(1, "Could not open group SNPReference/EMSParameters in the snp-reference-input-file " + strSnpReferenceFileName);
            }
            tsv5->close();
            delete tsv5;

            tsv5 = group5->openTsv("EMSNArray", affx::FILE5_OPEN);
            if( tsv5->open("EMSNArray", affx::FILE5_OPEN)!= affx::FILE5_OK)
            {
                bSuccessful=false;
                Verbose::warn(1, "Could not open group SNPReference/EMSNArray in the snp-reference-input-file " + strSnpReferenceFileName);
            }
            tsv5->close();
            delete tsv5;

            tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if( tsv5->open("Parameters", affx::FILE5_OPEN)!= affx::FILE5_OK)
            {
                bSuccessful=false;
                Verbose::warn(1, "Could not open group SNPReference/Parameters in the snp-reference-input-file " + strSnpReferenceFileName);
            }
            tsv5->close();
            delete tsv5;

            tsv5 = group5->openTsv("SNPReference", affx::FILE5_OPEN);
            if( tsv5->open("SNPReference", affx::FILE5_OPEN) != affx::FILE5_OK)
            {
                bSuccessful=false;
                Verbose::warn(1, "Could not open group SNPReference/SNPReference in the snp-reference-input-file " + strSnpReferenceFileName);
            }
            else
            {
                int iColumnCount = tsv5->getColumnCount(0);
                if( iColumnCount != 6){

                    bSuccessful=false;
                    Verbose::warn(1, "Missing columns in the SNPReference/SNPReference section of the file  " + strSnpReferenceFileName);
                }
                else
                {
                    bool boolFinal=true;
                    affx::File5_TsvColumn* pColumnPointer = tsv5->getColumnPtr(0, 0);

                    boolFinal = boolFinal && (pColumnPointer->getColumnName() == "probeset_id");
                    pColumnPointer = tsv5->getColumnPtr(0, 1);
                    boolFinal = boolFinal && (pColumnPointer->getColumnName() == "muAB");
                    pColumnPointer = tsv5->getColumnPtr(0, 2);
                    boolFinal = boolFinal && (pColumnPointer->getColumnName() == "muAA" );
                    pColumnPointer = tsv5->getColumnPtr(0, 3);
                    boolFinal = boolFinal && (pColumnPointer->getColumnName() == "muBB" );
                    pColumnPointer = tsv5->getColumnPtr(0, 4);
                    boolFinal = boolFinal && (pColumnPointer->getColumnName() == "information" );
                    pColumnPointer = tsv5->getColumnPtr(0, 5);
                    boolFinal = boolFinal && (pColumnPointer->getColumnName() == "useInEMAlgorithm" );

                    if(! boolFinal)
                    {
                        bSuccessful=false;
                        Verbose::warn(1, "Incorrectly named columns in the SNPReference/SNPReference section of the file  " + strSnpReferenceFileName);
                    }
                }
            }
            tsv5->close();
            delete tsv5;
        }


        group5->close();
        delete group5;
        file5.close();
   }
   return bSuccessful;
}


bool CNCytoEngine::checkGenotypeCallOverrideFile(       const std::string &strGenotypeCallOverrideFileName,
                                                        const std::vector<std::string> &vCelFiles,
                                                        const std::string  &strSnpReferenceInputFileName )
{
    bool bSuccessful=true;

    // Load the probeset names from the snp-reference file that have non-zero information content.
    std::vector<std::string> vReferenceFileProbeSetNames;
    affx::File5_File file5;
    file5.open(strSnpReferenceInputFileName, affx::FILE5_OPEN_RO);
    affx::File5_Group* group5 = file5.openGroup("SNPReference", affx::FILE5_OPEN);
    affx::File5_Tsv* tsv5 = group5->openTsv("SNPReference", affx::FILE5_OPEN);
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        double dInformation;
        tsv5->get(0, 4, &dInformation);
        if(dInformation >0.0)
        {
            std::string strProbeSetName;
            tsv5->get(0, 0, &strProbeSetName);
            vReferenceFileProbeSetNames.push_back(strProbeSetName);
         }
    }
    tsv5->close();
    delete tsv5;
    group5->close();
    delete group5;
    file5.close();
    // Done loading the probeset names.



    // Load data from genotype-call-override file. There are 3 types of data to load and check:
    //  1)  The first row must have the form:   "probeset_id"  followed by a CEL file name for
    //      each column for which we have external genotype calls. We check for the existence of
    //      the "probeset_id".
    //  2)  Then we load the CEL file names found in the remainder of the first row.
    //  3)  Load the probeset names found in the first column.


    affx::TsvFile *tsv = new affx::TsvFile;
    tsv->open(strGenotypeCallOverrideFileName);

    //  1) Do a syntax check on the first row.  i.e check for the existence of "probeset_id".
    std::string strSyntaxCheck =  tsv->getColumnName(0,0);
    if(strSyntaxCheck != "probeset_id")
    {
        bSuccessful = false;
        Verbose::warn(1, "The first row of the first column of the --genotype-call-override file must be the word:  \'probeset_id\' . The remainder of the first row must be a tab separated list of the CEL files for which external genotype information is being provided.");
    }


    //  2) Load the names of the CEL files for which we have external genotype call data.
    std::vector<std::string> vGenoCelFileNames;
    std::string celFileName;
    int iNumberOfCelFiles = tsv->getColumnCount(0);
    for(int iIndex=1; iIndex < iNumberOfCelFiles; iIndex++)
    {
        celFileName = tsv->getColumnName(0, iIndex);
        vGenoCelFileNames.push_back(celFileName);
    }


    //  3) Load the probeset names contained in the --genotype-call-override file.
    std::vector<std::string> vGenoFileProbeSetNames;
    while (tsv->nextLevel(0)==affx::TSV_OK)
    {
        std::string strProbeSetName;
        tsv->get(0, 0, strProbeSetName);
        vGenoFileProbeSetNames.push_back(strProbeSetName);
    }

    //  Determine whether the set of probesets in the snp-reference file with non-zero information
    //  all have entries in the the genotype-call-override file.
    sort(vGenoFileProbeSetNames.begin(), vGenoFileProbeSetNames.end(), StringLessThan());
    for(int iIndex = 0;  iIndex < vReferenceFileProbeSetNames.size(); iIndex++)
    {
        bool stringFound = false;
        stringFound = binary_search(    vGenoFileProbeSetNames.begin(),
                                        vGenoFileProbeSetNames.end(),
                                        vReferenceFileProbeSetNames[iIndex],
                                        StringLessThan());

        if(!stringFound)
        {

            bSuccessful=false;
            Verbose::warn(1, "No external genotype-call information was found in the genotype-call-override file for the probeset " + vReferenceFileProbeSetNames[iIndex] + ".");

        }
    }

    //  Determine whether the set of CEL files for which we have external genotype-call information
    //  contains data for all the CEL files defined in the --cel-files command.
    sort(vGenoCelFileNames.begin(), vGenoCelFileNames.end(), StringLessThan());
    for(int iIndex = 0;  iIndex < vCelFiles.size(); iIndex++)
    {
        bool stringFound = false;
        stringFound = binary_search(    vGenoCelFileNames.begin(),
                                        vGenoCelFileNames.end(),
                                        Fs::basename(vCelFiles[iIndex]),
                                        StringLessThan());

        if(!stringFound)
        {

            bSuccessful=false;
            Verbose::warn(1, "No external genotype-call information was found in the genotype-call-override file for the CEL file " + vCelFiles[iIndex] + ".");

        }
    }

    char cGenotype;
    int iGetValue;
    tsv->rewind();
    while (tsv->nextLevel(0)==affx::TSV_OK)
    {
        for(int iIndex=1; iIndex < iNumberOfCelFiles; iIndex++)
        {
            tsv->get(0, iIndex, iGetValue);
            cGenotype = (char)(iGetValue);
            if( cGenotype != -1 &&  cGenotype != 0 &&  cGenotype != 1 && cGenotype != 2)
            {
                bSuccessful=false;
                Verbose::warn(1, "An invalid genotype call of " + getInt(iGetValue)  + " was found in the --genotype-call-override file.");
            }
        }
    }
    delete tsv;
  return bSuccessful;


}



bool CNCytoEngine::checkGenderOverrideFile(     const std::string &strGenderOverrideFileName,
                                                const std::vector<std::string> &vCelFiles)
{
    bool bSuccessful=true;

    affx::TsvFile *tsv = new affx::TsvFile;
    tsv->open(strGenderOverrideFileName);

    //  1) Do a syntax check on the first row.  i.e check for the existence of the words: "experiment" and "gender".
    std::string strSyntaxCheck1 =  tsv->getColumnName(0,0);
    std::string strSyntaxCheck2 =  tsv->getColumnName(0,1);
    if(strSyntaxCheck1 != "experiment" || strSyntaxCheck2 != "gender")
    {
        bSuccessful = false;
        Verbose::warn(1, "The first row of the gender-override-file must be the words:  \'experiment\' and \' gender\', tab separated.");
    }


    //  2) Load the names of the experiments i.e. CEL files for which we have external gender data.
    std::vector<std::string> vGenderCelFileNames;
    std::string strCelFileName;
    while (tsv->nextLevel(0)==affx::TSV_OK)
    {
        tsv->get(0, 0, strCelFileName);
        vGenderCelFileNames.push_back(strCelFileName);
    }


    //  Determine whether the set of CEL files for which we have external gender information
    //  contains data for all the CEL files defined in the --cel-files command.
    sort(vGenderCelFileNames.begin(), vGenderCelFileNames.end(), StringLessThan());
    for(int iIndex = 0;  iIndex < vCelFiles.size(); iIndex++)
    {
        bool stringFound = false;
        stringFound = binary_search(    vGenderCelFileNames.begin(),
                                        vGenderCelFileNames.end(),
                                        Fs::basename(vCelFiles[iIndex]),
                                        StringLessThan());

        if(!stringFound)
        {

            bSuccessful=false;
            Verbose::warn(1, "No external gender information was found in the gender-override file for the CEL file " + vCelFiles[iIndex] + ".");

        }
    }

    std::string strGetValue;
    tsv->rewind();
    while (tsv->nextLevel(0)==affx::TSV_OK)
    {
        tsv->get(0, 1, strGetValue);
        if( strGetValue != "F" && strGetValue != "M")
        {
            bSuccessful=false;
            Verbose::warn(1, "An invalid gender call of " + strGetValue  + " was found in the --gender-override-file.");
        }
    }
    delete tsv;

    return bSuccessful;

}


/**
 * @brief Run the specified analysis
 */
void CNCytoEngine::runImp() {
  // make this all one big try catch to make sure destructors are
  // called if exceptions thrown.
  try { // outer try for handling the exception.
    try { // inner try for memory clean up.
      Fs::ensureWriteableDirPath(getOpt("out-dir"));
      Fs::ensureWriteableDirPath(getOpt("temp-dir"));
      //
      if (getOptBool("keep-intermediate-data")) {
        string strAnalysisFilePath=Fs::join(getOpt("out-dir"),"analysis");
        Fs::ensureWriteableDirPath(strAnalysisFilePath);
      }
        Verbose::out(1, "MAJOR PROGRESS UPDATE: Running Copy Number Cyto Engine.");
        Verbose::out(1, "*");
        if (getOptVector("reference-cels").size() > 0)
        {
            std::vector<std::string> vCels = getOptVector("cels");
            setOpt("cels", getOptVector("reference-cels"));

            setOpt("copynumber-reference-file", getOpt("reference-file"));
            useReference(false);

            setOpt("cels", vCels);
            useReference(true);

        }
        else
        {
            if (getOptBool("create-reference"))
            {
                createReference();
            }
            else
            {
                setOpt("copynumber-reference-file", getOpt("reference-file"));
                useReference(true);
            }
        }
        Verbose::out(1, "MAJOR PROGRESS UPDATE: Done Running Copy Number Cyto Engine.");   
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

void CNCytoEngine::createReference()
{
    CNAnalysisEngine* engine = new CNAnalysisEngine;
    engine->setEngine(this);
    engine->createAnalysis();
    delete engine; engine = NULL;

    CNLog2RatioData* pdata = new CNLog2RatioData();
    pdata->setEngine(this);
    int iExperimentCount = getOptVector("cels").size();
    pdata->getExperiments()->deleteAll();
    for (int iExperimentIndex = 0; (iExperimentIndex < iExperimentCount); iExperimentIndex++)
    {
        CNExperiment* p = new CNExperiment;
        p->setIndex(iExperimentIndex);        
        p->setExperimentName(Fs::basename(getOptVector("cels")[iExperimentIndex]));
        pdata->getExperiments()->add(p);
    }
    if (iExperimentCount == 0) {throw(Except("No experiments to process."));}
    setOpt("sample-count", ::getInt(iExperimentCount));

    pdata->loadAnnotation(getOptBool("create-reference"));	
    calculateGenders(*pdata);
    if (getOptBool("run-geno-qc")) {callGenoQCEngine();}
    removeTempProbeFiles();
    setOpt("temp-reference-file", runReferencePart1(*pdata));

//    if (getOptBool("cyto2"))
    {
        if ((getOpt("local-gc-background-correction-reference-method") != "none") && (getOpt("local-gc-background-correction-reference-method") != ""))
        {
            // pdnn takes a lot of memory so we want to clear as much as we can before running it.
            CNExperimentArray v = *pdata->getExperiments();
            pdata->getExperiments()->nullAll();
            delete pdata; pdata = NULL;

            Verbose::out(1, "*");
            CNAnalysisMethodFactory amFactory;
            CNAnalysisMethod* am = NULL;
            am = amFactory.CNAnalysisMethodForString(getOpt("local-gc-background-correction-reference-method"));
            am->setEngine(this);
            am->run();
            delete am;
            Verbose::out(1, "*");

            pdata = new CNLog2RatioData();
            pdata->setEngine(this);
            pdata->setExperiments(v);
            v.nullAll();
            pdata->loadAnnotation(getOptBool("create-reference"));
        }
    }

    if (getOpt("snp-reference-input-file") != "") {pdata->loadSnpReferenceFile(getOpt("snp-reference-input-file"));}
    if (getOpt("gender-override-file") != "") {pdata->loadGenderOverrideFile(getOpt("gender-override-file"));}
    runReferencePart2(*pdata, getOpt("temp-reference-file"));
    CNExperimentArray arExperiments;
    CNExperimentArray::copyExperiments(*pdata->getExperiments(), arExperiments);
    delete pdata; pdata = NULL;
    if ((getOpt("wave-correction-reference-method") != "none") && (getOpt("wave-correction-reference-method") != ""))
    {
        CNProbeSetArray vDummy;
        CNAnalysisMethodFactory amFactory;
        CNAnalysisMethod* am = NULL;
        am = amFactory.CNAnalysisMethodForString(getOpt("wave-correction-reference-method"));
        am->setEngine(this);
        am->setup(arExperiments, 0, vDummy);
        am->run();
        delete am;
    }
    arExperiments.deleteAll();
    if(!getOptBool("keep-temp-reference-data"))
    {
      Fs::rm(getOpt("temp-reference-file"),false);
    }
}

AffxString CNCytoEngine::runReferencePart1(CNLog2RatioData& data)
{
    CNAnalysisMethodReference reference;
    reference.setEngine(this);
    reference.setup(*data.getExperiments(), *data.getProbeSets());
    reference.runPart1();
    return reference.getTempFileName();
}

void CNCytoEngine::runReferencePart2(CNLog2RatioData& data, const AffxString& strTempFileName)
{
    CNAnalysisMethodReference reference;
    reference.setEngine(this);
    reference.setup(*data.getExperiments(), *data.getProbeSets());
    reference.setTempFileName(strTempFileName);
    reference.runPart2();
}

void CNCytoEngine::useReference(bool bAnalysis)
{
    CNLog2RatioData data;
    data.setEngine(this);
    CNAnalysisEngine engine;
    engine.setEngine(this);
    engine.createAnalysis();
    data.loadAnnotation(getOptBool("create-reference"));
    int iExperimentCount = getOptVector("cels").size();
    data.getExperiments()->deleteAll();
    loadReferenceHeader(getOpt("reference-input"));

    for (int iExperimentIndex = 0; (iExperimentIndex < iExperimentCount); iExperimentIndex++)
    {
        CNExperiment* p = new CNExperiment;
        p->setIndex(iExperimentIndex);
        p->setExperimentName(Fs::basename(getOptVector("cels")[iExperimentIndex]));
        p->setCNReferenceHeader(m_vCNReferenceHeader);
//MG temporary code for Bitao
                if(getOpt("normal-diploid-files-file")!="")
                {
                        p->setExperimentNormalDiploidFileName(getOptVector("normal-diploid-files")[iExperimentIndex]);
                }
//


        data.getExperiments()->add(p);
    }
    if (iExperimentCount == 0) {throw(Except("No experiments to process."));}	
    calculateGenders(data);

    if (bAnalysis)
    {
        if (getOptBool("run-geno-qc")) {callGenoQCEngine();}
    }
    removeTempProbeFiles();

    if ((isOptDefined("geno-qc-file")) && (getOpt("geno-qc-file") != ""))
    {
        data.loadQCMetrics(getOpt("geno-qc-file"));
    }
    if (getOptBool("cyto2"))
    {
        data.loadCyto2ModelFile(getOpt("copynumber-reference-file"));
    }
    else
    {
        data.loadModelFile(getOpt("copynumber-reference-file"));
    }
    if (getOpt("snp-reference-input-file") != "") {data.loadSnpReferenceFile(getOpt("snp-reference-input-file"));}
    int iProcessCount = data.getProbeSets()->getProcessCount();
    CNProbeSetArray arProbeSetsToProcess;
    arProbeSetsToProcess.reserve(iProcessCount);
    int iProbeSetCount = data.getProbeSetsAlternateSort()->getCount();
    for (int iIndex = 0; (iIndex < iProbeSetCount); iIndex++)
    {
        if (data.getProbeSetsAlternateSort()->getAt(iIndex)->isProcess())
        {
            arProbeSetsToProcess.add(data.getProbeSetsAlternateSort()->getAt(iIndex));
        }
    }
    arProbeSetsToProcess.quickSort(1);
    CNAnalysisMethodFactory amFactory;
    CNAnalysisMethod* pChipstream = NULL;
    pChipstream = amFactory.CNAnalysisMethodForString(getOpt("chipstream"));
    pChipstream ->setEngine(this);


    if (!bAnalysis)
    {
        data.getAAlleleEstimates()->initialize(iExperimentCount, iProbeSetCount);
        data.getBAlleleEstimates()->initialize(iExperimentCount, iProbeSetCount);
        data.getGenotypeCalls()->initialize(iExperimentCount, iProbeSetCount);
    }
    for (int iExperimentIndex = 0; (iExperimentIndex < iExperimentCount); iExperimentIndex++)
    {


        std::string strExperimentName = (*data.getExperiments()->getAt(iExperimentIndex)).getExperimentName();
        Verbose::out(1, "MAJOR PROGRESS UPDATE: Running  Analysis of Sample: " + strExperimentName);

        pChipstream->setup(*data.getExperiments(), iExperimentIndex, *data.getProbeSetsAlternateSort());
        pChipstream->run();
        if (bAnalysis)
        {
            engine.process(*data.getExperiments()->getAt(iExperimentIndex), arProbeSetsToProcess, pChipstream->getProbes());
            if ((isOptDefined("cancer")) && (getOptInt("cancer") == 2))
            {
                engine.process(*data.getExperiments()->getAt(iExperimentIndex), arProbeSetsToProcess);
            }
            if (getOptBool("log2ratio-text-output")) {writeOutputText(*data.getExperiments()->getAt(iExperimentIndex), arProbeSetsToProcess);}
            if (getOptBool("log2ratio-hdf5-output")) {writeOutputHdf5(*data.getExperiments()->getAt(iExperimentIndex), arProbeSetsToProcess);}
        }
        else
        {
            for (int iIndex = 0; (iIndex < iProbeSetCount); iIndex++)
            {
                CNProbeSet* pobjProbeSet = data.getProbeSetsAlternateSort()->getAt(iIndex);
                data.getAAlleleEstimates()->set(iExperimentIndex, iIndex, pobjProbeSet->getAAlleleSignal());
                data.getBAlleleEstimates()->set(iExperimentIndex, iIndex, pobjProbeSet->getBAlleleSignal());
                data.getGenotypeCalls()->set(iExperimentIndex, iIndex, pobjProbeSet->getGenotypeCall());
            }
        }
    }
    delete pChipstream;

    if (bAnalysis)
    {
        if ((isOptDefined("geno-qc-file")) && (getOpt("geno-qc-file") != ""))
        {
            data.dumpCyto2Report(getOpt("geno-qc-file"));
        }
        else
        {
            data.dumpCyto2Report("");
        }
    }
    else
    {
        data.calculateMedianSignal(0, iProbeSetCount);
        data.calculateXYMedianSignal(0, iProbeSetCount);
        data.calculateMedianSignalsPerGenotype(0, iProbeSetCount);

        TmpFileFactory tempFactory;
        tempFactory.setTmpdir(getOpt("temp-dir"));
        AffxString strReferenceFile = tempFactory.genFilename_basic( "CNReferenceNormals_", ".a5.tmp");
        data.writeCyto2ModelFile(strReferenceFile);
        setOpt("copynumber-reference-file", strReferenceFile);
    }
}

void CNCytoEngine:: loadReferenceHeader(std::string strReferenceFileName)
{
        try
        {
            AffxString strHeaderLine;

            affx::File5_File file5;
            file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = NULL;
            if (getOptBool("cyto2"))
            {
                group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            }
            else
            {
                group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            }
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);

            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &strHeaderLine);
                m_vCNReferenceHeader.push_back(strHeaderLine);
            }
            tsv5->close();
            delete tsv5;

            group5->close();
            delete group5;
            file5.close();
        }
        catch(...)
        {
                throw(Except("Cannot open file: " + strReferenceFileName));
        }
}




/**
* @brief Call the geno-qc engine
*/
void CNCytoEngine::callGenoQCEngine()
{
    GenoQC engine;
    engine.setOpt("force", getOpt("force"));
    engine.setOpt("exec-guid", getOpt("exec-guid"));
    engine.setOpt("command-line", getOpt("command-line"));
    engine.setOpt("program-cvs-id", getOpt("program-cvs-id"));
    engine.setOpt("program-name", getOpt("program-name"));
    engine.setOpt("program-company", getOpt("program-company"));
    engine.setOpt("program-version", getOpt("program-version"));
    engine.setOpt("version-to-report", getOpt("program-version") + " " + getOpt("program-name"));
    engine.setOpt("cdf-file", getOpt("cdf-file"));
    engine.setOpt("spf-file", getOpt("spf-file"));
    engine.setOpt("chrX-probes", getOpt("chrX-probes-filtered"));
    engine.setOpt("chrY-probes", getOpt("chrY-probes-filtered"));
    engine.setOpt("qca-file", getOpt("qca-file"));
    engine.setOpt("qcc-file", getOpt("qcc-file"));
    engine.setOpt("out-file", getOpt("geno-qc-file"));
    std::vector<std::string> cels = getOptVector("cels");
    engine.setOpt("cels", cels);
    engine.setOpt("female-thresh", getOpt("female-gender-ratio-cutoff"));
    engine.setOpt("male-thresh", getOpt("male-gender-ratio-cutoff"));
    engine.run();
}

// This method is called by ChAS without having called checkOptions.
// In short, best not to refer to Engine options/state in this method
AffxString CNCytoEngine::getAnnotationParameter(const AffxString& strFileName, const AffxString& strParameterName)
{
    AffxString strParameterValue;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        AffxString str;
        try
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5=NULL;

            if (isValidReferenceFile(strFileName))
            {
                if(isCyto2Reference(strFileName))
                {
                    group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
                }
                else
                {
                    group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
                }
            }
            else
            { 
                   Err::errAbort("An attempt was made to extract a parameter value from the CNReference file " + strFileName + ".  The file is not a valid reference file."); 
            }

            if (group5 != NULL)
            {
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if (tsv5 != NULL)
            {
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &str);
                int iIndex = str.indexOf("=");
                if (str.startsWith("#%" + strParameterName + "="))
                {
                        strParameterValue = str.substring(iIndex + 1);
                }
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

/**
 * @brief Return the annotation file names
 * @param const AffxString& - The name of the reference file
 * @param AffxString& - The name of the snp annotation file to return
 * @param AffxString& - The name of the cn annotation file to return
 * @return bool - true if successful
 */
AffxString CNCytoEngine::getAnnotationFileName(const AffxString& strFileName)
{
    AffxString strAnnotationFileName;
    AffxString strOldAnnotationFileName;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        AffxString str;
        try
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = NULL;
            if (isCyto2Reference(strFileName))
            {
                group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            }
            else
            {
                group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            }
            if (group5 != NULL)
            {
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if (tsv5 != NULL)
            {
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &str);
                int iIndex = str.indexOf("=");
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons, favoring state over option?
                if (str.startsWith("#%affymetrix-algorithm-param-option-annotation-file=")) {strAnnotationFileName = str.substring(iIndex + 1);}
                if (str.startsWith("#%affymetrix-algorithm-param-apt-opt-annotation-file=")) {strAnnotationFileName = str.substring(iIndex + 1);}
                if (str.startsWith("#%affymetrix-algorithm-param-apt-opt-netaffx-snp-annotation-file=")) {strOldAnnotationFileName = str.substring(iIndex + 1);}
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
    if ((strAnnotationFileName == "") && (strOldAnnotationFileName != ""))
    {
        strAnnotationFileName = strOldAnnotationFileName;
    }
    return strAnnotationFileName;
}

/**
 * @brief Return the normalization types
 * @param const AffxString& - The name of the reference file
 * @param int& - The normalization type
 * @param bool& - The adapter normalization
 * @param bool& - The marker normalization
 * @param int& - The number of samples
 * @return bool - true if successful
 */
bool CNCytoEngine::getNormalizationTypes(const AffxString& strFileName, int& iNormalizationType, bool& bAdapterNormalization, int& iSampleCount)
{
    bool bSuccessful = false;
    AffxByteArray ba;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        AffxString str;
        try
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = NULL;
            if (getOptBool("cyto2"))
            {
                group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            }
            else
            {
                group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            }
            if (group5 != NULL)
            {
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if (tsv5 != NULL)
            {
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &str);
                int iIndex = str.indexOf("=");

                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons, favoring state over option?
                if (str.startsWith("#%affymetrix-algorithm-param-option-normalization-type="))
                {
                    iNormalizationType = ::getInt(str.substring(iIndex + 1));
                }
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons, favoring state over option?
                if (str.startsWith("#%affymetrix-algorithm-param-state-sample-count="))
                {
                    iSampleCount = ::getInt(str.substring(iIndex + 1));
                }
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons, favoring state over option?
                if (str.startsWith("#%affymetrix-algorithm-param-option-adapter-type-normalization="))
                {
                    ba.assign(str.substring(iIndex + 1));
                    bAdapterNormalization = ba.parsebool();
                }

                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if (str.startsWith("#%affymetrix-algorithm-param-apt-opt-normalization-type="))
                {
                    iNormalizationType = ::getInt(str.substring(iIndex + 1));
                }
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if ((str.startsWith("#%affymetrix-algorithm-param-apt-opt-celFileCount=")) ||
                    (str.startsWith("#%affymetrix-algorithm-param-apt-state-sample-count=")) )
                {
                    iSampleCount = ::getInt(str.substring(iIndex + 1));
                }
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if ((str.startsWith("#%affymetrix-algorithm-param-apt-opt-adapter-normalization=")) ||
                    (str.startsWith("#%affymetrix-algorithm-param-apt-opt-adapter-type-normalization=")))
                {
                    ba.assign(str.substring(iIndex + 1));
                    bAdapterNormalization = ba.parsebool();
                }
            }
            tsv5->close();
            delete tsv5;
            }
            group5->close();
            delete group5;
            }
            file5.close();
            bSuccessful = true;
        } catch(...) {return false;}
    }
    return bSuccessful;
}

/**
 * @brief Return the reference gender counts
 * @param const AffxString& - The name of the reference file
 * @param int& - The reference male count
 * @param int& - The reference female count
 * @param int& - The reference unknown gender count
 * @return bool - true if successful
 */
bool CNCytoEngine::getReferenceGenderCounts(const AffxString& strFileName, int& iReferenceMaleCount, int& iReferenceFemaleCount, int& iReferenceUnknownCount)
{
    iReferenceMaleCount = -1;
    iReferenceFemaleCount = -1;
    iReferenceUnknownCount = -1;
    bool bSuccessful = false;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        AffxString str;
        try
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = NULL;
            if (getOptBool("cyto2"))
            {
                group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            }
            else
            {
                group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            }
            if (group5 != NULL)
            {
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if (tsv5 != NULL)
            {
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &str);
                int iIndex = str.indexOf("=");

                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons, favoring state over option?
                if (str.startsWith("#%affymetrix-algorithm-param-apt-state-reference-male-sample-count=")) {iReferenceMaleCount = ::getInt(str.substring(iIndex + 1));}
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons, favoring state over option?
                if (str.startsWith("#%affymetrix-algorithm-param-apt-state-reference-female-sample-count=")) {iReferenceFemaleCount = ::getInt(str.substring(iIndex + 1));}
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons, favoring state over option?
                if (str.startsWith("#%affymetrix-algorithm-param-apt-state-reference-unknown-sample-count=")) {iReferenceUnknownCount = ::getInt(str.substring(iIndex + 1));}


                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if ((str.startsWith("#%affymetrix-algorithm-param-apt-opt-ReferenceMaleCount=")) || (str.startsWith("#%affymetrix-algorithm-param-apt-state-reference-male-sample-count="))) {iReferenceMaleCount = ::getInt(str.substring(iIndex + 1));}
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if ((str.startsWith("#%affymetrix-algorithm-param-apt-opt-ReferenceFemaleCount=")) || (str.startsWith("#%affymetrix-algorithm-param-apt-state-reference-female-sample-count="))) {iReferenceFemaleCount = ::getInt(str.substring(iIndex + 1));}
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if ((str.startsWith("#%affymetrix-algorithm-param-apt-opt-ReferenceUnknownCount=")) || (str.startsWith("#%affymetrix-algorithm-param-apt-state-reference-unknown-sample-count="))) {iReferenceUnknownCount = ::getInt(str.substring(iIndex + 1));}
            }
            tsv5->close();
            delete tsv5;
            }
            group5->close();
            delete group5;
            }
            file5.close();
            if (iReferenceMaleCount != -1) {bSuccessful = true;}
        } catch(...) {return false;}
    }
    return bSuccessful;
}

/**
 * @brief Calculate the genders of the specified samples (cel files)
 */
void CNCytoEngine::calculateGenders(CNLog2RatioData& data)
{
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Performing a preliminary gender check.");
    Verbose::out(1, "*");
    Verbose::out(1, "female-gender-ratio-cutoff\t" + getOpt("female-gender-ratio-cutoff"));
    Verbose::out(1, "male-gender-ratio-cutoff\t" + getOpt("male-gender-ratio-cutoff"));
    Verbose::out(1, "*");
    std::vector<std::string> celFiles = getOptVector("cels");

    std::vector< std::vector<probeid_t> > chrXProbes;
    std::vector< std::vector<probeid_t> > chrYProbes;
    loadProbes(data, getOpt("chrX-probes"), getOpt("chrY-probes"), chrXProbes, chrYProbes);
    setOpt("chrX-probes-filtered", saveProbes("chrX-probes", chrXProbes));
    setOpt("chrY-probes-filtered", saveProbes("chrY-probes", chrYProbes));
    CnProbeGenderCelListener objGenderCaller(chrXProbes, chrYProbes, getOptDouble("female-gender-ratio-cutoff"), getOptDouble("male-gender-ratio-cutoff"));
    affymetrix_fusion_io::FusionCELData cel;
    int iUnknownCount = 0;
    int iMaleCount = 0;
    int iFemaleCount = 0;
    for (unsigned int uiCelFileIndex = 0; (uiCelFileIndex < celFiles.size()); uiCelFileIndex++)
    {
        cel.SetFileName(celFiles[uiCelFileIndex].c_str());
        if (!cel.Read()) {error("Cannot read cel file for gender calling: " + celFiles[uiCelFileIndex]);}

        // Check for bad cel files.
        unsigned int uiRowCount = cel.GetRows();
        unsigned int uiColCount = cel.GetCols();
        setOpt("probe-rows-count", ::getInt(uiRowCount));
        setOpt("probe-columns-count", ::getInt(uiColCount));
        unsigned int uiCount = uiRowCount * uiColCount;
        setOpt("probe-count", ::getInt(uiCount));
        float fMin = 1000000;
        float fMax = 0;
                int maxWarnings = 3;
        for (int iProbeIndex = 0; iProbeIndex < uiCount; iProbeIndex++)
        {
            float fIntensity = cel.GetIntensity(iProbeIndex);
            if (fIntensity == 0)
            {
                          if(maxWarnings > 0)
                            {
                              Verbose::out(1,"CEL file contains zero intensity values:\t" + celFiles[uiCelFileIndex]);
                              maxWarnings--;
                            }
                          continue;
            }
            else if (fIntensity < 0)
            {
                throw(Except("CEL file contains negative intensity values:\t" + celFiles[uiCelFileIndex]));
            }
            else if ((double)fIntensity != (double)fIntensity)
            {
                throw(Except("CEL file contains NaN intensity values:\t" + celFiles[uiCelFileIndex]));
            }
            else if (!Util::isFinite(fIntensity))
            {
                throw(Except("CEL file contains infinite intensity values:\t" + celFiles[uiCelFileIndex]));
            }
                     fMin = Min(fMin, fIntensity);
            fMax = Max(fMax, fIntensity);
        }
//        Verbose::out(1, ::getDouble(fMin, 6) + "\t" + ::getDouble(fMax, 6));

        objGenderCaller.newChip(&cel);
        cel.Close();
		
	data.getExperiments()->getAt(uiCelFileIndex)->setRawIntensityRatioGender(objGenderCaller.getLastGender());
	data.getExperiments()->getAt(uiCelFileIndex)->setRawIntensityRatio(objGenderCaller.getLastRatio());
        if ((isOptDefined("cyto2")) && (getOptBool("cyto2")))
        {
            data.getExperiments()->getAt(uiCelFileIndex)->setCNCallGender(objGenderCaller.getLastGender());
            data.getExperiments()->getAt(uiCelFileIndex)->setCNCallGenderRatio(objGenderCaller.getLastRatio());
        }
	Verbose::out(1, "CNCytoEngine::useReference(...)\t" + data.getExperiments()->getAt(uiCelFileIndex)->getExperimentName() + "\tRawIntensityRatio\t" + ::getDouble(data.getExperiments()->getAt(uiCelFileIndex)->getRawIntensityRatio(), 6) + "\t" + ((data.getExperiments()->getAt(uiCelFileIndex)->getRawIntensityRatioGender() == "male") ? "XY" : "XX"));

        float fRatio = objGenderCaller.getLastRatio();
        switch(objGenderCaller.getLastGender())
        {
        case affx::Male: iMaleCount++; Verbose::out(1, "Sample " + ::getUnsignedInt(uiCelFileIndex + 1) + ":\t" + Fs::basename(celFiles[uiCelFileIndex]) + "\tMale\t" + ::getDouble(fRatio, 6)); break;
        case affx::Female: iFemaleCount++; Verbose::out(1, "Sample " + ::getUnsignedInt(uiCelFileIndex + 1) + ":\t" + Fs::basename(celFiles[uiCelFileIndex]) + "\tFemale\t" + ::getDouble(fRatio, 6)); break;
        default: iUnknownCount++; Verbose::out(1, "Sample " + ::getUnsignedInt(uiCelFileIndex + 1) + ":\t" + Fs::basename(celFiles[uiCelFileIndex]) + "\tUnknown\t" + ::getDouble(fRatio, 6));
        }
    }
    Verbose::out(1, "*");
    Verbose::out(1, "MAJOR PROGRESS UPDATE: The cel file input contains " + ::getInt(iFemaleCount) + " female, " + ::getInt(iMaleCount) + " male, and " + ::getInt(iUnknownCount) + " unknown genders.");
    if (getOptBool("create-reference"))
    {
        setOpt("reference-male-sample-count", ::getInt(iMaleCount));
        setOpt("reference-female-sample-count", ::getInt(iFemaleCount));
        setOpt("reference-unknown-sample-count", ::getInt(iUnknownCount));
        int iSampleCount = iMaleCount + iFemaleCount + iUnknownCount;
        setOpt("reference-sample-count", ::getInt(iSampleCount));

        if (iSampleCount < getOptInt("recommended-reference-sample-count")) {Verbose::out(1, "WARNING: Reference will be based on less than " + getOpt("recommended-reference-sample-count") + " samples.");}
        if (iMaleCount < getOptInt("recommended-reference-male-sample-count")) {Verbose::out(1, "WARNING: Reference will be based on less than " + getOpt("recommended-reference-male-sample-count") + " male samples.");}
        if (iFemaleCount < getOptInt("recommended-reference-female-sample-count")) {Verbose::out(1, "WARNING: Reference will be based on less than " + getOpt("recommended-reference-female-sample-count") + " female samples.");}
    }
    else
    {
        int iReferenceMaleCount = 0;
        int iReferenceFemaleCount = 0;
        int iReferenceUnknownCount = 0;
        if (getReferenceGenderCounts(getOpt("reference-file"), iReferenceMaleCount, iReferenceFemaleCount, iReferenceUnknownCount))
        {
            setOpt("reference-male-sample-count", ::getInt(iReferenceMaleCount));
            setOpt("reference-female-sample-count", ::getInt(iReferenceFemaleCount));
            setOpt("reference-unknown-sample-count", ::getInt(iReferenceUnknownCount));
            int iSampleCount = iReferenceMaleCount + iReferenceFemaleCount + iReferenceUnknownCount;
            setOpt("reference-sample-count", ::getInt(iSampleCount));
            Verbose::out(1, "MAJOR PROGRESS UPDATE: The reference input is based on " + ::getInt(iReferenceFemaleCount) + " female, " + ::getInt(iReferenceMaleCount) + " male, and " + ::getInt(iReferenceUnknownCount) + " unknown genders.");
            if (iSampleCount < getOptInt("recommended-reference-sample-count")) {Verbose::out(1, "WARNING: Reference is based on less than " + getOpt("recommended-reference-sample-count") + " samples.");}
            if (iReferenceMaleCount < getOptInt("recommended-reference-male-sample-count")) {Verbose::out(1, "WARNING: Reference is based on less than " + getOpt("recommended-reference-male-sample-count") + " male samples.");}
            if (iReferenceFemaleCount < getOptInt("recommended-reference-female-sample-count")) {Verbose::out(1, "WARNING: Reference is based on less than " + getOpt("recommended-reference-female-sample-count") + " female samples.");}
            if (!getOptBool("cyto2"))
            {
                if ((iMaleCount > 0) && (iReferenceMaleCount == 0))
                {
                    error("Reference input is not based on any male samples.");
                }
                if ((iFemaleCount > 0) && (iReferenceFemaleCount == 0))
                {
                    error("Reference input is not based on any female samples.");
                }
            }
        }
    }
}

/**
 * @brief Ouptut the SNP data in ASCII text format.
 * @param int - The experiment index to process.
 * @return bool - true if successful
 */
bool CNCytoEngine::writeOutputText(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets)
{
    Verbose::out(3, "CNCytoEngine::writeOutputText");
    AffxString strFileName;
    if (getOpt("set-analysis-name") != "")
    {
      strFileName = Fs::join(getOpt("out-dir"),getOpt("set-analysis-name")+"."+objExperiment.getExperimentName()+".txt");
    }
    else
    {
      strFileName = Fs::join(getOpt("out-dir"),objExperiment.getExperimentName()+".txt");
    }
    affx::TsvFile tsv;
    tsv.addHeader("guid", affxutil::Guid::GenerateNewGuid());


    // dump options.
    std::vector<std::string> optionNames;
    getOptionNames(optionNames,1);
    for(int i=0; i< optionNames.size(); i++) {
      std::vector<std::string> vals = getOptVector(optionNames[i],1);
      if (vals.size() > 1) {
        for(int j=0; j< vals.size(); j++) {
          tsv.addHeader("affymetrix-algorithm-param-apt-" + optionNames[i] + "-" + ::getInt(j+1), vals[j]);
        }
      }
      else {
        tsv.addHeader("affymetrix-algorithm-param-apt-" + optionNames[i],  getOpt(optionNames[i],1));
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
    lineSStr << objExperiment.getExperimentName() + "\t" + ::getDouble(objExperiment.getMadDiffCN(), 6) + "\t" + ::getDouble(objExperiment.getIqr(), 6) + "\t" + ::getDouble(objExperiment.getMeanAbsRle(), 6) + "\t" + ::getDouble(objExperiment.getMedianAutosomeMedian(), 6) + "\t" + objExperiment.getCNCallGender() + "\t" + ::getDouble(objExperiment.getGCCorrectionMetric(), 6) << std::endl;
    tsv.write_str(lineSStr.str());
    lineSStr.str("");
    lineSStr << "probeset_id\tChromosome\tPosition\tLog2Ratio\tAllelicDifference\tGenotypeCall\tGenotypeConfidence\tChrXPar1\tChrXPar2" << std::endl;
    tsv.write_str(lineSStr.str());

    for (int iProbeSetIndex = 0, iCidx = 0; (iProbeSetIndex < vProbeSets.getCount()); iProbeSetIndex++, iCidx=0) {
      CNProbeSet* pobjProbeSet = vProbeSets.getAt(iProbeSetIndex);
      tsv.set(0, iCidx++, pobjProbeSet->getProbeSetName());
      tsv.set(0, iCidx++, ::getInt(pobjProbeSet->getChromosome()));
      tsv.set(0, iCidx++, ::getInt((int)pobjProbeSet->getPosition()));
      tsv.set(0, iCidx++, ::getDouble(pobjProbeSet->getLog2Ratio(), 10));
      tsv.set(0, iCidx++, ::getDouble(pobjProbeSet->getAllelicDifference(), 10));
      tsv.set(0, iCidx++, ::getInt((int)pobjProbeSet->getGenotypeCall()));
      tsv.set(0, iCidx++, ::getDouble(pobjProbeSet->getGenotypeConfidence(), 10));
      tsv.set(0, iCidx++, ::getInt((int)pobjProbeSet->isPseudoAutosomalRegion()));
      tsv.set(0, iCidx++, "0");
      tsv.set(0, iCidx++, ::getInt(pobjProbeSet->getGCBinIndex()));
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
bool CNCytoEngine::writeOutputHdf5(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets)
{
    Verbose::out(3, "CNCytoEngine::writeOutputHdf5");
    AffxString strFileName;
    if (getOpt("set-analysis-name") != "")
    {
      strFileName = Fs::join(getOpt("out-dir"),getOpt("set-analysis-name")+"."+objExperiment.getExperimentName()+".a5");
    }
    else
    {
      strFileName = Fs::join(getOpt("out-dir"),objExperiment.getExperimentName()+".a5");
    }
    try
    {
        affx::File5_File file5;
        file5.open(strFileName, affx::FILE5_REPLACE);

        affx::File5_Tsv* tsv5 = file5.openTsv("CNLog2Ratio_Parameters", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "Parameter", affx::FILE5_DTYPE_STRING, 1024);
        tsv5->set_string(0, 0, "#%guid=" + affxutil::Guid::GenerateNewGuid()); tsv5->writeLevel(0);

        // dump options.
        vector<string> optionNames;
        getOptionNames(optionNames,1);
        for(int i=0; i< optionNames.size(); i++) {
            std::vector<std::string> vals = getOptVector(optionNames[i],1);
            if (vals.size() > 1)
            {
                for(int j=0; j< vals.size(); j++)
                {
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-" + optionNames[i] + "-" + ::getInt(j+1) + "=" + vals[j]); tsv5->writeLevel(0);
                }
            }
            else
            {
                tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-" + optionNames[i] + "=" + getOpt(optionNames[i],1)); tsv5->writeLevel(0);
            }
        }

        ///@todo aw: should we dump state too?

        tsv5->close();
        delete tsv5;

        tsv5 = file5.openTsv("CNLog2Ratio_Experiments", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "experiment", affx::FILE5_DTYPE_STRING, objExperiment.getExperimentName().length());
        tsv5->defineColumn(0, 1, "MadDiffCN", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 2, "Iqr", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 3, "MeanAbsRle", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 4, "MedianAutosomeMedian", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 5, "Gender", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 6, "GCCorrectionSize", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->set_string(0, 0, objExperiment.getExperimentName());
        tsv5->set_d(0, 1, objExperiment.getMadDiffCN());
        tsv5->set_d(0, 2, objExperiment.getIqr());
        tsv5->set_d(0, 3, objExperiment.getMeanAbsRle());
        tsv5->set_d(0, 4, objExperiment.getMedianAutosomeMedian());
        tsv5->set_string(0, 5, objExperiment.getCNCallGender());
        tsv5->set_d(0, 6, objExperiment.getGCCorrectionMetric());
        tsv5->writeLevel(0);
        tsv5->close();
        delete tsv5;

        int iMaxProbeSetNameLength = 0;
        for (int iProbeSetIndex = 0; (iProbeSetIndex < vProbeSets.getCount()); iProbeSetIndex++)
        {
            CNProbeSet* pobjProbeSet = vProbeSets.getAt(iProbeSetIndex);
            iMaxProbeSetNameLength = Max(iMaxProbeSetNameLength, (int)pobjProbeSet->getProbeSetName().length());
        }
        tsv5 = file5.openTsv("CNLog2Ratio_" + objExperiment.getExperimentName(), affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "probeset_id", affx::FILE5_DTYPE_STRING, iMaxProbeSetNameLength);
        tsv5->defineColumn(0, 1, "Chromosome", affx::FILE5_DTYPE_INT, 0);
        tsv5->defineColumn(0, 2, "Position", affx::FILE5_DTYPE_INT, 0);
        tsv5->defineColumn(0, 3, "Log2Ratio", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 4, "AllelicDifference", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 5, "GenotypeCall", affx::FILE5_DTYPE_INT, 0);
        tsv5->defineColumn(0, 6, "GenotypeConfidence", affx::FILE5_DTYPE_DOUBLE, 0);
        tsv5->defineColumn(0, 7, "ChrXPar1", affx::FILE5_DTYPE_INT, 0);
        tsv5->defineColumn(0, 8, "ChrXPar2", affx::FILE5_DTYPE_INT, 0);
        for (int iProbeSetIndex = 0; (iProbeSetIndex < vProbeSets.getCount()); iProbeSetIndex++)
        {
            CNProbeSet* pobjProbeSet = vProbeSets.getAt(iProbeSetIndex);
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
bool CNCytoEngine::isDualNormalizeReference(const std::string& strReferenceFileName)
{
    int i = 0;
    int j = 0;
    try    {i = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "CopyNumber", "SketchSNP");}
    catch (...) {i = 0;}
    try    {j = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "Cyto2", "SketchSNP");}
    catch (...) {j = 0;}
    return (i > 0 || j > 0);
}

bool CNCytoEngine::isCN5Reference(const AffxString& strReferenceFileName)
{
    int i = 0;
    try    {i = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "CopyNumber", "Reference");}
    catch (...) {i = 0;}
    return (i > 0);
}

// This method is called by ChAS without having called checkOptions.
// In short, best not to refer to Engine options/state in this method
bool CNCytoEngine::isCyto2Reference(const AffxString& strReferenceFileName)
{
    int i = 0;
    try {i = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "Cyto2", "MedianSignals");}
    catch (...) {i = 0;}
    return (i > 0);
}

bool CNCytoEngine::isValidReferenceFile(const AffxString& strReferenceFileName)
{
    return isCyto2Reference(strReferenceFileName) || isCN5Reference(strReferenceFileName);
}

// Check the chip type
void CNCytoEngine::checkChipType()
{
    colrow_t numRows, numCols;
    int probeCount,  probeSetCount;
    vector<string> chipTypes;

    std::string strCdfFile = getOpt("cdf-file");
    std::string strSpfFile = getOpt("spf-file");
    if(strCdfFile != "") {
          EngineUtil::getCdfChipType(chipTypes, numRows, numCols, probeCount, probeSetCount, strCdfFile);
        }
        else if(strSpfFile != "") {
          EngineUtil::getSpfChipType(chipTypes, numRows, numCols, probeCount, probeSetCount, strSpfFile);
        }

    setOpt("num-rows", ToStr(numRows));
    setOpt("num-cols", ToStr(numCols));
    setOpt("probe-count", ToStr(probeCount));
    setOpt("probeset-count", ToStr(probeSetCount));

	if (!getOptBool("force")) {
		std::vector<std::string> cels = getOptVector("cels");
		EngineUtil::checkCelChipTypes(chipTypes, probeCount, cels, numRows, numCols);
		if (getOpt("chrX-probes") != "") {EngineUtil::checkTsvFileChipType(getOpt("chrX-probes"), chipTypes);}
		if (getOpt("chrY-probes") != "") {EngineUtil::checkTsvFileChipType(getOpt("chrY-probes"), chipTypes);}
	}
}

bool CNCytoEngine::isCyto2SnpReference(const AffxString& strReferenceFileName)
{
    int i = 0;
    try {i = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "SNPReference", "SNPReference");}
    catch (...) {i = 0;}
    return (i > 0);
    try {i = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "SNPReference", "EMSNParameters");}
    catch (...) {i = 0;}
    return (i > 0);
    try {i = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "SNPReference", "EMSNArray");}
    catch (...) {i = 0;}
    return (i > 0);
}

void CNCytoEngine::explain() {
    CNAnalysisMethodFactory factory;
    std::vector<SelfDoc>docs = factory.getDocs();

    bool found = false;
    for(unsigned int i = 0; i < docs.size(); i++)
    {
        if(docs[i].getDocName() == getOpt("explain"))
        {
            found = true;
            SelfDoc::printExplanation(docs[i],cout);
            break;
        }
    }
    if (!found) {
        cout << "'" << getOpt("explain") << "' is not a valid analysis method name" << endl;
    }
    docs.clear();
}

void CNCytoEngine::extraHelp() {
    CNAnalysisMethodFactory factory;
    printStandardMethods(cout);
    EngineUtil::printSelfDocs("Data transformations:", factory.getDocs());
}

bool CNCytoEngine::predictedIntensitiesExist(const string strInputReferenceFileName){

	bool bExists = true;
    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;
    if(  file5.open(strInputReferenceFileName, affx::FILE5_OPEN_RO) != affx::FILE5_OK)
    {
        Verbose::warn(1, "Could not open the --input-reference-file " + strInputReferenceFileName);
		return false;
    }
    else
    {

        if ((isOptDefined("cyto2")) && (getOptBool("cyto2")))
        {
            group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            if( group5->open("Cyto2", affx::FILE5_OPEN)!= affx::FILE5_OK)
            {
                Verbose::warn(1, "Could not open group Cyto2 the reference input file" + strInputReferenceFileName + ". Possibly an invalid input reference file." );
                delete group5;
				file5.close();
				return false;
            } 
        }
        else
        {
            group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
            if( group5->open("CN5", affx::FILE5_OPEN)!= affx::FILE5_OK)
            {
                Verbose::warn(1, "Could not open group CN5 in the input reference file " + strInputReferenceFileName + ". Possibly an invalid input reference file." );
                delete group5;
				file5.close();
                return false;
            } 
        }

        if ((isOptDefined("cyto2")) && (getOptBool("cyto2")))
        {
            tsv5 = group5->openTsv("ProbeEffects", affx::FILE5_OPEN);
            if( tsv5->open("ProbeEffects", affx::FILE5_OPEN)!= affx::FILE5_OK)
            {
                Verbose::warn(1, "Could not open the ProbeEffects section in the input reference file " + strInputReferenceFileName + " . Possibly an invalid inut reference file.");
				delete tsv5;
				tsv5 = NULL;
				bExists = false;
            }
        }
        else
        {
            tsv5 = group5->openTsv("CN5.plier.feature-response", affx::FILE5_OPEN);
            if( tsv5->open("CN5.plier.feature-response", affx::FILE5_OPEN)!= affx::FILE5_OK)
            {
                Verbose::warn(1, "Could not open the CN.plier.feature-response section in the input reference file " + strInputReferenceFileName + " . Possibly an invalid inut reference file.");
				delete tsv5;
				tsv5 = NULL;               
				bExists = false;
            }
        }
		if (tsv5 != NULL)
		{
			if (tsv5->getColumnCount(0) < 6) 
			{
				Verbose::warn(1, "The input reference file " + strInputReferenceFileName  + " does not contain pdnn predicted intensities.");
				 bExists = false;
			}
			tsv5->close();
			delete tsv5;
		}
		group5->close();
		delete group5;
		file5.close();
    }
    return bExists;
}

void CNCytoEngine::checkConsistency(const vector<string>& analysis) {
    // Disallow cn-state and cn-cyto2 specified together
    bool cn_state_present = false;
    bool cn_cyto2_present = false;
    const string cnState = "cn-state";
    const string cnCyto2 = "cn-cyto2";
    for (int i = 0; i < analysis.size(); i++) {
        if (analysis[i].find(cnState) != string::npos) {
            cn_state_present = true;
        }
        if (analysis[i].find(cnCyto2) != string::npos) {
            cn_cyto2_present = true;
        }
    }
    if (cn_state_present && cn_cyto2_present) {
        Err::errAbort("Analysis strings 'cn-state' and 'cn-cyto2' cannot be specified together.");
    }

    // Disallow allele-peaks when cyto2=false and
    // allelic-difference/allelic-difference-CytoScan when cyto2=true.
    string wrongStr;
    const bool cyto2flag = getOptBool("cyto2");
    if (cyto2flag) {
        wrongStr = "allelic-difference";
    } else {
        wrongStr = "allele-peaks";
    }
    for (int i = 0; i < analysis.size(); i++) {
        if (analysis[i].find(wrongStr) != string::npos) {
            error("'" + wrongStr + "' analysis cannot be run when the cyto2 input parameter is set to " + ToStr(cyto2flag) + ".");
        }
    }

    // Disallow allelic-difference and allelic-difference-CytoScan together
    const string adStr = "allelic-difference";
    if (!cyto2flag) {
        int iCount = 0;
        for (int i = 0; i < analysis.size(); i++) {
            if (analysis[i].find(adStr) != string::npos) {
                iCount++;
            }
        }
        if (iCount > 1) {
            error("'" + adStr + "' analysis cannot be specified more than once.");
        }
    }
}

void CNCytoEngine::loadProbes(CNLog2RatioData& data, const string& chrXProbeFile, const string& chrYProbeFile,
                  std::vector< std::vector<probeid_t> >& chrXProbes, 
                  std::vector< std::vector<probeid_t> >& chrYProbes)
{
    if (chrXProbeFile.empty())
        Err::errAbort("No chrX probe id file specified. Unable to filter chrXProbeFile probes.");

    if (chrYProbeFile.empty())
        Err::errAbort("No chrY probe id file specified. Unable to filter chrYProbeFile probes.");

    EngineUtil::readProbeFile(chrXProbes, chrXProbeFile);
    EngineUtil::readProbeFile(chrYProbes, chrYProbeFile);

	if (getOptBool("cytoscan-hd"))
	{
		// Filter the chrX and chrY probes to those that are members of probe sets where the processFlag > 0		
		std::map<unsigned int, CNProbeSet*> mapProbeSets;
		CNProbeSet* p = NULL;
		ChipLayout layout;		
		string strCdfFile = getOpt("cdf-file");
		string strSpfFile = getOpt("spf-file");
		if (!strCdfFile.empty()) {layout.openCdfAll(strCdfFile);}
		else if (!strSpfFile.empty()) {layout.openSpfAll(strSpfFile);}
		else {Err::errAbort("No cdf-file or spf-file specified");}
		unsigned int uiProbeID = 0;
		for (int iProbeSetIndex = 0; (iProbeSetIndex <data.getProbeSets()->getCount()); iProbeSetIndex++) 
		{
			p = data.getProbeSets()->getAt(iProbeSetIndex);
			ProbeListPacked pList = layout.getProbeListByName(p->getProbeSetName());
			if (!pList.isNull())
			{
				ProbeSet* ps = ProbeListFactory::asProbeSet(pList);
				for (unsigned int iAtomIndex = 0; (iAtomIndex < ps->atoms.size()); iAtomIndex++) 
				{
					for (unsigned int iProbeIndex = 0; (iProbeIndex < ps->atoms[iAtomIndex]->probes.size()); iProbeIndex++) 
					{
						uiProbeID = (ps->atoms[iAtomIndex]->probes[iProbeIndex]->id + 1);
						mapProbeSets.insert(std::pair <unsigned int, CNProbeSet*>(uiProbeID, p));
					}
				}
				delete ps;
			}
		}

		for (int i = 0; i < chrXProbes.size(); ++i)
		{
			vector<probeid_t>::iterator it = chrXProbes[i].begin();
			while (it != chrXProbes[i].end()) 
			{
				if (mapProbeSets.find(*it + 1) == mapProbeSets.end()) {chrXProbes[i].erase(it); continue;}
				it++;
			}
		}
		for (int i = 0; i < chrYProbes.size(); ++i)
		{
			vector<probeid_t>::iterator it = chrYProbes[i].begin();
			while (it != chrYProbes[i].end()) 
			{
				if (mapProbeSets.find(*it + 1) == mapProbeSets.end()) {chrYProbes[i].erase(it); continue;}
				it++;
			}
		}
	}
}

string CNCytoEngine::saveProbes(const string& basename, vector< vector<probeid_t> >& chrProbes)
{
    // Create temp name
    TmpFileFactory tempFactory;
    tempFactory.setTmpdir(getOpt("temp-dir"));
    AffxString strProbeFile = tempFactory.genFilename_basic(basename, ".tmp");

    affx::TsvFile *tsv = new affx::TsvFile;
    tsv->defineColumn(0,0,"probe_id", affx::TSV_TYPE_INT);

    bool channelCol = false;
    if (chrProbes.size() > 1)
    {
        tsv->defineColumn(0,1,"channel", affx::TSV_TYPE_INT);
        channelCol = true;
    }

    tsv->writeTsv(strProbeFile);

    for (int channelIdx = 0; channelIdx < chrProbes.size(); ++channelIdx)
    {
        for (int iIndex = 0; iIndex < chrProbes[channelIdx].size(); ++iIndex)
        {
            tsv->set(0,0,FROM_ZEROBASED(chrProbes[channelIdx][iIndex]));
            if (channelCol == true)
            {
                tsv->set(0,1,channelIdx);
            }
            tsv->writeLevel(0);
        }
    }
    delete tsv;
    return strProbeFile;
}

void CNCytoEngine::removeTempProbeFiles()
{
    Fs::rm(getOpt("chrX-probes-filtered"), false);
    Fs::rm(getOpt("chrY-probes-filtered"), false);
    setOpt("chrX-probes-filtered", "");
    setOpt("chrY-probes-filtered", "");
}













