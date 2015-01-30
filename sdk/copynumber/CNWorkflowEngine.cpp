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
 * @file CNWorkflowEngine.cpp
 *
 * @brief This file contains the CNWorkflowEngine class members.
 */

//
#include "copynumber/CNWorkflowEngine.h"
//
#include "copynumber/Annotation.h"
#include "copynumber/CNAnalysisEngine.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNLog2RatioEngine.h"
#include "copynumber/CNReferenceEngine.h"
//
#include "broadutil/BroadException.h"
#include "calvin_files/exception/src/ExceptionBase.h"
#include "chipstream/AdapterTypeNormTran.h"
#include "chipstream/CnProbeGenderCelListener.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/apt-probeset-genotype/ProbesetGenotypeEngine.h"
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "util/AffxBinaryFile.h"
#include "util/Fs.h"
#include "util/TmpFileFactory.h"
//
#include "../external/newmat/myexcept.h"
//
CNWorkflowEngine::Reg CNWorkflowEngine::reg;

CNWorkflowEngine * CNWorkflowEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == CNWorkflowEngine::EngineName())
        return (CNWorkflowEngine *)engine;
    return NULL;
}

/**
 * @brief Constructor
 */
CNWorkflowEngine::CNWorkflowEngine()
{
    m_apgEngine = NULL;
    m_pstdMethods = NULL;
    defineStdMethods();
    defineOptions();
}

/**
 * @brief Destructor
 */
CNWorkflowEngine::~CNWorkflowEngine()
{
    clear();
}

/**
 * @brief Initialize the data and free any memory allocated.
 */
void CNWorkflowEngine::clear()
{
    if (m_pstdMethods != NULL) {delete m_pstdMethods; m_pstdMethods = NULL;}
    if (m_apgEngine != NULL) {delete m_apgEngine; m_apgEngine = NULL;}
    // Clean up static data.
    CNExperiment::getQCMetricColumnNames()->deleteAll();
    CNAnalysisMethod::getCelFileParams()->clear();
    Annotation::getParams()->clear();
    CNAnalysisMethod::getParams()->clear();
}

/**
 * @brief Define any standard methods used by the engine
 */
void CNWorkflowEngine::defineStdMethods() {
//    m_stdMethods["brlmm"] = "quant-norm.sketch=50000,pm-only,brlmm.transform=ccs.K=4";
}

/**
 * @brief Print any standard methods used by the engine
 * @param std::ostream& - The specified output stream
 */
void CNWorkflowEngine::printStandardMethods(std::ostream &out) {
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
void CNWorkflowEngine::defineOptions() {

  defineOptionSection("Input Options");
  defineOption("", "config-file", PgOpt::STRING_OPT,
                    "The configuration file name as passed from GTC or the Cyto Browser.",
                    "");
  defineOption("", "reference-input", PgOpt::STRING_OPT, "Input reference file name.", "");
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
                     "Text file specifying cel files to process, one per line with the first line being 'cel_files'.",
                     "");
  defineOption("", "special-snps",PgOpt::STRING_OPT,
                    "File containing all snps of unusual copy (chrX,mito,Y)",
                    "");
  defineOption("", "chrX-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrX. \
                     Used for copy number probe chrX/Y ratio gender calling.",
                     "");
  defineOption("", "chrY-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrY. \
                     Used for copy number probe chrX/Y ratio gender calling.",
                     "");
  defineOption("", "target-sketch", PgOpt::STRING_OPT,
                    "File specifying a target distribution to use for quantile normalization.",
                    "");
  defineOption("", "use-feat-eff", PgOpt::STRING_OPT,
                    "File defining a feature effect for each probe. "
                     "Note that precomputed effects should only be used for an appropriately similar analysis \
                     (i.e. feature effects for pm-only may be different than for pm-mm).",
                     "");
  defineOption("", "read-models-brlmmp", PgOpt::STRING_OPT,
                     "File to read precomputed BRLMM-P snp specific models from.",
                     "");
  defineOption("", "minSegSeparation", PgOpt::INT_OPT,
                    "Value used to skip over centromere in LOH.",  "1000000000");


  defineOptionSection("Output Options");
  defineOption("", "reference-output", PgOpt::STRING_OPT, "Output reference file name.", "");

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
  defineOption("", "adapter-type-normalization", PgOpt::BOOL_OPT,
                    "Adapter Type Normalization option. true = perform adapter type normalization.",
                    "true");
  defineOption("", "normalization-type", PgOpt::INT_OPT,
                    "Normalization option. 0 = none, 1 = 'quant-norm', 2 = 'med-norm.target=1000'",
                    "1");
  defineOption("", "adapter-parameters", PgOpt::STRING_OPT,
                     "Parameters to use when running adapter type normalization.",
                     "");
  defineOption("", "brlmmp-parameters", PgOpt::STRING_OPT,
                     "Parameters to use when running brlmmp.",
                     "");
  defineOption("", "allele-peaks-reporter-method", PgOpt::STRING_OPT, "String representing allele peaks reporter pathway desired.", "allele-peaks-reporter-method");
  defineOption("", "gc-correction-bin-count", PgOpt::INT_OPT, "The number of bins to use for GC content.", "25");
  defineOption("", "allele-peaks-kernel", PgOpt::STRING_OPT, "Allele Peaks Kernel", "Gaussian Kernel");
  defineOption("", "prior-size", PgOpt::INT_OPT,
                    "How many probesets to use for determining prior.",
                    "10000");


  defineOptionSection("Misc Options");
  defineOption("", "explain", PgOpt::STRING_OPT,
                    "Explain a particular operation (i.e. --explain cn-state or --explain loh).",
                    "");

  defineOptionSection("Execution Control Options");
  defineOption("", "mem-usage", PgOpt::INT_OPT,
                    "How many MB of memory to use for this run.",
                    "0");
  defineOption("", "use-disk", PgOpt::BOOL_OPT,
               "Store CEL intensities to be analyzed on disk.", "true");
  defineOption("", "disk-cache", PgOpt::INT_OPT,
               "Size of intensity memory cache in millions of intensities (when --use-disk=true).",
               "50");

  defineOptionSection("Advanced Options");
  defineOption("", "run-geno-qc", PgOpt::BOOL_OPT,
      "Run the GenoQC engine.",
                     "false");
  defineOption("", "run-probeset-genotype", PgOpt::BOOL_OPT,
      "Run the Probeset Genotype engine. (For testing purposes only.)",
                     "true");
  defineOption("", "wave-correction-reference-method", PgOpt::STRING_OPT,
                    "String representing wave correction pathway desired.",
                    "none"); // wave-correction-reference-method

  defineOption("", "keep-intermediate-data", PgOpt::BOOL_OPT,
                    "Set to true, this option will keep all, intensity values computed "
                    "while invoking any intensity adjustment method.", "false" );

  defineOption("", "reference-chromosome", PgOpt::INT_OPT,
                    "Reference chromosome", "2");
  defineOption("", "xx-cutoff", PgOpt::DOUBLE_OPT,
                    "XX cutoff", "0.8");
  defineOption("", "xx-cutoff-high", PgOpt::DOUBLE_OPT,
                    "XX cutoff high", "1.07");
  defineOption("", "y-cutoff", PgOpt::DOUBLE_OPT, "Y cutoff", "0.65");

  defineOption("", "wave-correction-log2ratio-adjustment-method", PgOpt::STRING_OPT,
                    "String representing wave correction log2ratio adjustment pathway desired.",
                    "none"); // wave-correction-log2ratio-adjustment-method

  defineOption("", "waviness-block-size", PgOpt::INT_OPT,
                    "marker count", "50");
  defineOption("", "waviness-genomic-span", PgOpt::INT_OPT,
                    "genomic segment length", "0");
  defineOption("", "cn-calibrate-parameters", PgOpt::STRING_OPT,
                    "SmoothSignal calibration parameters", "");

  defineOptionSection("Engine Options (Not used on command line)");
  defOptMult("", "cels", PgOpt::STRING_OPT,
                     "CEL files to process.",
                     "");
  defOptMult("", "arrs", PgOpt::STRING_OPT,
                     "ARR files to process. Must be paired with cels.",
                     "");
  defOptMult("", "result-files", PgOpt::STRING_OPT,
                     "Names (with path) of CHP files to output. Must be paired with cels.",
                     "");
  defineOption("", "male-gender-ratio-cutoff", PgOpt::DOUBLE_OPT,
                    "Male gender ratio cutoff", "0.71");
  defineOption("", "female-gender-ratio-cutoff", PgOpt::DOUBLE_OPT,
                    "Female gender ratio cutoff", "0.48");


  CNReferenceEngine objCNReferenceEngine;
  objCNReferenceEngine.setEngine(this);
  objCNReferenceEngine.defineSharedOptions();

  CNLog2RatioEngine objCNLog2RatioEngine;
  objCNLog2RatioEngine.setEngine(this);
  objCNLog2RatioEngine.defineSharedOptions();

  // Override default option from sub engine
  setOpt("log2-input", "false");
}

/**
 * @brief Add a cel file specification to the options
 * @param const std::string& - The cel file name
 * @param const std::string& - The sample (arr) file name
 * @param const std::string& - The result (cychp) file name
 */
void CNWorkflowEngine::addCel(const std::string& strCelFileName, const std::string& strSampleFileName, const std::string& strResultFileName)
{
    pushOpt("cels", strCelFileName);
    pushOpt("arrs", strSampleFileName);
    pushOpt("result-files", strResultFileName);
}

/**
 * @brief Define the states used by this engine
 */
void CNWorkflowEngine::defineStates() {

  defineOption("", "cytoscan-hd", PgOpt::BOOL_OPT, "The array is the CytoScanHD product.", "false");
  defineOption("", "reference-wave-count", PgOpt::INT_OPT, "wave-count-used", "-1");
  defineOption("", "gender-method", PgOpt::STRING_OPT,
                     "The method used for gender calling.",
                     "cn-probe-chrXY-ratio");
  defineOption("", "genotype-analysis", PgOpt::STRING_OPT,
                     "The analysis string used for probeset-genotype.",
                     "");
  defineOption("", "qmethod-spec", PgOpt::STRING_OPT,
                     "The analysis string used to calculate probe set signal values.",
                     "");

  defineOption("", "recommended-reference-sample-count", PgOpt::INT_OPT,
                    "The recommended number of samples for a Copy Number reference.", "44");
  defineOption("", "recommended-reference-female-sample-count", PgOpt::INT_OPT,
                    "The recommended number of female samples for a Copy Number reference.", "15");
  defineOption("", "recommended-reference-male-sample-count", PgOpt::INT_OPT,
                    "The recommended number of male samples for a Copy Number reference.", "15");
  defineOption("", "reference-sample-count", PgOpt::INT_OPT,
                    "The actual number of samples for the Copy Number reference.", "0");
  defineOption("", "reference-female-sample-count", PgOpt::INT_OPT,
                    "The actual number of female samples for the Copy Number reference.", "0");
  defineOption("", "reference-male-sample-count", PgOpt::INT_OPT,
                    "The actual number of male samples for the Copy Number reference.", "0");
  defineOption("", "reference-unknown-sample-count", PgOpt::INT_OPT,
                    "The actual number of unknown gender samples for the Copy Number reference.", "0");
  defineOption("", "yTarget", PgOpt::DOUBLE_OPT, "Store yTarget for sharing between modules.", "-0.56747");
  defineOption("", "temp-reference-file", PgOpt::STRING_OPT, "Temporary reference file name.", "");
  defineOption("", "reference-file", PgOpt::STRING_OPT, "Reference file name.", "");
  defineOption("", "num-rows", PgOpt::INT_OPT, "The number of rows on the chip.", "-1");
  defineOption("", "num-cols", PgOpt::INT_OPT, "The number of cols on the chip.", "-1");
  defineOption("", "probe-count", PgOpt::INT_OPT, "The number of probes on the chip.", "-1");
  defineOption("", "probeset-count", PgOpt::INT_OPT, "The number of probesets on the chip.", "-1");
  defineOption("", "channel-count", PgOpt::INT_OPT, "The number of data channels in the CEL file.", "-1");
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

  CNReferenceEngine objCNReferenceEngine;
  objCNReferenceEngine.setEngine(this);
  objCNReferenceEngine.defineSharedStates();

  CNLog2RatioEngine objCNLog2RatioEngine;
  objCNLog2RatioEngine.setEngine(this);
  objCNLog2RatioEngine.defineSharedStates();
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
            Err::errAbort("CNWorkflowEngine/Analysis: ComparePred instantiated with an invalid compare code = " + ToStr(k));
            return false;
        }
    };
};

template<> struct Analysis::ComparePred<0> {
    bool operator()(const Analysis* lhs, const Analysis* rhs) const {
        return AffxArray<int>::compare(lhs->Order, rhs->Order) < 0;
    }
};

/**
 * @brief Make sure that our options are sane. Call Err::errAbort if not.
 */
void CNWorkflowEngine::checkOptionsImp() {

  defineStates();

  // not used by the engine
  //setLibFileOpt("config-file");
  setLibFileOpt("reference-input");
  setLibFileOpt("cdf-file");
  setLibFileOpt("spf-file");
  setLibFileOpt("qcc-file");
  setLibFileOpt("qca-file");
  setLibFileOpt("special-snps");
  setLibFileOpt("chrX-probes");
  setLibFileOpt("chrY-probes");
  setLibFileOpt("target-sketch");
  setLibFileOpt("use-feat-eff");
  setLibFileOpt("read-models-brlmmp");

  if(getOpt("explain") != "") { explain(); exit(0); }

  if (getOpt("out-dir") == "") {error("Must specify an output directory.");}
  if (getOpt("temp-dir") == "") {
    setOpt("temp-dir", Fs::join(getOpt("out-dir"),"temp"));
    Fs::ensureWriteableDirPath(getOpt("temp-dir"), false);
  }

  std::vector<std::string>& analysis = getOptVector("analysis");

    if(analysis.size() == 0) {
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
//        analysis.push_back("mosaicism");
    }

    // Insure proper analysis order.
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
    }
    ar.deleteAll();

        if (getOptBool("keep-intermediate-data") )
        {
        error("Cannot use the --keep-intermediate-data option in apt-copynumber-workflow");
        }


    if ((getOpt("reference-input") != "") && (getOpt("reference-output") != ""))
    {
        error("Please specify either a reference-input file or a reference-output file, but not both.");
    }
    if ((getOpt("reference-input") == "") && (getOpt("reference-output") == ""))
    {
        error("Please specify either a reference-input file or a reference-output file.");
    }
    if (getOpt("reference-input") != "")
    {
        if (!isCN5Reference(getOpt("reference-input"))) {error("The reference-input specified either does not exist, or is not the correct type for the CNWorkflowEngine.");}
        setOpt("create-reference", "false");
        setOpt("reference-file", getOpt("reference-input"));
    }
    if (getOpt("reference-output") != "")
    {
        setOpt("create-reference", "true");
        setOpt("reference-file", getOpt("reference-output"));
        vector<string> celFiles;
        EngineUtil::getCelFiles(celFiles, this);
        if(celFiles.size() < 5)
        {
             Err::errAbort("At least 5 cels files must be specified when creating a reference file.");
        }

    }

  if(analysis.empty())
     error("Must specify at least one analysis to perform.");

    /* Read in cel file list from other file if specified. */
    vector<string> celFiles;
    EngineUtil::getCelFiles(celFiles, this);
    if(celFiles.size() == 0)
        Err::errAbort("No cel files specified.");
    setOpt("cels",celFiles);
    int celfile_count = celFiles.size();

    if (getOpt("brlmmp-parameters") == "")
    {
        // brlmm-p
//        if (m_optsCn.m_bCreateReference) {m_optsCn.m_strBrlmmpParameters = "CM=1.mix=1.bins=100.K=2.SB=0.003.MS=0.05";}
//        else {m_optsCn.m_strBrlmmpParameters = "CM=2.mix=1.bins=100.K=2.SB=0.003.MS=0.05";}
        // brlmm-p-plus-force
        if (getOptBool("create-reference")) {setOpt("brlmmp-parameters", "copynumber=true.CM=1.bins=100.mix=1.bic=2.HARD=3.SB=0.45.KX=1.KH=1.5.KXX=0.5.KAH=-0.6.KHB=-0.6.transform=MVA.AAM=2.0.BBM=-2.0.AAV=0.06.BBV=0.06.ABV=0.06.copyqc=0.wobble=0.05.MS=1");}
        else {setOpt("brlmmp-parameters", "CM=2.bins=100.mix=1.bic=2.HARD=3.SB=0.45.KX=1.KH=1.5.KXX=0.5.KAH=-0.6.KHB=-0.6.transform=MVA.AAM=2.0.BBM=-2.0.AAV=0.06.BBV=0.06.ABV=0.06.copyqc=0.wobble=0.05.MS=1");}
    }
    if ((getOptBool("create-reference")) && (getOpt("brlmmp-parameters").find(".MS=1") == std::string::npos))
    {
        Verbose::out(1, "WARNING: The Score Threshold is recommended to be set to 1 when creating a reference model file.");
    }

    if (getOpt("set-analysis-name") == "") {setOpt("set-analysis-name", "CN5");}
    if (getOptBool("run-geno-qc"))
    {
        setOpt("geno-qc-file", Fs::join(getOpt("out-dir"), getOpt("set-analysis-name") + ".GenoQC.txt"));
    }
//    if ((popts->m_iCnChpVersion < 1) || (popts->m_iCnChpVersion > 2)) {popts->errAbort("Invalid CNCHP version: " + ::getInt(popts->m_iCnChpVersion));}
    if (getOpt("annotation-file") == "") {error("Must specify an annotation-file.");}
    if (getOpt("cdf-file") == "" && getOpt("spf-file") == "") {error("Must specify a cdf-file or spf-file.");}
    if (getOpt("out-dir") == "") {error("Must specify an output directory.");}
    if (getOpt("reference-file") == "") {error("Must specify a reference-file.");}
    if ( !Fs::dirExists(Fs::dirname(getOpt("reference-file"))) ) {
      Fs::mkdirPath(Fs::dirname(getOpt("reference-file")), false);
    }
    if (getOpt("special-snps") == "") {error("Must specify special-snps.");}
    if (getOpt("chrX-probes") == "") {error("Must specify chrX-probes.");}
    if (getOpt("chrY-probes") == "") {error("Must specify chrY-probes.");}
    if (!getOptBool("create-reference"))
    {
        if (getOpt("target-sketch") != "") {error("Must not specify target-sketch. It will be read from the Reference.");}
        if (getOpt("use-feat-eff") != "") {error("Must not specify use-feat-eff. It will be read from the Reference.");}
        if (getOpt("read-models-brlmmp") != "") {error("Must not specify read-model-brlmmp. It will be read from the Reference.");}
    }
    if (celfile_count == 0)
    {
        error("Must specify at least one cel-file.");
    }

    if (Fs::isReadable(getOpt("reference-file"))) {
            if (affx::File5_File::isHdf5file(getOpt("reference-file"))) {Verbose::out(1, "Input reference is in HDF5 format.");}
            if (!getOptBool("create-reference"))
              {
                if (!affx::File5_File::isHdf5file(getOpt("reference-file"))){error("Text format reference-file is not supported in apt-copynumber-workflow.");}
                AffxString strAnnotationFileName;
                if (getAnnotationFileName(getOpt("reference-file"), strAnnotationFileName))
                  {
                    if (strAnnotationFileName.endsWith(".csv")) {strAnnotationFileName = strAnnotationFileName.substring(0, strAnnotationFileName.length() - 4) + ".db";}
                    if (Fs::basename(strAnnotationFileName) != Fs::basename(getOpt("annotation-file")))
                      {
                        error("Annotation file name does not match the reference file.");
                      }
                    int iNormalizationType = -1;
                    int iSampleCount = -1;
                    bool bAdapterTypeNormalization = false;
                    if (getNormalizationTypes(getOpt("reference-file"), iNormalizationType, bAdapterTypeNormalization, iSampleCount))
                      {
                        if(iSampleCount != -1)
                          {
                            setOpt("sample-count",ToStr(iSampleCount));
                          }
                        if (iNormalizationType != getOptInt("normalization-type"))
                          {
                            error("Selected Normalization Type does not match reference.");
                          }
                        if (bAdapterTypeNormalization != getOptBool("adapter-type-normalization"))
                          {
                            error("Selected Adapter Type Normalization does not match reference.");
                          }
                      } else {error("Cannot open reference-file: " + getOpt("reference-file"));}
                  }
                else
                  {
                    error("Cannot retreive the annotation file name from the reference file.");
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
        if (getOptVector("arrs").size() != celfile_count)
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
        if (getOptVector("result-files").size() != celfile_count)
        {
          error("When specifying CHP files, the number of chp file names specified with result-files must match the number of CEL files specified .");
        }
        for (unsigned int uiIndex = 0; (uiIndex < getOptVector("arrs").size()); uiIndex++)
        {
            AffxString strChpFileName = getOptVector("result-files")[uiIndex];
            if (strChpFileName != "")
            {
              if ( !Fs::dirExists(Fs::dirname(strChpFileName)) ) {
                Fs::mkdirPath(Fs::dirname(strChpFileName), false);
              }
            }
        }
    }

    CNReferenceEngine objCNReferenceEngine;
    objCNReferenceEngine.setEngine(this);
    objCNReferenceEngine.checkOptions();
    objCNReferenceEngine.clear();

    CNLog2RatioEngine objCNLog2RatioEngine;
    objCNLog2RatioEngine.setEngine(this);
    objCNLog2RatioEngine.checkOptions();
    objCNLog2RatioEngine.clear();

    if ((!getOptBool("cychp-output")) && (!getOptBool("cnchp-output")))
    {
        throw(Except("Must specify either cychp-ouptut or cnchp-output."));
    }
    checkChipType();

  Fs::ensureWriteableDirPath(getOpt("out-dir"), false);
  Fs::ensureWriteableDirPath(getOpt("temp-dir"), false);

    if ((getOpt("wave-correction-reference-method") != "none") && (getOpt("wave-correction-reference-method") != ""))
    {
        CNAnalysisMethodFactory amFactory;
        CNAnalysisMethod* am = NULL;
        am = amFactory.CNAnalysisMethodForString(getOpt("wave-correction-reference-method"));
        delete am;
    }
    if ((getOpt("wave-correction-log2ratio-adjustment-method") != "none") && (getOpt("wave-correction-log2ratio-adjustment-method") != ""))
    {
        CNAnalysisMethodFactory amFactory;
        CNAnalysisMethod* am = NULL;
        am = amFactory.CNAnalysisMethodForString(getOpt("wave-correction-log2ratio-adjustment-method"));
        delete am;
    }
}
/**
 * @brief checks if there is enough disk space for diskmarts in the
 * temp-dir and for chp files in out-dir
 */
void CNWorkflowEngine::checkDiskSpaceImp() {
  int64_t out_disk_available = -1, scratch_disk_available = -1;

  int cel_count = getOptVector("cels").size();
  int probeset_count = getOptInt("probeset-count");

  std::vector<std::string> result_files = getOptVector("result-files");
  std::string out_dir;
  if (!result_files.empty()) {
    //int dir_end = result_files[0].rfind(PATH_SEPARATOR);
    //if (dir_end != result_files[0].npos) {
    //  out_dir = string(result_files[0],0,dir_end);
    //}
    //else {
    //  out_dir = ".";
    //}
    out_dir=Fs::dirname(result_files[0]);
    // @todo iterate through all of the result dirs to check available
    // disk space
  }

  if (out_dir.empty()) {
    out_dir = getOpt("out-dir");
  }

  out_disk_available = Fs::getFreeDiskSpace(out_dir, false);

  uint64_t out_disk_needed = 0;
  if (getOptBool("cnchp-output")) {
    out_disk_needed += static_cast<uint64_t>(probeset_count) * cel_count * 42.08;
    if (getOptBool("text-output")) {
      out_disk_needed += static_cast<uint64_t>(probeset_count) * cel_count * 52.64;
    }
  }

  if (getOptBool("cychp-output")) {
    out_disk_needed += static_cast<uint64_t>(probeset_count) * cel_count * 42.08;
    if (getOptBool("text-output")) {
      out_disk_needed += static_cast<uint64_t>(probeset_count) * cel_count * 66.04;
    }
  }

  uint64_t scratch_disk_needed = static_cast<uint64_t>(cel_count) * probeset_count * 16.04;

  std::string temp_dir = getOpt("temp-dir");
  AptErr_t aptErr = APT_OK;
  bool same_disk = Fs::isSameVolume(temp_dir, out_dir, aptErr, false);
  if (temp_dir.empty() ||
      same_disk||
      ((aptErr != APT_OK ) && (temp_dir == out_dir))
      ) {
    out_disk_needed += scratch_disk_needed;
  }
  else {
    scratch_disk_available = Fs::getFreeDiskSpace(temp_dir, false);
    if (scratch_disk_needed > 0 &&
        scratch_disk_available >= 0 &&
        scratch_disk_needed >= scratch_disk_available) {
      // format in kb
      scratch_disk_needed = (scratch_disk_needed + 512) / 1024;
      scratch_disk_available = (scratch_disk_available + 512) / 1024;
      Err::errAbort("In " + temp_dir + ", need " + ToStr(scratch_disk_needed) + "kb of scratch diskspace.  Only " + ToStr(scratch_disk_available) + "kb available.");
    }
  }

  if (out_disk_needed > 0 &&
      out_disk_available >= 0 &&
      out_disk_needed >= out_disk_available) {
    // format in kb
    out_disk_needed = (out_disk_needed + 512) / 1024;
    out_disk_available = (out_disk_available + 512) / 1024;
    Err::errAbort("In " + out_dir + ", need " + ToStr(out_disk_needed) + "kb of diskspace.  Only " + ToStr(out_disk_available) + "kb available.");
  }
}

/**
 * @brief Run the specified analysis
 */
void CNWorkflowEngine::runImp() {
    // make this all one big try catch to make sure destructors are
  // called if exceptions thrown.
  try { // outer try for handling the exception.
    try { // inner try for memory clean up.
        Verbose::out(1, "*");
        CNAnalysisMethodFactory objFactory;
        std::vector<std::string> analysis = getOptVector("analysis");
        for (int i = 0; (i < (int)analysis.size()); i++)
        {
            Verbose::out(1, analysis[i]);
            // Verify analysis method parameters.
            CNAnalysisMethod* p = objFactory.CNAnalysisMethodForString(analysis[i]);
            if (p != NULL) {delete p;}
        }
        Verbose::out(1, "*");
        calculateGenders();
        Verbose::out(1, "*");
        if (m_apgEngine == NULL) {m_apgEngine = new ProbesetGenotypeEngine;}
        if (getOptBool("create-reference"))
        {
            // Reference Creation Workflow
            setupReferenceRunParameters();
            Verbose::out(1, m_apgEngine->getOpt("analysis"));
            Verbose::out(1, m_apgEngine->getOpt("qmethod-spec"));
            Verbose::out(1, "*");
            if (getOptBool("run-geno-qc")) {callGenoQCEngine();}
            if (getOptBool("run-probeset-genotype")) {callProbesetGenotypeEngine();}
            callCNReferenceEngine();
            callCNLog2RatioEngine();
        }
        else
        {
            // CopyNumber Workflow
            setupCopyNumberRunParameters();
            Verbose::out(1, m_apgEngine->getOpt("analysis"));
            Verbose::out(1, m_apgEngine->getOpt("qmethod-spec"));
            Verbose::out(1, "*");
            if (getOptBool("run-geno-qc")) {callGenoQCEngine();}
            if (getOptBool("run-probeset-genotype")) {callProbesetGenotypeEngine();}
            callCNLog2RatioEngine();
        }
    } // inner try end
    catch (...) {
      // Free up any allocated memory here...
      clear();
      removeTempDir(getOpt("temp-dir"));
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
  removeTempDir(getOpt("temp-dir"));
}

/**
 * @brief Setup the parameters for the apg runs.
 */
void CNWorkflowEngine::setupParameters()
{
    if (getOpt("temp-dir") == "") {
    setOpt("temp-dir",Fs::join(getOpt("out-dir"),"temp"));
  }
    m_apgEngine->setOpt("male-thresh", getOpt("male-gender-ratio-cutoff"));
    m_apgEngine->setOpt("female-thresh", getOpt("female-gender-ratio-cutoff"));
    m_apgEngine->setOpt("use-disk", ToStr(getOptBool("use-disk")));
    m_apgEngine->setOpt("disk-cache", ToStr(getOptInt("disk-cache")));
    m_apgEngine->setOpt("temp-dir", getOpt("temp-dir"));
    m_apgEngine->setOpt("prior-size", ToStr(getOptInt("prior-size")));
    m_apgEngine->setOpt("verbose",  ToStr(getOptInt("verbose")));
    m_apgEngine->setOpt("exec-guid", getOpt("exec-guid"));
    m_apgEngine->setOpt("command-line", getOpt("command-line"));
    m_apgEngine->setOpt("program-cvs-id", getOpt("program-cvs-id"));
    m_apgEngine->setOpt("program-name", getOpt("program-name"));
    m_apgEngine->setOpt("program-company", getOpt("program-company"));
    m_apgEngine->setOpt("program-version", getOpt("program-version"));
    m_apgEngine->setOpt("version-to-report", getOpt("program-version") + " " + getOpt("program-name"));
    m_apgEngine->setOpt("probeset-ids", getOpt("probeset-ids"));
    m_apgEngine->setOpt("set-analysis-name", getOpt("set-analysis-name"));
    m_apgEngine->setOpt("cdf-file", getOpt("cdf-file"));
    m_apgEngine->setOpt("spf-file", getOpt("spf-file"));
    m_apgEngine->setOpt("set-gender-method", getOpt("gender-method"));
    m_apgEngine->setOpt("chrX-probes", getOpt("chrX-probes"));
    m_apgEngine->setOpt("chrY-probes", getOpt("chrY-probes"));
    m_apgEngine->setOpt("special-snps", getOpt("special-snps"));
    m_apgEngine->setOpt("annotation-file", getOpt("annotation-file"));

    AffxString strAnalysisString;
    if (getOptBool("adapter-type-normalization"))
    {
        if (getOpt("adapter-parameters") == "") {setOpt("adapter-parameters", AdapterTypeNormTran::getDefaultOptions());}
        strAnalysisString = "adapter-type-norm";
        if (getOpt("adapter-parameters") != "")
        {
            strAnalysisString += ".";
            strAnalysisString += getOpt("adapter-parameters");
        }
        strAnalysisString += ",";
    }
    if (getOptInt("cancer") > 0)
    {
        strAnalysisString += "intensity-reporter.file='" + Fs::join(getOpt("out-dir"),getOpt("set-analysis-name")) + ".Intensities.a5',";
    }
    switch (getOptInt("normalization-type"))
    {
    case 0: break;
    case 1: strAnalysisString = strAnalysisString + "quant-norm.bioc=false,"; break;
    case 2: strAnalysisString = strAnalysisString + "med-norm.target=1000,"; break;
    default: throw(Except("Normalization type must be one of 1 = quant-norm, or 2 = med-norm.target=1000."));
    }
    strAnalysisString += "pm-only,brlmm-p." + getOpt("brlmmp-parameters");
    m_apgEngine->setOpt("analysis",strAnalysisString);
    if (getOptBool("create-reference"))
    {
        m_apgEngine->setOpt("qmethod-spec", "plier.optmethod=1.FixFeatureEffect=true");
        setOpt("log2-input", "false");
    }
    else
    {
        m_apgEngine->setOpt("qmethod-spec", "plier.optmethod=1.FixFeatureEffect=true.SafetyZero=0");
        setOpt("log2-input", "false");
    }
    setOpt("genotype-analysis", m_apgEngine->getOpt("analysis"));
    setOpt("qmethod-spec", m_apgEngine->getOpt("qmethod-spec"));

    m_apgEngine->setOpt("all-types", "true" );
    m_apgEngine->setOpt("force", getOpt("force"));
    m_apgEngine->setOpt("table-output", "false" );
    m_apgEngine->setOpt("summaries", "false" );
    m_apgEngine->setOpt("xda-chp-output", "false" );
    m_apgEngine->setOpt("cc-chp-output", "false" );
    m_apgEngine->setOpt("feat-details", "false" );
    m_apgEngine->setOpt("set-analysis-name", getOpt("set-analysis-name") );
    m_apgEngine->setOpt("feat-details", "false" );
    m_apgEngine->setOpt("include-quant-in-report-file-name", "true" );
}

/**
 * @brief Setup the parameters (specific to the reference run) for the aps and apg runs.
 */
void CNWorkflowEngine::setupReferenceRunParameters()
{
    setupParameters();
    m_apgEngine->setOpt("use-feat-eff", getOpt("use-feat-eff"));
    m_apgEngine->setOpt("target-sketch", getOpt("target-sketch"));
    m_apgEngine->setOpt("read-models-brlmmp", getOpt("read-models-brlmmp"));
    m_apgEngine->setOpt("a5-global-file", getOpt("reference-file"));
    m_apgEngine->setOpt("a5-calls", "true" );
    m_apgEngine->setOpt("a5-calls-use-global", "false" );
    m_apgEngine->setOpt("a5-feature-effects", "true" );
    m_apgEngine->setOpt("a5-feature-effects-use-global", "true" );
    m_apgEngine->setOpt("a5-feature-details", "false" );
    m_apgEngine->setOpt("a5-feature-details-use-global", "false" );
    m_apgEngine->setOpt("a5-sketch", "true" );
    m_apgEngine->setOpt("a5-sketch-use-global", "true" );
    m_apgEngine->setOpt("a5-summaries", "true" );
    m_apgEngine->setOpt("a5-summaries-use-global", "false" );
    m_apgEngine->setOpt("a5-write-models", "true" );
    m_apgEngine->setOpt("a5-write-models-use-global", "true" );
    m_apgEngine->setOpt("feat-effects", "false" );
    m_apgEngine->setOpt("write-sketch", "false" );
    m_apgEngine->setOpt("write-models", "false" );
    m_apgEngine->setOpt("out-dir", getOpt("out-dir"));
  m_apgEngine->setOpt("writeOldStyleFeatureEffectsFile", "true");

  string outDir = m_apgEngine->getOpt("out-dir");
  string setAnalysisName = m_apgEngine->getOpt("set-analysis-name");
  if ( !Fs::dirExists(outDir)) {
    Fs::mkdirPath(outDir, false);
  }

  setOpt("genotype-report-file", Fs::join(outDir, setAnalysisName + ".report.txt"));
  setOpt("expr-summary-file",    Fs::join(outDir, setAnalysisName + ".plier.summary.a5"));
    setOpt("genotype-calls-file",  Fs::join(outDir, setAnalysisName + ".calls.a5"));
    setOpt("genotype-confidences-file", Fs::join(outDir,setAnalysisName + ".confidences.a5"));
}

/**
 * @brief Setup the parameters (specific to the copynumber run) for the aps and apg runs.
 */
void CNWorkflowEngine::setupCopyNumberRunParameters()
{
    setupParameters();
    m_apgEngine->setOpt("a5-global-file", "" );
    m_apgEngine->setOpt("a5-calls", "true" );
    m_apgEngine->setOpt("a5-calls-use-global", "false" );
    m_apgEngine->setOpt("a5-summaries", "true" );
    m_apgEngine->setOpt("a5-summaries-use-global", "false" );
    m_apgEngine->setOpt("a5-feature-effects", "true" );
    m_apgEngine->setOpt("a5-feature-effects-use-global", "false" );
    m_apgEngine->setOpt("a5-feature-details", "false" );
    m_apgEngine->setOpt("a5-feature-details-use-global", "false" );
    if (getOptInt("normalization-type") == 1)
    {
        m_apgEngine->setOpt("a5-sketch", "true" );
    }
    else
    {
        m_apgEngine->setOpt("a5-sketch", "false" );
    }
    m_apgEngine->setOpt("a5-sketch-use-global", "false" );
    m_apgEngine->setOpt("a5-write-models", "true" );
    m_apgEngine->setOpt("a5-write-models-use-global", "false" );
    m_apgEngine->setOpt("feat-effects", "false" );
    m_apgEngine->setOpt("write-sketch", "false" );
    m_apgEngine->setOpt("write-models", "false" );
    m_apgEngine->setOpt("out-dir", getOpt("out-dir"));

  string outDir = m_apgEngine->getOpt("out-dir");
  string setAnalysisName = m_apgEngine->getOpt("set-analysis-name");

  setOpt("genotype-report-file", Fs::join(outDir, setAnalysisName + ".report.txt"));
    setOpt("expr-summary-file",    Fs::join(outDir, setAnalysisName + ".plier.summary.a5"));
    setOpt("genotype-calls-file",  Fs::join(outDir, setAnalysisName + ".calls.a5"));
    setOpt("genotype-confidences-file", Fs::join(outDir,setAnalysisName + ".confidences.a5"));

  m_apgEngine->setOpt("a5-global-input-file", getOpt("reference-file"));
  // grab "strGlobalPrefix" from the reference file (Cause it could be different.)
  // the report names below are relative to "m_a5_global_input_groupname"
    AffxString strGlobalPrefix;
    if (!getAnalysisName(getOpt("reference-file"), strGlobalPrefix)) {
        throw(Except("Cannot find global prefix in reference-file."));
    }
    else {
        m_apgEngine->setOpt("a5-input-group", strGlobalPrefix );
    }

    if (getOptInt("normalization-type") == 1)
    {
        // == SKETCH
        m_apgEngine->setOpt("a5-sketch-input-global", "true" );
        m_apgEngine->setOpt("a5-sketch-input-name", "target-sketch" );
    }
    else
    {
        m_apgEngine->setOpt("a5-sketch-input-global", "false" );
        m_apgEngine->setOpt("a5-sketch-input-name", "" );
    }
    // == FEATURE-EFFECTS
    m_apgEngine->setOpt("a5-feature-effects-input-global", "true" );
    m_apgEngine->setOpt("a5-feature-effects-input-name",
            strGlobalPrefix+".plier.feature-response" );
    // == MODELS
    m_apgEngine->setOpt("a5-models-input-global", "true" );
    m_apgEngine->setOpt("a5-models-input-name", strGlobalPrefix+".snp-posteriors" );
}

/**
* @brief Call the geno-qc engine
*/
void CNWorkflowEngine::callGenoQCEngine()
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
    engine.setOpt("chrX-probes", getOpt("chrX-probes"));
    engine.setOpt("chrY-probes", getOpt("chrY-probes"));
    engine.setOpt("qca-file", getOpt("qca-file"));
    engine.setOpt("qcc-file", getOpt("qcc-file"));
    engine.setOpt("out-file", getOpt("geno-qc-file"));
    std::vector<std::string> cels = getOptVector("cels");
    engine.setOpt("cels", cels);
    engine.run();
}

/**
 * @brief Call the apg engine.
 */
void CNWorkflowEngine::callProbesetGenotypeEngine()
{
    std::vector<std::string> cels = getOptVector("cels");
    m_apgEngine->setOpt("cels", cels);
    m_apgEngine->run();
}

/**
 * @brief Call the CNReference engine.
 */
void CNWorkflowEngine::callCNReferenceEngine()
{
    // Create temp file
    TmpFileFactory tempFactory;
    tempFactory.setTmpdir(getOpt("temp-dir"));
    setOpt("temp-reference-file", tempFactory.genFilename_basic( "CNReference_", ".a5.tmp"));

    CNReferenceEngine objCNReferenceEngine;
    objCNReferenceEngine.setEngine(this);
    objCNReferenceEngine.run();
}

/**
 * Call the CNLog2Ratio/CopyNumber engines. Turn off intermediate file generation.
 * @return - bool value. (true if successful)
 */
void CNWorkflowEngine::callCNLog2RatioEngine()
{
    CNLog2RatioEngine objCNLog2RatioEngine;
    objCNLog2RatioEngine.setEngine(this);
    objCNLog2RatioEngine.setOpt("call-copynumber-engine", "true");
    objCNLog2RatioEngine.setOpt("log2ratio-hdf5-output", "false");
    objCNLog2RatioEngine.setOpt("log2ratio-text-output", "false");
    objCNLog2RatioEngine.run();
}

/**
 * @brief Return the annotation file names
 * @param const AffxString& - The name of the reference file
 * @param AffxString& - The name of the snp annotation file to return
 * @param AffxString& - The name of the cn annotation file to return
 * @return bool - true if successful
 */
bool CNWorkflowEngine::getAnnotationFileName(const AffxString& strFileName, AffxString& strAnnotationFileName)
{
    bool bSuccessful = false;
    strAnnotationFileName = "";
    AffxString strOldAnnotationFileName;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        AffxString str;
        try
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            if (group5 != NULL)
            {
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if (tsv5 != NULL)
            {
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &str);
                int iIndex = str.indexOf("=");
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if (str.startsWith("#%affymetrix-algorithm-param-apt-opt-annotation-file=")) {strAnnotationFileName = str.substring(iIndex + 1);}
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if (str.startsWith("#%affymetrix-algorithm-param-apt-opt-netaffx-snp-annotation-file=")) {strOldAnnotationFileName = str.substring(iIndex + 1);}
            }
            tsv5->close();
            delete tsv5;
            }
            group5->close();
            delete group5;
            }
            file5.close();
            if (strAnnotationFileName != "") {bSuccessful = true;}
        } catch(...) {return false;}
    }
    if ((strAnnotationFileName == "") && (strOldAnnotationFileName != ""))
    {
        strAnnotationFileName = strOldAnnotationFileName;
        bSuccessful = true;
    }
    return bSuccessful;
}

/**
 * @brief Return the normalization types
 * @param const AffxString& - The name of the reference file
 * @param int& - The normalization type
 * @param bool& - The adapter normalization
 * @param bool& - The marker normalization
 * @return bool - true if successful
 */
// This method is called by GTC without having called checkOptions.
// In short, best not to refer to Engine options/state in this method
bool CNWorkflowEngine::getNormalizationTypes(const AffxString& strFileName, int& iNormalizationType, bool& bAdapterNormalization, int& iSampleCount)
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
            affx::File5_Group* group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            if (group5 != NULL)
            {
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if (tsv5 != NULL)
            {
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &str);
                int iIndex = str.indexOf("=");
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
 * @brief Return the analysis name
 * @param const AffxString& - The name of the reference file
 * @param AffxString& - The analysis name
 * @return bool - true if successful
 */
bool CNWorkflowEngine::getAnalysisName(const AffxString& strFileName, AffxString& strAnalysisName)
{
    bool bSuccessful = false;
    strAnalysisName = "";
    if (affx::File5_File::isHdf5file(strFileName))
    {
        AffxString str;
        try
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            if (group5 != NULL)
            {
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if (tsv5 != NULL)
            {
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &str);
                int iIndex = str.indexOf("=");
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if ((str.startsWith("#%affymetrix-algorithm-param-apt-set-analysis-name="))  ||
                    (str.startsWith("#%affymetrix-algorithm-param-apt-analysis-name="))  ||
                    (str.startsWith("#%affymetrix-algorithm-param-apt-opt-set-analysis-name=")) ||
                    (str.startsWith("#%affymetrix-algorithm-param-apt-opt-analysis-name=")))
                {
                    strAnalysisName = str.substring(iIndex + 1);
                }
            }
            tsv5->close();
            delete tsv5;
            }
            group5->close();
            delete group5;
            }
            file5.close();
            if (strAnalysisName != "") {bSuccessful = true;}
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
bool CNWorkflowEngine::getReferenceGenderCounts(const AffxString& strFileName, int& iReferenceMaleCount, int& iReferenceFemaleCount, int& iReferenceUnknownCount)
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
            affx::File5_Group* group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            if (group5 != NULL)
            {
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            if (tsv5 != NULL)
            {
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                tsv5->get(0, 0, &str);
                int iIndex = str.indexOf("=");
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
void CNWorkflowEngine::calculateGenders()
{
    CnProbeGenderCelListener objGenderCaller(getOpt("chrX-probes"), getOpt("chrY-probes"), getOptDouble("female-gender-ratio-cutoff"), getOptDouble("male-gender-ratio-cutoff"));
    affymetrix_fusion_io::FusionCELData cel;
    int iUnknownCount = 0;
    int iMaleCount = 0;
    int iFemaleCount = 0;
    std::vector<std::string> celFiles = getOptVector("cels");
    for (unsigned int uiCelFileIndex = 0; (uiCelFileIndex < celFiles.size()); uiCelFileIndex++)
    {
        cel.SetFileName(celFiles[uiCelFileIndex].c_str());
        if (!cel.Read()) {error("Cannot read cel file for gender calling: " + celFiles[uiCelFileIndex]);}
        objGenderCaller.newChip(&cel);
        cel.Close();
        switch(objGenderCaller.getLastGender())
        {
        case affx::Male: iMaleCount++; Verbose::out(1, "Sample " + ::getUnsignedInt(uiCelFileIndex + 1) + ":\t" + Fs::basename(celFiles[uiCelFileIndex]) + "\tMale"); break;
        case affx::Female: iFemaleCount++; Verbose::out(1, "Sample " + ::getUnsignedInt(uiCelFileIndex + 1) + ":\t" + Fs::basename(celFiles[uiCelFileIndex]) + "\tFemale"); break;
        default: iUnknownCount++; Verbose::out(1, "Sample " + ::getUnsignedInt(uiCelFileIndex + 1) + ":\t" + Fs::basename(celFiles[uiCelFileIndex]) + "\tUnknown");
        }
    }
    Verbose::out(1, "*");
    Verbose::out(1, "The cel file input contains " + ::getInt(iFemaleCount) + " female, " + ::getInt(iMaleCount) + " male, and " + ::getInt(iUnknownCount) + " unknown genders.");
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
            Verbose::out(1, "The reference input is based on " + ::getInt(iReferenceFemaleCount) + " female, " + ::getInt(iReferenceMaleCount) + " male, and " + ::getInt(iReferenceUnknownCount) + " unknown genders.");
            if (iSampleCount < getOptInt("recommended-reference-sample-count")) {Verbose::out(1, "WARNING: Reference is based on less than " + getOpt("recommended-reference-sample-count") + " samples.");}
            if (iReferenceMaleCount < getOptInt("recommended-reference-male-sample-count")) {Verbose::out(1, "WARNING: Reference is based on less than " + getOpt("recommended-reference-male-sample-count") + " male samples.");}
            if (iReferenceFemaleCount < getOptInt("recommended-reference-female-sample-count")) {Verbose::out(1, "WARNING: Reference is based on less than " + getOpt("recommended-reference-female-sample-count") + " female samples.");}
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

// This method is called by GTC without having called checkOptions.
// In short, best not to refer to Engine options/state in this method
bool CNWorkflowEngine::isCN5Reference(const AffxString& strReferenceFileName)
{
    return (affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "CopyNumber", "Reference") > 0);
}

bool CNWorkflowEngine::isCyto2Reference(const AffxString& strReferenceFileName)
{
    return (affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "Cyto2", "MedianSignals") > 0);
}

// Check the chip type
void CNWorkflowEngine::checkChipType()
{
    colrow_t numRows = 0, numCols = 0;
    int probeCount = 0, probeSetCount = 0, channelCount = 1;
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
    setOpt("channel-count", ToStr(channelCount));

	if (!getOptBool("force")) {
		std::vector<std::string> cels = getOptVector("cels");
		EngineUtil::checkCelChipTypes(chipTypes, probeCount, cels, numRows, numCols);
		if (getOpt("chrX-probes") != "") {EngineUtil::checkTsvFileChipType(getOpt("chrX-probes"), chipTypes);}
		if (getOpt("chrY-probes") != "") {EngineUtil::checkTsvFileChipType(getOpt("chrY-probes"), chipTypes);}
	}
}

void CNWorkflowEngine::explain() {
    CNAnalysisMethodFactory factory;
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

void CNWorkflowEngine::extraHelp() {
    CNAnalysisMethodFactory factory;
    printStandardMethods(cout);
    EngineUtil::printSelfDocs("Data transformations:", factory.getDocs());
}

// This method is called by GTC without having called checkOptions.
// In short, best not to refer to Engine options/state in this method
AffxString CNWorkflowEngine::getAnnotationParameter(const AffxString& strFileName, const AffxString& strParameterName)
{
    AffxString strParameterValue;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        AffxString str;
        try
        {
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
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

