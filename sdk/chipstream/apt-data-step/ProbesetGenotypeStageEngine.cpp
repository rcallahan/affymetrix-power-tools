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

//
#include "chipstream/apt-data-step/ProbesetGenotypeStageEngine.h"
//
#include "chipstream/AnalysisStream.h"
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/CelReader.h"
#include "chipstream/CelStatListener.h"
#include "chipstream/ChipStream.h"
#include "chipstream/CnProbeGenderCelListener.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/DiskIntensityMart.h"
#include "chipstream/SparseMart.h"
#include "chipstream/DmListener.h"
#include "chipstream/DTFactory.h"
#include "chipstream/EmGenderCelListener.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/FileGenders.h"
#include "chipstream/FileInbred.h"
#include "chipstream/GenoSeed.h"
#include "chipstream/GenoSeedTxt.h"
#include "chipstream/GenotypeInfo.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantBRLMM.h"
#include "chipstream/QuantBirdseedDev.h"
#include "chipstream/QuantBirdseedLegacy.h"
#include "chipstream/QuantBirdseedv1.h"
#include "chipstream/QuantBirdseedv2.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantLabelZ.h"
#include "chipstream/QuantLabelZMulti.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodGTypeCHPReport.h"
#include "chipstream/QuantMethodGTypeChipSummary.h"
#include "chipstream/QuantMethodGTypeExprChipSummary.h"
#include "chipstream/QuantMethodGTypeReport.h"
#include "chipstream/apt-data-step/AnalysisStage.h"

/* Avoid compiling these into windows builds until GTC team has been 
   informed and given time to update their builds as new libraries are requred. */
#ifndef WIN32 
#include "chipstream/QuantMethodSqlGTypeReport.h"
#include "chipstream/QuantMethodSqlExprReport.h"
#endif

#include "chipstream/SnpModelConverter.h"
#include "chipstream/QuantMethodMultiDataCCCHPReport.h"
#include "chipstream/QuantMethodRunReport.h"
#include "chipstream/SpecialSnps.h"
//
#include "broadutil/BroadException.h"
#include "broadutil/CelUtil.h"
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/fusion/src/FusionCDFData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "file/TsvFile/TsvFile.h"

#include "util/BaseEngine.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/PgOptions.h"
#include "util/Util.h"
#include "util/Verbose.h"
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

using namespace std;
using namespace affx;
using namespace affymetrix_calvin_exceptions;

ProbesetGenotypeStageEngine::Reg ProbesetGenotypeStageEngine::reg;

ProbesetGenotypeStageEngine * ProbesetGenotypeStageEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == ProbesetGenotypeStageEngine::EngineName())
        return (ProbesetGenotypeStageEngine *)engine;
    return NULL;
}

ProbesetGenotypeStageEngine::ProbesetGenotypeStageEngine(
        ChipLayout *chipLayout,
        std::vector<const char *> *probesetNames) {
    m_a5_global_input_file = NULL;
    m_a5_global_input_group = NULL;
    m_a5_global_output_file = NULL;
    m_a5_global_output_group = NULL;
    defineStdMethods();
    defineOptions();

    if(chipLayout != NULL) {
      if(probesetNames == NULL)
          Err::errAbort("Must specify probesetNames when providing a chipLayout to ProbesetGenotypeStageEngine constructor.");
      m_FreeChipLayout = false;
      m_ChipLayout = chipLayout;
      m_ProbesetNames = *probesetNames;
    } else {
      m_FreeChipLayout = true;
      m_ChipLayout = NULL;
    }
}

ProbesetGenotypeStageEngine::~ProbesetGenotypeStageEngine() {
}

void ProbesetGenotypeStageEngine::defineStdMethods() {
    m_stdMethods["brlmm"] = "quant-norm.sketch=50000,pm-only,brlmm.transform=ccs.K=4";
    m_stdMethods["brlmm-p"] = "quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100.K=2.SB=0.003.MS=0.05";
    m_stdMethods["brlmm-p.force"] = "quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100.K=2.SB=0.003.MS=1";
    m_stdMethods["brlmm-p-plus"] = "quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.HARD=3.SB=0.45.KX=1.KH=1.5.KXX=0.5.KAH=-0.6.KHB=-0.6.transform=MVA.AAM=2.0.BBM=-2.0.AAV=0.06.BBV=0.06.ABV=0.06.copyqc=0.000001.wobble=0.05.MS=0.05";
    m_stdMethods["brlmm-p-plus.force"] = "quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.HARD=3.SB=0.45.KX=1.KH=1.5.KXX=0.5.KAH=-0.6.KHB=-0.6.transform=MVA.AAM=2.0.BBM=-2.0.AAV=0.06.BBV=0.06.ABV=0.06.copyqc=0.000001.wobble=0.05.MS=1";
    m_stdMethods["birdseed"] = "quant-norm.sketch=50000,pm-only,birdseed";
    m_stdMethods["birdseed.force"] = "quant-norm.sketch=50000,pm-only,birdseed.conf-threshold=1";
    m_stdMethods["birdseed-v1"] = "quant-norm.sketch=50000,pm-only,birdseed-v1";
    m_stdMethods["birdseed-v1.force"] = "quant-norm.sketch=50000,pm-only,birdseed-v1.conf-threshold=1";
    m_stdMethods["birdseed-v2"] = "quant-norm.sketch=50000.target=1000,pm-only,birdseed-v2";
    m_stdMethods["birdseed-v2.force"] = "quant-norm.sketch=50000.target=1000,pm-only,birdseed-v2.conf-threshold=1";
    m_stdMethods["birdseed-dev"] = "quant-norm.sketch=50000.target=1000,pm-only,birdseed-dev";
    m_stdMethods["birdseed-dev.force"] = "quant-norm.sketch=50000.target=1000,pm-only,birdseed-dev.conf-threshold=1";
}

void ProbesetGenotypeStageEngine::defineOptions() {

  defineOptionSection("Input Options");
  defineOption("", "cel-files", PgOpt::STRING_OPT,
                     "Text file specifying cel files to process, one per line with the first line being 'cel_files'.",
                     "");
  defineOption("c", "cdf-file", PgOpt::STRING_OPT,
                     "File defining probe sets. Use either --cdf-file or --spf-file. ",
                     "");
  defineOption("", "spf-file", PgOpt::STRING_OPT,
                     "File defining probe sets in spf (simple probe format) which is like a text cdf file.",
                     "");
  defineOption("", "chrX-snps", PgOpt::STRING_OPT,
                    "File containing snps on chrX (non-pseudoautosomal region).",
                    "");
  defineOption("", "special-snps",PgOpt::STRING_OPT,
                    "File containing all snps of unusual copy (chrX,mito,Y)",
                    "");
  defineOption("", "chrX-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrX. "
                     "Used for copy number probe chrX/Y ratio gender calling. "
                     "[Experimental]",
                     "");
  defineOption("", "chrY-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrY. "
                     "Used for copy number probe chrX/Y ratio gender calling. "
                     "[Experimental]",
                     "");
  defineOption("s", "probeset-ids", PgOpt::STRING_OPT,
                     "Tab delimited file with column 'probeset_id' specifying probesets to genotype.",
                     "");
  defineOption("", "probeset-ids-reported", PgOpt::STRING_OPT,
                     "Tab delimited file with column 'probeset_id' specifying probesets to report. "
                     "This should be a subset of those specified with --probeset-ids if that option is used.",
                     "");
  defineOption("", "probe-class-file", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes and a 'class' designation. "
                     "Used to compute mean probe intensity by class for report file.",
                     "");
  defOptMult("", "chip-type", PgOpt::STRING_OPT,
                     "Chip types to check library and CEL files against. "
                     "Can be specified multiple times. "
                     "The first one is propigated as the chip type in the output files. "
                     "Warning, use of this option will override the usual check between chip types "
                     "found in the library files and cel files. You should use this option instead "
                     "of --force when possible. ",
                     "");
  defineOption("","annotation-file", PgOpt::STRING_OPT,
                     "Annotation file.",       
                     "");
  defineOption("", "genotype-markers-cn-file", PgOpt::STRING_OPT, 
                     "Tab delimited file with copy number calls for genotype probesets within copy number regions", 
                     "");
  defineOption("", "writeOldStyleFeatureEffectsFile" , PgOpt::BOOL_OPT, 
                   "Boolean value to determine whether or not old style feature effects files are written.",
                   "false"); 
  defineOption("", "file5-compact" , PgOpt::BOOL_OPT, 
                   "Should we output results in a compact file5 output.",
                   "false"); 
/* Avoid compiling these into windows builds until GTC team has been 
   informed and given time to update their builds as new libraries are requred. */
#ifndef WIN32
  defineOption("", "sqlite-output" , PgOpt::BOOL_OPT, 
                   "Shoul output some results in sqlite3 format?",
                   "false"); 
#endif
  defineOptionSection("Output Options");
  defineOption("", "table-output", PgOpt::BOOL_OPT,
                    "Output matching matrices of tab delimited genotype calls and confidences.",
                    "true");
  defineOption("", "output-forced-calls", PgOpt::BOOL_OPT,
                    "Output a separate file with forced calls.",
                    "false");
  defineOption("", "output-context", PgOpt::BOOL_OPT,
                    "Output a separate file with the allele context used. This is only relevant for marker type probesets which have multiple groups of probes for each allele based on the context of nearby SNPs.",
                    "false");
  defineOption("", "cc-chp-output", PgOpt::BOOL_OPT,
                    "Output resulting calls in directory called 'cc-chp' under out-dir. "
                    "This makes one AGCC Multi Data CHP file per cel file analyzed.",
                    "false");
  defineOption("", "xda-chp-output", PgOpt::BOOL_OPT,
                    "Output resulting calls in directory called 'chp' under out-dir. "
                    "This makes one GCOS XDA CHP file per cel file analyzed. "
                    "Note that this format is not supported beyond the Mapping500K chips, "
                    "for subsequent chips look at the CC CHP format instead.",
                    "false");
  defineOption("", "cc-chp-out-dir", PgOpt::STRING_OPT,
                     "Over-ride the default location for chp output.",
                     "");
  defineOption("", "xda-chp-out-dir", PgOpt::STRING_OPT,
                     "Over-ride the default location for chp output.",
                     "");
  defineOption("", "summaries", PgOpt::BOOL_OPT,
                     "Output the summary values from the quantifcation method for each allele. "
                     "For brlmm-p this will also write a file of transformed summary values "
                     "in contrast space used in the clustering.",
                    "false");
  defineOption("", "report-file", PgOpt::STRING_OPT,
                     "Over-ride the default report file name.",
                     "");

  defineOptionSection("Analysis Options");

  defOptMult("a", "analysis", PgOpt::STRING_OPT,
                   "String representing analysis pathway desired. "
                   "For example: 'quant-norm.sketch=50000,pm-only,brlmm'.",
                   "brlmm");
  defineOption("", "qmethod-spec", PgOpt::STRING_OPT,
                     "Quantification Method to use for summarizing alleles.",
                     "plier.optmethod=1");
  defineOption("", "read-models-brlmm", PgOpt::STRING_OPT,
                     "File to read precomputed BRLMM snp specific models from.",
                     "");
  defineOption("", "read-models-brlmmp", PgOpt::STRING_OPT,
                     "File to read precomputed BRLMM-P snp specific models from.",
                     "");
  defineOption("", "brlmmp-models-non-sequential", PgOpt::BOOL_OPT,
                     "File to read precomputed BRLMM-P snp specific models from is not in sequential order.",
                     "false");
  defineOption("", "read-models-birdseed", PgOpt::STRING_OPT,
                     "File to read precomputed birdseed snp specific models from.",
                     "");
  defineOption("", "write-models", PgOpt::BOOL_OPT,
                     "Should we write snp specific models out for analysis? [experimental]",
                     "false");
  defineOption("", "db-from-prior-models", PgOpt::STRING_OPT,
                     "File to write prior snp models to for random access.",
                     "");
  defineOption("", "db-from-posterior-models", PgOpt::STRING_OPT,
                     "File to write posterior snp models to for random access.",
                     "");
  defineOption("", "feat-effects", PgOpt::BOOL_OPT,
                     "Output feature effects when available.",
                     "false");
  defineOption("", "use-feat-eff", PgOpt::STRING_OPT,
                    "File defining a plier feature effect for each probe. "
                     "Note that precomputed effects should only be used for an appropriately similar analysis "
                     "(i.e. feature effects for pm-only may be different than for pm-mm).",
                     "");
  defineOption("", "feat-details", PgOpt::BOOL_OPT,
                     "Output the feature details (usually residuals) from the quantification method if available.",
                     "false");
  defineOption("", "target-sketch", PgOpt::STRING_OPT,
                    "File specifying a target distribution to use for quantile normalization.",
                    "");
  defineOption("", "write-sketch", PgOpt::BOOL_OPT,
                    "Write the quantile normalization distribution (or sketch) to a file for reuse with target-sketch option.",
                    "false");
  defineOption("", "dm-thresh", PgOpt::DOUBLE_OPT,
                    "Minimum DM p-value to seed clusters with.",
                     ".17");
  defineOption("","reference-profile", PgOpt::STRING_OPT,
                    "File specifying reference chip profile.",
                    "");
  defineOption("","write-profile", PgOpt::BOOL_OPT,
                    "Write the reference chip profile to a file for reuse.",
                    "false");
  defineOption("", "dm-hetmult", PgOpt::DOUBLE_OPT,
                     "DM hetmultiplier to balance het/hom calls, additive to log likelihood.",
                     "0");
  defineOption("", "prior-size", PgOpt::INT_OPT,
                    "How many probesets to use for determining prior.",
                    "0");
  defineOption("", "list-sample", PgOpt::BOOL_OPT,
                    "Only sample for prior from list specified via --probeset-ids, not entire chip.",
                    "false");
  defineOption("", "read-priors-brlmm", PgOpt::STRING_OPT,
                    "File to load BRLMM priors from. Prior format is tab separated id, center, var, and center.var.",
                    "");
  defineOption("", "write-prior", PgOpt::BOOL_OPT,
                    "Write prior out to file in output-dir.",
                    "false");
  defineOption("", "norm-size", PgOpt::INT_OPT,
                    "Do contrast normalization using a sample of this many snps (brlmm-p)",
                    "0");
  defineOption("", "write-norm", PgOpt::BOOL_OPT,
                    "Write covariate norm fcns to file",
                    "false");
  defineOption("", "set-analysis-name", PgOpt::STRING_OPT,
                  "Explicitly set the analysis name. This affects output file names (ie prefix) and various meta info.",
                    "");
  defineOption("", "include-quant-in-report-file-name", PgOpt::STRING_OPT,
                  "Include the quant method name in the expression report files. ",
                    "false");

  defineOptionSection("Gender Options");

  defineOption("", "set-gender-method", PgOpt::STRING_OPT,
                    "Explicitly force the use of a particular gender method for genotype calling. "
                    "Valid values include: cn-probe-chrXY-ratio, dm-chrX-het-rate, "
                    "em-cluster-chrX-het-contrast, user-supplied, and none. "
                    "If you are supplied seed genotype calls, you can also use "
                    "supplied-genotypes-chrX-het-rate. When not set, the default behavior depends on the analysis.",
                    "");
  defineOption("", "read-genders", PgOpt::STRING_OPT,
                    "Explicitly read genders from a file.",
                    "");
  defineOption("", "read-inbred", PgOpt::STRING_OPT,
                    "Read penalty for hets by level of inbreeding per sample.",
                    "");
  defineOption("", "no-gender-force", PgOpt::BOOL_OPT,
                    "Perform analysis even without a suitable gender method for genotype calling.",
                    "false");
  defineOption("", "em-gender", PgOpt::BOOL_OPT,
                    "Enable EM Gender calling if special-snps or chrX-snp file is provided.",
                    "true");
  defineOption("", "female-thresh", PgOpt::DOUBLE_OPT,
                    "Threshold for calling females when using cn-probe-chrXY-ratio method.",
                    "0.48");
  defineOption("", "male-thresh", PgOpt::DOUBLE_OPT,
                    "Threshold for calling females when using cn-probe-chrXY-ratio method.",
                    "0.71");

  defineOptionSection("Misc Options");
  defineOption("", "explain", PgOpt::STRING_OPT,
                    "Explain a particular operation (i.e. --explain brlmm or --explain brlmm-p).",
                    "");

  defineOptionSection("Advanced Options");
  defineOption("", "kill-list", PgOpt::STRING_OPT,
                    "Do not use the probes specified in file for computing results. [experimental]",
                    "");
  defineOption("", "dm-out", PgOpt::BOOL_OPT,
                    "Output any initial seed calls used by BRLMM (seed default is DM calls). Only relevant for BRLMM.",
                    "false");
  defineOption("", "all-types", PgOpt::BOOL_OPT,
                    "Try and analyze all probeset types rather than just genotyping. "
                     "[Experimental]",
                    "false");
  defineOption("", "genotypes", PgOpt::STRING_OPT,
                     "File to read seed genotypes from instead of using DM to generate. [experimental]",
                     "");
  defineOption("", "select-probes", PgOpt::BOOL_OPT,
                    "Output estimates of which probes are most accurate",
                    "false");
  defineOption("", "call-coder-max-alleles", PgOpt::INT_OPT,
                    "For encoding/decoding calls, the max number of alleles per marker to allow.",
                    "6");
  defineOption("", "call-coder-type", PgOpt::STRING_OPT,
                    "The data size used to encode the call.",
                    "UCHAR");
  defineOption("", "call-coder-version", PgOpt::STRING_OPT,
                    "The version of the encoder/decoder to use",
                    "1.0");

  defineOptionSection("Execution Control Options");
  defineOption("", "use-disk", PgOpt::BOOL_OPT,
                     "Store CEL intensities to be analyzed on disk.", "true");
  defineOption("", "disk-cache", PgOpt::INT_OPT,
                     "Size of memory cache when working off disk in megabytes.",
                     "50");

  defineOptionSection("A5 output options");

  defineOption("","a5-global-file",PgOpt::STRING_OPT,
                     "Filename for the A5 global output file. "
                     "[Experimental]",
                     "");
  defineOption("","a5-global-file-no-replace",PgOpt::BOOL_OPT,
                     "Append or create rather than replace. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-group",PgOpt::STRING_OPT,
                     "Group name where to put results in the A5 output files. "
                     "[Experimental]",
                     "");
  defineOption("","a5-calls",PgOpt::BOOL_OPT,
                     "Output the genotype calls and confidences in A5 format. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-calls-use-global",PgOpt::BOOL_OPT,
                     "Use the global A5 file for calls and confidences."
                     "[Experimental]",
                     "false");
  defineOption("", "a5-summaries", PgOpt::BOOL_OPT,
                     "Output the summary values from the quantifcation method for each allele in A5 format. "
                     "[Experimental]",
                    "false");
  defineOption("","a5-summaries-use-global",PgOpt::BOOL_OPT,
                     "Use the global A5 file for summaries. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-feature-effects",PgOpt::BOOL_OPT,
                     "Output feature effects in A5 format. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-feature-effects-use-global",PgOpt::BOOL_OPT,
                     "Use the global A5 file for feature effects."
                     "[Experimental]",
                     "false");
  defineOption("","a5-feature-details",PgOpt::BOOL_OPT,
                     "Output feature level residuals in A5 format. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-feature-details-use-global",PgOpt::BOOL_OPT,
                     "Use the global A5 file for residuals. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-sketch",PgOpt::BOOL_OPT,
                     "Output normalization sketch in A5 format. "
                     "--write-sketch option will override this option. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-sketch-use-global",PgOpt::BOOL_OPT,
                     "Put the sketch in the global A5 output file. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-write-models",PgOpt::BOOL_OPT,
                     "Output genotype models/posteriors in A5 format. "
                     "--write-models option will override this option. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-write-models-use-global",PgOpt::BOOL_OPT,
                     "Put the models in the global A5 output file. "
                     "[Experimental]",
                     "false");
  
  defineOptionSection("A5 input options");

  defineOption("","a5-global-input-file",PgOpt::STRING_OPT,
                     "Filename for the group in the global input file."
                     "[Experimental]",
                     "");
  defineOption("","a5-input-group",PgOpt::STRING_OPT,
                     "Group name for input. Defaults to --a5-group or if that "
                     "is not set, then '/'. "
                     "[Experimental]",
                     "");

  defineOption("","a5-sketch-input-global",PgOpt::BOOL_OPT,
                     "Read the sketch from the global A5 input file. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-sketch-input-file",PgOpt::STRING_OPT,
                     "Read the sketch from the an A5 input file. "
                     "[Experimental]",
                     "");
  defineOption("","a5-sketch-input-group",PgOpt::STRING_OPT,
                     "Group name to read the sketch from. Defaults to --a5-input-group. "
                     "[Experimental]",
                     "");
  defineOption("","a5-sketch-input-name",PgOpt::STRING_OPT,
                     "The name of the data section. Defaults to 'target-sketch'. "
                     "[Experimental]",
                     "");

  defineOption("","a5-feature-effects-input-global",PgOpt::BOOL_OPT,
                     "Read the feature effects global A5 input file. "
                     "[Experimental]",
                     "false");
  defineOption("","a5-feature-effects-input-file",PgOpt::STRING_OPT,
                     "Read the feature effects from the an A5 input file. "
                     "[Experimental]",
                     "");
  defineOption("","a5-feature-effects-input-group",PgOpt::STRING_OPT,
                     "Group name to read the feature effects from. Defaults to --a5-input-group. "
                     "[Experimental]",
                     "");
  defineOption("","a5-feature-effects-input-name",PgOpt::STRING_OPT,
                     "The name of the data section. Defaults to XXX.feature-response "
                     "where XXX is the analysis name and quant method. IE 'brlmm-p.plier'. "
                     "[Experimental]",
                     "");

  defineOption("","a5-models-input-global",PgOpt::BOOL_OPT,
                     "Read the Models from the global A5 input file. "
                     "The tsv5 name must be 'XXX.snp-posteriors'. " // jhg
                     "[Experimental]",
                     "false");
  defineOption("","a5-models-input-file",PgOpt::STRING_OPT,
                     "Read the models from the an A5 input file. "
                     "[Experimental]",
                     "");
  defineOption("","a5-models-input-group",PgOpt::STRING_OPT,
                     "The group name where the models are located. Defaults to "
                     "the analysis name. "
                     "[Experimental]",
                     "");
  defineOption("","a5-models-input-name",PgOpt::STRING_OPT,
                     "The name of the data section. Defaults to XXX.snp-posteriors "
                     "where XXX is the analysis name. IE 'brlmm-p'. "
                     "[Experimental]",
                     "");

  defineOptionSection("Engine Options (Not used on command line)");
  defOptMult("", "cels", PgOpt::STRING_OPT,
                     "Cel files to process.",
                     "");
  defOptMult("", "result-files", PgOpt::STRING_OPT,
                     "CHP file names to output. Must be paired with cels.",
                     "");
}

void ProbesetGenotypeStageEngine::defineStates() {
  defineOption("", "analysis-count", PgOpt::INT_OPT,
                     "The number of user-specified analysis streams.",
                     "-1");
  defineOption("", "num-rows", PgOpt::INT_OPT,
                     "The number of rows on the chip.",
                     "-1");
  defineOption("", "num-cols", PgOpt::INT_OPT,
                     "The number of cols on the chip.",
                     "-1");
  defineOption("", "probe-count", PgOpt::INT_OPT,
                     "The number of probes on the chip.",
                     "-1");
  defineOption("", "probeset-count", PgOpt::INT_OPT,
                     "The number of probesets on the chip.",
                     "-1");
  defineOption("", "channel-count", PgOpt::INT_OPT,
                     "The number of data channels in the CEL file.",
                     "-1");
  defineOption("", "need-gc", PgOpt::BOOL_OPT,
                     "Do we need GC probes for the analysis.",
                     "false");
  defineOption("", "need-mismatch", PgOpt::BOOL_OPT,
                     "Do we need mismatch probes for the analysis.",
                     "false");
  defineOption("", "need-allelematch", PgOpt::BOOL_OPT,
                     "Do we need matching allele probes for the analysis.",
                     "false");
}

/**
 * Check the special-snp-version field in birdseed prior file. Call
 * Err::errAbort() if mismatch.
 *
 * @param path           - The birdseed prior file name
 * @param specialSnpType - The chiptype expected
 */
void ProbesetGenotypeStageEngine::checkBirdseedTsvPriorFile(const string &path, const string &specialSnpType) {
  affx::TsvFile tsv;
  string chip;
  int specialSeen = 0;
  if(tsv.open(path) != TSV_OK)
    Err::errAbort("Couldn't open file: '" + path + "' to read.");
  bool specialSnpMatch = specialSnpType.empty(); // treat empty as wildcard
  while(tsv.headersFindNext("special-snp-version", chip) == TSV_OK) {
    specialSeen++;
    if(chip == specialSnpType)
      specialSnpMatch = true;
  }
  if(!specialSnpMatch && specialSeen > 0) {
    Err::errAbort( "Special Snp File: '" + specialSnpType + "' not supported by prior file: '" + path + "'.");
  }
}

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void ProbesetGenotypeStageEngine::checkOptionsImp() {

  defineStates();

  setLibFileOpt("cdf-file");
  setLibFileOpt("spf-file");
  setLibFileOpt("chrX-snps");
  setLibFileOpt("special-snps");
  setLibFileOpt("chrX-probes");
  setLibFileOpt("chrY-probes");
  setLibFileOpt("probeset-ids");
  setLibFileOpt("probeset-ids-reported");
  setLibFileOpt("probe-class-file");
  setLibFileOpt("annotation-file");
  setLibFileOpt("genotype-markers-cn-file");
  setLibFileOpt("read-models-brlmm");
  setLibFileOpt("read-models-brlmmp");
  setLibFileOpt("read-models-birdseed");
  setLibFileOpt("use-feat-eff");
  setLibFileOpt("target-sketch");
  setLibFileOpt("reference-profile");
  setLibFileOpt("read-priors-brlmm");
  setLibFileOpt("kill-list");
  setLibFileOpt("a5-global-input-file");
  setLibFileOpt("a5-sketch-input-file");
  setLibFileOpt("a5-feature-effects-input-file");
  setLibFileOpt("a5-models-input-file");

  if(getOpt("explain") != "") { explain(); exit(0); }

  if (getOpt("out-dir") == "") {Err::errAbort("Must specify an output directory.");}
  if (getOpt("temp-dir") == "") { setOpt("temp-dir", Fs::join(getOpt("out-dir"), "temp")); }

  string cdfFile = getOpt("cdf-file");
  string spfFile = getOpt("spf-file");
  string killList = getOpt("kill-list");
  string specialSnps = getOpt("special-snps");
  string chrXSnps = getOpt("chrX-snps");
  string chrXProbes = getOpt("chrX-probes");
  string chrYProbes = getOpt("chrY-probes");
  vector<string> analysis = getOptVector("analysis");

  if(analysis.size() == 0) {
      analysis.push_back("brlmm");
      setOpt("analysis",analysis);
  }
  if(analysis.size() != 1) {
      Err::errAbort("Only one analysis is currently permitted.");
  }
  setOpt("analysis-count", ToStr(analysis.size()));

  // Check that a global file was given if the user wants global a5 output
  if(getOpt("a5-global-file") == "") {
      vector<string> optionNames;
      getOptionNames(optionNames);
      for(size_t i=0; i<optionNames.size(); i++) {
          string name = optionNames[i];
          if( name.find("use-global") != name.npos ) {
              if( getOptBool(name) ) {
                  Err::errAbort("--" + name + " option given, but no global file. Must specify --a5-global-file.");
              }
          }
      }
  }
  // Check that a global input file was given if the user wants global a5 input
  if(getOpt("a5-global-input-file") == "") {
      vector<string> optionNames;
      getOptionNames(optionNames);
      for(size_t i=0; i<optionNames.size(); i++) {
          string name = optionNames[i];
          if( name.find("input-global") != name.npos ) {
              if( getOptBool(name) ){
                  Err::errAbort("--" + name + " option given, but no global input file. Must specify --a5-global-input-file.");
              }
          }
      }
  }
  // make sure global input and output files are different
  if (getOpt("a5-global-input-file") == getOpt("a5-global-file"))
      if(getOpt("a5-global-input-file") != "")
          Err::errAbort("Unable to use the same A5 file for global input and output.");

  /* Read in cel file list from other file if specified. */
  vector<string> celFiles;
  EngineUtil::getCelFiles(celFiles, this);
  if(celFiles.size() == 0)
      Err::errAbort("No cel files specified.");
  setOpt("cels",celFiles);


  vector<string> resultFiles = getOptVector("result-files");
  if(resultFiles.size() > 0) {
      if(celFiles.size() != resultFiles.size())
        Err::errAbort("result-files option used but is not the same size as cel file listing");
      bool ccOut = getOptBool("cc-chp-output");
      bool xdaOut = getOptBool("xda-chp-output");
      if(ccOut && xdaOut) 
          Err::errAbort("cannot use multiple CHP output formats with the result-files option");
  }

  // If we're doing birdseed or brlmm-p we can save space by not loading the mismatch and gc mappings.
  bool needGc = false, needMisMatch = false, needAlleleMatch = false;
  for(size_t aIx = 0; aIx < analysis.size(); aIx++) {
      ///@todo would be nice to query AnalysisStreamFactory for this info
      string analysisName = analysis[aIx];
      map<string,string>::iterator iter;
      for(iter = m_stdMethods.begin(); iter != m_stdMethods.end(); iter++) {
          if(iter->first == analysisName) {
              analysisName = iter->second;
          }
      }
      if(analysis[aIx].find("pm-mm") != string::npos)
          needMisMatch = true;
      if(analysis[aIx].find("pm-gcbg") != string::npos)
          needGc = true;
      if(analysis[aIx].find("gc-bg") != string::npos)
          needGc = true;
      if(analysis[aIx].find("pm-sum") != string::npos)
          needAlleleMatch = true;

      ///@todo need to deal with DM listener, setting needMisMatch appropriately
      ///      this is a crude hack to get this. There may be cases where DM 
      ///      listener is run even though brlmm was not selected (and cases where
      ///      brlmm is run but DM seeds are not used).
      if(analysis[aIx].find(",brlmm.") != string::npos)
          needMisMatch = true;
  }
  setOpt("need-gc", ToStr(needGc));
  setOpt("need-mismatch",ToStr(needMisMatch));
  setOpt("need-allelematch",ToStr(needAlleleMatch));

  if(cdfFile == "" && spfFile == "") {
    Err::errAbort("Must specify either a cdf file or spf (simple probe format) file.");
  }
  if(killList != "" && cdfFile == "") {
    Err::errAbort("'--kill-list-file' option only supported on cdf files.");
   }
  if((chrXSnps  == "") && (specialSnps == "") && ((chrXProbes == "") || (chrYProbes == "")) && (getOpt("read-genders") == "") && !getOptBool("no-gender-force")) {
    Err::errAbort("Have to specify chrX snps, haploid snps, chrX/Y probe file, or a gender file in order to generate genders for genotype calling. Use the --no-gender-force option to force genotype calling.");
  }
  if(killList != "" && getOptBool("xda-chp-output")) {
    Err::errAbort("Can't specify a kill list and XDA CHP output at same time.");
  }
  if(killList != "" && getOptBool("cc-chp-output")) {
    Err::errAbort("Can't specify a kill list and CC CHP output at same time.");
  }
  if (chrXProbes != "" && chrYProbes == "") {
    Err::errAbort("Must provide a chrY Probe File when providing a chrX Probe File.");
  }
  if (chrXProbes == "" && chrYProbes != "") {
    Err::errAbort("Must provide a chrX Probe File when providing a chrY Probe File.");
  }
  if(getOpt("probeset-ids") != "" && getOptBool("xda-chp-output")) {
    Err::errAbort("Can't specify a subset and XDA CHP output at same time. Must use all probesets for CHP file output.");
  }
  if(getOpt("probeset-ids-reported") != "" && getOptBool("xda-chp-output")) {
    Err::errAbort("Can't specify a reporting subset and XDA CHP output at same time. Must use all probesets for CHP file output.");
  }
  if(getOpt("probeset-ids-reported") != "" && getOptBool("xda-chp-output")) {
    Err::errAbort("Can't specify a reporting subset and AGCC CHP output at same time. Must use all probesets for CHP file output.");
  }

  if(getOptBool("list-sample") && getOpt("probeset-ids") == "") {
    Err::errAbort("Can't set listSample without specifying a probeset-ids file.");
  }

  /* Some sanity checks. */
  if(analysis.empty())
      Err::errAbort("Must specify at least one analysis to preform.");

  if ((getOptInt("prior-size") == 0) &&
      (getOpt("read-priors-brlmm") == "")) {
    setOpt("prior-size", ToStr(10000)); // our default value.
  } else {
    setOpt("prior-size", getOpt("prior-size")); // our default value.
  }
  int numProbesets = 0;
  if(cdfFile != "") {
    /* Check to make sure that we aren't using more probesets for prior than are available. */
    ///@todo should have a similar check for spf input
    affymetrix_fusion_io::FusionCDFData cdf;
    cdf.SetFileName(cdfFile.c_str());
    if(!cdf.ReadHeader()) {
      Err::errAbort("Can't read header for " + cdfFile + ". Description: " + cdf.GetError());
    }
    affymetrix_fusion_io::FusionCDFFileHeader &cdfHeader = cdf.GetHeader();
    numProbesets = cdfHeader.GetNumProbeSets();
    if ((numProbesets < getOptInt("prior-size")) &&
        getOpt("read-priors-brlmm") == "") {
      Err::errAbort("Can't have size of prior: " + ToStr(getOptInt("prior-size")) +
                    " larger than number of probesets: " + ToStr(numProbesets));
    }
    cdf.Close();
  }
  /// @todo fix hack for discriminating snp5 vs snp6
  ///       aw: propose fix where we only allow xda chp output for
  ///           a specific set of chip types. Anything else is not
  ///           allowed.
  /* cws: If probesets are greater than 300K it has to be snp5
     or snp6 which shouldn't be used with xda format, has to be a
     better way to do this though... */
  if(numProbesets > 300000 && getOptBool("xda-chp-output")) {
    Err::errAbort("Can't make XDA style probesets with newer chips, use --cc-chp-output instead to get new format chp files.");
  }

  // special SNPs 
  bool needSpecialSnps = false;
  for(size_t aIx = 0; aIx < analysis.size(); aIx++) {
    if(analysis[aIx].find("birdseed") != string::npos) {
      needSpecialSnps = true;
      if ((getOpt("read-models-birdseed") == "")) {
        Err::errAbort("Birdseed requires models with '--read-models-birdseed'.");
      }
    }
  }
  if ((chrXSnps != "") && needSpecialSnps)
     Err::errAbort("Must use --special-snps rather than --chrX-snps when using birdseed "
                   "as birdseed expects chrY snps.");
  if ((specialSnps == "") && needSpecialSnps)
    Err::errAbort("Must specify '--special-snps' when special SNPS are needed.");

  // Check chip types
  vector<string> chipTypesInLayout;

  /* Get the intial info about the chip and check cel files to make sure
     they match. */
  colrow_t numRows = 0, numCols = 0;
  int probeCount = 0, probeSetCount=0, channelCount = 1;

  // Get channel count from first CEL file.
  affymetrix_fusion_io::FusionCELData cel;  
  try {
    if(!(Fs::fileExists(celFiles[0].c_str()))) {
      Err::errAbort("File: " + celFiles[0] + " doesn't appear to exist.");
      
    }
    cel.SetFileName(celFiles[0].c_str());
    if(!cel.Read()) 
      Err::errAbort("\nCan't read cel file: " + cel.GetFileName() + 
                    "\n>>> Error reported: " + StringUtils::ConvertWCSToMBS(cel.GetError()));
  }
  catch(const Except &e) {
    Err::errAbort("Error opening cel file: " + celFiles[0] + " error message is: " + e.what());
  }
  catch(...) {
    Err::errAbort("Error opening cel file: " + celFiles[0]);
  }

  // If this is a AGCC CEL file, then get channels.  If this is a GCOS
  // CEL file, then an empty vector is returned.
  std::vector<std::wstring> data_channels = cel.GetChannels();
  if (!data_channels.empty()) {
    channelCount = data_channels.size();
  }

  if (cdfFile != "") {
    EngineUtil::getCdfChipType(chipTypesInLayout, numRows, numCols, probeCount, probeSetCount, cdfFile);
  }
  else if (spfFile != "") {
    EngineUtil::getSpfChipType(chipTypesInLayout, numRows, numCols, probeCount, probeSetCount, spfFile);
  }
  else {
    Err::errAbort("Must specify either a cdf file or spf (simple probe format) file.");
  }

  setOpt("num-rows", ToStr(numRows));
  setOpt("num-cols", ToStr(numCols));
  setOpt("probe-count", ToStr(probeCount));
  setOpt("probeset-count", ToStr(probeSetCount));
  setOpt("channel-count", ToStr(channelCount));

  if(chipTypesInLayout.empty() || chipTypesInLayout[0] == "" || probeCount == 0) 
      Err::errAbort("Problem determining ChipType in file: " + 
              ( cdfFile != "" ? cdfFile : spfFile));

  /* Did the user "force" a set of chip types via options? */
  vector<string> chipTypesSupplied = getOptVector("chip-type");

  /* Figure out what chip type to report */
  if(chipTypesSupplied.size() > 0) {
      setOpt("chip-type", chipTypesSupplied[0]);
  } else if(chipTypesInLayout.size() > 0) {
      setOpt("chip-type", chipTypesInLayout[0]);
  } else {
      Err::errAbort("Unable to figure out a chip type.");
  }

  /* Do Chip Type Check */
  if(!getOptBool("force")) {
    vector<string> chipTypesToCheck;

    if(chipTypesSupplied.size() > 0) {
        chipTypesToCheck = chipTypesSupplied;
        EngineUtil::checkChipTypeVectors(chipTypesSupplied, chipTypesInLayout);
    } else {
        chipTypesToCheck = chipTypesInLayout;
    }

    EngineUtil::checkCelChipTypes(chipTypesToCheck, probeCount, celFiles, numRows, numCols);

    // Check model files
    // We only check the first (primary) chip type.
    // This ensures that model files for foo.cdf
    // but not foo.bar.cdf will fail the check
    if (getOpt("read-models-birdseed") != "") {
        EngineUtil::checkTsvFileChipType(getOpt("read-models-birdseed"), chipTypesToCheck);
        checkBirdseedTsvPriorFile(getOpt("read-models-birdseed"), "");
    }
    if (getOpt("read-models-brlmm") != "") {
        EngineUtil::checkTsvFileChipType(getOpt("read-models-brlmm"), chipTypesToCheck);
    }
    if (getOpt("read-models-brlmmp") != "") {
        EngineUtil::checkTsvFileChipType(getOpt("read-models-brlmmp"), chipTypesToCheck);
    }
    
    // Check prior files
    ///@todo check prior files chip type -- not tested
    //if (getOpt("read-priors-brlmm") != "") {
    //    EngineUtil:: checkTsvFileChipType(getOpt("read-priors-brlmm"), chipTypesToCheck);
    //}
    
    // Check special SNPs files
    if (specialSnps != "") {
        EngineUtil::checkTsvFileChipType(specialSnps, chipTypesToCheck);
    }
    if (chrXSnps != "") {
        EngineUtil::checkTsvFileChipType(chrXSnps, chipTypesToCheck);
    }
    
    // And other files
    if (getOpt("probe-class-file") != "") {
        EngineUtil::checkTsvFileChipType(getOpt("probe-class-file"), chipTypesToCheck);
    }
    if (chrXProbes != "") {
        EngineUtil::checkTsvFileChipType(chrXProbes, chipTypesToCheck);
    }
    if (chrYProbes != "") {
        EngineUtil::checkTsvFileChipType(chrYProbes, chipTypesToCheck);
    }
  }

  if(!Fs::isWriteableDir(getOpt("out-dir")))
    if(Fs::mkdirPath(getOpt("out-dir"), false) != APT_OK)
      APT_ERR_ABORT("Can't make or write to directory: " + getOpt("out-dir"));
  
  makeTempDir(getOpt("temp-dir"));
}

/**
 * @brief checks if there is enough disk space for diskmarts in the
 * temp-dir and for chp files in out-dir
 */
void ProbesetGenotypeStageEngine::checkDiskSpaceImp() {
  int64_t out_disk_available = -1, scratch_disk_available = -1;
  int cel_count = getOptVector("cels").size();
  int probeset_count = getOptInt("probeset-count");
  //  int row_count = getOptInt("num-rows");
  //  int col_count = getOptInt("num-cols");
  //  int channel_count = getOptInt("channel-count");
  //  int analysis_stream_count = getOptInt("analysis-count");
  
  std::vector<std::string> result_files = getOptVector("result-files");

  FsPath out_path;
  if (!result_files.empty()) {
    FsPath tmp_path(result_files[0]);
    if (tmp_path.getDirName()=="") {
      out_path.setDirName(tmp_path.getDirName());
    }
    // @todo iterate through all of the result dirs to check available disk space
  }

  if (out_path.empty()) {
    out_path.setDirName(getOpt("cc-chp-out-dir"));
  }
  if (out_path.empty()) {
    out_path.setDirName(getOpt("out-dir"));
  }

  out_disk_available = out_path.getFreeDiskSpace();

  uint64_t out_disk_needed = 0;  
  if (getOptBool("cc-chp-output")) {
    // @todo - why 35.49? CS
    out_disk_needed = static_cast<uint64_t>(probeset_count) * cel_count * 35.49;
  }
  
  // (1 + analysis_stream_count): one initial raw diskmart + one diskmart for
  // each analysis stream
  uint64_t scratch_disk_needed = 0;

  FsPath tempdir_path;
  tempdir_path.setDirName(getOpt("temp-dir"));

  bool same_vol = tempdir_path.isSameVolume(out_path);
  if (tempdir_path.empty() ||
      (same_vol==true) ||
      (tempdir_path == out_path)) {
    out_disk_needed += scratch_disk_needed;
  }
  else {
    scratch_disk_available = tempdir_path.getFreeDiskSpace();
    if (scratch_disk_needed > 0 && 
        scratch_disk_available >= 0 &&
        scratch_disk_needed >= scratch_disk_available) {
      // format in kb
      scratch_disk_needed = (scratch_disk_needed + 512) / 1024;
      scratch_disk_available = (scratch_disk_available + 512) / 1024;
      Err::errAbort("In " + tempdir_path.asString() + ", need " +
                    ToStr(scratch_disk_needed) + "kb of scratch diskspace. "
                    "Only " + ToStr(scratch_disk_available) + "kb available.");
    }
  }

  if (out_disk_needed > 0 &&
      out_disk_available >= 0 &&
      out_disk_needed >= out_disk_available) {
    // format in kb
    out_disk_needed = (out_disk_needed + 512) / 1024;
    out_disk_available = (out_disk_available + 512) / 1024;
    Err::errAbort("In " + out_path.asString() + ", need " + ToStr(out_disk_needed) + "kb of diskspace. "
                  "Only " + ToStr(out_disk_available) + "kb available.");
  }
}
  
/*   For copynumber chp files, the formula is */
/* 42.08 * num_probesets_processed */

/**
 * @brief Get a list from a text file to select a subset of probe sets for
 * analysis.
 *
 * @param fileName - Path to file for reading probe subset from.
 * @param vNames - vector of probeset names
 * @param sNames - set of probeset names
 */
void ProbesetGenotypeStageEngine::loadPSFile(const std::string &fileName,
                        std::vector<const char *> *vNames,
                        std::set<const char *, Util::ltstr> *sNames) {
  affx::TsvFile tsv;
  std::string probeset;
  tsv.bind(0, "probeset_id", &probeset, TSV_BIND_REQUIRED);
  if(tsv.open(fileName) != TSV_OK)
    Err::errAbort("Couldn't open file: " + fileName + " to read.");
  int psCount = 0;
  while(tsv.nextLevel(0) == TSV_OK) {
    if(vNames != NULL)
        vNames->push_back(Util::cloneString(probeset.c_str()));
    if(sNames != NULL) {
        if(vNames != NULL)
            sNames->insert(vNames->at(psCount));
        else
            sNames->insert(Util::cloneString(probeset.c_str()));
    }
    psCount++;
  }
  tsv.close();
}

/**
 * Read the design of the chip from the cdf or spf file.
 *
 * @param layout        - Object to be filled in.
 * @param probesetNames - Vector of probeset names to be filled in
 * @param killList      - Vector of probes to ignore
 * @param justStats     - Boolean indicating whether we want to keep probe layout info
 * @param  psTypesToLoad - What probeset types to load
 */
void ProbesetGenotypeStageEngine::loadLayout(ChipLayout &layout,
                   std::vector<const char *> &probesetNames,
                   probeidmap_t &killList,
                   bool justStats,
                   std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad) {

  std::set<const char *, Util::ltstr> probeSetsToLoad;
  vector<string> requiredCols;

  string cdfFile = getOpt("cdf-file");
  string spfFile = getOpt("spf-file");
  int computePrior = getOptInt("prior-size");
  bool listSample = getOptBool("list-sample");

  ///@todo populate probeSetsToLoad -- right now we are simply loading everything

  /* If we are sampling for a prior we need to load up all of the
     probesets, so empty the set which indicates this to openCdf(). */
  if(computePrior > 0 && !listSample) {
    probeSetsToLoad.clear();
  }

  vector<bool> probeSubset;  //leave empty as we simply load everything

  try {
    if(cdfFile != "") {
      Verbose::out(1, ToStr("Opening layout file: ") + cdfFile);
      if(!layout.openCdf(       cdfFile, 
                                probeSetsToLoad, 
                                &probesetNames,
                                probeSubset, 
                                "", 
                                killList, 
                                justStats, 
                                psTypesToLoad)) {
        Err::errAbort("Couldn't open layout file: " + cdfFile);
      }
    }
    else if(spfFile != "") {
      Verbose::out(1, ToStr("Opening layout file: ") + spfFile);
      layout.openSpf(spfFile, probeSetsToLoad, &probesetNames, probeSubset,
                           "", justStats, psTypesToLoad);
    }
    else {
      Err::errAbort("Must have either a cdf file or spf file.");
    }
  }
  catch(const Except &e) {
    Err::errAbort(e.what());
  }
  catch(...) {
    Err::errAbort("Unknown exception caught - terminating.");
  }
  Verbose::out(2, "Loaded " + ToStr(probesetNames.size()) + " probesets.");
}

/**
 * Open a text file and read in snps that are possibly haploid. Currently
 * only chrX snps are supported.
 *
 * @param SpecialSnps -
 * @param layout      -
 * @param permissive  -
 */
void ProbesetGenotypeStageEngine::fillInSpecialSnps(SpecialSnpMap &SpecialSnps,
                       ChipLayout &layout, bool permissive) {
  std::string specialType, chip;

  string specialSnps = getOpt("special-snps");

  readSpecialSnpsFile(&SpecialSnps, &chip, &specialType, specialSnps);

  // check the chip type
  affx::TsvFile tsv;
  string snpName;
  tsv.bind(0, "probeset_id", &snpName, TSV_BIND_REQUIRED);
  /* sanity checks */
  if(tsv.open(specialSnps) != TSV_OK)
    Err::errAbort("Couldn't open file: " + specialSnps + " to read.");
  tsv.close();

  // make sure the SNP in the specialSnps file is also in CDF
  if (!permissive) {
      for (SpecialSnpMap::const_iterator it = SpecialSnps.begin();
           it != SpecialSnps.end(); ++it) {
          const string &snpName = (*it).first;
          if (!layout.containsProbeSet(snpName)) {
              Err::errAbort("Can't find special snp: " + snpName + " in cdf file. Wrong special snp file?");
          }
      }
  }
}



/**
 * Read the prior from a file or generate one from the data.
 *
 * @param outFile - File to output prior if desired.
 * @param sample - Vector of probesets to generate prior from, if using data.
 * @param layout - Grouping and location of probesets.
 * @param iMart - Intensity mart of data for generating prior.
 * @param chipStream - Vector of modification objects to pass raw data through.
 * @param pmAdjuster - Probe specific adjustment.
 * @param qBRLMM - Quantification method for calling snps.
 */
void ProbesetGenotypeStageEngine::setPrior(string &outFile, vector<string> &sample,
              ChipLayout &layout, IntensityMart &iMart,
              vector<ChipStream *> &chipStream, PmAdjuster &pmAdjuster,
              QuantBRLMM *qBRLMM) {
  if(getOpt("read-priors-brlmm") != "") {
    Verbose::out(1, "Loading priors from: " + getOpt("read-priors-brlmm") );
    ClusterPrior prior = QuantBRLMM::loadPrior(getOpt("read-priors-brlmm"));
    if(getOptBool("write-prior")) {
      QuantBRLMM::writePrior(prior, outFile);
    }
    qBRLMM->setPrior(prior);
  }
  else if(getOptInt("prior-size") > 0) {
    Verbose::out(1, "Computing prior with: " + ToStr(sample.size()) + " initial snps.");
    ClusterPrior prior = qBRLMM->makePrior(sample, layout, iMart, chipStream, pmAdjuster);
    if(prior.snpCount < 10) {
      Err::errAbort(ToStr("Only ") + ToStr(prior.snpCount) +
                    " SNP clusters had enough information to generate prior. Please try using more cel files.");
    }
    Verbose::out(1, ToStr(prior.snpCount) + " SNP clusters used to generate prior.");
    if(getOptBool("write-prior")) {
      QuantBRLMM::writePrior(prior, outFile);
    }
    qBRLMM->setPrior(prior);
  }
}

/**
 * Make covariate normalizations from a file or generate one from the data.
 *
 * @param outFile - File to output normalization parameters if desired.
 * @param sample - Vector of probesets to generate normalizations from, if using data.
 * @param layout - Grouping and location of probesets.
 * @param iMart - Intensity mart of data for generating normalizers.
 * @param chipStream - Vector of modification objects to pass raw data through.
 * @param pmAdjuster - Probe specific adjustment.
 * @param qLabelZ - Quantification method for calling snps.
 */
void ProbesetGenotypeStageEngine::setNormLabelZ(string &outFile, vector<std::string> sample,
                   ChipLayout &layout, IntensityMart &iMart,
                   vector<ChipStream *> &chipStream, PmAdjuster &pmAdjuster,
                   QuantLabelZ *qLabelZ) {
  if(getOptInt("norm-size") > 0) {
    Verbose::out(1, "Computing normalizer with: " + ToStr(sample.size()) + " initial snps.");
    // here's where we do the magic
    // take the sample of snps
    // build a normalizer map of normalization across samples
    // 1 normalization predictor per sample
    vector<NormalizationPredictor> normMap = qLabelZ->makeNorm(sample, layout, iMart,
                                                               chipStream, pmAdjuster);
    if(getOptBool("write-norm")) {
      QuantLabelZ::writeNormMap(normMap, outFile);
    }
    qLabelZ->setNorm(normMap);
  }
}

/**
 * Construct an object for output that contains the program an
 * algorithm information.
 * @param as - AnalysisStream to get information from.
 * @param layout - Layout to get information from.
 * @param probesetNames - List of probeset names to output
 *
 * @return - Information for making CHP file.
 */
AnalysisInfo ProbesetGenotypeStageEngine::makeAnalysisInfo(AnalysisStream *as,
                              ChipLayout &layout,
                              const std::vector<const char*> &probesetNames,
                              GenderCalls *genders) {

  AnalysisInfo info;
  info.m_MaxPsNameLength = layout.getMaxNameLength();
  info.m_NumExpression = layout.getNumExpressionPSets();
  info.m_NumGenotyping = layout.getNumGenotypeingPSets();
  info.m_NumReporting = layout.getProbeSetCount();
  info.m_NumRows = layout.getXCount();
  info.m_NumCols = layout.getYCount();
  info.m_NumProbeSets = probesetNames.size();
  info.m_ProbeSetType = affxcdf::GenotypingProbeSetType;
  fillInAnalysisInfo(info, as);

  QuantGTypeMethod *method =
      dynamic_cast<QuantGTypeMethod *>(as->getQuantMethod());
  float min = method->getMinThresh();
  info.addParam("min-thresh",ToStr(min));
  float max = method->getMaxThresh();
  info.addParam("max-thresh",ToStr(max));

  info.addParam("gender_method_used",genders->getGenderName());

  // this is a "magic" value that works with GCOS/GTYPE according to Richard
  info.m_ProgID = "GDMTAnalysis.BRMBaseCall.1";

  // If CHP Output, then we need list of probeset names to output
  // Otherwise do not waste the memory
  if(getOptBool("cc-chp-output") || getOptBool("xda-chp-output")) 
    info.m_ProbesetNames = probesetNames;

  // Sanity check
  Err::check(info.m_ParamValues.size() == info.m_ParamNames.size(),
             "AnalysisInfo - Names and values out of sync.");
  return info;
}

void ProbesetGenotypeStageEngine::fillInAnalysisInfo(AnalysisInfo &info, QuantMethod *qMethod) {
    info.m_AlgVersion = qMethod->getVersion();
    //    info.m_AlgName = as->getName();
    //    info.addParam(prefix + "analysis-guid", as->getGuid());
    info.addParam("quantification-name", qMethod->getType());
    info.addParam("quantification-version", qMethod->getVersion());
    //    info.addParam(prefix + subprefix + "analysis-name", as->getName());

    info.addParam("quantification-scale", QuantMethod::scaleToTxt(qMethod->getScale()));
    info.addParam("quantification-type", QuantMethod::quantTypeToTxt(qMethod->getQuantType()));

}

void ProbesetGenotypeStageEngine::fillInAnalysisInfo(AnalysisInfo &info, std::string prefix) {
    ///@todo need to refactor to ensure consistency between meta info tracked in probeset-summarize
    ///      versus probeset-genotype
    ///@todo automatically pull state and options info?
    ///@todo sync up option name with label in chp header
    // FOR NOW -- PLEASE KEEP THIS IN SYNC WITH probeset-genotype
    vector<string> celFiles = getOptVector("cels");

    info.m_ProgramName = getOpt("program-name");
    info.m_ProgramVersion = getOpt("version-to-report");
    info.m_ProgramCompany = getOpt("program-company");
    info.m_ChipType = getOpt("chip-type");
    info.m_ProgID = "";
    info.m_ExecGuid = getOpt("exec-guid");
    info.m_AnalysisGuid = affxutil::Guid::GenerateNewGuid();

    // State and Execution info
    info.addParam("apt-engine", "ProbesetGenotypeStageEngine");
    info.addParam(prefix + "program-name", getOpt("program-name"));
    info.addParam(prefix + "command-line", getOpt("command-line"));
    info.addParam(prefix + "exec-guid", getOpt("exec-guid"));

    info.addParam(prefix + "time-str", getOpt("time-start"));
    info.addParam(prefix + "version", getOpt("version-to-report"));
    info.addParam(prefix + "cvs-id", getOpt("program-cvs-id"));
    info.addParam(prefix + "free-mem", getOpt("free-mem-at-start"));
    ///@todo add other state:
    ///@todo time-end state

    // Engine Options
    std::string subprefix = "opt-";
    // no do-ps-names
    info.addParam(prefix + subprefix + "chip-type", getOpt("chip-type"));
    info.addParam(prefix + subprefix + "probe-count", getOpt("probe-count"));
    info.addParam(prefix + subprefix + "force", getOpt("force"));
    // no precision
    info.addParam(prefix + subprefix + "out-dir", getOpt("out-dir"));
    // cc-expr-chp-out-dir
    info.addParam(prefix + subprefix + "cc-md-chp-out-dir", getOpt("cc-chp-out-dir"));
    info.addParam(prefix + subprefix + "xda-chp-out-dir", getOpt("xda-chp-out-dir"));
    info.addParam(prefix + subprefix + "cdf-file",       Fs::basename(getOpt("cdf-file")));
    // no spf-file 
    // no pgf-file
    // no clf-file
    // no bgp-file
    info.addParam(prefix + subprefix + "ps-list-file",   Fs::basename(getOpt("probeset-ids")));
    // no meta-ps-file
    // no-qc-groups-file
    info.addParam(prefix + subprefix + "kill-list",      Fs::basename(getOpt("kill-list")));
    info.addParam(prefix + subprefix + "temp-dir", getOpt("temp-dir"));
    info.addParam(prefix + subprefix + "use-disk", getOpt("use-disk"));
    info.addParam(prefix + subprefix + "diskCache", getOpt("disk-cache"));
    info.addParam(prefix + subprefix + "cel-count", ToStr(celFiles.size()));
    for(uint32_t i = 0; i < celFiles.size(); i++) {
      std::string paramName = prefix + subprefix + "cel-" + ToStr(i+1);
      info.addParam(paramName, Fs::basename(celFiles[i]));
    }

    // Algorithm Options
    info.addParam(prefix + subprefix + "analysis-name", info.m_AlgName);
    info.addParam(prefix + subprefix + "set-analysis-name", getOpt("set-analysis-name"));
    //    info.addParam(prefix + subprefix + "set-analysis-name", getOpt("set-analysis-name"));

    info.addParam(prefix + subprefix + "feat-effect-file", Fs::basename(getOpt("use-feat-eff")));
    info.addParam(prefix + subprefix + "target-sketch-file", Fs::basename(getOpt("target-sketch")));
    info.addParam(prefix + subprefix + "reference-profile-file",Fs::basename(getOpt("reference-profile")));
    info.addParam(prefix + subprefix + "do-residuals", getOpt("feat-details"));
    info.addParam(prefix + subprefix + "do-feature-effects", getOpt("feat-effects"));
    info.addParam(prefix + subprefix + "write-sketch", getOpt("write-sketch"));
    info.addParam(prefix + subprefix + "write-profile",getOpt("write-profile"));

    info.addParam(prefix + subprefix + "a5-write-sketch", getOpt("a5-sketch"));
    info.addParam(prefix + subprefix + "set-gender-method", getOpt("set-gender-method"));
    info.addParam(prefix + subprefix + "chrXProbeFile", Fs::basename(getOpt("chrX-probes")));
    info.addParam(prefix + subprefix + "chrYProbeFile", Fs::basename(getOpt("chrY-probes")));
    info.addParam(prefix + subprefix + "read-priors-brlmm", Fs::basename(getOpt("read-priors-brlmm")));
    info.addParam(prefix + subprefix + "initial-calls", ToStr("false"));
    info.addParam(prefix + subprefix + "output-summaries", getOpt("summaries"));
    info.addParam(prefix + subprefix + "list-sample", getOpt("list-sample"));
    info.addParam(prefix + subprefix + "write-prior", getOpt("write-prior"));
    info.addParam(prefix + subprefix + "chrX-file", Fs::basename(getOpt("chrX-snps")));
    info.addParam(prefix + subprefix + "chrX-snps", Fs::basename(getOpt("chrX-snps")));
    info.addParam(prefix + subprefix + "no-gender-force", getOpt("no-gender-force"));
    info.addParam(prefix + subprefix + "special-snps", Fs::basename(getOpt("special-snps")));
    info.addParam(prefix + subprefix + "compute-prior", ToStr(getOptInt("prior-size")));
    info.addParam(prefix + subprefix + "norm-size", getOpt("norm-size"));
    info.addParam(prefix + subprefix + "dm-thresh", getOpt("dm-thresh"));
    info.addParam(prefix + subprefix + "dm-het-mult", getOpt("dm-hetmult"));
    info.addParam(prefix + subprefix + "read-genders",        Fs::basename(getOpt("read-genders")));
    info.addParam(prefix + subprefix + "read-inbred",         Fs::basename(getOpt("read-inbred")));
    info.addParam(prefix + subprefix + "model-file-brlmm"   , Fs::basename(getOpt("read-models-brlmm")));
    info.addParam(prefix + subprefix + "model-file-brlmmp"  , Fs::basename(getOpt("read-models-brlmmp")));
    info.addParam(prefix + subprefix + "model-file-birdseed", Fs::basename(getOpt("read-models-birdseed")));
    info.addParam(prefix + subprefix + "genotypes-file"     , Fs::basename(getOpt("genotypes")));
    info.addParam(prefix + subprefix + "dm-calls-out", getOpt("dm-out"));
    info.addParam(prefix + subprefix + "qmethod-spec", getOpt("qmethod-spec"));
    info.addParam(prefix + subprefix + "genotype-only", ToStr(!getOptBool("all-types")));

    // Quantification Related Meta Info
    // Add user defined meta data parameters
    vector<pair<string, string> > metaData = getMetaDataDescription();
    for (int i = 0; i < metaData.size(); i++) {
      pair<string,string> p = metaData[i];
      info.addClientMetaInfo(p.first, p.second);
    }
}


/**
 * Construct an object for output that contains the program an
 * algorithm information.
 * @param as - AnalysisStream to get information from.
 * @param layout - Layout to get information from.
 * @param probesetNames - List of probeset names to output
 *
 * @return - Information for making CHP file.
 */
AnalysisInfo ProbesetGenotypeStageEngine::makeAnalysisInfo(
                                                           ChipLayout &layout,
                                                           const std::vector<const char*> &probesetNames,
                                                           GenderCalls *genders,
                                                           const std::string &analysisName) {

  AnalysisInfo info;
  info.m_MaxPsNameLength = layout.getMaxNameLength();
  info.m_NumExpression = layout.getNumExpressionPSets();
  info.m_NumGenotyping = layout.getNumGenotypeingPSets();
  info.m_NumReporting = layout.getNumGenotypeingPSets();
  info.m_NumRows = layout.getXCount();
  info.m_NumCols = layout.getYCount();
  info.m_AnalysisName = analysisName;
  info.m_AlgName = analysisName;
  info.m_NumProbeSets = layout.getProbeSetCount(); //probesetNames.size();
  info.m_ProbeSetType = affxcdf::GenotypingProbeSetType;
  // @todo refactor - what is output prefix?
  fillInAnalysisInfo(info, "apt-");

  // @todo refactor - make this work
  // QuantGTypeMethod *method =
  //     dynamic_cast<QuantGTypeMethod *>(as->getQuantMethod());
  // float min = method->getMinThresh();
  // info.addParam("min-thresh",ToStr(min));
  // float max = method->getMaxThresh();
  // info.addParam("max-thresh",ToStr(max));

  info.addParam("gender_method_used",genders->getGenderName());

  // this is a "magic" value that works with GCOS/GTYPE according to Richard
  info.m_ProgID = "GDMTAnalysis.BRMBaseCall.1";

  // If CHP Output, then we need list of probeset names to output
  // Otherwise do not waste the memory
  if(getOptBool("cc-chp-output") || getOptBool("xda-chp-output")) 
    info.m_ProbesetNames = probesetNames;

  // Sanity check
  Err::check(info.m_ParamValues.size() == info.m_ParamNames.size(),
             "AnalysisInfo - Names and values out of sync.");
  return info;
}


/**
 * Use a string specification to create an analysis stream.
 *
 * @param analysis - vector to be filled in with analysis streams.
 * @param asFactory - factory for converting strings to streams.
 * @param layout - Probeset specifications.
 * @param probesetNames - Names probesets in order for chp files
 * @param numProbeSets -
 */
void ProbesetGenotypeStageEngine::createAnalysisFromSpec(
                            std::vector<AnalysisStream *> &analysisStreams,
                            AnalysisStreamFactory &asFactory,
                            ChipLayout &layout,
                            std::vector<const char*> &probesetNames,
                            unsigned int numProbeSets) {
  // *Warning* People really did a cut and paste number on this
  // function and then things have evolved since then, be careful.

  string outDir = getOpt("out-dir");
  bool writeModels = getOptBool("write-models");
  string precompFeatureEff = getOpt("use-feat-eff");
  string qMethodSpec = getOpt("qmethod-spec");
  vector<string> celFiles = getOptVector("cels");

  // Create all the analysis paths requested.
  vector<string> analysisStrings = getOptVector("analysis");
  for(int i = 0; i < analysisStrings.size(); i++) {
    AnalysisStream *as = NULL;
    if(analysisStrings[i] == "")
        Err::errAbort("Invalid (blank) analysis specification");
    as = asFactory.constructGTypeAnalysisStream(analysisStrings[i],
                                                layout,
                                                m_stdMethods,
                                                getOpt("set-analysis-name"));
    QuantMethod *method = as->getQuantMethod();

    // BRLMM method detected, do stuff
    if(InstanceOf(method, QuantBRLMM)) {
      QuantBRLMM *qBRLMM = NULL;
      qBRLMM =  static_cast<QuantBRLMM *>(method);

      if(getOpt("read-models-brlmm") != "") {
        qBRLMM->readSnpModels(getOpt("read-models-brlmm"));
      }
      //
      {
        string outfile = Fs::join(outDir,as->getName() + ".snp-models");
        if (writeModels) {
          qBRLMM->setSnpModelOut(outfile,affx::TsvReport::FMT_TSV);
        }
        else if (getOptBool("a5-write-models")) {
            if (getOptBool("a5-write-models-use-global")) {
                if(m_a5_global_output_group == NULL)
                    Err::errAbort("--a5-write-models-use-global option given, but no global file. Must specify --a5-global-file.");
                qBRLMM->m_SnpModelsTsv.setA5SharedGroup(m_a5_global_output_group);
            }
            qBRLMM->setSnpModelOut(outfile,affx::TsvReport::FMT_A5);
        }
      }

      // Add quantification method.
      QuantMethodFactory factory(QuantMethodFactory::Expression);
      if(precompFeatureEff != "")
        factory.readPrecompFeatureEffectsFromFile(layout.getProbeCount(), precompFeatureEff, layout);
      QuantExprMethod *eMethod = factory.quantExprMethodForString(qMethodSpec, layout,
                                                                  QuantMethodFactory::Expression);
      // If we want to report allele specific summaries add a reporter for that.
      qBRLMM->setQuantExprMethod(eMethod);
    } // end of stuff for BRLMM detected

    if(InstanceOf(method, QuantBirdseed)) {
      QuantBirdseed *qBirdseed = NULL;
      qBirdseed =  static_cast<QuantBirdseed *>(method);

      // Add quantification method.
      QuantMethodFactory factory(QuantMethodFactory::Expression);
      if(precompFeatureEff != "")
        factory.readPrecompFeatureEffectsFromFile(layout.getProbeCount(), precompFeatureEff, layout);
      QuantExprMethod *eMethod = factory.quantExprMethodForString(qMethodSpec, layout,
                                                                  QuantMethodFactory::Expression);
      qBirdseed->setQuantExprMethod(eMethod);
      qBirdseed->setVerbosity(getOptInt("verbose"));
      //
      if(writeModels) {
        string outfile = Fs::join(outDir,as->getName()+".snp-models.txt");
        qBirdseed->setClusterOutFile(outfile);
      }
      else if (getOptBool("a5-write-models")) {
          ///@todo implement models in A5 for birdseed
          Err::errAbort("--a5-write-models for birdseed not implemented");
      }
    } // end of stuff for Birdseed detected

    // LabelZ || LabelZMulti
    if(InstanceOf(method, QuantLabelZ)) {
      QuantLabelZ *qLabelZ = NULL;
      qLabelZ =  static_cast<QuantLabelZ *>(method);
      
      //
      if (getOpt("read-models-brlmmp") != "") {
        qLabelZ->readSnpPriorMap(getOpt("read-models-brlmmp"));
      }
      else if (getOptBool("a5-models-input-global") || getOpt("a5-models-input-file") != "") {
        string groupName = "/";
        string dataName = as->getName() + ".snp-posteriors";
        if(getOpt("a5-models-input-name") != "") 
            dataName = getOpt("a5-models-input-name");
        if(getOpt("a5-models-input-group") != "") 
            groupName = getOpt("a5-models-input-group");
        else if(getOpt("a5-input-group") != "") 
            groupName = getOpt("a5-input-group");
        else if(getOpt("a5-group") != "") 
            groupName = getOpt("a5-group");
    
        if (getOptBool("a5-models-input-global")) {
            Verbose::out(1,"Loading models from global A5 file, group '" + groupName + "', data '" + dataName + "'");
            if(m_a5_global_input_file == NULL)
                Err::errAbort("--a5-models-input-global option given, but no global input file. Must specify --a5-global-file.");
            affx::File5_Group *a5group = m_a5_global_input_file->openGroup(groupName,affx::FILE5_OPEN);
            affx::File5_Tsv* tsv5 = a5group->openTsv(dataName,affx::FILE5_OPEN);
            qLabelZ->readSnpPriorMap_tsv5(tsv5);
            tsv5->close();
            Freez(tsv5);
            a5group->close();
            Freez(a5group);
        } else if (getOpt("a5-models-input-file")!="") {
            Verbose::out(1,"Loading models from '" + getOpt("a5-models-input-file") + "' A5 file, group '" + groupName + "', data '" + dataName + "'");
            affx::File5_File  *a5file = new affx::File5_File();
            a5file->open(getOpt("a5-models-input-file"),affx::FILE5_OPEN|affx::FILE5_RO);
            affx::File5_Group *a5group = a5file->openGroup(groupName,affx::FILE5_OPEN);
            affx::File5_Tsv* tsv5 = a5group->openTsv(dataName,affx::FILE5_OPEN);
            qLabelZ->readSnpPriorMap_tsv5(tsv5);
            tsv5->close();
            Freez(tsv5);
            a5group->close();
            Freez(a5group);
            a5file->close();
            Freez(a5file);
        }
      }
      // block for scoping outfile
      {
        string outfile = Fs::join(outDir,as->getName() + ".snp-posteriors");
        if(writeModels) {
          qLabelZ->writeSnpPosteriorTsv(outfile,affx::TsvReport::FMT_TSV);
        }
        else if (getOptBool("a5-write-models")) {
            if (getOptBool("a5-write-models-use-global")) {
                if(m_a5_global_output_group == NULL)
                    Err::errAbort("--a5-write-models-use-global option given, but no global file. Must specify --a5-global-file.");
                qLabelZ->m_SnpPosteriorTsv.setA5SharedGroup(m_a5_global_output_group);
            }
            qLabelZ->writeSnpPosteriorTsv(outfile,affx::TsvReport::FMT_A5);
        }
      }
      // block for scoping outfile
      {
        string outfile = Fs::join(outDir,as->getName() + ".normalized-summary");
        if (getOptBool("summaries")) {
          qLabelZ->writeNormSummaryTsv(outfile, celFiles,affx::TsvReport::FMT_TSV);
        }
        else if (getOptBool("a5-summaries")) {
            if (getOptBool("a5-summaries-use-global")) {
                if(m_a5_global_output_group == NULL) 
                    Err::errAbort("--a5-summaries-use-global option given, but no global file. Must specify --a5-global-file.");
                qLabelZ->m_NormSummaryTsv.setA5SharedGroup(m_a5_global_output_group);
            } 
            qLabelZ->writeNormSummaryTsv(outfile, celFiles,affx::TsvReport::FMT_A5);
        }
      }
      if (getOptBool("select-probes"))
        {
          string outfile =Fs::join(outDir,as->getName() + ".select-probes");
          qLabelZ->writeSnpProbeTsv(outfile,affx::TsvReport::FMT_TSV);
        }

      // Add quantification method.
      QuantMethodFactory factory(QuantMethodFactory::Expression);
      QuantMethodFactory::setupQuantMethodFactory(
            factory,
            m_a5_global_input_file,
            m_a5_global_input_group,
            m_a5_global_output_file,
            m_a5_global_output_group,
            getOptInt("probe-count"),
            getOpt("out-dir"),
            getOpt("use-feat-eff"),
            getOptBool("a5-feature-effects-input-global"),
            getOpt("a5-feature-effects-input-file"),
            getOpt("a5-feature-effects-input-name"),
            getOpt("a5-feature-effects-input-group"),
            getOpt("a5-input-group"),
            getOpt("a5-group"),
            getOpt("set-analysis-name"),
            getOpt("qmethod-spec"),
            layout
      );
      QuantExprMethod *eMethod = factory.quantExprMethodForString(qMethodSpec, layout, QuantMethodFactory::Expression);
      qLabelZ->setQuantExprMethod(eMethod);
    } // end of stuff for LabelZ detected

    analysisStreams.push_back(as);
  }
}

/**
 * Get the intial calls from the DmListener and give them to analysis
 * methods that need them.
 *
 * @param layout - Probeset definitions.
 * @param dmCaller - Object that does the DM calls.
 * @param knownGenotypes - Known genotype object to fill in.
 * @param haploidSnps - Snps on chrX.
 * @param analysis - Analysis methods to be done.
 * @param precompGTypes -
 */
void ProbesetGenotypeStageEngine::setInitialCalls(ChipLayout &layout, GenoSeed &seeds,
                     std::map<const char *,std::vector<affx::GType>, Util::ltstr > *knownGenotypes,
                     std::map<std::string,bool> &haploidSnps,
                     std::vector<AnalysisStream *> &analysis,
                     std::vector<const char *> toRunProbesets) {

  QuantMethod *qMethod = NULL;
  int genotypable_probeset_count = 0;

  // Set the initial calls in DM
  for(unsigned int i = 0; i < analysis.size(); i++) {
    AnalysisStream *as = analysis[i];
    qMethod = as->getQuantMethod();
    if(InstanceOf(qMethod, QuantBRLMM)) {
      // Get the initial genotype calls
      for(uint32_t i = 0; i < toRunProbesets.size(); i++) {
          const char *name = toRunProbesets[i];
        ProbeListPacked pList = layout.getProbeListByName(name);
        // 
        if (pList.isNull()) {
          Verbose::warn(1,std::string("ProbesetGenotypeStageEngine::setInitialCalls(): missing probeset named: '")+name+"'");
          continue;
        }
        //
        if((pList.get_type() == ProbeSet::GenoType && (pList.block_cnt() == 2 || pList.block_cnt() == 4)) || pList.get_type() == ProbeSet::Marker || pList.get_type() == ProbeSet::MultichannelMarker) {
            genotypable_probeset_count++;
            vector<GType> callVec = seeds.getGenoCalls(name);
            (*knownGenotypes)[name] = callVec;
        }
      }
      // Set the known genotypes here.
      QuantBRLMM *qBRLMM = NULL;
      qBRLMM =  static_cast<QuantBRLMM *>(qMethod);
      qBRLMM->setHaploidSnps(haploidSnps);
      qBRLMM->setKnownGenoTypes(*knownGenotypes);
      Verbose::out(1, "Setting " + ToStr(knownGenotypes->size()) + " seeds for " + ToStr(toRunProbesets.size()) + " probesets.");
    }
    else if(InstanceOf(qMethod,QuantLabelZ)) {
      // Get the initial genotype calls
      for(uint32_t i = 0; i < toRunProbesets.size(); i++) {
        const char *name = toRunProbesets[i];
        ProbeListPacked pList = layout.getProbeListByName(name);
        // didnt find it.
        if (pList.isNull()) {
          Verbose::warn(1,std::string("ProbesetGenotypeStageEngine::setInitialCalls(): missing probeset named: '")+name+"'");
          continue;
        }
        //
        if((pList.get_type() == ProbeSet::GenoType && (pList.block_cnt() == 2 || pList.block_cnt() == 4)) || pList.get_type() == ProbeSet::Marker || pList.get_type() == ProbeSet::MultichannelMarker) {
        //if ((pList.get_type() == ProbeSet::GenoType || pList.get_type() == ProbeSet::Marker) && pList.block_cnt() == 2) {
            genotypable_probeset_count++;
            if (seeds.checkGenoCallsName(name)) {
                vector<GType> callVec = seeds.getGenoCalls(name);
                (*knownGenotypes)[name] = callVec;
            }
        }
      }
      // Set the known genotypes here.
      QuantLabelZ *qLabelZ = NULL;
      qLabelZ = static_cast<QuantLabelZ *>(qMethod);
      qLabelZ->setKnownGenoTypes(knownGenotypes);
      Verbose::out(1, "Setting " + ToStr(knownGenotypes->size()) + " seeds for " + ToStr(genotypable_probeset_count) + " genotypable probesets.");
    }
  }
}


void ProbesetGenotypeStageEngine::checkSnpLists(const map<string,bool> &haploidSnps, ChipLayout *layout) {
  map<string, bool>::const_iterator hIter, hEnd = haploidSnps.end();
  for(hIter = haploidSnps.begin(); hIter != hEnd; hIter++) {
    if(!layout->containsProbeSet(hIter->first)) {
      Err::errAbort("Can't find snp: '" + hIter->first + "' in cdf file. Wrong chrx file?");
    }
  }
}

/**
 * When selecting probes to sample from we have a few different things to consider:
 * - Are we selecting from a certain type of probeset (i.e. genotyping).
 * - Are we selecting from all probesets.
 * - are we selecting from a specific list (i.e. --list-sample)
 *
 * If we're selecting from a specific list then life is easy, just use the list.
 * otherwise we have to check to make sure that we're selecting from the right type.
 *
 * @param toSample - We'll return the list to sample in here.
 * @param layout - Layout has information about probesets types.
 * @param probesetNames - Names has names of all the probesets.
 * @param toRunProbesetNames - This is the list of probesets that will be run.
 * @param psTypesToSample - What types of probesets should we be sampling from.
 */
void ProbesetGenotypeStageEngine::makeToSampleVector(vector<const char *> &toSample, 
                        map<std::string, int> &sampleMap,
                        ChipLayout *layout, const vector<const char *> &probesetNames,
                        const vector<const char *> &toRunProbesetNames,
                        std::set<affxcdf::GeneChipProbeSetType> &psTypesToSample) {
  int i;
  toSample.clear();
  sampleMap.clear();
  /* if the user specified --list-sample life is easy, just do what they say. */
  if(getOptBool("list-sample")) {
    toSample = toRunProbesetNames;
    for (i = 0; i < toRunProbesetNames.size(); i++) {
      sampleMap.insert(make_pair(string(toRunProbesetNames[i]), i));
    }
  }
  else {
    std::set<affxcdf::GeneChipProbeSetType>::const_iterator end = psTypesToSample.end();
    int pos = 0;
    for(i = 0; i < probesetNames.size(); i++) {
      bool keep =  psTypesToSample.empty() ||
        psTypesToSample.find((affxcdf::GeneChipProbeSetType)layout->getProbeSetType(i)) != end;
      if(keep) {
        toSample.push_back(probesetNames[i]);
        sampleMap.insert(make_pair(string(probesetNames[i]), pos++));
      }
    }
  }
}

/**
 * Figure out which snps are going to be used for a representation of all snps
 *
 * @param sampleSize -
 * @param layout - Description of chip format.
 * @param selectFrom -
 * @param SnpSampleNames - Vector to be filled in with pseudo random snp probeset names.
 */
void ProbesetGenotypeStageEngine::ObtainSnpSample(
                     ChipLayout *layout,
                     std::vector<const char *> &selectFrom,
                     std::vector<std::string> &SnpSampleNames) {
  int sampleSize = getOptInt("prior-size");
  vector<const char *> selected; // temporary for selection, no new memory allocated outside of stl
  if(sampleSize > selectFrom.size()) {
    Verbose::out(2, "Sample size larger than number of probesets. Resetting to: " + ToStr(selectFrom.size()));
    sampleSize = selectFrom.size();
    setOpt("prior-size", ToStr(sampleSize));
  }

  if(sampleSize < 100) {
    Verbose::out(1, "Warning!: " + ToStr(sampleSize) +
                 " seems like a small number of SNPs to represent all!");
  }

  SnpSampleNames.reserve(sampleSize);
  selected.resize(sampleSize);
  // time. Use the meaning of life, the universe and
  // everything...
  reproducible_random_sample(selectFrom.begin(), selectFrom.end(),
                             selected.begin(), sampleSize, 42);
  for(int i = 0; i < selected.size(); i++) {
    SnpSampleNames.push_back(selected[i]);
  }
}

/**
 * Sort probeset names according to given map
 * 
 * @param probesetsToSort - vector of probeset names
 * @param mapOrder - map of probeset names to relative positions
 */
void ProbesetGenotypeStageEngine::sortProbesetsByMapOrder(vector<std::string>::iterator begin, vector<std::string>::iterator end, map<std::string, int> &mapOrder) {
  if (begin < end - 1) {
    std::vector<std::string>::iterator right = end - 1;
    std::vector<std::string>::iterator left = begin;
    std::string temp;
    int pivot = mapOrder[*right];
    
    while (left < right) {
      while (left < right && mapOrder[*left] <= pivot) {
        left++;
      }
      while (left < right && pivot < mapOrder[*right]) {
        right--;
      }
      
      if (left < right) {
        temp = *left;
        *left = *right;
        *right = temp;
      }
    }
    
    sortProbesetsByMapOrder(begin, right, mapOrder);
    sortProbesetsByMapOrder(left, end, mapOrder);
  }  
}

/**
 * Add required reporters to analysis.
 *
 * @param qMethod - QuantGTypeMethod
 * @param as  - AnalysisStream
 * @param gender - genders
 */
void ProbesetGenotypeStageEngine::addReporters(
                  QuantGTypeMethod *qMethod,
                  AnalysisStream *as,
                  vector<ChipSummary *> &chipSummaries,
                  GenderCalls *gender,
                  std::set<const char *, Util::ltstr> *probeSetsToReport,
                  SpecialSnpMap &SpecialSnps) {

  AnalysisInfo info = as->getInfo();
  bool doCompact = getOptBool("file5-compact");
  string outDir = getOpt("out-dir");
  vector<string> celFiles = getOptVector("cels");

  // Add Quant Listener to Get Summary Stats
  QuantMethodGTypeChipSummary *runSummary = new QuantMethodGTypeChipSummary(qMethod->getType(), gender, SpecialSnps);
  as->addReporter(runSummary);

  // Add Expression Summary Stat
  QuantMethodGTypeExprChipSummary *exprSummary = new QuantMethodGTypeExprChipSummary(qMethod->getType());
  qMethod->addExprReporter(exprSummary);

  //// Add Chip Summary Report
  // TSV
  {
    QuantMethodRunReport *runReport = new QuantMethodRunReport(celFiles);
    if(getOpt("report-file") != "") {
        string dir = Fs::dirname(getOpt("report-file"));
        if(dir == "")
            dir = outDir;
        string name = Fs::basename(getOpt("report-file"));
        runReport->setDirPath(dir);
        runReport->setFilename(name);
    } else {
        runReport->setDirPath(outDir);
        runReport->setFilename(as->getName()+".report");
    }
    runReport->setFormat(TsvReport::FMT_TSV); // <-- TSV (!)
    //
    runReport->registerChipSummary(runSummary);
    runReport->registerChipSummary(exprSummary);
    for(int i = 0; i<chipSummaries.size(); i++)
      runReport->registerChipSummary(chipSummaries[i]);
    as->addReporter(runReport);
  }

  // Add Text Reporter  (calls & confs)
  if (getOptBool("table-output")) {
    QuantMethodGTypeReport *reporter = new QuantMethodGTypeReport(getOptBool("output-forced-calls"),getOptBool("output-context"));
    // @todo: hack -jhg
    reporter->setDirPath(outDir);
    reporter->setFilename("QuantMethodGTypeReport-TXT"); // not used -- for debugging.
    reporter->setFileprefix(as->getName());
    reporter->setFormat(TsvReport::FMT_TSV);
    //
    reporter->setPrecision(5);
    //
    reporter->registerProbeSetsToReport(probeSetsToReport);
    //
    as->addReporter(reporter);
  }

  //
  if (getOptBool("a5-calls")) {
    QuantMethodGTypeReport *reporter = new QuantMethodGTypeReport(getOptBool("output-forced-calls"),getOptBool("output-context"));
    reporter->setDirPath(outDir);
    reporter->setFileprefix(as->getName());
    reporter->setFilename("QuantMethodGTypeReport-A5"); // not used -- for debugging.
    if(getOptBool("a5-calls-use-global")) {
        if(m_a5_global_output_group == NULL)
            Err::errAbort("--a5-calls-use-global option given, but no global file. Must specify --a5-global-file.");
        reporter->setA5SharedGroup(m_a5_global_output_group);
    }
    reporter->setFormat(affx::TsvReport::FMT_A5);
    if(doCompact) {
      reporter->setCompactFile5Format(info.m_MaxPsNameLength);
    }
    //
    reporter->registerProbeSetsToReport(probeSetsToReport);
    //
    as->addReporter(reporter);
  }

  // Add AGCC CHP Reporter
  if(getOptBool("cc-chp-output")) {
      string ccchpDir;
      /* If chp directory not specified use  output/cc-chp/ by default. */
      if(getOpt("cc-chp-out-dir") == "") {
        ccchpDir = Fs::join(outDir,"cc-chp");
      }
      else {
        ccchpDir = getOpt("cc-chp-out-dir");
      }
      Util::chompLastIfSep(ccchpDir);
      QuantMethodMultiDataCCCHPReport *ccchpReporter = new QuantMethodMultiDataCCCHPReport(ccchpDir);
      vector<string> resultFiles = getOptVector("result-files");
      if(resultFiles.size() > 0)
          ccchpReporter->setChpFileNames(resultFiles);
      ccchpReporter->registerChipSummary(runSummary);
      ccchpReporter->registerChipSummary(exprSummary);
      for(int i = 0; i<chipSummaries.size(); i++)
        ccchpReporter->registerChipSummary(chipSummaries[i]);
      as->addReporter(ccchpReporter);
  }

  // Add XDA CHP Reporter
  if(getOptBool("xda-chp-output")) {
      string chpDir;
      /* If chp directory not specified use  output/chp/ by default. */
      if (getOpt("xda-chp-out-dir") == "") {
        chpDir = Fs::join(outDir,"chp");
      }
      else {
        chpDir = getOpt("xda-chp-out-dir");
      }
      Util::chompLastIfSep(chpDir);
      QuantMethodGTypeCHPReport *chpReporter = new QuantMethodGTypeCHPReport(
              info, chpDir, qMethod->getType());
      vector<string> resultFiles = getOptVector("result-files");
      if(resultFiles.size() > 0)
          chpReporter->setChpFileNames(resultFiles);
      as->addReporter(chpReporter);
  }
/* Avoid compiling these into windows builds until GTC team has been 
   informed and given time to update their builds as new libraries are requred. */
#ifndef WIN32
  if (getOptBool("sqlite-output")) {
    string databaseFile = Fs::join(outDir,as->getName() + ".sqlite");
    QuantMethodSqlGTypeReport *reporter = new QuantMethodSqlGTypeReport(databaseFile, "call", "confidence", info.m_MaxPsNameLength);
    as->addReporter(reporter);
    QuantMethodSqlExprReport *sumReporter = new QuantMethodSqlExprReport(databaseFile + ".summary.db", "summary", info.m_MaxPsNameLength);
    sumReporter->addStdHeaders(sumReporter, getOpt("exec-guid"), as->getGuid(), getOpt("time-start"),
                               getOpt("command-line"), getOpt("version-to-report"), info);
    if(InstanceOf(qMethod,QuantLabelZ)) {
      static_cast<QuantLabelZ *>(qMethod)->addExprReporter(sumReporter);
    } 
    else if (InstanceOf(qMethod,QuantBRLMM)) {
      static_cast<QuantBRLMM *>(qMethod)->addExprReporter(sumReporter);
    } 
    else if (InstanceOf(qMethod,QuantBirdseed)) {
      static_cast<QuantBirdseed *>(qMethod)->addExprReporter(sumReporter);
    }
                                                   
  }
#endif
  // Add txt expression reporter
  if (getOptBool("summaries") || getOptBool("feat-effects") || getOptBool("feat-details")) {
    addReporters_Expression(affx::TsvReport::FMT_TSV,info, qMethod,as,chipSummaries,gender,probeSetsToReport, doCompact);
  }
  if (getOptBool("a5-summaries") || getOptBool("a5-feature-effects") || getOptBool("a5-feature-details")) {
    addReporters_Expression(affx::TsvReport::FMT_A5,info, qMethod,as,chipSummaries,gender,probeSetsToReport, doCompact);
  }
}

/**
 * Add required reporters to analysis.
 *
 * @param qMethod - QuantGTypeMethod
 * @param as  - AnalysisStream
 * @param gender - genders
 */
void ProbesetGenotypeStageEngine::makeReporters(
                                                const std::string &qMethodName,
                                                std::vector<QuantMethodReport *> &reporters,
                                                std::vector<QuantMethodReport *> &exprReporters,
                                                AnalysisInfo &info,
                                                vector<ChipSummary *> &chipSummaries,
                                                GenderCalls *gender,
                                                std::set<const char *, Util::ltstr> *probeSetsToReport,
                                                SpecialSnpMap &SpecialSnps) {

  bool doCompact = getOptBool("file5-compact");
  string outDir = getOpt("out-dir");
  vector<string> celFiles = getOptVector("cels");

  // Add Quant Listener to Get Summary Stats
  QuantMethodGTypeChipSummary *runSummary = new QuantMethodGTypeChipSummary(qMethodName, gender, SpecialSnps);
  reporters.push_back(runSummary);

  // Add Expression Summary Stat
  QuantMethodGTypeExprChipSummary *exprSummary = new QuantMethodGTypeExprChipSummary(qMethodName);
  exprReporters.push_back(exprSummary);

  //// Add Chip Summary Report
  // TSV
  {
    QuantMethodRunReport *runReport = new QuantMethodRunReport(celFiles);
    if(getOpt("report-file") != "") {
        string dir = Fs::dirname(getOpt("report-file"));
        if(dir == "")
            dir = outDir;
        string name = Fs::basename(getOpt("report-file"));
        runReport->setDirPath(dir);
        runReport->setFilename(name);
    } else {
        runReport->setDirPath(outDir);
        runReport->setFilename(info.m_AnalysisName+".report");
    }
    runReport->setFormat(TsvReport::FMT_TSV); // <-- TSV (!)
    //
    runReport->registerChipSummary(runSummary);
    runReport->registerChipSummary(exprSummary);
    for(int i = 0; i<chipSummaries.size(); i++)
      runReport->registerChipSummary(chipSummaries[i]);
    reporters.push_back(runReport);
  }

  // Add Text Reporter  (calls & confs)
  if (getOptBool("table-output")) {
    QuantMethodGTypeReport *reporter = new QuantMethodGTypeReport(getOptBool("output-forced-calls"),getOptBool("output-context"));
    // @todo: hack -jhg
    reporter->setDirPath(outDir);
    reporter->setFilename("QuantMethodGTypeReport-TXT"); // not used -- for debugging.
    reporter->setFileprefix(info.m_AnalysisName);
    reporter->setFormat(TsvReport::FMT_TSV);
    //
    reporter->setPrecision(5);
    //
    reporter->registerProbeSetsToReport(probeSetsToReport);
    //
    reporters.push_back(reporter);
  }

  //
  if (getOptBool("a5-calls")) {
    QuantMethodGTypeReport *reporter = new QuantMethodGTypeReport(getOptBool("output-forced-calls"),getOptBool("output-context"));
    reporter->setDirPath(outDir);
    reporter->setFileprefix(info.m_AnalysisName);
    reporter->setFilename("QuantMethodGTypeReport-A5"); // not used -- for debugging.
    if(getOptBool("a5-calls-use-global")) {
        if(m_a5_global_output_group == NULL)
            Err::errAbort("--a5-calls-use-global option given, but no global file. Must specify --a5-global-file.");
        reporter->setA5SharedGroup(m_a5_global_output_group);
    }
    reporter->setFormat(affx::TsvReport::FMT_A5);
    if(doCompact) {
      reporter->setCompactFile5Format(info.m_MaxPsNameLength);
    }
    //
    reporter->registerProbeSetsToReport(probeSetsToReport);
    //
    reporters.push_back(reporter);
  }

  // Add AGCC CHP Reporter
  if(getOptBool("cc-chp-output")) {
      string ccchpDir;
      /* If chp directory not specified use  output/cc-chp/ by default. */
      if (getOpt("cc-chp-out-dir") == "") {
        ccchpDir = Fs::join(outDir,"cc-chp");
      }
      else {
        ccchpDir = getOpt("cc-chp-out-dir");
      }
      Util::chompLastIfSep(ccchpDir);
      QuantMethodMultiDataCCCHPReport *ccchpReporter = new QuantMethodMultiDataCCCHPReport(ccchpDir);
      vector<string> resultFiles = getOptVector("result-files");
      if(resultFiles.size() > 0)
          ccchpReporter->setChpFileNames(resultFiles);
      ccchpReporter->registerChipSummary(runSummary);
      ccchpReporter->registerChipSummary(exprSummary);
      for(int i = 0; i<chipSummaries.size(); i++)
        ccchpReporter->registerChipSummary(chipSummaries[i]);
      reporters.push_back(ccchpReporter);
  }

  // Add XDA CHP Reporter
  if(getOptBool("xda-chp-output")) {
      string chpDir;
      /* If chp directory not specified use  output/chp/ by default. */
      if(getOpt("xda-chp-out-dir") == "")
        chpDir = Fs::join(outDir,"chp");
      else
        chpDir = getOpt("xda-chp-out-dir");
      Util::chompLastIfSep(chpDir);
      QuantMethodGTypeCHPReport *chpReporter = new QuantMethodGTypeCHPReport(
              info, chpDir, qMethodName);
      vector<string> resultFiles = getOptVector("result-files");
      if(resultFiles.size() > 0)
          chpReporter->setChpFileNames(resultFiles);
      reporters.push_back(chpReporter);
  }
/* Avoid compiling these into windows builds until GTC team has been 
   informed and given time to update their builds as new libraries are requred. */
#ifndef WIN32
  if(getOptBool("sqlite-output")) {
    string databaseFile = Fs::join(outDir,info.m_AnalysisName+".sqlite");
    QuantMethodSqlGTypeReport *reporter = new QuantMethodSqlGTypeReport(databaseFile, "call", "confidence", info.m_MaxPsNameLength);
    reporters.push_back(reporter);
    QuantMethodSqlExprReport *sumReporter = new QuantMethodSqlExprReport(databaseFile + ".summary.db", "summary", info.m_MaxPsNameLength);
    sumReporter->addStdHeaders(sumReporter, getOpt("exec-guid"), info.m_AnalysisGuid, getOpt("time-start"),
                               getOpt("command-line"), getOpt("version-to-report"), info);
    exprReporters.push_back(sumReporter);
  }
#endif
  // Add txt expression reporter
  if (getOptBool("summaries") || getOptBool("feat-effects") || getOptBool("feat-details")) {
    addReporters_Expression(affx::TsvReport::FMT_TSV,info, exprReporters,chipSummaries,gender,probeSetsToReport, doCompact);
  }
  if (getOptBool("a5-summaries") || getOptBool("a5-feature-effects") || getOptBool("a5-feature-details")) {
    addReporters_Expression(affx::TsvReport::FMT_A5,info, exprReporters,chipSummaries,gender,probeSetsToReport, doCompact);
  }
}


/**
 * Any methods that must have seed genotype calls?  BRLMM
 *
 * @param analysis - which methods are called
 */
int ProbesetGenotypeStageEngine::DetectNeedForSeeds(vector<AnalysisStream *> analysis)
{
  int needSeeds=0;
  QuantMethod *qMethod = NULL;
  for(unsigned int asIx = 0; asIx < analysis.size(); asIx++) {
    AnalysisStream *as = analysis[asIx];
    qMethod = as->getQuantMethod();
    if(InstanceOf(qMethod, QuantBRLMM))
      needSeeds = 1;
  }
  return(needSeeds);
}


/**
 *  Guess what: does setup for analysis methods first time through the loop
 *  gets samples of snps for various purposes
 *  calls genders
 *  builds priors, normalization functions
 *
 * @param analysis - what methods need setting up
 * @param layout - this particular subset of the chip
 * @param iMart - intensities from cel files
 * @param priorSnpNames - names of snp sample for prior (BRLMM)
 * @param normSnpNames - names of snps sampled for normalization (LabelZ = BRLMM-p)
 * @param fileNames - names of celfiles
 * @param genders - genders to pass to quant methods
 * @param inbred - sample covariate to pass to quant methods, like genders
 * @param haploidSnps - chrX snps which get special treatment
 * @param specialSnps -
 */
void ProbesetGenotypeStageEngine::DoSetup(std::vector<AnalysisStream *> &analysis,
             ChipLayout &layout,
             IntensityMart &iMart,
//              IntensityMart &priorIntensityMart,
//              IntensityMart &normIntensityMart,
             vector<string> &priorSnpNames,
             vector<string> &normSnpNames,
             vector<string> &fileNames,
             GenderCalls *genders,
             InbredStatus *inbred,
             map<string,bool> &haploidSnps,
             SpecialSnpMap &SpecialSnps,
             CopyNumberMap &copyNumberMap,
             vector<ChipSummary *> &chipSummaries,
             const vector<const char *> &probesetNames,
             std::set<const char *, Util::ltstr> *probeSetsToReport) {

  QuantMethod *qMethod = NULL;

  string outDir = getOpt("out-dir");
  vector<string> celFiles = getOptVector("cels");

  // refactor out the hack that sets up everything
  // now deal with the setup for streams
  // including learning about samples/snps

  for(unsigned int asIx = 0; asIx < analysis.size(); asIx++) {
    AnalysisStream *as = analysis[asIx];
    qMethod = as->getQuantMethod();

    AnalysisInfo info = makeAnalysisInfo(as, layout, probesetNames, genders);
    as->setInfo(info);

    // If we are learning a prior this is the time to do it.
    vector<affx::Gender> genderVec = genders->getGenders();
    vector<double> inbredVec = inbred->getInbredStatus();

    // if(InstanceOf(qMethod, QuantBRLMM)) {
    //   QuantBRLMM *qBRLMM = NULL;
    //   qBRLMM =  static_cast<QuantBRLMM *>(qMethod);
    //   qBRLMM->setGenders(genderVec);
    //   qBRLMM->setHaploidSnps(haploidSnps);
    //   /// @todo refactor dmCallsOut out of QuantBRLMM and into engine
    //   if(getOptBool("dm-out")) {
    //     qBRLMM->setDmCallsOut(outDir + PATH_SEPARATOR + "dm.calls", celFiles);
    //   }
    //   addReporters(qBRLMM, as, chipSummaries, genders, probeSetsToReport);
    //   string outPrefix = outDir + PATH_SEPARATOR + ToStr(as->getName()) + ".prior.txt";
    //   setPrior(outPrefix, priorSnpNames, layout, iMart,
    //            *(as->getChipStream()), *(as->getPmAdjuster()), qBRLMM);

    // } // end BRLMM detected

    // if(InstanceOf(qMethod, QuantBirdseed)) {
    //   QuantBirdseed *qBirdseed = NULL;
    //   qBirdseed =  static_cast<QuantBirdseed *>(qMethod);
    //   qBirdseed->setGenders(genderVec);
    //   addReporters(qBirdseed, as, chipSummaries, genders, probeSetsToReport);
    // } // end Birdseed detected
    
    //  Note that QuantLabelZMulti will also satisfy the "InstanceOf" check on next line.
    if (InstanceOf(qMethod,QuantLabelZ)){
      QuantLabelZ *qLabelZ = NULL;
      qLabelZ = static_cast<QuantLabelZ *>(qMethod);
      // just report QC output
      if (!SpecialSnps.empty())
          qLabelZ->setSpecialSnps(SpecialSnps);
      if(!haploidSnps.empty())
        qLabelZ->setHaploidSnps(haploidSnps);
      qLabelZ->setGenders(genderVec);
      // add sample covariate inbreeding penalty here
      // this is highly indirect at this point, but perhaps there's a general model for sample covariates
      qLabelZ->setInbredHetPenalty(inbredVec);
      // send genders to the summary reporter for all chip files
      addReporters(qLabelZ, as, chipSummaries, genders, probeSetsToReport, SpecialSnps);
      // set normalization if required
      // that is, precompute the normalization parameters based on
      // covariates for each sample
      string outPrefix = Fs::join(outDir,as->getName()+".norm.txt");
      if (getOptInt("norm-size")>0)
        setNormLabelZ(outPrefix, normSnpNames, layout, iMart,
                      *(as->getChipStream()), *(as->getPmAdjuster()), qLabelZ);

    } // end labelZ detected

    // QuantLabelZMulti already partially set as it is also a QuantLabelZ
    if (InstanceOf(qMethod,QuantLabelZMulti)){
      QuantLabelZMulti *qLabelZMulti = NULL;
      qLabelZMulti = dynamic_cast<QuantLabelZMulti *>(qMethod);
      ///@todo we should really make this generic to all QuantMethods
      qLabelZMulti->registerProbeSetsToReport(probeSetsToReport);
      if(!copyNumberMap.empty()) 
        qLabelZMulti->setCopyNumberMap(copyNumberMap);
    }

    as->addStdHeaders(getOpt("exec-guid"),
                      getOpt("time-start"),
                      getOpt("command-line"),
                      getOpt("version-to-report"),
                      info);
    // this opens the file
    as->prepare(iMart);
  }
}

void ProbesetGenotypeStageEngine::setIterationData(ChipLayout &layout, 
                      vector<AnalysisStream *> &analysis,
                      vector<string> &priorSnpNames, vector<string> &normSnpNames,
                      vector<const char *> &toRunProbesets,
                      SpecialSnpMap &specialSnps) {
  QuantMethod *qMethod = NULL;
  // Set the initial data that doesn't need seeds
  for(unsigned int i = 0; i < analysis.size(); i++) {
    AnalysisStream *as = analysis[i];
    qMethod = as->getQuantMethod();
    if(InstanceOf(qMethod, QuantBirdseed)) {
      QuantBirdseed *qBirdseed = NULL;
      qBirdseed =  static_cast<QuantBirdseed *>(qMethod);
      // Load up just the probeset priors we need for this iteration.
      if(getOpt("read-models-birdseed") != "") {
        set<const char *, Util::ltstr> toLoad;
        for(int i = 0; i < priorSnpNames.size(); i++) {
          toLoad.insert(priorSnpNames[i].c_str());
        }
        for(int i = 0; i < normSnpNames.size(); i++) {
          toLoad.insert(normSnpNames[i].c_str());
        }
        // Always need to load up the snps to be processed this time through.
        for(int i = 0; i < toRunProbesets.size(); i++) {
          toLoad.insert(toRunProbesets[i]);
        }
        qBirdseed->setPriorFile(getOpt("read-models-birdseed"), getOpt("chip-type"), &toLoad, specialSnps);
      }
    }
  }
}

void ProbesetGenotypeStageEngine::fillInDesiredOrder(std::vector<int> &desiredOrder) {
    vector<bool> probeSeen(m_ChipLayout->getProbeCount(), false);
    int psCount = m_ChipLayout->getProbeSetCount();
    APT_ERR_ASSERT(psCount > 0, "Must have positive number of probesets to count.");
    // @todo - this should be a function rather than inline...
    for(int i = 0; i < psCount ; i++) {
        ProbeListPacked pl = m_ChipLayout->getProbeListAtIndex(i);
        if (!pl.isNull()) {
            for(int pIx = 0; pIx < pl.probe_cnt(); pIx++) {
                if(!probeSeen[pl.get_probeId(pIx)]) {
                    desiredOrder.push_back(pl.get_probeId(pIx));
                    probeSeen[pl.get_probeId(pIx)]  = true;
                }
            }
        }
    }
    
    for (int i = 0; i < probeSeen.size(); i++) {
        if (!probeSeen[i]) {
            desiredOrder.push_back(i);
        } 
    }

}

void  ProbesetGenotypeStageEngine::fillInSelectedGenders(GenderCalls **selectedGenders, 
                                                         std::string &analysisName, 
                                                         std::map<std::string, GenderCalls *> &genderCallsMap,
                                                         NoGenderCalls *noGenderCalls) {

    if(analysisName.find("birdseed-v2") != string::npos || analysisName.find("birdseed-dev") != string::npos) {
        if(genderCallsMap.find("cn-probe-chrXY-ratio")!=genderCallsMap.end()) {
            *selectedGenders = genderCallsMap["cn-probe-chrXY-ratio"];
        } else if(genderCallsMap.find("em-cluster-chrX-het-contrast")!=genderCallsMap.end()) {
            *selectedGenders = genderCallsMap["em-cluster-chrX-het-contrast"];
        }
    } 
    else if(analysisName.find("birdseed") != string::npos) {
        if(genderCallsMap.find("em-cluster-chrX-het-contrast")!=genderCallsMap.end()) {
            *selectedGenders = genderCallsMap["em-cluster-chrX-het-contrast"];
        }
    }  
    else if(analysisName.find("brlmm-p") != string::npos) {
        if(genderCallsMap.find("cn-probe-chrXY-ratio")!=genderCallsMap.end()) {
            *selectedGenders = genderCallsMap["cn-probe-chrXY-ratio"];
        } else if(genderCallsMap.find("em-cluster-chrX-het-contrast")!=genderCallsMap.end()) {
            *selectedGenders = genderCallsMap["em-cluster-chrX-het-contrast"];
        }
    }
    else if(analysisName.find("brlmm") != string::npos) {
        if(genderCallsMap.find("dm-chrX-het-rate")!=genderCallsMap.end()) {
            *selectedGenders = genderCallsMap["dm-chrX-het-rate"];
        } else if(genderCallsMap.find("supplied-genotypes-chrX-het-rate")!=genderCallsMap.end()) {
            *selectedGenders = genderCallsMap["supplied-genotypes-chrX-het-rate"];
        }
    }
    
    // user specified gender method trumps all
    string genderMethod = getOpt("set-gender-method");
    if(genderMethod != "") {
        if(genderMethod == "none") {
            *selectedGenders = noGenderCalls;
        } else if(genderCallsMap.find(genderMethod)!=genderCallsMap.end()) {
            *selectedGenders = genderCallsMap[genderMethod];
        } else {
            Err::errAbort("No " + genderMethod + " genders to use, but " + genderMethod + " genders requested with --gender-method.");
        }
    }

}

void ProbesetGenotypeStageEngine::fillInChipLayout() {
    std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad; // What types of probesets are we interested in analyzing.
    probeidmap_t killList;                 // Probes in cdf that we don't want to use.
    colrow_t numRows = Convert::toUnsignedInt(getOpt("num-rows"));
    colrow_t numCols = Convert::toUnsignedInt(getOpt("num-cols"));
    if(m_FreeChipLayout) {
        
        // Does the user want to toss out some probes?
        if(getOpt("kill-list") != "") {
            ChipLayout::fillInKillList(getOpt("kill-list"), killList, numRows, numCols);
        }
        if(!getOptBool("all-types")) {
            psTypesToLoad.insert(affxcdf::MarkerProbeSetType);
            psTypesToLoad.insert(affxcdf::GenotypingProbeSetType);
            psTypesToLoad.insert(affxcdf::MultichannelMarkerProbeSetType);
        }
        
        Freez(m_ChipLayout);
        m_ChipLayout = new ChipLayout;
        m_ChipLayout->setNeedMismatch(getOptBool("need-mismatch"));
        m_ChipLayout->setNeedGc(getOptBool("need-gc"));
        m_ChipLayout->setNeedPmAlleleMatch(getOptBool("need-allelematch"));
        
        loadLayout(  *m_ChipLayout, 
                     m_ProbesetNames,    // filled in and used later
                     killList, 
                     false, 
                     psTypesToLoad
                     );
    }
}

/**
   This is the "main()" equivalent of the engine. This function has
   two main phases: 1) Creation an initialization of AnalysisStream
   and associated objects. 2) Looping through the list of cdf
   probesets to do the genotyping calls.

   Phase 2 is definitely where most of the time and cpu cycles are
   spent. The list of things that happens for each iteration of phase
   2 is:
   - Read in cel files only keeping data needed for current iteration in RAM.
   - While reading cel files make DM calls.
   - Once cel files are read loop through probesets making BRLMM genotype
   calls.

   Special things that happen in the first iteration of phase 2:
   - Calculate prior if desired.
   - Determine gender using DM calls if desired.
   - Learn chipstream parameters necessary.
   - Attach any reporters to the analysis streams.

   Top of the call stack... This function is getting to be a real
   beast, time for some refactoring...
*/
void ProbesetGenotypeStageEngine::runImp() {
  // *Warning* This function has evolved into quite a behemoth with a
  // series of feature requests and no time to refactor it, be
  // careful with modifications.
  // probeset names act as a common pool for a number of different
  // things so we don't have to copy the relatively large (in RAM)
  // probe names over and over. They are known to be used for the
  // chp file writers, meta probesets, and prior selection but may
  // also be used for other things as time goes on. Change or delete
  // them prematurely at your own risk.
  std::vector<const char *> toRunProbesets; /// names of probesets to compute/process
  bool needPrior = true;                 // Are we computing a prior or is it being read from file?
  vector<string> priorSnpNames;          // Names of snps that will be used to estimate prior.
  //  bool needNorm = false;                 // are we computing norm? - this variable is not used?
  vector<string> normSnpNames;           // Names of snps that will be used to estimate normalizers.
  IntensityMart *iMart = NULL;
  map<string,bool> haploidSnps;          // chrx snps, used for gender and estimated with different model
  SpecialSnpMap SpecialSnps;             // all special snps
  CopyNumberMap copyNumberMap;           // A map of names of markers to a vector of copy Numbers, indexed by sample. 
  std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad; // What types of probesets are we interested in analyzing.
  std::set<affxcdf::GeneChipProbeSetType> psTypesToSample; // What types of probesets are we interested in sampling.
  CelStatListener *celStats = NULL;
  std::set<const char *, Util::ltstr> probeSetsToReport; // Probe sets to report on

  // seed and gender stuff
  //  GenoSeed *seeds = NULL;                // The seed genotypes to use for determining gender and seeding brlmm
  GenoSeedTxt *precompSeeds = NULL;      // Precomputed seeds that were read in.
  DmListener *dmCaller = NULL;

  // gender stuff
  map<string, GenderCalls *> genderCallsMap; /// Map of gender call methods available
  FileGenders *fileGenderReader = NULL;
  EmGenderCelListener *emCelListener = NULL;
  CnProbeGenderCelListener *cnProbeCelListener = NULL;
  NoGenderCalls *noGenderCalls = NULL;

  // inbred sample covariate stuff
  map<string, InbredStatus *> InbredMap; /// Map of Inbred status methods [only one currently, but who knows?]
  //  FileInbred *fileInbredReader = NULL;
  NoInbredStatus *noInbredReader = NULL;
  CelReader reader;
  PsBoard *board = new PsBoard(); 
  board->setOptions(this);
  // make this all one big try catch to make sure destructors are
  // called if exceptions thrown.
  try { // outer try for handling the exception.
    try { // inner try for memory clean up.
        vector<string> celFiles = getOptVector("cels");
        vector<string> analysisStrings = getOptVector("analysis");

        // set up null covariates in case nothing is supplied
        // sample covariates: gender, inbreeding status
        noGenderCalls = new NoGenderCalls(celFiles.size());
        noInbredReader = new NoInbredStatus(celFiles.size());

        // mem tracking stuff
        std::vector<ChipSummary *> chipSummaries(0);  // ChipSummaries to add to reporters

        Verbose::out(1, "Beginning analysis of " + ToStr(celFiles.size()) + " cel files.");

        /* Are we doing a subset of all the probe sets? */
        string psListFile = getOpt("probeset-ids");
        string readPriorsBrlmm = getOpt("read-priors-brlmm");
        //        int computePrior = getOptInt("prior-size");
        //        bool listSample = getOptBool("list-sample");
        /* Are we only reporting a subset of those processed */
        string psReportListFile = getOpt("probeset-ids-reported");

        // read in the layout of probesets on the chip.
        fillInChipLayout(); // Fills in m_ChipLayout and m_ProbesetNames

        vector<string> fileNames = celFiles;
        board->set("chiplayout", m_ChipLayout);
        // Load up marker by sample specific cn info
        string genotype_markers_cn_file = getOpt("genotype-markers-cn-file");
        if(genotype_markers_cn_file != "") {
            readCopyNumberFile(copyNumberMap, genotype_markers_cn_file);
        }

        // We need to handle chrX snps separately for some algorithms, read in those snps if specified.
        // special snps trumps chrX snps specification
        if (getOpt("special-snps") != ""){
            fillInSpecialSnps(SpecialSnps, *m_ChipLayout, true);
            haploidFromSpecial(SpecialSnps,haploidSnps);
        }        if(getOpt("chrX-probes") != "") { // @todo
            /// @todo should we require a minimum number of X and Y probes?
            cnProbeCelListener = new CnProbeGenderCelListener(getOpt("chrX-probes"), getOpt("chrY-probes"), 
                                                              getOptDouble("female-thresh"), getOptDouble("male-thresh"));
            reader.registerCelListener(cnProbeCelListener);
            chipSummaries.push_back(cnProbeCelListener);
            genderCallsMap[cnProbeCelListener->getGenderName()] = cnProbeCelListener;
        }

        needPrior = false;
        string diskDir = getOpt("temp-dir");

        // @todo refactor create DataStore here...
        std::map<const char *,std::vector<affx::GType>, Util::ltstr > knownGenotypes;
        map<string,bool>::iterator start, end = haploidSnps.end();
        vector<ProbeListPacked> haploidProbeSets;
        for(start = haploidSnps.begin(); start != end; start++) {
            ProbeListPacked ps = m_ChipLayout->getProbeListByName(start->first);
            if(!ps.isNull()) {
              haploidProbeSets.push_back(ps);
            }
            //else
                // Verbose::out(1,"Unable to find probeset " + start->first + " in layout for haploidSnps. Ignoring.");
        }
        if(!haploidProbeSets.empty()) {
            /// @todo is 1 enough? probably not, should we require a minimum
            /// number of haploid probe sets?

            // Hack to turn this off for DMET_Plus
            if(getOptBool("em-gender")) {
                emCelListener = new EmGenderCelListener(haploidProbeSets);
                reader.registerCelListener(emCelListener);
                chipSummaries.push_back(emCelListener);
                genderCallsMap[emCelListener->getGenderName()] = emCelListener;
            }
        }

        // @todo refactor 
        // 1) Add multiple channels for datastore
        // 2) Add feature effects for multiple channels
        // 3) SpfFile v4 reading and writing
        // 4) Freez layout and work off disk...
        // @todo refactor - setup data store here...
        Verbose::out(2, "Reading cel files.");
        DataStore *dStore = NULL;
        int stageIx = 0; 
        TmpFileFactory tmpFactory;
        tmpFactory.setTmpdir(getOpt("temp-dir"));
        // string stageTemp = tmpFactory.genFilename_basic(ToStr("stage-") + ToStr(stageIx), ".a5"); 
        string stageTemp = Fs::join(diskDir,ToStr("stage-") + ToStr(stageIx) + "-temp.a5");
        dStore = new DataStore(stageTemp);

        vector<int> desiredOrder;
        fillInDesiredOrder(desiredOrder);

        dStore->setProbeOrder(desiredOrder);
      //  string probeInfoTemp = tmpFactory.genFilename_basic(ToStr("stage-probe-info") + ToStr(stageIx), ".a5");
        string probeInfoTemp = Fs::join(diskDir,"stage-probe-info-temp.a5"); 
        DataStore probeInfo(probeInfoTemp);
        probeInfo.setProbeOrder(desiredOrder);
        probeInfo.setValues(*m_ChipLayout, *board);
        board->setProbeInfo(&probeInfo);
         
        reader.setFiles(celFiles);
        reader.registerCelListener(dStore);
        reader.readFiles();

        // Null covariate first
        InbredStatus *selectedInbred = noInbredReader;
        GenderCalls *selectedGenders = noGenderCalls;        

        
        string analysisString = getOpt("analysis");
        vector<AnalysisStage> stages = AnalysisStage::makeStages(analysisString, m_stdMethods);
        string analysisName =  getOpt("set-analysis-name");
        if (analysisName.empty()) {
            if(m_stdMethods.find(analysisString) != m_stdMethods.end()) {
                analysisName = analysisString;
            }
            else {
                analysisName = stages[0].name;
                for(int i = 1; i < stages.size(); i++) {
                    analysisName = analysisName + "." + stages[i].name;
                }
            }
        }
        
        fillInSelectedGenders(&selectedGenders, analysisName, genderCallsMap, noGenderCalls);

        Verbose::out(1,"Using gender method " + selectedGenders->getGenderName() + " for genotype calling.");
        Verbose::out(1,"Using inbred covariates " + selectedInbred->getInbredName() + " for genotype calling.");

        // @todo refactor - setup blackboard calls, genders, etc.
        GenotypeInfo gInfo;
        gInfo.m_Genders = selectedGenders;
        gInfo.m_Inbred = selectedInbred;
        gInfo.m_HaploidSnps = &haploidSnps;
        gInfo.m_SpecialSnps = &SpecialSnps;
        gInfo.m_CopyNumberMap = &copyNumberMap;
        board->setGenotypeInfo(&gInfo);

        AnalysisInfo info = makeAnalysisInfo(*m_ChipLayout, toRunProbesets, selectedGenders, analysisName);
        //unsigned int dotMod = max(int(toRunProbesets.size()/40), 1);

        board->setAnalysisInfo(&info);
        DTFactory dtFactory;

        std::vector<QuantMethodReport *> reporters;
        std::vector<QuantMethodReport *> exprReporters;
        makeReporters(stages[stages.size() - 1].name, reporters, exprReporters, info, chipSummaries, selectedGenders, NULL, SpecialSnps);
        board->setQuantGTypeReporters(&reporters);
        board->setQuantExpressionReporters(&exprReporters);
        Verbose::out(1,"Preprocessing cel file data");
        for(int i = 0; i < stages.size(); i++) {

          Verbose::out(1, "Doing stage: " + ToStr(i) + " " + stages[i].name);
          time_t startTime = time(NULL);
          //string nextStore = tmpFactory.genFilename_basic(ToStr("stage-") + ToStr(++stageIx), ".a5"); 
          string nextStore = Fs::join(diskDir,ToStr("stage-") + ToStr(++stageIx) + "-temp.a5");
          DataStore *out = new DataStore(nextStore);
          out->initFromDataStore(*dStore);
          //          runChipStreamStage(stages[i], *board, *dStore, *out);
          DataTransform *transform = dtFactory.dataTransformForString(stages[i].spec, *board);
          bool updated = transform->transformData(*board, *dStore, *out);
          Freez(transform);
          if (updated) {
            DataStore *toFree = dStore;
            dStore = out;
            Freez(toFree);
          }
          else {
            Freez(out);
          }
          time_t endTime = time(NULL);
          int t = int(  (float)(endTime - startTime) / 60.0 * 100); // convert to minutes
          Verbose::out(1, "Stage finished (" + ToStr((float)t/100) + " min)");
        }
        Verbose::out(1,"Done Preprocessing cel file data");
        if (!getOpt("db-from-prior-models").empty()) {
            SnpModelConverter conv;
            if (!getOpt("read-models-brlmmp").empty()) {
                conv.convertToDbModel(getOpt("read-models-brlmmp"), getOpt("db-from-prior-models"));
            }
            else if (!getOpt("read-models-birdseed").empty()) {
                conv.convertToDbModel(getOpt("read-models-birdseed"), getOpt("db-from-prior-models"));
            }
        }
        if (!getOpt("db-from-posterior-models").empty()) {
            SnpModelConverter conv;
            string input;
            board->get("model-file-written", &input);
            string file5Name = analysisName + ".snp-posteriors";
            conv.convertToDbModel(input, getOpt("db-from-posterior-models"), file5Name);
        }
        // Tell analyses that we are done computing.
        Freez(dStore);
    } // inner try end
    catch (...) {
      // cleanup -- KEEP IN SYNC WITH MEMORY FREEZ BELOW WHEN NO EXCEPTION
      ///@todo leaking probeSetToReport
      Freez(precompSeeds);
      Freez(dmCaller);
      Freez(emCelListener);
      Freez(cnProbeCelListener);
      Freez(fileGenderReader);
      Freez(noGenderCalls);
      Freez(noInbredReader);
      Freez(celStats);

      if(m_FreeChipLayout) {
        Freez(m_ChipLayout);
        for(uint32_t i = 0; i < m_ProbesetNames.size(); i++) {
            FreezArray(m_ProbesetNames[i]);
        }
        m_ProbesetNames.clear();
      }
      if(!Util::sameString(getOpt("probeset-ids"),"")) {
          for(uint32_t i = 0; i < toRunProbesets.size(); i++) {
              FreezArray(toRunProbesets[i]);
          }
      }
      toRunProbesets.clear();
      Freez(iMart);
      closeGlobalA5();
      removeTempDir(getOpt("temp-dir"));

      /// @todo reporters and main chip summary generater not freed

      // re-throw the unchanged exception.
      throw;
    }
  } // outer try end
  /* When things go wrong see if we can die gracefully here. */
  catch(Except &e) {
    Verbose::out(0,"");
    Err::errAbort(e.what());
  }
  catch(const std::bad_alloc &e) {
    Verbose::out(0,"");
    Err::errAbort("Ran out of memory. "
                  "Try quitting other applications.");
  }
  catch(CalvinException &ce) {
    Verbose::out(0,"");
    Err::errAbort("Affymetrix GeneChip Command Console library has thrown an exception. "
                  "Description: '" +StringUtils::ConvertWCSToMBS(ce.Description()) + "'");
  }
  catch(BroadException &e) {
    Verbose::out(0,"");
    Err::errAbort("Exception. msg: '" + ToStr(e.m_msg) + "' source file: '" + ToStr(e.m_sourcefile) +
                  ":" + ToStr(e.m_sourceline) + "'");
  }
  catch(const std::exception &e) {
    Verbose::out(0,"");
    Err::errAbort("Exception caught. "
                  "Most problems with this program are memory related. "
                  "Message is: " + ToStr(e.what()));
  }
  catch(const BaseException &e) {
    Verbose::out(0,"");
    Err::errAbort("newmat issue: " + ToStr(e.what()));
  }
  catch(...) {
    Verbose::out(0,"");
    Err::errAbort("Unknown exception caught "
                  "(most problems with this program are memory related).");
  }

  // cleanup -- KEEP IN SYNC WITH MEMORY FREEZ ABOVE IN CATCH
  ///@todo leaking probeSetToReport
  Freez(precompSeeds);
  Freez(dmCaller);
  Freez(emCelListener);
  Freez(cnProbeCelListener);
  Freez(fileGenderReader);
  Freez(noGenderCalls);
  Freez(noInbredReader);
  Freez(celStats);
  Freez(iMart);
  
  if(m_FreeChipLayout) {
     Freez(m_ChipLayout);
    for(uint32_t i = 0; i < m_ProbesetNames.size(); i++) {
        FreezArray(m_ProbesetNames[i]);
    }
    m_ProbesetNames.clear();
  }
  if(!Util::sameString(getOpt("probeset-ids"),"")) {
      for(uint32_t i = 0; i < toRunProbesets.size(); i++) {
          FreezArray(toRunProbesets[i]);
      }
  }
  toRunProbesets.clear();
  closeGlobalA5();
  removeTempDir(getOpt("temp-dir"));
  
  /// @todo reporters and main chip summary generater not freed

}

void ProbesetGenotypeStageEngine::closeGlobalA5() {
  // A5 close global groups
  if (m_a5_global_output_group!=NULL) {
    m_a5_global_output_group->close();
    delete m_a5_global_output_group;
    m_a5_global_output_group=NULL;
  }
  if (m_a5_global_input_group!=NULL) {
    m_a5_global_input_group->close();
    delete m_a5_global_input_group;
    m_a5_global_input_group=NULL;
  }

  // A5 close global files
  if(m_a5_global_output_file != m_a5_global_input_file) {
    affx::TsvReport::closeA5File(m_a5_global_output_file);
    affx::TsvReport::closeA5File(m_a5_global_input_file);
  } else {
    affx::TsvReport::closeA5File(m_a5_global_output_file);
    m_a5_global_input_file = NULL;
  }
}

void ProbesetGenotypeStageEngine::printStandardMethods(std::ostream &out) {
      map<string,string>::iterator iter;
      out << endl << "Standard Methods:" << endl;
      unsigned int maxLength = 0;
      for(iter = m_stdMethods.begin(); iter != m_stdMethods.end(); iter++) {
        if(iter->first.size() > maxLength)
          maxLength = iter->first.size();
      }
      for(iter = m_stdMethods.begin(); iter != m_stdMethods.end(); iter++) {
        unsigned int currentLength = 0;
        out << " '" << iter->first << "' ";
        currentLength = iter->first.size();
        while(currentLength < maxLength + 1) {
          out.put(' ');
          currentLength++;
        }
        out << iter->second << endl;
      }
}

void ProbesetGenotypeStageEngine::addReporters_Expression(affx::TsvReport::TsvReportFmt_t report_format,
                            AnalysisInfo& info,
                            QuantGTypeMethod *qMethod,
                            AnalysisStream *as,
                            vector<ChipSummary *> &chipSummaries,
                            GenderCalls *gender,
                            std::set<const char *, Util::ltstr> *probeSetsToReport,
                                                     bool doCompact) {
  //
  if(InstanceOf(qMethod,QuantLabelZ) || InstanceOf(qMethod,QuantBRLMM) || InstanceOf(qMethod,QuantBirdseed)) {
    QuantExprMethod *eMethod = NULL;
    if(InstanceOf(qMethod,QuantLabelZ)) {
      eMethod = static_cast<QuantLabelZ *>(qMethod)->getQuantExprMethod();
    } else if (InstanceOf(qMethod,QuantBRLMM)) {
      eMethod = static_cast<QuantBRLMM *>(qMethod)->getQuantExprMethod();
    } else if (InstanceOf(qMethod,QuantBirdseed)) {
      eMethod = static_cast<QuantBirdseed *>(qMethod)->getQuantExprMethod();
    }

    QuantMethodExprReport *report = new QuantMethodExprReport(info.getNumCols());
    report->m_qMethod=as->getQuantMethod();
    report->setIsHeaderBuffer(1);
    report->setDirPath(getOpt("out-dir"));
    report->m_qMethod=eMethod;
    report->setPrecision(5);
    report->setWriteOldStyleFeatureEffectsFile(getOptBool("writeOldStyleFeatureEffectsFile"));

    // We do not do this for the expression reporters as output
    // is controlled by the quant method. If you set this then
    // nothing will be reported as the probeset names reported
    // have the extra allele and context stuff appended.
    //report->registerProbeSetsToReport(probeSetsToReport);
    if(getOptBool("include-quant-in-report-file-name"))
        report->setFileprefix(as->getName()+"."+eMethod->getType());
    else
        report->setFileprefix(as->getName());
    report->setFilename("QuantMethodExprReport-DEBUG"); // filename not used -- for debugging

    if(report_format == affx::TsvReport::FMT_TSV) {
        report->m_DoSummary=getOptBool("summaries");
        report->m_DoResiduals=getOptBool("feat-details");
        report->m_DoFeatureEffects=getOptBool("feat-effects");
    } else {
        report->m_DoSummary=getOptBool("a5-summaries");
        report->m_DoResiduals=getOptBool("a5-feature-details");
        report->m_DoFeatureEffects=getOptBool("a5-feature-effects");
        
        if(getOptBool("a5-summaries-use-global")) {
            if(m_a5_global_output_group == NULL)
                Err::errAbort("--a5-summaries-use-global option given, but no global file. Must specify --a5-global-file.");
            report->m_summary_a5_group = m_a5_global_output_group;
            Verbose::out(1,"Using global A5 file for allele summary output.");
        }
    
        if(getOptBool("a5-feature-details-use-global")) {
            if(m_a5_global_output_group == NULL)
                Err::errAbort("--a5-feature-details-use-global option given, but no global file. Must specify --a5-global-file.");
            report->m_residuals_a5_group = m_a5_global_output_group;
            Verbose::out(1,"Using global A5 file for residuals output.");
        }
    
        if(getOptBool("a5-feature-effects-use-global")) {
            if(m_a5_global_output_group == NULL)
                Err::errAbort("--a5-feature-effects-use-global option given, but no global file. Must specify --a5-global-file.");
            report->m_feffects_a5_group = m_a5_global_output_group;
            Verbose::out(1,"Using global A5 file for feature effects output.");
        }
    }

    report->setFormat(report_format);
    if(doCompact && report_format == affx::TsvReport::FMT_A5) {
      report->setCompactFile5Format(info.m_MaxPsNameLength + 2);
    }
    report->addStdHeaders(report,
                          getOpt("exec-guid"),
                          as->getGuid(),
                          getOpt("time-start"),
                          getOpt("command-line"),
                          getOpt("version-to-report"),
                          info);

    if(InstanceOf(qMethod,QuantLabelZ)) {
      static_cast<QuantLabelZ *>(qMethod)->addExprReporter(report);
    } else if (InstanceOf(qMethod,QuantBRLMM)) {
      static_cast<QuantBRLMM *>(qMethod)->addExprReporter(report);
    } else if (InstanceOf(qMethod,QuantBirdseed)) {
      static_cast<QuantBirdseed *>(qMethod)->addExprReporter(report);
    }
  }
}

void ProbesetGenotypeStageEngine::addReporters_Expression(affx::TsvReport::TsvReportFmt_t report_format,
                                                          AnalysisInfo& info,
                                                          std::vector<QuantMethodReport *> &exprReporters,
                                                          vector<ChipSummary *> &chipSummaries,
                                                          GenderCalls *gender,
                                                          std::set<const char *, Util::ltstr> *probeSetsToReport,
                                                          bool doCompact) {
  //

    QuantMethodExprReport *report = new QuantMethodExprReport(info.getNumCols());
    report->setIsHeaderBuffer(1);
    report->setDirPath(getOpt("out-dir"));
    //    report->m_qMethod=eMethod;
    report->setPrecision(5);
    report->setWriteOldStyleFeatureEffectsFile(getOptBool("writeOldStyleFeatureEffectsFile"));

    // We do not do this for the expression reporters as output
    // is controlled by the quant method. If you set this then
    // nothing will be reported as the probeset names reported
    // have the extra allele and context stuff appended.
    //report->registerProbeSetsToReport(probeSetsToReport);
    if(getOptBool("include-quant-in-report-file-name"))
      report->setFileprefix(info.m_AnalysisName+".emethod");
    else
        report->setFileprefix(info.m_AnalysisName);
    report->setFilename("QuantMethodExprReport-DEBUG"); // filename not used -- for debugging

    if(report_format == affx::TsvReport::FMT_TSV) {
        report->m_DoSummary=getOptBool("summaries");
        report->m_DoResiduals=getOptBool("feat-details");
        report->m_DoFeatureEffects=getOptBool("feat-effects");
    } else {
        report->m_DoSummary=getOptBool("a5-summaries");
        report->m_DoResiduals=getOptBool("a5-feature-details");
        report->m_DoFeatureEffects=getOptBool("a5-feature-effects");
        
        if(getOptBool("a5-summaries-use-global")) {
            if(m_a5_global_output_group == NULL)
                Err::errAbort("--a5-summaries-use-global option given, but no global file. Must specify --a5-global-file.");
            report->m_summary_a5_group = m_a5_global_output_group;
            Verbose::out(1,"Using global A5 file for allele summary output.");
        }
    
        if(getOptBool("a5-feature-details-use-global")) {
            if(m_a5_global_output_group == NULL)
                Err::errAbort("--a5-feature-details-use-global option given, but no global file. Must specify --a5-global-file.");
            report->m_residuals_a5_group = m_a5_global_output_group;
            Verbose::out(1,"Using global A5 file for residuals output.");
        }
    
        if(getOptBool("a5-feature-effects-use-global")) {
            if(m_a5_global_output_group == NULL)
                Err::errAbort("--a5-feature-effects-use-global option given, but no global file. Must specify --a5-global-file.");
            report->m_feffects_a5_group = m_a5_global_output_group;
            Verbose::out(1,"Using global A5 file for feature effects output.");
        }
    }

    report->setFormat(report_format);
    if(doCompact && report_format == affx::TsvReport::FMT_A5) {
      report->setCompactFile5Format(info.m_MaxPsNameLength + 2);
    }
    report->addStdHeaders(report,
                          getOpt("exec-guid"),
                          info.m_AnalysisGuid,
                          getOpt("time-start"),
                          getOpt("command-line"),
                          getOpt("version-to-report"),
                          info);

  //   if(InstanceOf(qMethod,QuantLabelZ)) {
  //     static_cast<QuantLabelZ *>(qMethod)->addExprReporter(report);
  //   } else if (InstanceOf(qMethod,QuantBRLMM)) {
  //     static_cast<QuantBRLMM *>(qMethod)->addExprReporter(report);
  //   } else if (InstanceOf(qMethod,QuantBirdseed)) {
  //     static_cast<QuantBirdseed *>(qMethod)->addExprReporter(report);
  //   }
  // }
    exprReporters.push_back(report);
}

void ProbesetGenotypeStageEngine::fillInAnalysisInfo(AnalysisInfo &info, AnalysisStream *as, std::string prefix) {
    ///@todo need to refactor to ensure consistency between meta info tracked in probeset-summarize
    ///      versus probeset-genotype
    ///@todo automatically pull state and options info?
    ///@todo sync up option name with label in chp header
    // FOR NOW -- PLEASE KEEP THIS IN SYNC WITH probeset-genotype
    vector<string> celFiles = getOptVector("cels");
    assert(as);
    QuantMethod *qMethod = as->getQuantMethod();
    Err::check(qMethod != NULL, "Options::fillInAnalysisInfo() - Must have quantification method.");
    info.m_AlgVersion = qMethod->getVersion();
    info.m_AlgName = as->getName();
    info.m_ProgramName = getOpt("program-name");
    info.m_ProgramVersion = getOpt("version-to-report");
    info.m_ProgramCompany = getOpt("program-company");
    info.m_ChipType = getOpt("chip-type");
    info.m_ProgID = "";
    info.m_ExecGuid = getOpt("exec-guid");
    info.m_AnalysisGuid = as->getGuid();

    // State and Execution info
    info.addParam("apt-engine", "ProbesetGenotypeStageEngine");
    info.addParam(prefix + "program-name", getOpt("program-name"));
    info.addParam(prefix + "command-line", getOpt("command-line"));
    info.addParam(prefix + "exec-guid", getOpt("exec-guid"));
    info.addParam(prefix + "analysis-guid", as->getGuid());
    info.addParam(prefix + "time-str", getOpt("time-start"));
    info.addParam(prefix + "version", getOpt("version-to-report"));
    info.addParam(prefix + "cvs-id", getOpt("program-cvs-id"));
    info.addParam(prefix + "free-mem", getOpt("free-mem-at-start"));
    ///@todo add other state:
    ///@todo time-end state

    // Engine Options
    std::string subprefix = "opt-";
    // no do-ps-names
    info.addParam(prefix + subprefix + "chip-type", getOpt("chip-type"));
    info.addParam(prefix + subprefix + "probe-count", getOpt("probe-count"));
    info.addParam(prefix + subprefix + "force", getOpt("force"));
    // no precision
    info.addParam(prefix + subprefix + "out-dir", getOpt("out-dir"));
    // cc-expr-chp-out-dir
    info.addParam(prefix + subprefix + "cc-md-chp-out-dir", getOpt("cc-chp-out-dir"));
    info.addParam(prefix + subprefix + "xda-chp-out-dir", getOpt("xda-chp-out-dir"));
    info.addParam(prefix + subprefix + "cdf-file",       Fs::basename(getOpt("cdf-file")));
    // no spf-file 
    // no pgf-file
    // no clf-file
    // no bgp-file
    info.addParam(prefix + subprefix + "ps-list-file",   Fs::basename(getOpt("probeset-ids")));
    // no meta-ps-file
    // no-qc-groups-file
    info.addParam(prefix + subprefix + "kill-list",      Fs::basename(getOpt("kill-list")));
    info.addParam(prefix + subprefix + "temp-dir", getOpt("temp-dir"));
    info.addParam(prefix + subprefix + "use-disk", getOpt("use-disk"));
    info.addParam(prefix + subprefix + "diskCache", getOpt("disk-cache"));
    info.addParam(prefix + subprefix + "cel-count", ToStr(celFiles.size()));
    for(uint32_t i = 0; i < celFiles.size(); i++) {
      std::string paramName = prefix + subprefix + "cel-" + ToStr(i+1);
      info.addParam(paramName, Fs::basename(celFiles[i]));
    }

    // Algorithm Options
    info.addParam(prefix + subprefix + "analysis-name", as->getName());
    info.addParam(prefix + subprefix + "set-analysis-name", getOpt("set-analysis-name"));
    info.addParam(prefix + subprefix + "analysis-spec", as->getAnalysisSpec());
    info.addParam(prefix + subprefix + "feat-effect-file", Fs::basename(getOpt("use-feat-eff")));
    info.addParam(prefix + subprefix + "target-sketch-file", Fs::basename(getOpt("target-sketch")));
    info.addParam(prefix + subprefix + "reference-profile-file",Fs::basename(getOpt("reference-profile")));
    info.addParam(prefix + subprefix + "do-residuals", getOpt("feat-details"));
    info.addParam(prefix + subprefix + "do-feature-effects", getOpt("feat-effects"));
    info.addParam(prefix + subprefix + "write-sketch", getOpt("write-sketch"));
    info.addParam(prefix + subprefix + "write-profile",getOpt("write-profile"));

    info.addParam(prefix + subprefix + "a5-write-sketch", getOpt("a5-sketch"));
    info.addParam(prefix + subprefix + "set-gender-method", getOpt("set-gender-method"));
    info.addParam(prefix + subprefix + "chrXProbeFile", Fs::basename(getOpt("chrX-probes")));
    info.addParam(prefix + subprefix + "chrYProbeFile", Fs::basename(getOpt("chrY-probes")));
    info.addParam(prefix + subprefix + "read-priors-brlmm", Fs::basename(getOpt("read-priors-brlmm")));
    info.addParam(prefix + subprefix + "initial-calls", ToStr("false"));
    info.addParam(prefix + subprefix + "output-summaries", getOpt("summaries"));
    info.addParam(prefix + subprefix + "list-sample", getOpt("list-sample"));
    info.addParam(prefix + subprefix + "write-prior", getOpt("write-prior"));
    info.addParam(prefix + subprefix + "chrX-file", Fs::basename(getOpt("chrX-snps")));
    info.addParam(prefix + subprefix + "chrX-snps", Fs::basename(getOpt("chrX-snps")));
    info.addParam(prefix + subprefix + "no-gender-force", getOpt("no-gender-force"));
    info.addParam(prefix + subprefix + "special-snps", Fs::basename(getOpt("special-snps")));
    info.addParam(prefix + subprefix + "compute-prior", ToStr(getOptInt("prior-size")));
    info.addParam(prefix + subprefix + "norm-size", getOpt("norm-size"));
    info.addParam(prefix + subprefix + "dm-thresh", getOpt("dm-thresh"));
    info.addParam(prefix + subprefix + "dm-het-mult", getOpt("dm-hetmult"));
    info.addParam(prefix + subprefix + "read-genders",        Fs::basename(getOpt("read-genders")));
    info.addParam(prefix + subprefix + "read-inbred",         Fs::basename(getOpt("read-inbred")));
    info.addParam(prefix + subprefix + "model-file-brlmm"   , Fs::basename(getOpt("read-models-brlmm")));
    info.addParam(prefix + subprefix + "model-file-brlmmp"  , Fs::basename(getOpt("read-models-brlmmp")));
    info.addParam(prefix + subprefix + "model-file-birdseed", Fs::basename(getOpt("read-models-birdseed")));
    info.addParam(prefix + subprefix + "genotypes-file"     , Fs::basename(getOpt("genotypes")));
    info.addParam(prefix + subprefix + "dm-calls-out", getOpt("dm-out"));
    info.addParam(prefix + subprefix + "qmethod-spec", getOpt("qmethod-spec"));
    info.addParam(prefix + subprefix + "genotype-only", ToStr(!getOptBool("all-types")));

    // Quantification Related Meta Info
    info.addParam("quantification-name", qMethod->getType());
    info.addParam("quantification-version", qMethod->getVersion());
    info.addParam("quantification-scale", QuantMethod::scaleToTxt(qMethod->getScale()));
    info.addParam("quantification-type", QuantMethod::quantTypeToTxt(qMethod->getQuantType()));

    // Add user defined meta data parameters
    vector<pair<string, string> > metaData = getMetaDataDescription();
    for (int i = 0; i < metaData.size(); i++) {
      pair<string,string> p = metaData[i];
      info.addClientMetaInfo(p.first, p.second);
    }
}

void ProbesetGenotypeStageEngine::extraHelp() {
    ChipStreamFactory cFactory;
    PmAdjusterFactory aFactory;
    QuantMethodFactory qFactory(QuantMethodFactory::Expression);
    AnalysisStreamFactory asFactory;

    printStandardMethods(cout);
    EngineUtil::printSelfDocs("Data transformations:", cFactory.getDocs());
    EngineUtil::printSelfDocs("Pm Intensity Adjustments:", aFactory.getDocs());
    EngineUtil::printSelfDocs("Quantification Methods:", qFactory.getDocs());
    EngineUtil::printSelfDocs("Analysis Streams:", asFactory.getDocs());
}

void ProbesetGenotypeStageEngine::explain() {
    EngineUtil::explainParameter(getOpt("explain"));
}
