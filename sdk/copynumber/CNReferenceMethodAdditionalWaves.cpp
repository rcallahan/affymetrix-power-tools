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
#include "copynumber/CNReferenceMethodAdditionalWaves.h"
//
#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNCychp.h"
#include "copynumber/CNLog2RatioData.h"
//
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "portability/affy-base-types.h"

#include "util/AffxMultiDimensionalArray.h"
#include "util/TmpFileFactory.h"
#include "util/Verbose.h"
//
#include <algorithm>
#include <limits>
//


/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNReferenceMethodAdditionalWaves::explainSelf()
    {
    CNReferenceMethodAdditionalWaves obj;
    SelfDoc doc;
    doc.setDocName(obj.getType());
    doc.setDocDescription(obj.getDescription());
    doc.setDocOptions(obj.getDefaultDocOptions());
    return doc;
    }

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> CNReferenceMethodAdditionalWaves::getDefaultDocOptions()
    {
    std::vector<SelfDoc::Opt> opts;

    // SelfDoc::Opt(name, type, value, default, min, max, description)

    SelfDoc::Opt dTrim = {"trim", SelfDoc::Opt::Double, "2.0", "2.0", "NA", "NA", "Log2Ratio Trim value."};
    opts.push_back(dTrim);

    SelfDoc::Opt dPercentile = {"percentile", SelfDoc::Opt::Double, "0.75", "0.75", "NA", "NA", "High Percentile value."};
    opts.push_back(dPercentile);

    SelfDoc::Opt iWaveCount = {"additional-wave-count", SelfDoc::Opt::Integer, "1", "1", "0", "NA", "Number of waves to add to the reference."};
    opts.push_back(iWaveCount);

    SelfDoc::Opt bDemean = {"demean", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", "Demean the input to the SVD."};
    opts.push_back(bDemean);

    SelfDoc::Opt dCNQCCutoff = {"cn-qc-cutoff", SelfDoc::Opt::Double, "0.27", "0.27", "NA", "NA", "If the CN QC values if over this cutoff then the sample fails QC."};
    opts.push_back(dCNQCCutoff);

    SelfDoc::Opt dSNPQCCutoff = {"snp-qc-cutoff", SelfDoc::Opt::Double, "1.1", "1.1", "NA", "NA", "If the SNP QC values if below this cutoff then the sample fails QC."};
    opts.push_back(dSNPQCCutoff);

    SelfDoc::Opt iWaveSegCountCutoff = {"waviness-seg-count-cutoff", SelfDoc::Opt::Integer, "100", "100", "0", "NA", "The waviness seg count cutoff."};
    opts.push_back(iWaveSegCountCutoff);

    SelfDoc::Opt bUseHighWavinessSegCounts = {"use-high-waviness-seg-count", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA", "Use only those cychp files that have a waviness-seg-count > the cutoff if true, else use only those cychp files <= the cutoff."};
    opts.push_back(bUseHighWavinessSegCounts);

    SelfDoc::Opt bForce = {"force", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", "Force the job to run even if there is a mismatch between the cychp files and the input CN reference."};
    opts.push_back(bForce);

    SelfDoc::Opt bKeepTempData = {"keep-temp-data", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", "If true, then do not delete the temporary data files used by the module."};
    opts.push_back(bKeepTempData);

    SelfDoc::Opt selectedQC = {"selected-qc", SelfDoc::Opt::String, "snp-qc", "snp-qc", "", "", "Choose the QC option to use.  Available values: snp-qc, contrast-qc, and contrast-qc-nsp"};
    opts.push_back(selectedQC);

//    SelfDoc::Opt dCutoff = {"wave-cutoff", SelfDoc::Opt::Double, "0.000001", "0.000001", "0", "NA", "The cutoff for wave calculation."};
//    opts.push_back(dCutoff);

    return opts;
    }


/**
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 *
 * @param param - Map of key/value pairs to initialize the object.
 *
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate * CNReferenceMethodAdditionalWaves::newObject(
    std::map<std::string,std::string>& params)
    {
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNReferenceMethodAdditionalWaves* pMethod = new CNReferenceMethodAdditionalWaves();
    std::string strPrefix = getPrefix();

    pMethod->m_dTrim = setupDoubleParameter("trim", strPrefix, params, doc, opts);
    pMethod->m_dPercentile = setupDoubleParameter("percentile", strPrefix, params, doc, opts);
    pMethod->m_iWaveCount = setupIntParameter("additional-wave-count", strPrefix, params, doc, opts);
    pMethod->m_bDemean = setupBoolParameter("demean", strPrefix, params, doc, opts);
    pMethod->m_dCNQCCutoff = setupDoubleParameter("cn-qc-cutoff", strPrefix, params, doc, opts);
    pMethod->m_dSNPQCCutoff = setupDoubleParameter("snp-qc-cutoff", strPrefix, params, doc, opts);
    pMethod->m_iWavinessSegCountCutoff = setupIntParameter("waviness-seg-count-cutoff", strPrefix, params, doc, opts);
    pMethod->m_bUseHighWavinessSegCount = setupBoolParameter("use-high-waviness-seg-count", strPrefix, params, doc, opts);
    pMethod->m_bForce = setupBoolParameter("force", strPrefix, params, doc, opts);
    pMethod->m_bKeepTempData = setupBoolParameter("keep-temp-data", strPrefix, params, doc, opts);
    pMethod->m_strSelectedQC = setupStringParameter("selected-qc", strPrefix, params, doc, opts);
//    pMethod->m_dWaveCutoff = setupDoubleParameter("wave-cutoff", strPrefix, params, doc, opts);

    if ((pMethod->m_iWavinessSegCountCutoff == 0) && (pMethod->m_bUseHighWavinessSegCount == false))
    {
        Err::errAbort("Cannot specify waviness-seg-count-cutoff=0 and use-high-waviness-seg-count=false for the same run.");
    }

    return pMethod;
    }
