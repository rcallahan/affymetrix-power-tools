////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

// put this last so it gets "STRUCT_ALIGNMENT" from the CEL reader.
#include "chipstream/QuantBirdseedDev.h"
//
#include "chipstream/QuantBRLMM.h"
#include "chipstream/QuantMethodFactory.h"
//
#include "birdseed-dev/FitSNPGaussiansPriors3.h"
#include "birdseed-dev/GenotypeCaller.h"
#include "birdseed-dev/PriorsReader.h"
#include "broadutil/BroadUtil.h"
#include "broadutil/CelUtil.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <sstream>
//

using namespace affx;
using namespace birdseed::dev;

SelfDoc QuantBirdseedDev::explainSelf() { 
    SelfDoc doc;
    setupSelfDoc(doc);
    return doc;
}

std::vector<SelfDoc::Opt> QuantBirdseedDev::getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt> opts;
    SelfDoc::Opt confThreshold = {"conf-threshold", SelfDoc::Opt::Double, "0.1", "0.1", "0", "1.0",
                                "Confidence must be <= this value to be considered a call."};
    opts.push_back(confThreshold);
    SelfDoc::Opt correctionFactor = {"correction-factor", SelfDoc::Opt::Double, "-1.0", "-1.0", "-1.0", "1000000",
                                "Scaling factor for priors. If <= 0, correction factor is calculated internally."};
    opts.push_back(correctionFactor);
    birdseed::dev::getSelfDocOptions(&opts);
    return opts;
}

SelfCreate *QuantBirdseedDev::newObject(std::map<std::string,std::string> &param) {
    SelfDoc doc = explainSelf();
    double confThreshold;
    fillInValue(confThreshold, "conf-threshold", param, doc);
    double correctionFactor;
    fillInValue(correctionFactor, "correction-factor", param, doc);
    birdseed::dev::setSelfDocOptions(&doc, param);
    QuantBirdseedDev *obj = new QuantBirdseedDev(confThreshold, correctionFactor);

    /// @todo hack to get the self doc correctly set with user supplied values
    /// this should really go in the constructor
    birdseed::dev::setSelfDocOptions(obj, param);
    return obj;
}

// Index is birdseed::dev::CallAndConfidence::Call
// This appears to produce the right result, despite affx::BB == 0 while Call::BB == 0
static const GType kBirdseedCallMap[] = {affx::AA, affx::AB, affx::BB, affx::NN};

/** 
 * @brief Do the heavy lifting of estimation.
 */
void QuantBirdseedDev::computeEstimate()
{
    // Without doing this step, the values are the same as what is output by -a pm-only,sumz, so perhaps
    // this isn't necessary.  When this code was executed, some intensities were negative.
    //transformData(m_AValues, m_BValues);
    birdseed::dev::IntensityMatrix intensities(m_AValues.size());
    assert(m_AValues.size() == m_BValues.size());
    for (size_t i = 0; i < m_AValues.size(); ++i) {
        intensities[i][0] = m_AValues[i];
        intensities[i][1] = m_BValues[i];
    }

    int birdseedVerbosity = (verbosity <= 3? 0: 2);
    birdseed::dev::GenderAwareGenotypeCaller caller(intensities,
                                               m_Genders,
                                               m_PriorsReader.get(),
                                               m_ProbesetName.c_str(),
                                               m_CorrectionFactor,
                                               birdseedVerbosity,
                                               m_ClusterOstrm.get());

    m_Calls.reserve(intensities.numRows());
    m_Confidences.reserve(intensities.numRows());
    m_Distances.resize(intensities.numRows());
    
    for (size_t i = 0; i < intensities.numRows(); ++i) {
        birdseed::dev::CallAndConfidence callAndConfidence = caller.nextCall();
        assert(callAndConfidence.call >= 0 && callAndConfidence.call < ARRAY_LEN(kBirdseedCallMap));
        // convert from birdseed enum to affx::GType
        m_Calls.push_back(kBirdseedCallMap[callAndConfidence.call]);
        m_Confidences.push_back(callAndConfidence.confidence);
        m_Distances[i].resize(birdseed::dev::MAX_NUM_CLUSTERS);
        m_Distances[i][affx::AA] = callAndConfidence.pvalues[birdseed::dev::Priors::AA_INDEX];
        m_Distances[i][affx::AB] = callAndConfidence.pvalues[birdseed::dev::Priors::AB_INDEX];
        m_Distances[i][affx::BB] = callAndConfidence.pvalues[birdseed::dev::Priors::BB_INDEX];
    }
}

// We should test for the format of the file here.
void QuantBirdseedDev::setPriorFile(const std::string& priorsPath, const string &chipType, 

                                 std::set<const char *, Util::ltstr> *probeSetsToLoad, 
                                 const std::map<std::string, std::pair< int, int > > &specialSnps)
{
  // @todo In the interest of software testing, we will only read the tsv format file.
  m_PriorsReader.reset(new PriorsReader(specialSnps,
                                        new TsvPriorsReader(priorsPath, probeSetsToLoad)));
  Verbose::out(2, "Loaded: " + ToStr(m_PriorsReader->numPriors()) + " priors from: " + ToStr(priorsPath));
}

/******************************************************************/
/**************************[END OF QuantBirdseedDev.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
