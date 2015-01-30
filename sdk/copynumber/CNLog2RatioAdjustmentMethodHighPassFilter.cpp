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
#include "copynumber/CNLog2RatioAdjustmentMethodHighPassFilter.h"

#include "copynumber/DataBlock.h"
#include "copynumber/Indexer.h"

#include "util/Verbose.h"

#include <new>
#include <utility>
#include <algorithm>

using namespace std;


static bool pci_sortf(const ProbeCellInfo & A, const ProbeCellInfo & B) {
        return (A.probe_set_idx < B.probe_set_idx);
        }

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNLog2RatioAdjustmentMethodHighPassFilter::explainSelf() 
    {
    CNLog2RatioAdjustmentMethodHighPassFilter obj;
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
std::vector<SelfDoc::Opt> CNLog2RatioAdjustmentMethodHighPassFilter::
        getDefaultDocOptions() 
    {
    std::vector<SelfDoc::Opt> opts;

    SelfDoc::Opt MiniBlockRows = {"l2r-mini-block-rows",
        SelfDoc::Opt::Integer, "8", "8", "NA", "NA", "Mini block rows."};
    opts.push_back(MiniBlockRows);

    SelfDoc::Opt MiniBlockCols = {"l2r-mini-block-cols",
        SelfDoc::Opt::Integer, "8", "8", "NA", "NA", "Mini block columns."};
    opts.push_back(MiniBlockCols);

    SelfDoc::Opt GlobalSmoothWeight = {"l2r-global-smooth-weight",
        SelfDoc::Opt::Double, "256.0", "256.0", "NA", "NA",
        "Global smooth weight."};
    opts.push_back(GlobalSmoothWeight);

    SelfDoc::Opt LocalSmoothWeight = {"l2r-local-smooth-weight",
        SelfDoc::Opt::Double, "64.0", "64.0", "NA", "NA",
        "Local smooth weight."};
    opts.push_back(LocalSmoothWeight);

    SelfDoc::Opt Converged = {"l2r-converged", SelfDoc::Opt::Double,
        "0.0001", "0.0001", "NA", "NA", "Converged."};
    opts.push_back(Converged);

    SelfDoc::Opt Use = {"use", SelfDoc::Opt::Boolean, "true", "true",
        "NA", "NA", "Use the high pass filter"};
    opts.push_back(Use);

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
SelfCreate * CNLog2RatioAdjustmentMethodHighPassFilter::
        newObject(std::map<std::string,std::string>& params) 
    {
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNLog2RatioAdjustmentMethodHighPassFilter * pMethod =
        new CNLog2RatioAdjustmentMethodHighPassFilter();
    std::string strPrefix = getPrefix();

    pMethod->mini_block_rows = setupIntParameter("l2r-mini-block-rows",
        strPrefix, params, doc, opts);

    pMethod->mini_block_cols = setupIntParameter("l2r-mini-block-cols",
        strPrefix, params, doc, opts);

    pMethod->global_smooth_weight =
        setupDoubleParameter("l2r-global-smooth-weight",
        strPrefix, params, doc, opts);

    pMethod->local_smooth_weight =
        setupDoubleParameter("l2r-local-smooth-weight",
        strPrefix, params, doc, opts);

    pMethod->converged =
        setupDoubleParameter("l2r-converged", strPrefix, params, doc, opts);

    pMethod->use_hpf = setupBoolParameter("use", strPrefix, params, doc, opts);
    
    return pMethod;
    }


/**
 * @brief Constructor
 */
CNLog2RatioAdjustmentMethodHighPassFilter::
        CNLog2RatioAdjustmentMethodHighPassFilter()
    {
    mini_block_cols = 0;
    mini_block_rows = 0;
    global_smooth_weight = 0.0;
    local_smooth_weight = 0.0;
    converged = 0.0;
    }


/**
 * @brief Destructor
 */
CNLog2RatioAdjustmentMethodHighPassFilter::
    ~CNLog2RatioAdjustmentMethodHighPassFilter() { }

/**
 * @brief Run the calculations.
 */
void CNLog2RatioAdjustmentMethodHighPassFilter::run()
    {
    Verbose::out(1,
        "CNLog2RatioAdjustmentMethodHighPassFilter::run(...) start");

    if (!use_hpf) { Verbose::out(1,
         "Use of high pass filter set to false. Exiting HPF");
        return;
        }

    // Note that for chips so far, a dimension swap needs to occur to
    // use a right handed coordinate system.

    // for example 2680 rows in SNP6 and CytoSCanHD;
    int image_rows = m_pEngine->getOptInt("probe-columns-count");

    // for example 2572 columns in SNP6 and CytoSCanHD;
    int image_cols = m_pEngine->getOptInt("probe-rows-count");


    // fetch pointers to copy number probes and probe sets
    CNProbeArray * cn_probes = getProbes();
    CNProbeSetArray * cn_probe_sets = CNAnalysisMethod::getProbeSets();

    // this order produces a valid CNPobe::getProbeSetIndex()
    cn_probe_sets->quickSort(0);

    // For speed we need a fast lookup table - no maps
    vector<ProbeCellInfo>pc_info_vec;

    // go through all the probes used only for copy number and collect
    // their array positions, probe ids and indices into the array of
    // probe sets.
    for (int i=0; i<cn_probes->size(); i++) {
        // CNProbe index into the CNProbeSet array
        int probe_set_idx = cn_probes->at(i)->getProbeSetIndex();

        // Only non-polymorphic probes are high pass filtered
        if(!cn_probe_sets->at(probe_set_idx)->processAsCN()) continue;

        // It may turn out due to custom annotation files or product
        // development quirks that some copy number probes are out of
        // the workflow and might still carry raw intensities.  Ignore
        // them.
        if(!cn_probe_sets->at(probe_set_idx)->processAsCNNormalize()) continue;

        // The probe ID can be used to obtain probe cell array coordinates
        int probe_id = cn_probes->at(i)->getProbeID();

        // Here are the coordinates in a right hand coordinate system
        pair<int,int> row_col =
            probeid_to_rc(probe_id,image_rows,image_cols);

        // Assemble the above into a struct and push back
        ProbeCellInfo PCI;
        PCI.row = row_col.first;
        PCI.col = row_col.second;
        PCI.probe_id = probe_id;
        PCI.probe_set_idx = probe_set_idx;
        PCI.autosome = false;
        if (cn_probe_sets->at(probe_set_idx)->getChromosome() < 23) {
            PCI.autosome = true;
            }
        pc_info_vec.push_back(PCI);
        }

    //  If there is nothing to do just exit
    if (!pc_info_vec.size()) return;


    // Make an instance of the high pass filter.  Resizing or clearing is
    // not implemented.  A new one is always needed.
    DataBlock * DB = new DataBlock(image_rows, image_cols,
                        mini_block_rows, mini_block_cols);


    // Load up the DataBlock with whatever was observed.  The elements
    // in the Datablock with no entries are missing values.  Missing
    // values are okay. 
    for (vector<ProbeCellInfo>::iterator iter=pc_info_vec.begin();
            iter != pc_info_vec.end(); iter++) {
        if (!iter->autosome) continue;
        double log2r = cn_probe_sets->at(iter->probe_set_idx)->getLog2Ratio();
        DB->set_value(iter->row,iter->col,log2r);
        }


    // Load the full hierarchy of the data block and run to convergence.
    DB->load_array_level();
    DB->mrf_smooth(local_smooth_weight, global_smooth_weight, converged);

    // order by array index to avoid using a lookup map
    sort(pc_info_vec.begin(),pc_info_vec.end(),pci_sortf);

    // Subtract the background estimate from estimated DataBlock.
    int count=0;
    double log2r_backg = 0.0;
    vector<ProbeCellInfo>::iterator pci_iter = pc_info_vec.begin();

    double low_pass_norm = 0.0;
    int low_pass_count=0;

    while (pci_iter != pc_info_vec.end()) {
        int probe_set_idx = pci_iter->probe_set_idx;
        count += 1;

        double backg = DB->background(pci_iter->row,pci_iter->col);
        log2r_backg += backg;
        low_pass_norm += backg*backg;
        low_pass_count++;
        
        pci_iter++;

        // the backgrounds of probes in a probe set are averaged
        if (pci_iter != pc_info_vec.end()) {
            if (pci_iter->probe_set_idx == probe_set_idx) continue;
            }

        // subtract the low pass from the log2r to get the high pass
        double log2r_obs = cn_probe_sets->at(probe_set_idx)->getLog2Ratio();
        log2r_obs -= log2r_backg/count;
        cn_probe_sets->getAt(probe_set_idx)->setLog2Ratio(log2r_obs);

        log2r_backg = 0.0;
        count = 0;
        }

    // free memory that was allocated
    delete DB;

    low_pass_norm = 128.0*sqrt(low_pass_norm/low_pass_count);
    getExperiment()->setL2Gradient(low_pass_norm);

    // out of courtesy
    cn_probe_sets->quickSort(1);
    }
