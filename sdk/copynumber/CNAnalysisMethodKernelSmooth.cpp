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
 * @file CNAnalysisMethodKernelSmooth.cpp
 *
 * @brief This file contains the CNAnalysisMethodKernelSmooth class members.
 */
#include "copynumber/CNAnalysisMethodKernelSmooth.h"
//
#include "copynumber/GKernel.h"
#include "copynumber/kernel_smooth.h"
//


/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodKernelSmooth::explainSelf()
    {
    CNAnalysisMethodKernelSmooth obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodKernelSmooth::getDefaultDocOptions()
    {
    std::vector<SelfDoc::Opt> opts;

    // SelfDoc::Opt(name, type, value, default, min, max, description)

    SelfDoc::Opt sigma_span = {"sigma_span", SelfDoc::Opt::Double,
        "50.0", "50.0", "0.0", "NA", "Probes spanned by one stdev"};
    opts.push_back(sigma_span);

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
SelfCreate * CNAnalysisMethodKernelSmooth::newObject(
    std::map<std::string,std::string>& params)
    {
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodKernelSmooth* pMethod = new CNAnalysisMethodKernelSmooth();
    std::string strPrefix = getPrefix();

    pMethod->m_sigma_span = setupDoubleParameter("sigma_span",
        strPrefix, params, doc, opts);

    /*
    pMethod->m_iXChromosome = setupIntParameter("xChromosome",
        strPrefix, params, doc, opts);
    */

    return pMethod;
    }

//constructor
CNAnalysisMethodKernelSmooth::CNAnalysisMethodKernelSmooth()
{
        m_bLocalProbeSetsDetermined=false;
}



/**
 * @brief Run the analysis
 */
void CNAnalysisMethodKernelSmooth::run()
    {
    Verbose::out(1, "CNAnalysisMethodKernelSmooth::run(...) start");

    isSetup();

    determineLocalProbeSets();

    vector<int>chromosomes = getChromosomes(getProbeSets());

    GaussianKernel GK(m_sigma_span);

    for (unsigned int k=0; k<chromosomes.size(); k++)
        {
        int chr = chromosomes[k];
        pair<int,int>chr_span = getChrBounds(chr, getProbeSets());

        unsigned int start = chr_span.first;
        unsigned int span = (unsigned int)(chr_span.second - chr_span.first);
        vector<double>smooth_me(span);
        for (unsigned int i=0; i<span; i++)
            smooth_me[i] = getProbeSets()->getAt(i + start)->getLog2Ratio();

        kernel_smooth(smooth_me, GK);

        for (unsigned int i=0; i<span; i++)
            {
            getProbeSets()->getAt(i+start)->setSmoothedLog2Ratio(smooth_me[i]);
            }
        }
    if (m_pEngine->getOptBool("keep-intermediate-data"))
    {
        writeSmoothedLog2Ratios("kernelSmoothed", getProbeSets());
    }

//    Verbose::progressEnd(1, "Done");
    Verbose::out(1, "CNAnalysisMethodKernelSmooth::run(...) end");
    }
void CNAnalysisMethodKernelSmooth::determineLocalProbeSets()
{
        if(m_bLocalProbeSetsDetermined)
        {
                return;
        }
        const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
        int iNumberOfProbeSets=CNAnalysisMethod::getProbeSets()->getCount();
        for (int iIndex = 0; iIndex<iNumberOfProbeSets; iIndex++)
        {
                if( CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAll() )
                {
						if (isCytoScanHD && (!CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAsCN())) {continue;}
                        getProbeSets()->add( CNAnalysisMethod::getProbeSets()->getAt(iIndex));
                }
        }
        m_bLocalProbeSetsDetermined=true;
}

