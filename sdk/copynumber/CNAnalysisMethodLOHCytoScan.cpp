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
 * @file CNAnalysisMethodLOHCytoScan.cpp
 *
 * @brief This file contains the CNAnalysisMethodLOHCytoScan class members.
 */
#include "copynumber/CNAnalysisMethodLOHCytoScan.h"
#include "file5/File5_Tsv.h"
#include "file/TsvFile/TsvFile.h"
//
#include "util/AffxStatistics.h"
#include "util/Guid.h"
//



CNAnalysisMethodLOHCytoScan::CNAnalysisMethodLOHCytoScan()
{

        m_bLocalProbeSetsDetermined = false;    
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodLOHCytoScan::explainSelf()
{
    CNAnalysisMethodLOHCytoScan obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodLOHCytoScan::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)

  SelfDoc::Opt lohCS_errorrate = {"lohCS_errorrate", SelfDoc::Opt::Double, "0.05", "0.05", "NA", "NA", "LOH CS Error Rate"};
  opts.push_back(lohCS_errorrate);

  SelfDoc::Opt lohCS_beta = {"lohCS_beta", SelfDoc::Opt::Double, "0.001", "0.001", "NA", "NA", "LOH CS Beta"};
  opts.push_back(lohCS_beta);

  SelfDoc::Opt lohCS_alpha = {"lohCS_alpha", SelfDoc::Opt::Double, "0.01", "0.01", "NA", "NA", "LOH CS Alpha"};
  opts.push_back(lohCS_alpha);

  SelfDoc::Opt lohCS_separation = {"lohCS_separation", SelfDoc::Opt::Integer, "1000000", "1000000", "NA", "NA", "LOH CS Separation"};
  opts.push_back(lohCS_separation);

  SelfDoc::Opt lohCS_nMinMarkers = {"lohCS_nMinMarkers", SelfDoc::Opt::Integer, "10", "10", "NA", "NA", "LOH CS Minimum Marker Count"};
  opts.push_back(lohCS_nMinMarkers);

  SelfDoc::Opt lohCS_NoCallThreshold = {"lohCS_NoCallThreshold", SelfDoc::Opt::Double, "0.05", "0.05", "NA", "NA", "LOH CS No Call Threshold"};
  opts.push_back(lohCS_NoCallThreshold);

  SelfDoc::Opt lohCS_MinGenomicSpan = {"lohCS_minGenomicSpan", SelfDoc::Opt::Integer, "1000000", "1000000", "NA", "NA", "LOH CS Minimum Genomic Span"};
  opts.push_back(lohCS_MinGenomicSpan);

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
SelfCreate* CNAnalysisMethodLOHCytoScan::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodLOHCytoScan* pMethod = new CNAnalysisMethodLOHCytoScan();
    std::string strPrefix = getPrefix();

    pMethod->m_dLohCNErrorRate = setupDoubleParameter("lohCS_errorrate", strPrefix, params, doc, opts);
    pMethod->m_dLohCNBeta = setupDoubleParameter("lohCS_beta", strPrefix, params, doc, opts);
    pMethod->m_dLohCNAlpha = setupDoubleParameter("lohCS_alpha", strPrefix, params, doc, opts);
    pMethod->m_iLohCNSeparation = setupIntParameter("lohCS_separation", strPrefix, params, doc, opts);
    pMethod->m_iLohCNMinimumMarkerCount = setupIntParameter("lohCS_nMinMarkers", strPrefix, params, doc, opts);
    pMethod->m_dLohCNNoCallThreshold = setupDoubleParameter("lohCS_NoCallThreshold", strPrefix, params, doc, opts);
    pMethod->m_iLohCNMinGenomicSpan = setupIntParameter("lohCS_minGenomicSpan", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief run the analysis
 */
void CNAnalysisMethodLOHCytoScan::run()
{
    Verbose::out(1, "CNAnalysisMethodLOHCytoScan::run(...) start");
    isSetup();

    getEngine()->setOpt("minSegSeparation", ::getInt(m_iLohCNSeparation));

    determineLocalProbeSets();

    std::vector<char> vHomHet(getProbeSets()->size());
    int iMarkerCount = 0;
    int iHetCutoff = 0;
    if (!lohPreProcessing(vHomHet, iMarkerCount, iHetCutoff)) {throw(Except("LohPreProcessing failed."));}
    int iLastChromosome = getProbeSets()->at(getProbeSets()->size() - 1)->getChromosome();
//    Verbose::progressBegin(1, "CNAnalysisMethodLOH::run(...) ", iLastChromosome, 1, iLastChromosome);
    for (int iChromosome = 1; (iChromosome <= iLastChromosome); iChromosome++)
    {
        try
        {
//        Verbose::progressStep(1);
        int iChromosomeProbeSetCount = 0;
        for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
        {
            CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
            if ((vHomHet[iIndex] >= 0) && (pobjProbeSet->getChromosome() == iChromosome))
            {
                iChromosomeProbeSetCount++;
            }
        }
        std::vector<char> vChromosomeGenotypeCalls(iChromosomeProbeSetCount);
        std::vector<int> vChromosomePositions(iChromosomeProbeSetCount);
        std::vector<int> vChromosomeProbeSetIndexes(iChromosomeProbeSetCount);
        int iProbeSetIndex = 0;
        for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
        {
            CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
            if ((vHomHet[iIndex] >= 0) && (pobjProbeSet->getChromosome() == iChromosome))
            {
                vChromosomeGenotypeCalls[iProbeSetIndex] = vHomHet[iIndex];
                vChromosomePositions[iProbeSetIndex] = pobjProbeSet->getPosition();
                vChromosomeProbeSetIndexes[iProbeSetIndex] = iIndex;
                iProbeSetIndex++;
            }
        }
        int iStartIndex = 0;
        for (int iIndex = 1; (iIndex <= iChromosomeProbeSetCount); iIndex++)
        {
            if ((iIndex == iChromosomeProbeSetCount) || ((getProbeSets()->at(vChromosomeProbeSetIndexes[iIndex])->getPosition() - getProbeSets()->at(vChromosomeProbeSetIndexes[iIndex - 1])->getPosition()) > m_iLohCNSeparation))
            {
                int iSegmentProbeSetCount = ((iIndex - 1) - iStartIndex + 1);
                if (iSegmentProbeSetCount >= iMarkerCount)
                {
                    std::vector<char> vSegmentGenotypeCalls(iSegmentProbeSetCount);
                    std::vector<int> vSegmentPositions(iSegmentProbeSetCount);
                    std::vector<char> vLoh(iSegmentProbeSetCount);
                    for (int iSegmentIndex = 0; (iSegmentIndex < iSegmentProbeSetCount); iSegmentIndex++)
                    {
                        vSegmentGenotypeCalls[iSegmentIndex] = vChromosomeGenotypeCalls[iStartIndex + iSegmentIndex];
                        vSegmentPositions[iSegmentIndex] = vChromosomePositions[iStartIndex + iSegmentIndex];
                    }
                    lohFind(vSegmentGenotypeCalls, vSegmentPositions, iMarkerCount, iHetCutoff, m_iLohCNMinGenomicSpan, vLoh);
                    for (int iSegmentIndex = 0; (iSegmentIndex < iSegmentProbeSetCount); iSegmentIndex++)
                    {
                        getProbeSets()->at(vChromosomeProbeSetIndexes[iStartIndex + iSegmentIndex])->setLoh(vLoh[iSegmentIndex]);
                    }
                }
                iStartIndex = iIndex;
            }
        }
        }
        catch(...) {Verbose::out(1, "CNAnalysisMethodLOHCytoScan::run(...) Working on Chromosome " + ::getInt(iChromosome)); throw;}
    }
    imputeMissingLOHCalls();
   //    Verbose::progressEnd(1, "Done");
    Verbose::out(1, "CNAnalysisMethodLOHCytoScan::run(...) end");
}
   

void CNAnalysisMethodLOHCytoScan::imputeMissingLOHCalls()
{        
        std::vector<CNProbeSet*> vLOHProbeSets; 
        for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
        {
                CNProbeSet* pProbeSet = getProbeSets()->getAt(iIndex); 
                if( pProbeSet->processAsSNP() )
                {
                        vLOHProbeSets.push_back(pProbeSet);   
                } 
        }

        // This is code to dump out the loh status before the impute algorithm is invoked.
        bool testFlag = false;
        if(testFlag)
        {
                affx::TsvFile *tsv = new affx::TsvFile;
                tsv->writeTsv(getExperiment()->getExperimentName() + ".LohProbes.txt");
                tsv->defineColumn(    0,0,"position", affx::TSV_TYPE_DOUBLE);
                tsv->defineColumn(    0,1,"lohState", affx::TSV_TYPE_DOUBLE);
                tsv->defineColumn(    0,2,"Chromosome", affx::TSV_TYPE_DOUBLE);
                tsv->defineColumn(    0,3,"Genotype", affx::TSV_TYPE_DOUBLE);

            for (int iIndex = 0; iIndex < vLOHProbeSets.size(); iIndex++)
            {
                CNProbeSet* pProbeSet = vLOHProbeSets[iIndex];
                tsv->set(0,0,pProbeSet->getPosition());
                tsv->set(0,1,pProbeSet->getLoh());
                tsv->set(0,2,pProbeSet->getChromosome());
                tsv->set(0,3,pProbeSet->getGenotypeCall());
                tsv->writeLevel(0);
            }
            delete tsv;
        }


        int iLeftIndex=0;
        int iRightIndex=1;
                

        int iLeftLohStatus=0;  
        int iRightLohStatus=0;  
        CNProbeSet* pLeftProbeSet = NULL;      
        CNProbeSet* pRightProbeSet =  NULL;     

        while(iRightIndex < vLOHProbeSets.size())
        {


                pLeftProbeSet = vLOHProbeSets[iLeftIndex];      
                pRightProbeSet =  vLOHProbeSets[iRightIndex];     

                iLeftLohStatus=pLeftProbeSet->getLoh();  
                iRightLohStatus=pRightProbeSet->getLoh();  

                if(!pLeftProbeSet->lohCalled())
                {
                        if(!pRightProbeSet->lohCalled() ) 
                        {
                                // If we don't get an loh call on the right we advance.
                                iRightIndex++;
                                // If we have advanced beyond the collection of probesets we clean up since we never get back in the loop. 
                                if(iRightIndex == vLOHProbeSets.size())
                                {
                                        for(int iIndex = iLeftIndex; iIndex < iRightIndex; iIndex++)
                                        {
                                                vLOHProbeSets[iIndex]->setLoh(0);         
                                        }
                                }
                                continue;
                        } 
                        else
                        {
                        //  At this point we have just finished advancing the right pointer and thus 3 things can happen.  
                        //  RightIndex passes out of the chromosome or over the centromere.
                        //  Or hits a marker with valid loh call.
                        //  Note that we put the 3rd case last since we can get a valid loh call across a boundary but want to 
                        //  deal with it as in the first two cases.
                        //  Note that the right pointer cannot point outside the array since that is the condition on the while loop.
                                if( (vLOHProbeSets[iRightIndex]->getPosition() - vLOHProbeSets[iRightIndex-1]->getPosition()) > m_iLohCNSeparation)
                                {
                                        for(int iIndex = iLeftIndex; iIndex < iRightIndex; iIndex++)
                                        {
                                                vLOHProbeSets[iIndex]->setLoh(0);
                                        }
                                        // Note that iRightIndex points across a boundary so we do not set it's loh. ie. < not <= in above for loop. 
                                        iLeftIndex=iRightIndex;
                                        iRightIndex=iLeftIndex+1;
                                        // Note that here we can be pointing beyond the array bounds, but even though we never
                                        // get back in the loop we have cleaned up. i.e. set the last set of markers to have non-loh.
                                        continue;                             
                                }
                                if( vLOHProbeSets[iRightIndex]->getChromosome() != vLOHProbeSets[iRightIndex-1]->getChromosome() )
                                {
                                        for(int iIndex = iLeftIndex; iIndex < iRightIndex; iIndex++)
                                        {
                                                vLOHProbeSets[iIndex]->setLoh(0);
                                        }
                                        iLeftIndex=iRightIndex;
                                        iRightIndex=iLeftIndex+1;
                                        continue;                             
                                }
                                // We now have an loh call on the right but none on the left.
                                {
                                        for(int iIndex = iLeftIndex; iIndex <= iRightIndex; iIndex++)
                                        {
                                                vLOHProbeSets[iIndex]->setLoh(iRightLohStatus);
                                        }
                                        // We now move the left and right indices.  Note that how they move depends on whether or not we
                                        // are within a chromsome or if the right index is at the end of a chromsome.  In the first case
                                        // we advance the left to the right index since we need this points loh status to continue. In the
                                        // second case we advance both left and right to the first two positions in the next chromsome.
                                        if(iRightIndex+1 == vLOHProbeSets.size()) 
                                        {
                                                iRightIndex++;
                                                continue;
                                        }
                                        if
                                        (    
                                             ( vLOHProbeSets[iRightIndex+1]->getChromosome() != vLOHProbeSets[iRightIndex]->getChromosome() ) 
                                                 ||   
                                             ( (vLOHProbeSets[iRightIndex+1]->getPosition()- vLOHProbeSets[iRightIndex]->getPosition()) > m_iLohCNSeparation )  
                                        )
                                        {     
                                                iLeftIndex=iRightIndex+1;
                                                iRightIndex=iLeftIndex+1;
                                        }
                                        else
                                        { 
                                                iLeftIndex=iRightIndex;
                                                iRightIndex=iLeftIndex+1;
                                        }
                                        continue;                         
                                }

                        }   
                }
                else
                {
                        if(!pRightProbeSet->lohCalled() )
                        {
                                iRightIndex++;
                                if(iRightIndex == vLOHProbeSets.size())                                                  
                                {
                                        for(int iIndex = iLeftIndex; iIndex < iRightIndex; iIndex++)
                                        {
                                                vLOHProbeSets[iIndex]->setLoh(iLeftLohStatus);
                                        }  
                                        continue;                
                                }
                        }
                        else
                        {
                         // Again, at this point 3 things can happen.   
                         // RightIndex passes out of the chromosome or over the centromere.
                         // Or hits a marker with valid loh call.
                         // Note that we put the 3rd case last since we can get a valid loh call across a boundary but want to 
                         // deal with it as in the first two cases.
                              
                               if( (vLOHProbeSets[iRightIndex]->getPosition() - vLOHProbeSets[iRightIndex-1]->getPosition()) > m_iLohCNSeparation)
                                {
                                        for(int iIndex = iLeftIndex; iIndex < iRightIndex; iIndex++)
                                        {
                                                vLOHProbeSets[iIndex]->setLoh(iLeftLohStatus);
                                        }
                                         
                                        iLeftIndex=iRightIndex;
                                        iRightIndex=iLeftIndex+1;
                                        // Note that here we can be pointing beyond the array bounds, but even though we never
                                        // get back in the loop we have cleaned up. i.e. set the last set of markers to have non-loh.
                                        continue;
                                }
                                if( vLOHProbeSets[iRightIndex]->getChromosome() != vLOHProbeSets[iRightIndex-1]->getChromosome() )
                                {
                                        for(int iIndex = iLeftIndex; iIndex < iRightIndex; iIndex++)
                                        {
                                                vLOHProbeSets[iIndex]->setLoh(iLeftLohStatus);
                                        }
                                        iLeftIndex=iRightIndex;
                                        iRightIndex=iLeftIndex+1;
                                        continue;
                                }
                                // We now have an loh call on the left and the right.
                                {
                                        // At this point we have 2 cases.  The left and right probesets have the same loh state 
                                        // Or one has loh the other non-loh.
                                        if(iRightLohStatus==iLeftLohStatus)
                                        {
                                                for(int iIndex = iLeftIndex; iIndex <= iRightIndex; iIndex++)
                                                {
                                                        vLOHProbeSets[iIndex]->setLoh(iLeftLohStatus);
                                                }
                                        }
                                        else
                                        {
                                                if(iRightLohStatus==0)
                                                {
                                                        for(int iIndex = iLeftIndex+1; iIndex <= iRightIndex; iIndex++)
                                                        {
                                                                vLOHProbeSets[iIndex]->setLoh(0);
                                                        }
                                                }
                                                else
                                                {
                                                        for(int iIndex = iLeftIndex; iIndex < iRightIndex; iIndex++)
                                                        {
                                                                vLOHProbeSets[iIndex]->setLoh(0);
                                                        }
                                                }
                                        }
                                        if(iRightIndex+1 == vLOHProbeSets.size())
                                        {
                                                iRightIndex++;
                                                continue;
                                        }
                                       if
                                        (
                                             ( vLOHProbeSets[iRightIndex+1]->getChromosome() != vLOHProbeSets[iRightIndex]->getChromosome() )
                                                 ||
                                             ( (vLOHProbeSets[iRightIndex+1]->getPosition()- vLOHProbeSets[iRightIndex]->getPosition()) > m_iLohCNSeparation )
                                        )
                                        { 
                                                iLeftIndex=iRightIndex+1;
                                                iRightIndex=iLeftIndex+1;
                                        }
                                        else
                                        { 
                                                iLeftIndex=iRightIndex;
                                                iRightIndex=iLeftIndex+1;
                                        }
                                        continue;
                                }
                        }
               }
        }
}



/**
 * @brief Preprocess the genotype calls to get the marker count cutoff and the het call cutoff
 * @param std::vector<char>& - The vector of hom and/or het calls
 * @param int& - The marker count cutoff to calculate
 * @param int& - The het call cutoff to calculate
 * @return bool - Was the function successful or not (true = successful)
 */
bool CNAnalysisMethodLOHCytoScan::lohPreProcessing(std::vector<char>& vHomHet, int& iMarkerCount, int& iHetCutoff)
{
    int iHetCount = 0;
    int iCallCount = 0;
    int iNoCallCount = 0;
    double dMaximumGenotypeConfidence = 0;
    int iSnpCount = 0;
    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
        pobjProbeSet->setLoh(numeric_limits<float>::quiet_NaN());
        if (pobjProbeSet->processAsSNP()) {iSnpCount++;}
        vHomHet[iIndex] = pobjProbeSet->getGenotypeCall();
        if (vHomHet[iIndex] == (char)2) {vHomHet[iIndex] = (char)0;}
        if (pobjProbeSet->getGenotypeConfidence() > m_dLohCNNoCallThreshold) {vHomHet[iIndex] = (char)-1;}
        if (vHomHet[iIndex] == 1) {iHetCount++;}
        if (vHomHet[iIndex] >= 0)
        {
            iCallCount++;
            dMaximumGenotypeConfidence = Max((float)dMaximumGenotypeConfidence, pobjProbeSet->getGenotypeConfidence());
        }
        else
        {
            if (pobjProbeSet->processAsSNP()) {iNoCallCount++;}
        }
    }
    double dErrorRate = m_dLohCNErrorRate;
    if (((dMaximumGenotypeConfidence * 1.1) >= m_dLohCNNoCallThreshold) && (m_dLohCNNoCallThreshold == 0.05))
    {
        dErrorRate = ((double)iNoCallCount / (double)iSnpCount);
        dErrorRate = (0.002926 + 6.607106 * dErrorRate);
        dErrorRate = ceil(dErrorRate * 100) / 100.0;
        if (dErrorRate < m_dLohCNErrorRate) {dErrorRate = m_dLohCNErrorRate;}
    }
    double dHetRate = ((double)iHetCount / (double) iCallCount);
    if (dHetRate > 0.5) {dHetRate = 0.5;}
    if (dHetRate < 0.2) {dHetRate = 0.2;}
    int iIterations = 1;
    int iStart = 20;
    int n = 0;
    while (iIterations < 3)
    {
        iHetCutoff = AffxStatistics::binoinv((1 - m_dLohCNAlpha), iStart, dErrorRate);
        double a = AffxStatistics::norminv(m_dLohCNBeta);
        double sqrtn = (-a * 0.5 * sqrt(1 - dHetRate) / sqrt(dHetRate) + sqrt(a * a * (1 - dHetRate) / (4 * dHetRate) + iHetCutoff / dHetRate));
        n = (int)ceil(sqrtn * sqrtn);
        if (n == iStart) {break;}
        else {iStart = n;}
        iIterations++;
    }
    iMarkerCount = Max(n, m_iLohCNMinimumMarkerCount);
    if ((iMarkerCount % 2) == 0) {iMarkerCount++;}
    if (iMarkerCount < 3)
    {
        throw(Except("MarkerCount must be >= 3.")); return false;
    }



    // The hetRate is used in the CNAnalysisMethodSegment method to assign a confidence to the loh determination.  In order
    // that this value is accessible there we pass it to the CNExperiment object. Note that the het frequency is also 
    // computed and stored in a postprocessing function in the CNCytoEngine.  The two computations are the same.
    // We also need the dErrorRate when calculating confidences.  This value is also passed to the experiment object and then accessed
    // when segmenting.
    getExperiment()->setHetFrequency(dHetRate);
    getExperiment()->setPError(dErrorRate);

    return true;
}

/**
 * @brief Find the runs off LOH
 * @param std::vector<char>& - The genotype calls
 * @param std::vector<int>& - The chromosome positions
 * @param int - The marker count cutoff
 * @param int - The het call cutoff
 * @param int - The minimum genomic span for a run of LOH
 * @param std::vector<char>& - The resulting LOH calls (1 = LOH, 0 = !LOH)
 */
bool CNAnalysisMethodLOHCytoScan::lohFind(      std::vector<char>& vGenotypeCalls, 
                                                std::vector<int>& vPositions, 
                                                int iMarkerCount, 
                                                int iHetCutoff, 
                                                int iMinGenomicSpan, 
                                                std::vector<char>& vLoh)
{
    int iWindowCenter = (int)floor(iMarkerCount / 2.0);
    if ((iWindowCenter * 2) != iMarkerCount - 1) {throw(Except("Window size must be odd.")); return false;}
    int iCurrentState = 0;
    for (int iIndex = 0; (iIndex < (int)vLoh.size()); iIndex++)
    {
        vLoh[iIndex] = 0;
    }
    int iStateLoh = 1;
    int iStateNonLoh = 0;
    std::vector<int> vHetCount(vGenotypeCalls.size());
    for (int iIndex = 0; (iIndex < (int)vGenotypeCalls.size()); iIndex++)
    {
        vHetCount[iIndex] = 0;
        if ((vGenotypeCalls[iIndex] != 0) && (vGenotypeCalls[iIndex] != 1))
        {
            throw(Except("LOH Failed")); return false;
        }
    }
    int iCurrentCallIndex = iWindowCenter;
    int iWindowOffset = iWindowCenter - iHetCutoff - 1;
    if (iWindowOffset < 1) {throw(Except("LOH Failed")); return false;}
    int iLastCallIndex = (int)vGenotypeCalls.size() - iWindowCenter - 1;
    int iCount = 0;
    for (int iIndex = 0; (iIndex < iMarkerCount); iIndex++)
    {
        if (vGenotypeCalls[iIndex] == 1) {iCount++;}
    }
    vHetCount[iCurrentCallIndex] = iCount;

    while (iCurrentCallIndex < iLastCallIndex)
    {
        int iPrevCallIndex = iCurrentCallIndex;
        if (vHetCount[iPrevCallIndex] >= iHetCutoff)
        {
            if (iCurrentState == 1)
            {
                // in a LOH region transitioning back to non-LOH
                // Change back to non-LOH state
                // really now need to fill in from
                // iPrevCallIndex-1:iPrevCallIndex+iWindowOffset -1 as also LOH
                iCurrentState = 0;
                int iMinCallIndex = Max(0, iPrevCallIndex - 1);
                int iMaxCallIndex = iPrevCallIndex - 1;
                int iHetsSeen = 0;
                // starting from the window center go out as far as the
                // second het or the end of the window which ever occurs first
                // Take special care to stop at the last hom
                while ((iHetsSeen < 2) && (iMaxCallIndex < iPrevCallIndex + iWindowOffset - 1) && (iMaxCallIndex < iLastCallIndex))
                {
                    iMaxCallIndex = Min(iMaxCallIndex + 1, iPrevCallIndex + iWindowOffset - 1);
                    if (vGenotypeCalls[iMaxCallIndex] == 1)
                    {
                        iHetsSeen++;
                        if (iHetsSeen == 2)
                        {
                            iMaxCallIndex--;
                            if (vGenotypeCalls[iMaxCallIndex] == 1)
                            {
                                iMaxCallIndex--;
                            }
                        }
                    }
                }
                for (int iIndex = iMinCallIndex; (iIndex <= iMaxCallIndex); iIndex++)
                {
                    vLoh[iIndex] = iStateLoh;
                }
                if (iLastCallIndex >= (iPrevCallIndex + iWindowOffset))
                {
                    for (int iIndex = (iPrevCallIndex + 1); (iIndex <= (iPrevCallIndex + iWindowOffset)); iIndex++)
                    {
                        vHetCount[iIndex] = vHetCount[iIndex - 1] - vGenotypeCalls[iIndex - 1 - iWindowCenter] + vGenotypeCalls[iIndex + iWindowCenter];
                    }
                }
            }
            else
            {
                // in a non-LOH region and not transitioning, move to the next marker
                iCurrentCallIndex++;
                vHetCount[iCurrentCallIndex] = vHetCount[iPrevCallIndex] - vGenotypeCalls[iPrevCallIndex - iWindowCenter] + vGenotypeCalls[iCurrentCallIndex + iWindowCenter];
            }
        }
        else
        {
            if (iCurrentState == 0)
            {
                // if we are transitioning to the LOH state then we should
                // mark back to the start of the Window as LOH
                iCurrentState = 1;
                int iMinCallIndex = iPrevCallIndex;
                int iHetsSeen = 0;
                // starting from the window center go out as far as the
                // second het or the end of the window which ever occurs first
                // Take special care to stop at the last hom
                while ((iHetsSeen < 2) && (iMinCallIndex > (iPrevCallIndex - iWindowOffset)) && (iMinCallIndex > 0))
                {
                    iMinCallIndex = Max(iMinCallIndex - 1, 0);
                    if (vGenotypeCalls[iMinCallIndex] == 1)
                    {
                        iHetsSeen++;
                        if (iHetsSeen == 2)
                        {
                            iMinCallIndex++;
                            if (vGenotypeCalls[iMinCallIndex] == 1)
                            {
                                iMinCallIndex++;
                            }
                        }
                    }
                }
                for (int iIndex = iMinCallIndex; (iIndex <= iCurrentCallIndex); iIndex++)
                {
                    vLoh[iIndex] = iStateLoh; // we know definitively that these are LOH
                }
                // we could also make a guess about those that are
                // following but instead we will move forward
                // conservatively
            }
            vLoh[iCurrentCallIndex] = iStateLoh;
            iCurrentCallIndex++;
            vHetCount[iCurrentCallIndex] = vHetCount[iPrevCallIndex] - vGenotypeCalls[iPrevCallIndex - iWindowCenter] + vGenotypeCalls[iCurrentCallIndex + iWindowCenter];
        }
    }
    // deal with markers that don't have a full window around them
    for (int iIndex = 0; (iIndex < iWindowCenter); iIndex++)
    {
        vLoh[iIndex] = vLoh[iWindowCenter + 1];
    }
    for (int iIndex = iCurrentCallIndex; (iIndex < (int)vLoh.size()); iIndex++)
    {
        vLoh[iIndex] = vLoh[iCurrentCallIndex-1];
    }

    // We move up the segment until iBeginIndex points to the first marker having iStateLoh.
    int iBeginIndex =0;
    while(iBeginIndex < vLoh.size() && vLoh[iBeginIndex] != iStateLoh )
    {
        iBeginIndex++;
    } 
    if(iBeginIndex == vLoh.size())
    {
        return true;
    }


    int iSegmentStartIndex=iBeginIndex;
    int iStartPosition=vPositions[iBeginIndex];
    int iPreviousState = vLoh[iBeginIndex];
    int iPreviousPosition=vPositions[iBeginIndex];

    for(int iIndex=iBeginIndex+1; iIndex < vLoh.size(); iIndex++)
    {
    
        int iPresentPosition = vPositions[iIndex];
        int iPresentState = vLoh[iIndex];

        // Case 1:  We have hit the end of the segment.   
        if(iIndex == vPositions.size()-1) 
        {
            // SubCase a: In this case the state of the final marker is the same as the state of the penultimate marker
            // and this state is LOH but the length of the region is small.
            if(	    iPresentState == iPreviousState  && 
                    iPresentState == iStateLoh          && 
                    iPresentPosition-iStartPosition < iMinGenomicSpan) 
            {
                for(int iNewIndex = iSegmentStartIndex; iNewIndex < vPositions.size(); iNewIndex++)
                {
                    vLoh[iNewIndex] = iStateNonLoh;
                }
            }
            // SubCase b: In this case the state of the final marker is non-loh, but the penultimate is LOH.
            // Also the region of LOH is small we set it to non-loh
            if(     iPreviousState == iStateLoh          && 
                    iPresentState == iStateNonLoh          && 
                    iPreviousPosition-iStartPosition < iMinGenomicSpan) 
            {
                for(int iNewIndex = iSegmentStartIndex; iNewIndex < vPositions.size(); iNewIndex++)
                {
                    vLoh[iNewIndex] = iStateNonLoh;
                }
            }
            // SubCase c:  In this case the final marker is LOH and penultimate is non-Loh. 
            // We set the final marker to non-LOH.
            if(     iPreviousState == iStateNonLoh          && 
                    iPresentState == iStateLoh) 
            {
                vLoh[iPresentPosition] = iStateNonLoh;
            }        
            break;
        }

        if(	iPresentState != iPreviousState &&
                iPresentState == iStateNonLoh      &&
                iPreviousPosition-iStartPosition >  iMinGenomicSpan)
        { 
            int i=0;
            i++;
        } 

        // Case 2:  We have changed state at the end of an LOH region into a non-LOH region 
        if(	iPresentState != iPreviousState &&
                iPresentState == iStateNonLoh      &&
                iPreviousPosition-iStartPosition <  iMinGenomicSpan)
        { 
            for(int iNewIndex = iSegmentStartIndex; iNewIndex < iIndex; iNewIndex++)
            {
                vLoh[iNewIndex] = iStateNonLoh;
            }
        } 
       
        // Case 3:  We have changed state at the end of a non-Loh region into an LOH region. 
        // So we reset the iStartPosition and iSegmentStartIndex.
        if(	iPresentState != iPreviousState &&
                iPresentState == iStateLoh)
        {
            iStartPosition = iPresentPosition;
            iSegmentStartIndex = iIndex; 
        } 

        iPreviousPosition = iPresentPosition;
        iPreviousState = iPresentState;
    }

    return true;
}

void CNAnalysisMethodLOHCytoScan::determineLocalProbeSets()
{
    if (m_bLocalProbeSetsDetermined) {return;}
    int iNumberOfProbeSets = CNAnalysisMethod::getProbeSets()->getCount();
    for (int iIndex = 0; (iIndex < iNumberOfProbeSets); iIndex++)
    {
        if (CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAsSNP() )
        {
             getProbeSets()->add( CNAnalysisMethod::getProbeSets()->getAt(iIndex));
        }
    }
    m_bLocalProbeSetsDetermined = true;
}





