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
 * @file CNAnalysisMethodSegment.cpp
 *
 * @brief This file contains the CNAnalysisMethodSegment class members.
 */
#include "copynumber/CNAnalysisMethodSegment.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "util/AffxStatistics.h"
#include "util/Guid.h"
//


/**
 * @brief Constructor
 */
CNAnalysisMethodSegment::CNAnalysisMethodSegment()
{
}

/**
 * @brief Create new segments
 * @param int - The segment type to create
 * @param CNSegmentArray& - The segment vector to load
 * @param int - The affx sex code
 * @param int - The X chromosome numeric value
 * @param int - The Y chromosome numeric value
 * @return int - The number of new segments created
 */
int CNAnalysisMethodSegment::newSegments(    int iSegmentType,
                                              CNProbeSetArray* vProbeSets,
                                             int iMinSegSeparation   )
{
    AffxMultiDimensionalArray<int> vMarkerDistances(vProbeSets->getCount());
    int iPrevMarkerPosition = 0;
    if (vProbeSets->getCount() == 0) {return 0;}
    int iSegmentCount = 0;
    int iMarkerCount = 0;
    double dSumConfidence = 0.0;
    double dConfidence=0.0;
    int iHetsInSegment=0;
    int iHomsInSegment=0;
    CNSegment* pobjSegment = NULL;
    CNProbeSet* pobjProbeSet = NULL;
    int iStartIndex = 0;
    // If Segment type loh-segment, loh-cyto2, cn-neutral-log, normal-diploid
    if ((iSegmentType == 2) || (iSegmentType == 3) || (iSegmentType == 4) || (iSegmentType == 5))
    {
        for (;(iStartIndex < vProbeSets->getCount()); iStartIndex++)
        {
            pobjProbeSet = vProbeSets->getAt(iStartIndex);
            if ((double)pobjProbeSet->getLoh() == (double)pobjProbeSet->getLoh()) // Not = NaN
            {
                break;
            }
        }
    }
    int iCurrCNState = -1;
    bool cytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
    if (iStartIndex < vProbeSets->getCount())
    {
        pobjProbeSet = vProbeSets->getAt(iStartIndex);
        if (cytoScanHD)
        {
            // Case 1: SNP marker with CN=-1 sits at the beginning of the chromosome
            unsigned char currentChromosome = pobjProbeSet->getChromosome();
            for (int iIndex = iStartIndex; (iIndex < vProbeSets->getCount()); iIndex++)
            {
                if (currentChromosome != vProbeSets->getAt(iIndex)->getChromosome()) break;
                if (vProbeSets->getAt(iIndex)->getCNState() != -1)
                {
                    iCurrCNState = vProbeSets->getAt(iIndex)->getCNState();
                    break;
                }
            }
        }
        char cCall = pobjProbeSet->getSegmentTypeCall(iSegmentType, getExperiment()->getCNCallGenderAsInt(), m_iXChromosome, m_iYChromosome);
        float fConfidence = pobjProbeSet->getSegmentTypeConfidence(iSegmentType);
        bool bPar = pobjProbeSet->isPseudoAutosomalRegion();
        if (cytoScanHD)
        {
            // Override cCall in the two cases where CN state makes a difference.
		    if (iSegmentType == 4) cCall = (((iCurrCNState == 2) && (pobjProbeSet->getLoh() == 1)) ? (char)1 : (char)0); // CNNeutralLOH
		    if (iSegmentType == 5) cCall = (((iCurrCNState == 2) && (pobjProbeSet->getLoh() == 0)) ? (char)1 : (char)0); // Normal Diploid
        }

        pobjSegment = new CNSegment;
        pobjSegment->setSegmentType(iSegmentType);
        pobjSegment->setChromosome(pobjProbeSet->getChromosome());
        pobjSegment->setStartPosition(pobjProbeSet->getPosition());
        pobjSegment->setEndPosition(pobjProbeSet->getPosition());
        pobjSegment->setCall(cCall);
        pobjSegment->setConfidence((float)0.5);
        pobjSegment->setPseudoAutosomalRegion(bPar);
        iMarkerCount = 0;
        dSumConfidence = 0;
        iPrevMarkerPosition = pobjProbeSet->getPosition();

        for (int iIndex = iStartIndex; (iIndex < vProbeSets->getCount()); iIndex++)
        {
            pobjProbeSet = vProbeSets->getAt(iIndex);
            cCall = pobjProbeSet->getSegmentTypeCall(iSegmentType, getExperiment()->getCNCallGenderAsInt(), m_iXChromosome, m_iYChromosome);
            fConfidence = pobjProbeSet->getSegmentTypeConfidence(iSegmentType);
            bPar = pobjProbeSet->isPseudoAutosomalRegion();
			
            // Case 2: If a SNP marker falls within a CN segment (including CN or SNP), assign SNP marker with the same CN state
            // Case 3: If a SNP marker falls outside any CN segment, assign CN state for this SNP marker based on the CN segment just right before it.
			if (pobjProbeSet->getCNState() != -1) {iCurrCNState = pobjProbeSet->getCNState();}

			if (cytoScanHD)
			{
                // Override cCall in the two cases where CN state makes a difference.
		        if (iSegmentType == 4) cCall = (((iCurrCNState == 2) && (pobjProbeSet->getLoh() == 1)) ? (char)1 : (char)0); // CNNeutralLOH
		        if (iSegmentType == 5) cCall = (((iCurrCNState == 2) && (pobjProbeSet->getLoh() == 0)) ? (char)1 : (char)0); // Normal Diploid
			}
            if ((iSegmentType == 2) || (iSegmentType == 3) || (iSegmentType == 4) || (iSegmentType == 5))
            {
                if ((double)pobjProbeSet->getLoh() != (double)pobjProbeSet->getLoh()) // NaN
                {
                    continue;
                }
            }
            bool newSegmentTestFlag = false;
            if(iSegmentType == 2)
            { 
                newSegmentTestFlag = (	
                    (pobjProbeSet->getChromosome() != pobjSegment->getChromosome()) 
                 || (cCall != pobjSegment->getCall()) 
                 || ((getExperiment()->getCNCallGenderAsInt() == affx::Male) && (iSegmentType == 1) && (bPar != pobjSegment->isPseudoAutosomalRegion()))
                 || (fabs((float)iPrevMarkerPosition - (float) (pobjProbeSet->getPosition())) > (float) iMinSegSeparation) 
                                     );
            }
            else
            {
                newSegmentTestFlag = (
                    (pobjProbeSet->getChromosome() != pobjSegment->getChromosome()) 
                 || (cCall != pobjSegment->getCall()) 
                 || ((getExperiment()->getCNCallGenderAsInt() == affx::Male) && (iSegmentType == 1) && (bPar != pobjSegment->isPseudoAutosomalRegion()))
                                     ); 
            }

            if (newSegmentTestFlag)
            {
                pobjSegment->setMeanMarkerDistance((float)vMarkerDistances.mean(iMarkerCount - 1));
                pobjSegment->setSegmentName(affxutil::Guid::GenerateNewGuid());
                pobjSegment->setMarkerCount(iMarkerCount);
                pobjSegment->setConfidence((float)(dSumConfidence / (double)iMarkerCount));

                if(cytoScanHD)
                {
                    if(iSegmentType==2)
                    {
                        dConfidence = assignConfidenceSigmoidal(	iHetsInSegment,
                                                                        iHetsInSegment + iHomsInSegment,
                                                                        getExperiment()->getHetFrequency(),
                                                                        getExperiment()->getPError(),
                                                                        3,
                                                                        10);
                        pobjSegment->setConfidence((float) dConfidence); 
                    }
                    if(iSegmentType==4)
                    {
                        dConfidence = assignConfidenceSigmoidal(	iHetsInSegment,
                                                                        iHetsInSegment + iHomsInSegment,
                                                                        getExperiment()->getHetFrequency(),
                                                                        getExperiment()->getPError(),
                                                                        3,
                                                                        10);
                        pobjSegment->setConfidence((float) (dConfidence/2.0));
                    }

               }

                m_vSegments.add(pobjSegment);
                iSegmentCount++;

                if (cytoScanHD && pobjProbeSet->getChromosome() != pobjSegment->getChromosome())
                {
                    // Case 1: SNP marker with CN=-1 sits at the beginning of the chromosome
                    iCurrCNState = -1;
                    unsigned char currentChromosome = pobjProbeSet->getChromosome();
                    for (int iInner = iIndex; (iInner < vProbeSets->getCount()); iInner++)
                    {
                        if (currentChromosome != vProbeSets->getAt(iIndex)->getChromosome()) break;
                        if (vProbeSets->getAt(iInner)->getCNState() != -1)
                        {
                            iCurrCNState = vProbeSets->getAt(iInner)->getCNState();
                            break;
                        }
                    }
                     // Override cCall in the two cases where CN state makes a difference.
	                if (iSegmentType == 4) cCall = (((iCurrCNState == 2) && (pobjProbeSet->getLoh() == 1)) ? (char)1 : (char)0); // CNNeutralLOH
	                if (iSegmentType == 5) cCall = (((iCurrCNState == 2) && (pobjProbeSet->getLoh() == 0)) ? (char)1 : (char)0); // Normal Diploid
                }

                pobjSegment = new CNSegment;
                pobjSegment->setSegmentType(iSegmentType);
                pobjSegment->setChromosome(pobjProbeSet->getChromosome());
                pobjSegment->setStartPosition(pobjProbeSet->getPosition());
                pobjSegment->setCall(cCall);
                pobjSegment->setConfidence((float)0.5);
                pobjSegment->setPseudoAutosomalRegion(bPar);
                iHetsInSegment=0;
                iHomsInSegment=0;
                iMarkerCount = 0;
                dSumConfidence = 0;
                iPrevMarkerPosition = pobjProbeSet->getPosition();
            }
            iMarkerCount++;

            if(pobjProbeSet->getGenotypeCall() == (char) 2 || pobjProbeSet->getGenotypeCall() == (char) 0)
            {
                iHomsInSegment++;
            }

            if(pobjProbeSet->getGenotypeCall() == (char) 1 )
            {
                iHetsInSegment++;
            }


            dSumConfidence += fConfidence;
            pobjSegment->setEndPosition(pobjProbeSet->getPosition());
            if (iMarkerCount > 1) {vMarkerDistances.set((iMarkerCount - 2), (pobjProbeSet->getPosition() - iPrevMarkerPosition));}
            iPrevMarkerPosition = pobjProbeSet->getPosition();
        }
    }
    if (pobjSegment != NULL)
    {
        pobjSegment->setConfidence((float)(dSumConfidence / (double)iMarkerCount));

        if(cytoScanHD)
        { 
            if(iSegmentType==2)
            {
                dConfidence = assignConfidenceSigmoidal(	iHetsInSegment,
                                                                iHetsInSegment + iHomsInSegment,
                                                                getExperiment()->getHetFrequency(),
                                                                getExperiment()->getPError(),
                                                                3,
                                                                10);
                 pobjSegment->setConfidence((float) dConfidence);
            }
            if(iSegmentType==4)
            {
                dConfidence = assignConfidenceSigmoidal(	iHetsInSegment,
                                                                iHetsInSegment + iHomsInSegment,
                                                                getExperiment()->getHetFrequency(),
                                                                getExperiment()->getPError(),
                                                                3,
                                                                10);
                pobjSegment->setConfidence((float) (dConfidence/2.0));
            }
        }  

        pobjSegment->setMeanMarkerDistance((float)vMarkerDistances.mean(iMarkerCount - 1));
        pobjSegment->setMarkerCount(iMarkerCount);
        pobjSegment->setSegmentName(affxutil::Guid::GenerateNewGuid());
        m_vSegments.add(pobjSegment); iSegmentCount++;
    }
    return iSegmentCount;
}




