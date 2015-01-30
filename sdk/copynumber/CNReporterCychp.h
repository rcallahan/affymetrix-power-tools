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

#ifndef _CNReporterCychp_H_
#define _CNReporterCychp_H_
/**
 * @file CNReporterCychp.h
 *
 * @brief This header contains the CNReporterCychp class definition.
 */

#include "copynumber/CNReporter.h"

using namespace std;

/**
 * @brief  A class for storing chromosome data
 */
class CNChromosome
{
public:
    unsigned char m_ucChromosome;
    AffxString m_strDisplay;
    unsigned int m_uiStartIndex;
    unsigned int m_uiMarkerCount;
    float m_fMinSignal;
    float m_fMaxSignal;
    float m_fMedianSignal;
    float m_fMedianCnState;
    float m_fHomFrequency;
    float m_fHetFrequency;
    float m_fMosaicism;
    float m_fChromosomeLOH;
};

/**
 * @brief  The SNP data CN Reporter.
 */
class CNReporterCychp : public CNReporter
{
private:
    unsigned int m_uiSegmentID;

public:
    CNReporterCychp() {}
    virtual ~CNReporterCychp() {}

    virtual void defineOptions(BaseEngine& e);
    virtual void checkOptions(BaseEngine& e);
    virtual void run();

private:
    AffxString getFileName();
    AffxString getCelFileName();
    void setArrayName();
    void writeCychpFile(const AffxString& strFileName);
    int getMaxProbeSetNameLength();
    void addParentHeader(affymetrix_calvin_io::GenericDataHeader* pHeader, AffxString& strFileName);
    void addCNReferenceHeader(affymetrix_calvin_io::GenericDataHeader* pHeader);
    void addCNReferenceSelectedInfo(affymetrix_calvin_io::GenericDataHeader* pHeader);
    void addHeaderParams(affymetrix_calvin_io::GenericDataHeader* pHeader);
    void addHeaderClientMetaInfo(affymetrix_calvin_io::GenericDataHeader* pHeader);
    void addHeaderAlgorithmParams(affymetrix_calvin_io::GenericDataHeader* pHeader);
    void addHeaderChipSummaryParams(affymetrix_calvin_io::GenericDataHeader* pHeader);
    void loadChromosomes(AffxArray<CNChromosome>& arChromosomes);
    void addChromosomesSummaryDataSetHeader(affymetrix_calvin_io::DataGroupHeader& group, int iRowCount);
    void addProbeSetsCopyNumberDataSetHeader(affymetrix_calvin_io::DataGroupHeader& group, int iRowCount);
    void addProbeSetsAllelePeaksDataSetHeader(affymetrix_calvin_io::DataGroupHeader& group, int iRowCount);
    void addAlgorithmDataMarkerABSignalDataSetHeader(affymetrix_calvin_io::DataGroupHeader& group, int iRowCount);
    void addSegmentsDataSets(affymetrix_calvin_io::DataGroupHeader& group);
    void setupSegmentsDataSetHeader(affymetrix_calvin_io::DataGroupHeader& group, int iRowCount, short nSegmentType);
    void addGenotypingDataSetHeader(affymetrix_calvin_io::DataGroupHeader& group, int iRowCount);
    void writeDataGroupHeader(affymetrix_calvin_io::GenericFileWriter* pWriter, int iGroupIndex, std::vector<int>& vOffsets);
    void writeChromosomesSummaryDataSet(affymetrix_calvin_io::GenericFileWriter* pWriter, int iGroupIndex, int iSetIndex, std::vector<int>& vOffsets, AffxArray<CNChromosome>& arChromosomes);
    void writeProbeSetsCopyNumberDataSet(affymetrix_calvin_io::GenericFileWriter* pWriter, int iGroupIndex, int iSetIndex, std::vector<int>& vOffsets, CNProbeSetArray& arProbeSets);
    void writeProbeSetsAllelePeaksDataSet(affymetrix_calvin_io::GenericFileWriter* pWriter, int iGroupIndex, int iSetIndex, std::vector<int>& vOffsets, CNProbeSetArray& arProbeSets);
    void writeSegmentsDataSet(affymetrix_calvin_io::GenericFileWriter* pWriter, int iGroupIndex, int iSetIndex, std::vector<int>& vOffsets, CNSegmentArray& arSegments, short nSegmentType);
    void writeAlgorthmDataMarkerABSignalDataSet(affymetrix_calvin_io::GenericFileWriter* pWriter, int iGroupIndex, int iSetIndex, std::vector<int>& vOffsets, CNProbeSetArray& arProbeSets);
    void writeGenInfoDataSet(affymetrix_calvin_io::GenericFileWriter* pWriter, int iGroupIndex, int iSetIndex, std::vector<int>& vOffsets, CNProbeSetArray& arProbeSets);
};

#endif


