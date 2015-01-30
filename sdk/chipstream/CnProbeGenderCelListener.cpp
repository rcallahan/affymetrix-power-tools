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
#include "chipstream/CnProbeGenderCelListener.h"
//
#include "chipstream/EngineUtil.h"
//
#include "file/TsvFile/TsvFile.h"
//

using namespace affx;
using namespace std;

/** 
 * Process another cel files worth of data.
 * @param cel - Filehandle to open cel file.
 */
void CnProbeGenderCelListener::newChip(affymetrix_fusion_io::FusionCELData *cel) {
  std::vector<ChipSummary::Metric> metrics;
  m_SummaryStats.push_back(metrics);
  m_CelNames.push_back(cel->GetFileName());
  float fRatio = 0;

  affx::Gender gender = callGender(cel, m_ChrXProbes, m_ChrYProbes, fRatio);

  if (m_ZWGenderCalling) {
      gender = affx::flipMaleFemaleGender(gender);
  }

  m_Genders.push_back(gender);
  m_vRatios.push_back(fRatio);

  m_SummaryStats[m_SummaryStats.size()-1].push_back(ChipSummary::Metric(m_GenderName + "_gender", affx::getGenderString(gender)));
  setValid(true);
}

static inline double sumIntensities(affymetrix_fusion_io::FusionCELData *cel, 
                                    const vector< vector<probeid_t> > &probes)
{
    double ret = 0.0;

    vector< wstring > channelNames = cel->GetChannels();
    // If this is a GCOS CEL file, then make up any old name for a
    // single channel -- subsequent call to SetActiveDataGroup will be
    // a no-op in this case, and GetIntensity will proceed as normal.
    if (channelNames.empty()) {
      std::wstring temp = L"Default Group";
      channelNames.push_back(temp);
    }
    for ( unsigned int i = 0; i < probes.size(); i++ ) {
        if ( probes[i].size() ) {
            // sum intensities for channel i
            cel->SetActiveDataGroup( channelNames[i] );
            for (vector<probeid_t>::const_iterator it = probes[i].begin(); it != probes[i].end(); ++it) {
	      ret += cel->GetIntensity(*it);
            }
        }
    }
    return ret;
}

affx::Gender CnProbeGenderCelListener::callGender(affymetrix_fusion_io::FusionCELData *cel, 
                                                  const vector< vector<probeid_t> > &chrxProbes,
                                                  const vector< vector<probeid_t> > &chryProbes,
                                                  float& fRatio)
{
	///@todo check to see if chrxProbes can be a non-zero vector of empty channels
    unsigned int num_Xprobes = 0, num_Yprobes = 0;
    for ( unsigned int i = 0 ; i < chrxProbes.size(); i ++ )
        num_Xprobes += chrxProbes[i].size();
    for ( unsigned int i = 0 ; i < chryProbes.size(); i ++ )
        num_Yprobes += chryProbes[i].size();

    if (num_Xprobes == 0) {return affx::UnknownGender;}
    if (num_Yprobes == 0) {return affx::UnknownGender;}

    float meanX = sumIntensities(cel, chrxProbes)/num_Xprobes;
    float meanY = sumIntensities(cel, chryProbes)/num_Yprobes;

    m_SummaryStats[m_SummaryStats.size()-1].push_back(ChipSummary::Metric(m_GenderName + "_gender_mean"+m_DipSexChrLetter,meanX));
    m_SummaryStats[m_SummaryStats.size()-1].push_back(ChipSummary::Metric(m_GenderName + "_gender_mean"+m_HapSexChrLetter,meanY));

    affx::Gender call = affx::UnknownGender;
    if (meanX == 0.0) {
        fRatio = numeric_limits<double>::quiet_NaN();
        m_SummaryStats[m_CelNames.size()-1].push_back(ChipSummary::Metric(m_GenderName + "_gender_ratio",(double)0.0));
        //strm << "Gender call for sample " << (m_CelNames.size() + 1) << ": meanY=" << meanY << "; meanX=" << meanX << "; ratio=Undefined" ;
        //Verbose::out(2, strm.str());
        return call;
    }

    fRatio = meanY / meanX;

    m_SummaryStats[m_CelNames.size()-1].push_back(ChipSummary::Metric(m_GenderName + "_gender_ratio",fRatio));

    if (fRatio < m_FemaleThresh) {
        call = affx::Female;
    } else if (fRatio > m_MaleThresh) {
        call = affx::Male;
    }
    std::stringstream strm;
    strm << "Gender call for sample " << cel->GetFileName() << ": mean" << m_HapSexChrLetter << "=" << meanY << "; mean" << m_DipSexChrLetter << "=" << meanX << "; ratio=" << fRatio;
    Verbose::out(2, strm.str());
    return call;
}

CnProbeGenderCelListener::CnProbeGenderCelListener(const string &chrXProbeFile, 
                                                   const string &chrYProbeFile,
                                                   double femaleThresh, 
                                                   double maleThresh, 
                                                   bool chrZW) {
    m_ZWGenderCalling = chrZW;

    if (m_ZWGenderCalling) {
        m_DipSexChrLetter = "Z";
        m_HapSexChrLetter = "W";
    }
    else {
        m_DipSexChrLetter = "X";
        m_HapSexChrLetter = "Y";
    }

    if (chrXProbeFile.empty())
        Err::errAbort("No chr"+m_DipSexChrLetter+" probe id file specified. Unable to run CnProbeGenderCelListener.");

    if (chrYProbeFile.empty())
        Err::errAbort("No chr"+m_HapSexChrLetter+" probe id file specified. Unable to run CnProbeGenderCelListener.");

    EngineUtil::readProbeFile(m_ChrXProbes, chrXProbeFile);
    EngineUtil::readProbeFile(m_ChrYProbes, chrYProbeFile);

    initialize(femaleThresh, maleThresh);
}

CnProbeGenderCelListener::CnProbeGenderCelListener(const std::vector< std::vector<probeid_t> >& chrXProbes,
                                                   const std::vector< std::vector<probeid_t> >& chrYProbes,
                                                   double femaleThresh, 
                                                   double maleThresh, 
                                                   bool chrZW)
{
    m_ZWGenderCalling = chrZW;

    if (m_ZWGenderCalling) {
        m_DipSexChrLetter = "Z";
        m_HapSexChrLetter = "W";
    }
    else {
        m_DipSexChrLetter = "X";
        m_HapSexChrLetter = "Y";
    }
    m_ChrXProbes = chrXProbes;
    m_ChrYProbes = chrYProbes;
    initialize(femaleThresh, maleThresh);
}

void CnProbeGenderCelListener::initialize(double femaleThresh, double maleThresh)
{
    stringstream strm;
    for ( unsigned int i = 0 ; i < m_ChrXProbes . size(); i ++ )
        strm << "Read " << m_ChrXProbes[i].size() << " " << m_DipSexChrLetter << " probes on channel "<< i << ",";
    for ( unsigned int i = 0 ; i < m_ChrYProbes . size(); i ++ )
        strm << "Read " << m_ChrYProbes[i].size() << " " << m_HapSexChrLetter << " probes on channel "<< i << ",";
    Verbose::out(2, strm.str());

    if(m_ChrXProbes.empty())
        Verbose::out(1, "WARNING: No chr" + m_DipSexChrLetter + " probes loaded. Unable to run CnProbeGenderCelListener.");

    if(m_ChrYProbes.empty())
        Verbose::out(1, "WARNING: No chr" + m_HapSexChrLetter + " probes loaded. Unable to run CnProbeGenderCelListener.");

    m_GenderName="cn-probe-chr"+m_DipSexChrLetter+m_HapSexChrLetter+"-ratio";
    m_GenderDescription="copy number probe chr"+m_DipSexChrLetter+"/"+m_HapSexChrLetter+" ratio";

    m_FemaleThresh = femaleThresh;
    m_MaleThresh = maleThresh;

    declareMetric(m_GenderName + "_gender_mean"+m_DipSexChrLetter, ChipSummary::Metric::Double);
    declareMetric(m_GenderName + "_gender_mean"+m_HapSexChrLetter, ChipSummary::Metric::Double);
    declareMetric(m_GenderName + "_gender_ratio", ChipSummary::Metric::Double);
    declareMetric(m_GenderName + "_gender",       ChipSummary::Metric::String);
}



