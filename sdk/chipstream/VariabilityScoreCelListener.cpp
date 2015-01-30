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
#include "chipstream/VariabilityScoreCelListener.h"
//
#include "chipstream/GenoUtility.h"
#include "chipstream/apt-geno-qc/GenoQC.h"
//
#include "stats/stats.h"
//

using namespace std;
using namespace affx;

VariabilityScoreCelListener::VariabilityScoreCelListener(std::vector<ProbeListPacked>& probeSets, 
                                                         QCProbesetOptions& psOpts,
                                                         ChipLayout* layout, 
                                                         const std::string& label) :
  m_ProbeSets(probeSets),
  m_psOpts(psOpts),
  m_layout(layout),
  m_Label(label)
{
  Verbose::out(2, "Initializing Multichannel VariabilityScoreCelListener Okay!");

  declareMetric(m_Label+"_CV_GC",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_CV_AT",ChipSummary::Metric::Double);
}

/** 
 * Loop though an open cel file and fill in the statistics for each
 * probe subset.
 * 
 * @param masks - Different probe subsets to calculate stats for.
 * @param cel - Cel file to get data from.
 * @param chipStats - Vector to fill in with stats for each mask.
 */
void 
VariabilityScoreCelListener::calcStats(vector<ProbeListPacked>& probeSets, 
                                       affymetrix_fusion_io::FusionCELData* cel,
                                       CumulativeStats<double> &chipStats, 
                                       int channel)
{
  Verbose::out(4, "Total probesets to compute D_Score: " + ToStr(probeSets.size()));
  for (int i=0; i<probeSets.size(); i++) {
    //float Min, Max; //## turn off V-score
    float mean, std;
    ProbeSet *ps = ProbeListFactory::asProbeSet(probeSets[i]);
    //probesetMinMax(cel, ps, Min, Max); //## turn off V-score, min max cal;

    probesetMeanStd(cel, ps, mean, std);
    //For debug only
    //Verbose::out(4, "probesetname: " + ToStr(ps->name) + "; mean: " + ToStr(mean)
    //	 + ";std: " + ToStr(std));
    double cv_score = std/mean;
    // chipStats1.addData(d_score); //## turn off V-score
    chipStats.addData(cv_score);    


    // double d_score = (Max-Min)/(Max+Min); //##turn ooff V-score
    //Verbose::out(4, ToStr(ps->name)+":"+ToStr(channel)+","+ToStr(Min)+","+ToStr(Max)+","+ToStr(d_score));
    chipStats.addData(cv_score);
    Freez(ps);
  }
}

/** 
 * Loop though an open cel file and fill in the statistics for each
 * probe subset.
 * 
 * @param masks - Different probe subsets to calculate stats for.
 * @param cel - Cel file to get data from.
 * @param chipStats - Vector to fill in with stats for each mask.
 */
void VariabilityScoreCelListener::calcStats(vector<ProbeListPacked>& probeSets, 
                                            affymetrix_fusion_io::FusionCELData* cel,
                                            CumulativeStats<double>& chipStats1, 
                                            CumulativeStats<double>& chipStats2, 
                                            int channel)
{
  Verbose::out(4, "Total probesets to compute D_Score: " + ToStr(probeSets.size()));
  for (int i=0; i<probeSets.size(); i++){
    float Min, Max, mean, std;
    ProbeSet *ps = ProbeListFactory::asProbeSet(probeSets[i]);
    probesetMinMax(cel, ps, Min, Max);
    probesetMeanStd(cel, ps, mean, std);
    //For debug only
    //Verbose::out(4, "probesetname: " + ToStr(ps->name) + "; mean: " + ToStr(mean)
    //	 + ";std: " + ToStr(std));
    double d_score = (Max-Min)/(Max+Min);
    double cv_score = std/mean;
    chipStats1.addData(d_score);
    chipStats2.addData(cv_score);
    Freez(ps);
  }
}

/*
 * Cel file of data.
 *
 * @param cel - Handle to an open cel phone.
 */

void VariabilityScoreCelListener::newChip(affymetrix_fusion_io::FusionCELData *cel) {
  m_CelNames.push_back(cel->GetFileName());

  const int channel_cnt=2;
	std::vector<std::wstring> channel_names = cel->GetChannels();

  std::vector<double> medians; // medians[channel] V_score
  std::vector<double> medians2; // medians2[channel] CV_score
  medians.resize(channel_cnt);
  medians2.resize(channel_cnt);

  // std::vector<CumulativeStats<double> > chipStats;
  // chipStats.resize(channel_cnt);
  std::vector<CumulativeStats<double> > chipStats2;
  chipStats2.resize(channel_cnt);

  Verbose::out(4, "calculating medians for "+ ToStr(cel->GetFileName()));
 
  for (int channel = 0; channel<channel_cnt; channel++) {
    cel->SetActiveDataGroup(channel_names[channel]);

//    chipStats[channel].setFull();
//    calcStats(m_ProbeSets, cel, chipStats[channel], channel);
//    chipStats[channel].setUpFullMethod();
//    medians[channel] = chipStats[channel].getMedian();

    //chipStats[channel].setFull(); //turning off V-score method
    chipStats2[channel].setFull(); //turning on FULL method

    // old, turn off V-score
    // calcStats(m_ProbeSets, &CEL, chipStats[channel], channel);
    // dual cal, V and CV score
    // calcStats(m_ProbeSets, &CEL, chipStats[channel], chipStats2[channel],channel); 

    calcStats(m_ProbeSets, cel, chipStats2[channel],channel);  //CV_score only
    //chipStats[channel].setUpFullMethod();  //## turn off V_score
    chipStats2[channel].setUpFullMethod();   
    //medians[channel] = chipStats[channel].getMedian(); //## turn off V_score
    medians2[channel] = chipStats2[channel].getMedian();
  }

  int chipIx=getNextChipIdx();

  // setMetric(chipIx,m_Label+"_GC_V",medians[1]);
  // setMetric(chipIx,m_Label+"_AT_V",medians[0]);
  setMetric(chipIx,m_Label+"_CV_AT", medians2[1]);
  setMetric(chipIx,m_Label+"_CV_GC", medians2[0]);

  ///@todo option to report other metrics?
  setValid(true);
}
