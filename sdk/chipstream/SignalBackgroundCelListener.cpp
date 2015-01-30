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
#include "chipstream/SignalBackgroundCelListener.h"
//
#include "chipstream/GenoUtility.h"
#include "chipstream/apt-geno-qc/GenoQC.h"
//
#include "stats/stats.h"
//

using namespace std;
using namespace affx;

SignalBackgroundCelListener::SignalBackgroundCelListener(std::vector<ProbeListPacked>& probeSets,
                                                         QCProbesetOptions& psOpts,
                                                         ChipLayout* layout,
                                                         std::string label) :
  m_ProbeSets(probeSets),
  m_psOpts(psOpts),
  m_layout(layout),
  m_Label(label)
{
  Verbose::out(2, "Initializing Multichannel SignalBackgroundCelListener Okay!");

  //
  declareMetric(m_Label+"_AT_B_IQR",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_AT_B",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_AT_FLD",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_AT_SBR",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_AT_S_IQR",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_AT_S",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_A_signal_mean",ChipSummary::Metric::Double);
  //declareMetric(m_Label+"_AsignalToTsignal",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_C_signal_mean",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_GC_B_IQR",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_GC_B",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_GC_FLD",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_GC_SBR",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_GC_S_IQR",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_GC_S",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_G_signal_mean",ChipSummary::Metric::Double);
  //declareMetric(m_Label+"_GsignalToCsignal",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_T_signal_mean",ChipSummary::Metric::Double);
}

/**
 * Add a new probe subset via a mask associated with identifier name.
 *
 * @param name - string identifier associated with this mask (i.e. "all", "pm")
 * @param mask - bit mask with probes to calculate stats for set to true.
 */
void SignalBackgroundCelListener::addProbeMask(const std::string &name,
                                               const std::vector<bool> &mask) {
  m_MaskNames.push_back(name);
  m_ProbeMasks.push_back(mask);
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
SignalBackgroundCelListener::calcStats(const std::vector<std::vector<bool> > &masks,
                                       affymetrix_fusion_io::FusionCELData* cel,
                                       std::vector<CumulativeStats<double> > &chipStats)
{
  assert(chipStats.size() == masks.size());
  /*Turn on full method where medians and IQR could be computed */
  for (int sIx = 0; sIx < chipStats.size(); sIx++){
    chipStats[sIx].setFull();
  }

  int numCels = cel->GetNumCells();
  if (numCels == 0) {
    Err::errAbort("SignalBackgroundCelListener::calcStats -- No data found in one of the channels for CEL file: " + ToStr(cel->GetFileName()));
  }

  /* Get each entry in a cel file and then add it to different masks
     as appropriate. */

  for(int celIx = 0; celIx < numCels; celIx++) {
    std::vector<std::vector<bool> >::const_iterator maskIx;
    std::vector<CumulativeStats<double> >::iterator statIx;
    float intensity = cel->GetIntensity(celIx);
    /* Check all the masks to see if they include this data point. */
    for(maskIx = masks.begin(), statIx = chipStats.begin();
        maskIx != masks.end() && statIx != chipStats.end();
        maskIx++, statIx++) {
      if(*(maskIx->begin() + celIx) == true) {
        statIx->addData(intensity);
      }
    }
  }

  // calculate full method metrics like medians and IQR
  for (int sIx = 0; sIx < chipStats.size(); sIx++) {
    chipStats[sIx].setUpFullMethod();
  }

}

// /*
//  * Cel file of data.
//  *
//  * @param cel - Handle to an open cel phone.
//  */
// void SignalBackgroundCelListener::newChip(affymetrix_fusion_io::FusionCELData *cel) {
//   //std::vector<std::vector<bool> >::iterator maskIx;
//   std::vector<CumulativeStats<double> > chipStats(m_ProbeMasks.size());
//   calcStats(m_ProbeMasks, cel, chipStats);
//
//   std::vector<ChipSummary::Metric> metrics;
//   for(int i = 0; i<chipStats.size(); i++) {
//     metrics.push_back(ChipSummary::Metric(m_MaskNames[i] + "_mean",chipStats[i].getMean()));
//   }
//   m_SummaryStats.push_back(metrics);
//   setValid(true);
// }
//
// void SignalBackgroundCelListener::newMultiChannelChip(int celIx) {

void SignalBackgroundCelListener::newChip(affymetrix_fusion_io::FusionCELData* cel)
{
  double A_signal=0.0;
  double T_signal=0.0;
  double G_signal=0.0;
  double C_signal=0.0;

  const int channel_cnt=2;

  std::vector<double> signal; //signs[channel]
  std::vector<double> signal_median; //signal_median[channel]
  std::vector<double> signal_IQR; //signal_IQR[channel] inter-quartile range
  signal.resize(channel_cnt);
  signal_median.resize(channel_cnt);
  signal_IQR.resize(channel_cnt);

  std::vector<double> bg;   //bg[channel]
  std::vector<double> bg_median; //bg_median[channel]
  std::vector<double> bg_IQR; //bg_IQR[channel]
  bg.resize(channel_cnt);
  bg_median.resize(channel_cnt);
  bg_IQR.resize(channel_cnt);

  std::vector<double> signalbackground; //signalbackground[channel]
  std::vector<double> signalbackgroundFLD; //T score type of metric take two medians and avg of IQR
  signalbackground.resize(channel_cnt); //resize to fit two channels
  signalbackgroundFLD.resize(channel_cnt); //resize to fit two channels

  Verbose::out(4, "calculating signal-background statistics for "+ ToStr(cel->GetFileName()));

	std::vector<std::wstring> channel_names = cel->GetChannels();

  for (int channel = 0; channel< channel_cnt; channel++){
    // @todo how do we know it has two?
    cel->SetActiveDataGroup(channel_names[channel]);

    std::vector<std::vector<bool> >::iterator maskIx;
    std::vector<CumulativeStats<double> > chipStats(m_ProbeMasks.size());

    calcStats(m_ProbeMasks, cel, chipStats);

#define PRINTIT(_x) { printf("   %-30s %10.6f\n",#_x,_x); }

//    for (int j=0;j<m_ProbeMasks.size();j++) {
//      printf("JHG: channel: %2d j: %2d --------------------\n",channel,j);
//      PRINTIT(chipStats[j].getMean());
//      PRINTIT(chipStats[j].getMedian());
//      PRINTIT(chipStats[j].getIQR());
//    }

    // GC channel
    if (channel == 0) {
      for(int j=0; j<m_ProbeMasks.size(); j++){
        if (m_MaskNames[j].compare("GC") == 0){
          //chipStats[j].setUpFullMethod();
          //Verbose::out(4, "GC median is " + ToStr(chipStats[j].getMedian()));
          //Verbose::out(4, "GC IQR is " + ToStr(chipStats[j].getIQR()));
          signal[channel] = chipStats[j].getMean();
          signal_median[channel] = chipStats[j].getMedian();
          signal_IQR[channel] = chipStats[j].getIQR();
        }
        else if (m_MaskNames[j].compare("AT") == 0){
          bg[channel] = chipStats[j].getMean();
          bg_median[channel] = chipStats[j].getMedian();
          bg_IQR[channel] = chipStats[j].getIQR();
        }
        else if (m_MaskNames[j].compare("G") == 0){ //getting A signal in AT channel
          G_signal = chipStats[j].getMean();
        }
        else if (m_MaskNames[j].compare("C") == 0){ //getting T signal in AT channel
          C_signal = chipStats[j].getMean();
        }
        else{
          // A and T signal are the backgrounds
          // we can pull them off if we need it
        }
      }
    }
    // AT channel
    else if (channel == 1) {
      for(int j=0; j<m_ProbeMasks.size(); j++){
        if (m_MaskNames[j].compare("AT") == 0){
          signal[channel] = chipStats[j].getMean();
          signal_median[channel] = chipStats[j].getMedian();
          signal_IQR[channel] = chipStats[j].getIQR();
        }
        else if (m_MaskNames[j].compare("GC") == 0){
          bg[channel] = chipStats[j].getMean();
          bg_median[channel] = chipStats[j].getMedian();
          bg_IQR[channel] = chipStats[j].getIQR();
        }
        else if (m_MaskNames[j].compare("A") == 0){ //getting A signal in AT channel
          A_signal = chipStats[j].getMean();
        }
        else if (m_MaskNames[j].compare("T") == 0){ //getting T signal in AT channel
          T_signal = chipStats[j].getMean();
        }
        else{
          // G and C signals are the background in this channel
          // we can pull them off if need it
        }
      }
    }

//    PRINTIT(signal[channel]);
//    PRINTIT(signal_median[channel]);
//    PRINTIT(signal_IQR[channel]);
//    PRINTIT(bg[channel]);
//    PRINTIT(bg_median[channel]);
//    PRINTIT(bg_IQR[channel]);
//    PRINTIT(A_signal);
//    PRINTIT(C_signal);
//    PRINTIT(G_signal);
//    PRINTIT(T_signal);

    //Verbose::out(4, "signal for channel " + ToStr(channel) + ": "
    //	 + ToStr(signal[channel]));
    //Verbose::out(4, "background for channel " + ToStr(channel) + ": "
    //+ ToStr(bg[channel]));

    // what?
    signalbackground[channel] = (signal[channel] / bg[channel]);
    signalbackgroundFLD[channel] =
      ((signal_median[channel] - bg_median[channel]) * (signal_median[channel] - bg_median[channel])) /
      (0.5 * ((signal_IQR[channel] * signal_IQR[channel]) + (bg_IQR[channel] * bg_IQR[channel])));
  }

  //double A_T_signalratio = A_signal / T_signal;
  //double G_C_signalratio = G_signal / C_signal;

//  std::string AT_B_IQR   = m_Label+"_AT_B_IQR";
//  std::string AT_B_label = m_Label+"_AT_B";
//  std::string AT_FLD     = m_Label+"_AT_FLD";
//  std::string AT_SBR     = m_Label+"_AT_SBR";
//  std::string AT_S_IQR   = m_Label+"_AT_S_IQR";
//  std::string AT_S_label = m_Label+"_AT_S";
//  std::string A_sig_mean = m_Label+"_A_signal_mean";
//  //std::string Asig_Tsig  = m_Label+"_AsignalToTsignal";
//  std::string C_sig_mean = m_Label+"_C_signal_mean";
//  std::string GC_B_IQR   = m_Label+"_GC_B_IQR";
//  std::string GC_B_label = m_Label+"_GC_B";
//  std::string GC_FLD     = m_Label+"_GC_FLD";
//  std::string GC_SBR     = m_Label+"_GC_SBR";
//  std::string GC_S_IQR   = m_Label+"_GC_S_IQR";
//  std::string GC_S_label = m_Label+"_GC_S";
//  std::string G_sig_mean = m_Label+"_G_signal_mean";
//  //std::string Gsig_Csig  = m_Label+"_GsignalToCsignal";
//  std::string T_sig_mean = m_Label+"_T_signal_mean";
//
  // what are we looking at?
  int chipIx=getNextChipIdx();

  // Populate the Summary Metrics

  // these appear to be the the correct order.
  // but the output is swapped around
  //setMetric(chipIx , m_Label+"Asig_Tsig      , A_T_signalratio);
  //setMetric(chipIx , m_Label+"Gsig_Csig      , G_C_signalratio);
  setMetric(chipIx   , m_Label+"_AT_B"          , bg[1]);
  setMetric(chipIx   , m_Label+"_AT_B_IQR"      , bg_IQR[1]);
  setMetric(chipIx   , m_Label+"_AT_FLD"        , signalbackgroundFLD[1]);
  setMetric(chipIx   , m_Label+"_AT_S"          , signal[1]);
  setMetric(chipIx   , m_Label+"_AT_SBR"        , signalbackground[1]);
  setMetric(chipIx   , m_Label+"_AT_S_IQR"      , signal_IQR[1]);
  setMetric(chipIx   , m_Label+"_A_signal_mean" , A_signal);
  setMetric(chipIx   , m_Label+"_C_signal_mean" , C_signal);
  setMetric(chipIx   , m_Label+"_GC_B"          , bg[0]);
  setMetric(chipIx   , m_Label+"_GC_B_IQR"      , bg_IQR[0]);
  setMetric(chipIx   , m_Label+"_GC_FLD"        , signalbackgroundFLD[0]);
  setMetric(chipIx   , m_Label+"_GC_S"          , signal[0]);
  setMetric(chipIx   , m_Label+"_GC_SBR"        , signalbackground[0]);
  setMetric(chipIx   , m_Label+"_GC_S_IQR"      , signal_IQR[0]);
  setMetric(chipIx   , m_Label+"_G_signal_mean" , G_signal);
  setMetric(chipIx   , m_Label+"_T_signal_mean" , T_signal);

  ///@todo option to report other metrics?
  setValid(true);
}
