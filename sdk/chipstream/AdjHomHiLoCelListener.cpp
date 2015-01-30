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
#include "chipstream/AdjHomHiLoCelListener.h"
//
#include "chipstream/GenoUtility.h"
#include "chipstream/QuantMethodReport.h"
#include "chipstream/apt-geno-qc/GenoQC.h"
//
#include "stats/stats.h"

using namespace std;
using namespace affx;

/** 
 * Process another cel files worth of data.
 */
void AdjHomHiLoCelListener::newChip(affymetrix_fusion_io::FusionCELData *cel) {
  std::vector<ChipSummary::Metric> metrics;
  m_SummaryStats.push_back(metrics);
  m_CelNames.push_back(cel->GetFileName());

  // Pass along CEL file to underlying cel listeners
  m_RandHiLo->newChip(cel);
  m_Set1HiLo->newChip(cel);
  m_Set2HiLo->newChip(cel);

  // Get Metrics out
  vector<ChipSummary::Metric> randMetrics = m_RandHiLo->getMetrics(m_CelNames.size()-1);
  vector<ChipSummary::Metric> set1Metrics = m_Set1HiLo->getMetrics(m_CelNames.size()-1);
  vector<ChipSummary::Metric> set2Metrics = m_Set2HiLo->getMetrics(m_CelNames.size()-1);

  // Compute Adjusted Hi Lo
  double adjhilo = 0;
  ///@todo all for "null" metric
  Metric rand("foo",2.2), set1("foo",2.2), set2("foo",2.2);

  if(!m_RandHiLo->getMetric(m_CelNames.size()-1, m_LabelRand, rand))
      Err::errAbort("Unable to find HomHiLo for " + m_LabelRand);
  if(!m_Set1HiLo->getMetric(m_CelNames.size()-1, m_LabelSet1, set1))
      Err::errAbort("Unable to find HomHiLo for " + m_LabelSet1);
  if(!m_Set2HiLo->getMetric(m_CelNames.size()-1, m_LabelSet2, set2))
      Err::errAbort("Unable to find HomHiLo for " + m_LabelSet2);

  if(fabs(set1.m_Double-set2.m_Double) > m_HiLoDiffCutoff)
      adjhilo = m_HiLoFailedValue;
  else
      adjhilo = rand.m_Double;

  // Populate the Summary Metrics
  m_SummaryStats[m_CelNames.size()-1].push_back(ChipSummary::Metric(m_Label, adjhilo));

  // Pass along underlying metrics
  for(int i = 0; i< randMetrics.size(); i++)
    m_SummaryStats[m_CelNames.size()-1].push_back(randMetrics[i]);
  for(int i = 0; i< set1Metrics.size(); i++)
    m_SummaryStats[m_CelNames.size()-1].push_back(set1Metrics[i]);
  for(int i = 0; i< set2Metrics.size(); i++)
    m_SummaryStats[m_CelNames.size()-1].push_back(set2Metrics[i]);

  setValid(true);
}
