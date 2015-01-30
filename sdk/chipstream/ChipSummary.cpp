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

#include "chipstream/ChipSummary.h"
//
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <utility>
#include <vector>
//



ChipSummary::Metric::Metric() :
  m_Name(""),
  m_Type(Unset),
  m_Double(0.0),
  m_Integer(0),
  m_String("")
{
  // nothing
}

ChipSummary::Metric::Metric(std::string name, int value) : m_Double(0.0), m_Integer(0), m_String("")
{
  m_Name = name;
  m_Type = Integer;
  m_Integer = value;
}

ChipSummary::Metric::Metric(std::string name, double value) : m_Double(0.0), m_Integer(0), m_String("")
{
  m_Name = name;
  m_Type = Double;
  m_Double = value;
}

ChipSummary::Metric::Metric(std::string name, std::string value) : m_Double(0.0), m_Integer(0), m_String("")
{
  m_Name = name;
  m_Type = String;
  m_String = value;
}

//ChipSummary::Metric& ChipSummary::Metric::operator=(const ChipSummary::Metric& from)
//{
//  m_Name=from.m_Name;
//  m_Type=from.m_Type;
//  m_Double=from.m_Double;
//  m_Integer=from.m_Integer;
//  m_String=from.m_String;
//}

/**
 * Constructor
 */
ChipSummary::ChipSummary() {
    m_Valid = false;
    m_SummaryChipCount=0;
    m_nextChipIdx=0;
}

/**
 * DeConstructor
 */
ChipSummary::~ChipSummary()
{
  // nothing
} 

int ChipSummary::getSummaryStatsSize() {
  return m_SummaryStats.size();
}

//
int ChipSummary::setSummaryChipCount(int cnt) {
  m_SummaryChipCount=cnt;
  return m_SummaryChipCount;
};
int ChipSummary::getSummaryChipCount() {
  return m_SummaryChipCount;
};

int ChipSummary::setNextChipIdx(int idx) {
  m_nextChipIdx=idx;
  return m_nextChipIdx;
}
int ChipSummary::getNextChipIdx() {
  int idx=m_nextChipIdx;
  m_nextChipIdx++;
  return idx;
}

/**
 * Get vector of metrics
 */
ChipSummary::metricVec_t ChipSummary::getMetrics(int chip)
{
    if(!isValid())
        Err::errAbort("ChipSummary::getMetrics called when not valid");
    if(m_SummaryStats.size() <= chip)
        Err::errAbort("ChipSummary::getMetrics requested chip (" + 
                  ToStr(chip) + ") out of range.");
    /// @todo why dont we fscking sort them before whining?
    for(int i=0; i<m_SummaryStats[chip].size(); i++)
        checkStat(m_SummaryStats[chip][i], m_MetricDefs[i]);
    return m_SummaryStats[chip];
}

/**
 * Get named metric
 */
bool ChipSummary::getMetric(int chip, std::string name, ChipSummary::Metric &metric){
    if(!isValid())
        Err::errAbort("ChipSummary::getMetrics called when not valid");
    if(m_SummaryStats.size() <= chip)
        Err::errAbort("ChipSummary::getMetrics requested chip (" + 
                  ToStr(chip) + ") out of range.");
    for(int i=0; i<m_SummaryStats[chip].size(); i++) {
        ///@todo not dealing with non-unique names
        if(m_SummaryStats[chip][i].m_Name == name) {
            checkStat(m_SummaryStats[chip][i], m_MetricDefs[i]);
            metric = m_SummaryStats[chip][i];
            return true;
        }
    }
    return false;
}

/**
 * Get vector of metric names -- This will exactly match the number and order
 * of metrics returned by getMetrics().
 */
ChipSummary::metricDefVec_t ChipSummary::getMetricDefs(){
    return m_MetricDefs;
}

/** 
 * Is the summary information valid/set
 */
bool ChipSummary::isValid() { return m_Valid; }

/** 
 * Set Valid State
 */
bool ChipSummary::setValid(bool setTo) { m_Valid = setTo; return m_Valid; }

/**
 * Check that the metric available matches the predefined list
 */
void ChipSummary::checkStat(ChipSummary::Metric &metric,const ChipSummary::MetricDef& metricdef)
{
  if (metric.m_Name != metricdef.m_name) {
    Err::errAbort("Chip summary metric name mismatch. Found: '"+metric.m_Name+"' Expected: '"+metricdef.m_name+"'");
  }
  if (metric.m_Type != metricdef.m_type) {
    Err::errAbort("Chip summary metric type mismatch for '"+metric.m_Name+"'.");
  }
}

/// @todo change all the fscking "m_MetricDefs.push_back(...)" to use this.
/// And since we are declaring them, why not allocate the metric and fill it in later by name?!?
void
ChipSummary::declareMetric(const ChipSummary::MetricDef& metricdef)
{
  m_MetricDefs.push_back(metricdef);
}

void
ChipSummary::declareMetric(const std::string& metric_name,Metric::Type metric_type)
{
  m_MetricDefs.push_back(ChipSummary::MetricDef(metric_name,metric_type));
}

void
ChipSummary::declareMetrics(const ChipSummary::metricDefVec_t& metricDefVec)
{
  m_MetricDefs.insert(m_MetricDefs.end(),metricDefVec.begin(),metricDefVec.end());
}

void
ChipSummary::declareMetric(const std::string& metric_name,Metric::Type metric_type,int precision)
{
  m_MetricDefs.push_back(ChipSummary::MetricDef(metric_name,metric_type,precision));
}

int 
ChipSummary::metricNameToIndex(const std::string& metric_name)
{
  for (int i=0;i<m_MetricDefs.size();i++) {
    if (m_MetricDefs[i].m_name==metric_name) {
      return i;
    }
  }
  //
  APT_ERR_ABORT("Metric name "+metric_name+" not found.");
  return -1;
}

ChipSummary::Metric* 
ChipSummary::getMetricPtr(int chipIdx,const std::string& metric_name)
{
  if (m_SummaryStats.size()<=chipIdx) {
    m_SummaryStats.resize(chipIdx+1);
  }

  if (m_SummaryStats[chipIdx].size()!=m_MetricDefs.size()) {
    m_SummaryStats[chipIdx].resize(m_MetricDefs.size());
  }

  int idx=metricNameToIndex(metric_name);  
  APT_ERR_ASSERT(idx<m_SummaryStats[chipIdx].size(),"internal error.");

  return &m_SummaryStats[chipIdx][idx];
}

void
ChipSummary::setMetric(int chipIdx,const ChipSummary::Metric& val)
{
  ChipSummary::Metric* mptr=getMetricPtr(chipIdx,val.m_Name);
  *mptr=val;
}

void
ChipSummary::setMetric(int chipIdx,const std::string& metric_name,int val)
{
  //printf("setMetric(%d,'%s',%d)\n",chipIdx,metric_name.c_str(),val);

  ChipSummary::Metric* mptr=getMetricPtr(chipIdx,metric_name);
  //
  mptr->m_Name=metric_name;
  mptr->m_Type=Metric::Integer;
  mptr->m_Integer=val;
}

void
ChipSummary::setMetric(int chipIdx,const std::string& metric_name,double val)
{
  //printf("setMetric(%d,'%s',%f)\n",chipIdx,metric_name.c_str(),val);

  ChipSummary::Metric* mptr=getMetricPtr(chipIdx,metric_name);
  //
  mptr->m_Name=metric_name;
  mptr->m_Type=Metric::Double;
  mptr->m_Double=val;
}

void
ChipSummary::setMetric(int chipIdx,const std::string& metric_name,const std::string& val)
{
  //printf("setMetric(%d,'%s','%s')\n",chipIdx,metric_name.c_str(),val.c_str());

  ChipSummary::Metric* mptr=getMetricPtr(chipIdx,metric_name);
  //
  mptr->m_Name=metric_name;
  mptr->m_Type=Metric::String;
  mptr->m_String=val;
}

void
ChipSummary::setMetrics(int chipIdx,const ChipSummary::metricVec_t& metrics)
{
  for (int i=0;i<metrics.size();i++) {
    setMetric(chipIdx,metrics[i]);
  }
}

void ChipSummary::addHeaderComment(const std::string& comment)
{
  m_headercomments.push_back(comment);
}

std::vector<std::string> ChipSummary::getHeaderComments()
{
  return m_headercomments;
}

