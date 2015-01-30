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

#ifndef CHIPSUMMARY_H
#define CHIPSUMMARY_H

#include <cstring>
#include <string>
#include <utility>
#include <vector>
//

///@todo AW: This is a bit of a mess because we had to add in after the fact
///          the need to know up front what metrics will be available. Originally
///          the list of available metrics (eg getMetricNames) was only available
///          after the analysis was run.

/** 
 * Class for objects that provide chip summary metrics. This is basically
 * an interface with some helper methods for exposing per-chip metrics
 * to various reporters.
 */
class ChipSummary {
  public:

    /**
    * Class for storing metrics
    */
    class Metric {
        public:
      /// @todo Could we have a more generic name?
      enum Type {Unset, Integer, Double, String};

      Metric();
      Metric(std::string name, int value);
      Metric(std::string name, double value);
      Metric(std::string name, std::string value);
      //Metric& operator=(const Metric& metric);

      std::string m_Name;
      Type m_Type;
      double      m_Double;
      int         m_Integer;
      std::string m_String;
      // @todo add docstring to header of report file.
      // std::string m_docstring
    };

  class MetricDef {
  public:
    std::string  m_name;
    Metric::Type m_type;
    int          m_precision;

    MetricDef(const std::string& name,Metric::Type type) {
      m_name=name;
      m_type=type;
      m_precision=-1;
    };
    MetricDef(const std::string& name,Metric::Type type,int precision) {
      m_name=name;
      m_type=type;
      m_precision=precision;
    };
  };

  typedef std::vector<Metric> metricVec_t;
  typedef std::vector<MetricDef> metricDefVec_t;

    /**
    * Constructor
    */
    ChipSummary();

    /**
    * DeConstructor
    */
    virtual ~ChipSummary();

    /**
    * Get vector of metrics
    */
  metricVec_t getMetrics(int chip);

    /**
    * Get named metric
    */
    bool getMetric(int chip, std::string name, Metric &metric);

    /**
    * Get vector of metric names -- This will exactly match the number and order
    * of metrics returned by getMetrics().
    */
    metricDefVec_t getMetricDefs();

    /** 
    * Is the summary information valid/set
    */
    bool isValid();
    
  protected:
    /** 
    * Set Valid State
    */
    bool setValid(bool setTo);
    
  private:
    /**
    * Check that the metric available matches the predefined list
    */
  void checkStat(Metric &metric, const ChipSummary::MetricDef& metricdef);

public:  
  void declareMetric(const ChipSummary::MetricDef& metricdef);
  void declareMetric(const std::string& metric_name,Metric::Type metric_type);
  void declareMetric(const std::string& metric_name,Metric::Type metric_type,int precision);

  void declareMetrics(const ChipSummary::metricDefVec_t& metricDefVec);

  int metricNameToIndex(const std::string& metric_name);
  ChipSummary::Metric* getMetricPtr(int chipIdx,const std::string& metric_name);

  // Why enforce an order? Have the program reorder them.
  void setMetric(int chipIdx,const ChipSummary::Metric& metric);
  void setMetric(int chipIdx,const std::string& metric_name,int val);
  void setMetric(int chipIdx,const std::string& metric_name,double val);
  void setMetric(int chipIdx,const std::string& metric_name,const std::string& val);
  //
  void setMetrics(int chipIdx,const ChipSummary::metricVec_t& metrics);

  //
  int getSummaryStatsSize();
  int setSummaryChipCount(int cnt);
  int getSummaryChipCount();
  //
  int setNextChipIdx(int idx);
  int getNextChipIdx();

  // a comment to 
  void addHeaderComment(const std::string& comment);
  std::vector<std::string> getHeaderComments();

private:
  /// Vector of the metrics to be reported
  /// use "declareMetric".
  ChipSummary::metricDefVec_t m_MetricDefs;

  std::vector<std::string> m_headercomments;

protected:
  /// Stats for ChipSummary
  std::vector<ChipSummary::metricVec_t> m_SummaryStats;
  int m_SummaryChipCount;
  
  int m_nextChipIdx;

private:
    /// Are the summary states ready for consumption
    bool m_Valid;

};

#endif /* CHIPSUMMARY_H */
