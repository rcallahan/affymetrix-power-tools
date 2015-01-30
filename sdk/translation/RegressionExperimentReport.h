////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
 * @file   RegressionExperimentReport.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  Class for generating the regression results as a report file.
 */

#ifndef TRANSLATION_REGRESSION_EXPERIMENT_REPORT_H
#define TRANSLATION_REGRESSION_EXPERIMENT_REPORT_H

#include "translation/ExperimentReport.h"
#include "translation/RunTimeEnvironment.h"
//
#include "file/TsvFile/TsvFile.h"
//

const std::string REGRESSION_MARKER_TSV_COLUMNS[] = { "Experiment", "Gene", "ProbeSet", "A1", "A2", "Ref", "Var", "Call" };
const size_t REGRESSION_NUM_MARKER_TSV_COLUMNS = 8;

const std::string REGRESSION_HAPLOTYPE_TSV_COLUMNS[] = { "Experiment", "Gene", "Call", "Call_Count", "Known_Count" };
const size_t REGRESSION_NUM_HAPLOTYPE_TSV_COLUMNS = 5;


class RegressionMarkerRow
{
public:
  std::string m_experiment;
  std::string m_gene;
  std::string m_probeSet;
  std::string m_a1;
  std::string m_a2;
  std::string m_ref;
  std::string m_var;
  std::string m_call;
};

class RegressionMarkerRowComp
{
public:
  bool operator()(const RegressionMarkerRow & a,
                  const RegressionMarkerRow & b) const;

};


class RegressionHaplotypeRow
{
public:
  std::string m_experiment;
  std::string m_gene;
  std::string m_call;
  std::string m_callCount;
  std::string m_knownCount;
};


class RegressionExperimentReport : public ExperimentReport
{
public:

  RegressionExperimentReport() {
    m_open_called = false;
  }

  virtual bool       close(bool abort = false);
  virtual bool       generate(class RunTimeEnvironment & rte,
                              class TranslationTableModel & ttm,
                         std::map<std::string, ExperimentResults*> & er);
  std::string        getHaplotypeReportName(const RunTimeEnvironment & rte);
  std::string        getMarkerReportName(const RunTimeEnvironment & rte);
  static std::string getDmet3MarkerFileNameExt() {
    return ToStr("dmet3_marker.reg");
  }
  static std::string getDmet3HaplotypeFileNameExt() {
    return ToStr("dmet3_haplotype.reg");
  }

  virtual std::string name() {
    return std::string("RegressionExperimentReport");
  }
  virtual ExperimentReportTypeEnum report_type() {
    return REGRESSION;
  }



private:
  affx::TsvFile m_haplotype_tsv;
  std::string   m_haplotypeReportName;
  affx::TsvFile m_marker_tsv;
  std::string   m_markerReportName;
  bool          m_open_called;

  void         _close_haplotype_file();
  void         _close_marker_file();
  bool         _generateHaplotypeCallResults(RunTimeEnvironment & rte,
                                       class ExperimentGeneResults & egr);
  bool         _generateMarkerCallResults(RunTimeEnvironment & rte,
                                    class ExperimentGeneResults & egr);
  std::string _getBaseReportName(const RunTimeEnvironment & rte) const;
  void        _open_haplotype_file(const RunTimeEnvironment & rte);
  void        _open_marker_file(const RunTimeEnvironment & rte);
  void        _report_haplotype(const RegressionHaplotypeRow & rhr);
  void        _report_marker(const RegressionMarkerRow & rmr);

};

#endif /* TRANSLATION_REGRESSION_EXPERIMENT_REPORT_H */
