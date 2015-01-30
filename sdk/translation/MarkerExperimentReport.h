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
 * @file   MarkerExperimentReport.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  DEPRECATED DMET2: Class for the marker report as report files.
 */

#ifndef TRANSLATION_MARKER_EXPERIMENT_REPORT_H
#define TRANSLATION_MARKER_EXPERIMENT_REPORT_H

#include "translation/ExperimentReport.h"
//
#include "file/TsvFile/TsvFile.h"
//


// Marker

const std::string MARKER_FILE_EXT = "marker.rpt";
const std::string MARKER_TSV_COLUMNS[] = {
  "Experiment", "Gene", "Sample", "Functional Change", "External ID",
  "ProbeSet ID", "Basecall", "Reference Base", "Variant Base",  "Call",
  "Haplotype Marker", "Change for Variant", "Variant cDNA Change",
  "Variant DNA Change", "dbSNP ID", "Validated", "Allele Defining Marker",
  "Relevant Alleles"
};

const size_t NUM_MARKER_TSV_COLUMNS = 18;

class MarkerRow
{
public:
  // DO NOT ALPHABETIZE, these are in report order.
  std::string m_experiment;
  std::string m_gene;
  std::string m_sample;
  std::string m_functionalChange;
  std::string m_externalId;
  std::string m_probeSet;
  std::string m_baseCall;
  std::string m_referenceBase;
  std::string m_variantBase;
  std::string m_call;
  std::string m_haplotypeMarker;
  std::string m_changeForVariant;
  std::string m_variantCDNAChange;
  std::string m_variantDNAChange;
  std::string m_dbSNPId;
  std::string m_validated;
  std::string m_alleleDefiningMarker;
  std::string m_relevantAlleles;
};

class MarkerRowComp
{
public:
  bool operator()(const MarkerRow & a,
                  const MarkerRow & b) const;

};


class MarkerExperimentReport : public ExperimentReport
{
public:


  virtual bool         close(bool abort = false);
  virtual bool         generate(class RunTimeEnvironment & rte,
                                class TranslationTableModel & ttm,
                                std::map<std::string, ExperimentResults*> &er);

  std::string          getMarkerReportName(const class RunTimeEnvironment & rte) ;
  virtual std::string  name() {
    return std::string("MarkerExperimentReport");
  }
  void                 open(const class RunTimeEnvironment & rte);
  virtual ExperimentReportTypeEnum report_type() {
    return MARKER;
  }

  MarkerExperimentReport() {
    m_open_called = false;
  }

private:
  std::string         m_guid;
  class affx::TsvFile m_marker_tsv;
  bool                m_open_called;
  std::string         m_reportName;

};


#endif /* TRANSLATION_MARKER_EXPERIMENT_REPORT_H */
