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
 * @file   HaplotypeExperimentReport.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  DEPRECATED DMET2 Class for generating the haplotype report as a report file. 
 */

#ifndef TRANSLATION_HAPLOTYPEEXPERIMENTREPORT_H
#define TRANSLATION_HAPLOTYPEEXPERIMENTREPORT_H

#include "translation/ExperimentGeneResults.h"
#include "translation/ExperimentReport.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h" // includes Verbose.h
//

using namespace std;

class HaplotypeExperimentReport : public ExperimentReport
{
public:
  std::string m_guid;
  affx::TsvFile m_tsv;

  void open(const std::string& fileName);
  void report(const ExperimentGeneResults& result);
  void close();

  /*  
  //string reportHeader;

  // ROW DATA
  // A std::vector index represents a table row.
  // The collection of std::vectors represents the set of experiment results
  // All data is represented as std::strings. Any processing of numbers
  // will translate too and from std::strings in the calculation methods.


  std::vector<std::string> m_experiment;
  std::vector<std::string> m_geneId;
  std::vector<std::string> m_sample;
  std::vector<std::string> m_call;
  std::vector<std::string> m_callCount;
  std::vector<std::string> m_knownCount;
  std::vector<std::string> m_UNKExists;
  std::vector<std::string> m_basecallRate;
  std::vector<std::string> m_basecallCount;
  std::vector<std::string> m_noCallCount;
  std::vector<std::string> m_possibleRareAlleleCount;
  std::vector<std::string> m_notAvailableCount;

  
  bool generate();

  bool input_translated_experiment_results(  const ExperimentGeneResults & );

  ExperimentReportTypeEnum report_type() {
    return HAPLOTYPE;
  }
  */
};

#endif /* TRANSLATION_HAPLOTYPEEXPERIMENTREPORT_H */ 
