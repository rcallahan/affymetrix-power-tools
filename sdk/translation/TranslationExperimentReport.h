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
 * @file   TranslationExperimentReport.h
 * @author Mybrid Spalding
 * @date   Mon Jul 14 10:02:58 PDT 2008
 * @brief  Single view object class for the three DMET3 reports (comprehensive, summary and uncalled) given their similarity. 
 */

#ifndef TRANSLATION_TRANSLATION_EXPERIMENT_REPORT_H
#define TRANSLATION_TRANSLATION_EXPERIMENT_REPORT_H

#include "translation/ExperimentGeneResults.h"
#include "translation/ExperimentReport.h"
#include "translation/RunTimeEnvironment.h"
#include "translation/TranslationTableModel.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h" // includes Verbose.h
#include "util/Guid.h"
//
#include <map>
//


extern const std::string TRANSLATION_REPORT_FILE_EXT;
extern const std::string SUMMARY_REPORT_FILE_EXT;
extern const std::string TRANSLATION_TSV_COLUMNS[];


class TranslationInterpretationCodeMap
{
public:
  enum TICM_ENUM {
    NoHAP, UNIQ, NC_PRA_NA, UNIQ_UNK, MULT, MULT_UNK, UNDH,
  };
  
  std::string  m_reportString;
};

class GeneHaplotypeCallMapElement
{
public:
  std::string m_gene;
  int         m_copyNumber;
  std::vector< std::string > m_knownCall;
  std::vector< std::string > m_unknownCall;
  GeneHaplotypeCallMapElement() {
    m_copyNumber = 2;
  }
};

class TranslationRow
{
public:

  // CONDITION AND STATE VARIABLES
  unsigned int m_copyNumber;
  bool         m_geneCopyNumberIndicator;
  bool         m_isSummaryRow;
  bool         m_isGeneSummary;
  bool         m_isUncalledRow;

  // DATA IN COLUMN ORDER, please do not alphabetize
  std::string m_index;
  std::string m_chpFileName;
  std::string m_gene;
  std::string m_knownCall;
  std::string m_unknownCall;
  std::string m_interpretationCode;
  std::string m_summaryFlag;
  std::string m_relevantAlleles;
  std::string m_markerName; //external id
  std::string m_probeSet;
  std::string m_baseCall;
  std::string m_referenceBase;
  std::string m_referenceReportBase;
  std::string m_variantBase;
  std::string m_variantReportBase;
  std::string m_variantUncalledBase;
  std::string m_call;
  std::string m_haplotypeMarker;
  std::string m_changeForVariant;
  std::string m_CDNAChange; //variantCDNAChange
  std::string m_genomicPosition;
  std::string m_dbSNPId;
  std::string m_validated;
  std::string m_originalBasecall;
  std::string m_overrideComment;

  TranslationRow() {
    m_isSummaryRow            = false;
    m_isGeneSummary           = false;
    m_isUncalledRow           = false;
    m_copyNumber              = 2;
    m_geneCopyNumberIndicator = false;
    m_call                    = "";
  }
  void dump();

};

class SummaryRowIndex {
public:
  int         m_copyNumber;
  std::string m_gene;
  bool        m_isGeneSummary;
  bool        m_isSummaryRow;
  std::string m_knownCall;
  std::string m_probeSet;
  std::string m_markerName;
  std::string m_summaryFlag;
  bool operator<( const SummaryRowIndex &b) const;
  
};


class TranslationExperimentReport : public ExperimentReport
{
public:

  TranslationExperimentReport(class SampleInfoTableModel *sitm = NULL,
                              class GenotypeOverrideTableModel *gotm = NULL);

  virtual bool  close(bool abort = false);
  virtual bool  generate(class RunTimeEnvironment & rte,
                         TranslationTableModel & ttm,
                         std::map<std::string, ExperimentResults*> &er);

  std::string  getComprehensiveReportName(const RunTimeEnvironment & rte);
  std::string  getSummaryReportName(const RunTimeEnvironment & rte);
  std::string  getUncalledReportName(const RunTimeEnvironment & rte);

  virtual std::string  name() { return std::string("TranslationExperimentReport");  }
  virtual ExperimentReportTypeEnum report_type() { return TRANSLATION; }

  std::string getSummaryFileName();
  std::string getComprehensiveFileName();

private:

  // Model pointers cached for convienance
  class GenotypeOverrideTableModel  *m_gotm;
  class SampleInfoTableModel        *m_sitm;

  affx::TsvFile m_comprehensiveTsv;
  std::string   m_comprehensiveReportName;
  int           m_experimentCount;
  int           m_experimentFileCount;
  std::string   m_guid;
  bool          m_openCalled;
  int           m_numFixedColumns;
  std::string   m_summaryReportName;
  std::map< SummaryRowIndex, TranslationRow >
                m_summaryRows;
  affx::TsvFile m_summaryTsv;
  int           m_totalColumnsWithSampleInfo;
  std::string   m_uncalledReportName;
  affx::TsvFile m_uncalledTsv;
  

  // HELPER functions, one time use 
  bool    _closeComprehensiveFile(bool abort);
  bool    _closeSummaryFile(bool abort);
  bool    _closeUncalledFile(bool abort);

  void    _generateBaseCall(ExperimentGeneResults *egr,
                            int ttRow,
                            int markerCallResultsIndex,
                            TranslationRow * cRow,
                            TranslationTableModel & ttm,
                            bool isMultiAllelicRow,
                            bool isProbeSetPresent);



  void   _generateGeneHaplotypeCallMap(RunTimeEnvironment & rte,
                                       std::map<std::string, ExperimentResults *>&er,
                                       std::map< std::string, GeneHaplotypeCallMapElement > & geneHaplotypeMap);
  void   _generateInterpretationCode(const ExperimentGeneResults & egr,
                                     TranslationTableModel & ttm,
                                     TranslationRow * cRow,
                                     int ttmRow,
                                     int numKnownCalls,
                                     int numUnknownCalls,
                                     std::map< std::string, std::string> & geneIntrepretationCode,
                                     int geneCopyNumber);

  void   _generateRelevantAlleles(TranslationTableModel & ttm,
                                  TranslationRow * cRow,
                                  int ttmRow,
                                  bool isMultAllelicRow);

  void   _generateReferenceBase(int ttRow,
                                TranslationRow * cRow,
                                TranslationTableModel & ttm);


  void   _generateVariantBase(int ttRow,
                              TranslationRow * cRow,
                              TranslationTableModel & ttm,
                              bool isMultiAllelicRow);

  int    _getNumColumns(bool includeSampleColumns = false);

  std::string _getFirstCHPFileRoot(const RunTimeEnvironment & rte);



  void   _open(const RunTimeEnvironment & rte);
  void   _openComprehensiveFile(const RunTimeEnvironment & rte,
                                std::vector<std::string > &sampleInfoHeaders);
  void   _openSummaryFile(const RunTimeEnvironment & rte,
                          std::vector<std::string > &sampleInfoHeaders);
  void   _openUncalledFile(const RunTimeEnvironment & rte);

  void   _report(const RunTimeEnvironment & rte,
                 const TranslationRow & cRow,
                 std::vector< std::string > & sampleInfo);
  void   _reportSummary(const TranslationRow & cRow,
                        std::vector< std::string > & sampleInfo);
  void   _reportSummarySorted(std::vector< std::string > & sampleInfo);
  
  void   _setIsUncalledRow(const RunTimeEnvironment & rte, TranslationRow & cRow);

};


#endif /* TRANSLATION_TRANSLATION_EXPERIMENT_REPORT_H */
