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
 * @file   TranslationExperimentReport.cpp
 * @author Mybrid Spalding
 * @date   Mon Jul 14 10:02:58 PDT 2008
 * @brief  Single view object class for the three DMET3 reports 
 * (comprehensive, summary and uncalled) given their similarity. 
 */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//
#include "translation/TranslationExperimentReport.h"
//
#include "translation/ExperimentResults.h"
#include "translation/GenotypeOverrideTableModel.h"
#include "translation/SampleInfoTableModel.h"
#include "translation/TranslationTableModel.h"
//
#include "util/Fs.h"
#include "util/Util.h"
//
#include "pcrecpp.h"
//
#include <stdio.h>

using namespace std;
using namespace affx;

std::map< std::string, std::string > geneInterpretationCode;
pcrecpp::RE RE_VAR("Var");
pcrecpp::RE RE_WILDCARD("NA|NoCall|NC|PRA|PossibleRareAllele|NotAvailable");

const std::string COMPREHENSIVE_REPORT_FILE_EXT = "_comprehensive";
const std::string SUMMARY_REPORT_FILE_EXT       = "_summary";
const std::string UNCALLED_REPORT_FILE_EXT      = "_uncalled";
const std::string COMPLETE_REPORT_FILE_EXT      = ".rpt";
const std::string TEMP_REPORT_FILE_EXT          = ".tmp";

const std::string COMPREHENSIVE_TSV_COLUMNS[]   = { "Index", "CHP Filename", "Gene", "Known Call", "Unknown Call",
    "Interpretation Code", "Summary Flag", "Relevant Alleles",
    "Common Name", "Probe Set ID", "Basecall", "Reference Base",
    "Variant Base",  "Call", "Haplotype Marker", "Change for Variant",
    "cDNA Change", "Genome Position", "dbSNP RS ID", "Original Basecall",
     "Override Comment" };


const std::string UNCALLED_TSV_COLUMNS[] = { "CHP Filename", "Gene", "Common Name", "Probe Set ID", "Basecall",
                                        "Override Comment", "Reference Allele", "Variant Allele"
                                      };

const TranslationInterpretationCodeMap INTERPRETATION_CODES[] = { {std::string("NoHAP")}, {std::string("UNIQ")}, {std::string("NC/PRA/NA")}, {std::string("UNIQ+UNK")}, {std::string("MULT")}, {std::string("MULT+UNK")}, {std::string("UNDH") } };

/*****************************************************************************/
/**
 * TranslationRow::dump
 * Synopsis:
 *  Debug routine only.
 */
/*****************************************************************************/
void TranslationRow::dump() {
  
      cerr << "m_index:               " << m_index << endl;
      cerr << "m_chpFileName:         " << m_chpFileName << endl; 
      cerr << "m_gene:                " << m_gene << endl;
      cerr << "m_knownCall:           " << m_knownCall << endl;
      cerr << "m_unknownCall:         " << m_unknownCall << endl;
      cerr << "m_interpretationCode:  " << m_interpretationCode << endl;
      cerr << "m_summaryFlag :        " << m_summaryFlag << endl;
      cerr << "m_relevantAlleles:     " << m_relevantAlleles << endl;
      cerr << "m_markerName:          " << m_markerName << endl;
      cerr << "m_probeSet:            " << m_probeSet << endl;
      cerr << "m_baseCall:            " << m_baseCall << endl;
      cerr << "m_referenceBase:       " << m_referenceBase << endl;
      cerr << "m_referenceReportBase: " << m_referenceReportBase << endl;
      cerr << "m_variantBase:         " << m_variantBase << endl;
      cerr << "m_variantReportBase:   " << m_variantReportBase << endl;
      cerr << "m_variantUncalledBase: " << m_variantUncalledBase << endl;
      cerr << "m_call :               " << m_call << endl;
      cerr << "m_haplotypeMarker :    " << m_haplotypeMarker << endl;
      cerr << "m_changeForVariant :   " << m_changeForVariant << endl;
      cerr << "m_CDNAChange :         " << m_CDNAChange << endl;
      cerr << "m_genomicPosition :    " << m_genomicPosition << endl;
      cerr << "m_dbSNPId :            " << m_dbSNPId << endl;
      cerr << "m_validated :          " << m_validated << endl;
      cerr << "m_originalBasecall :   " << m_originalBasecall << endl;
      cerr << "m_overrideComment :    " << m_overrideComment << endl;

}
// end TranslationRow::dump
/*****************************************************************************/
/*****************************************************************************/
/**
 * SummaryRowIndexSort::operator():
 * Synopsis:
 *  An STL class operator used to insert std::map records for the summary report in
 *  sorted order. 
 *
 *  Sort by business rules.
 *  1.) Copy number designation goes at the very top under all circumstances.
 *  2.) Gene haplotype or non-haplotype markers with only references
 *      are at the bottom.
 *  3.) Gene haplotype or non-haplotype markers with variants are reported
 *      second.
 *
 *  Rows within the same group as determined by the rules above are to be
 *  sorted in gene order,  and if where genes are also equivalent
 *  then further sort by probe set id. 
 * 
 * @param b - The object being compared with. 
 * 
 * @return - true if this object is less than b. 
 *
 */
/*****************************************************************************/
bool SummaryRowIndex::operator<( const SummaryRowIndex &b) const
{


  // COPY NUMBER FIRST
  if ( (m_copyNumber < 2)  &&  ( b.m_copyNumber >= 2 ) ) {
     return true;
  }
  else if ((m_copyNumber >= 2) && ( b.m_copyNumber < 2) ) {
     return false;
  }
  else if (  m_isGeneSummary || b.m_isGeneSummary ) {
    // REFERENCE LAST
    
    if (!( m_isGeneSummary && b.m_isGeneSummary)) {

      if (m_isGeneSummary) {
        return false;
      }
      return true;
    }
    
  } else if ( (m_knownCall.empty() || b.m_knownCall.empty()) &&
              !(m_knownCall.empty() && b.m_knownCall.empty())) {

    // VARIANT CALLS, the filling
    if ( m_knownCall.empty() ) {
      return false;
    }
    return true;
    
  }

  // Objects are of the same copy number, reference or variant group
  // so sort by gene/probe set id.

  if ( m_gene == b.m_gene ) {
    return(  m_probeSet < b.m_probeSet );
  }

  return( m_gene < b.m_gene);

  
}
// end SummaryRowIndex::operator<
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::TranslationExperimentReport:
 * Synopsis:
 *
 *  Default constructor.
 *
 * 
 * @param sitm  - sample info single instance data model cached for conviencance
 * @param gotm  - genotype override single instance data model cached for convienance
 *
 */
/*****************************************************************************/
TranslationExperimentReport::TranslationExperimentReport(SampleInfoTableModel *sitm, GenotypeOverrideTableModel *gotm) : m_gotm(gotm), m_sitm(sitm)
{

  m_openCalled = false;
  m_experimentFileCount = 0;
  m_experimentCount = 0;

  m_numFixedColumns = (sizeof(COMPREHENSIVE_TSV_COLUMNS) / sizeof(COMPREHENSIVE_TSV_COLUMNS[0])) ;

  if ((sitm != NULL)) {
    m_totalColumnsWithSampleInfo = _getNumColumns() + sitm->m_rows.size();
  } else {
    m_totalColumnsWithSampleInfo = 0;
  }

}
// end TranslationExperimentReport::TranslationExperimentReport
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::close:
 * Synopsis:
 *
 *   Close the report files and set the open state.
 *
 * 
 * @param abort - if the close is happening as a resulf of an abort exception then do not message via Verbose due to a bug in Windows. 
 *
 */
/*****************************************************************************/
bool TranslationExperimentReport::close(bool abort)
{

  if (! abort) {
    Verbose::out(ADT_VERBOSE_NORMAL, "TranslationExperimentReport::generate called: ");
    Verbose::out(ADT_VERBOSE_EXCEPTION, "dmet3-report-guid=" + m_guid);
    Verbose::out(ADT_VERBOSE_NORMAL,  m_comprehensiveReportName + COMPLETE_REPORT_FILE_EXT);
    Verbose::out(ADT_VERBOSE_NORMAL,  m_summaryReportName + COMPLETE_REPORT_FILE_EXT);
    Verbose::out(ADT_VERBOSE_NORMAL,  m_uncalledReportName + COMPLETE_REPORT_FILE_EXT);
  }

  bool okClose = true;
  if (m_openCalled) {
    okClose = _closeComprehensiveFile(abort) && okClose;
    okClose = _closeSummaryFile(abort) && okClose;
    okClose = _closeUncalledFile(abort) && okClose;
    m_openCalled = false;
    m_experimentFileCount = 0;
    m_experimentCount = 0;
    m_guid.clear();
  }

  return okClose;
}
// end TranslationExperimentReport::close
/*****************************************************************************/
/*****************************************************************************/
/**
 * _generateOriginalBaseCall
 * Synopsis:
 *
 *   Helper function to generate that take the report alleles stored
 * in the translation table and applies them to the original basecall as
 * input by the override file. 
 *
 * 
 * @param probeSet - The unique probe set in the translation table.
 * @param originalBasecall - the original basecall to get the report allele for
 * @param ttm - the single instance TranslationTableModel
 *
 * @return - the input originalBasecall or a translated originalBasecall if
 * report alleles exist. 
 */
/*****************************************************************************/
static std::string _generateOriginalBasecall(const std::string probeSet, std::string originalBasecall, TranslationTableModel & ttm  ) {

  std::string first, second, translatedBasecall;

  if ( originalBasecall.empty() ) {
    return originalBasecall;
  }

  if (!pcrecpp::RE("^([^/]+)/?$").FullMatch(originalBasecall, &first)) {

    if (! pcrecpp::RE("^([^/]+)/([^/]+)$").FullMatch(originalBasecall, &first, &second)) {
      APT_ERR_ABORT(probeSet + ": unknown call? " + originalBasecall);
    }
  }

  
  if ( second.empty() ) {
    if ( first == "0" ) {
      translatedBasecall = "ZeroCopyNumber";
    }
    else {
      translatedBasecall = ttm.getReportAllele( probeSet, first );
    }
  }
  else {
    translatedBasecall = ttm.getReportAllele( probeSet, first ) + '/' +
      ttm.getReportAllele( probeSet, second );
  }

  return translatedBasecall;
  
}
// end _generateOriginalBaseCall
/*****************************************************************************/
/**
 * TranslationExperimentReport::generate
 * Synopsis:
 *
 *   Creates the following reports in tandem due to the similarity in
 * in requirements:
 * 1.) Comprehensive
 * 2.) Summary
 * 3.) Uncalled
 * 
 * Reporting is accomplished on a per experiment basis in order
 * to avoid using excessive memory. There is no reporting done
 * across experiments. All data collected during the generate
 * process is expunged before the next call to generate. 
 * 
 * The translation table (file) contains annotation data that is
 * not used during process and hence not found in the ExperimentResults
 * object passed in. Therefore the translation table is passed in so
 * that the annotation data can be "passed through" to the reports.
 * 
 * Some pass through data from the translation table to the the reports
 * does in fact have some processng or business logic swirling the data. The
 * relevant alleles are one such example.
 *
 * 
 * @param rte - the single instance run time environment
 * @param ttm - the single instance translation table with pass through data.
 * @param er  - the ExperimentResuls container of ExperimentGeneResults for per gene results. 
 * 
 * @return true - the report was generated without error. 
 *
 */
/*****************************************************************************/
bool TranslationExperimentReport::generate(RunTimeEnvironment & rte,
    TranslationTableModel & ttm,
    std::map<std::string, ExperimentResults *>&er)
{

  bool okReport = true;

  pcrecpp::RE reUnknown("UNK");
  if (! m_openCalled) {
    _open(rte);
  }

  std::map< std::string, GeneHaplotypeCallMapElement> ghcm;
  std::map< std::string, int > markerIndex;
  std::map< std::string, bool > geneSummaryWritten;

  _generateGeneHaplotypeCallMap(rte, er, ghcm);

  int count = 0;

  std::map<std::string, ExperimentResults*>::const_iterator itSER;
  std::vector< std::string > sampleInfo;

  // For each experiment, recapitulate the translation table.



  for (itSER = er.begin(), geneInterpretationCode.clear();
       (itSER != er.end()) && okReport ;  itSER++) {

    m_experimentCount++;

    // IGNORE Experiments with no marker calls for DMET2.
    if (rte.m_adtOpts.m_dmet2Calling && (itSER->second->m_markerCallCount == 0)) {
      continue;
    }

    if (m_sitm != NULL) {
      sampleInfo = m_sitm->getSampleInfo(itSER->second->m_experimentGuid.empty() ? itSER->second->m_experiment : itSER->second->m_experimentGuid);
    }


    std::vector<TranslationRow> cRows;
    std::map<std::string, int> probeSetRow;
    int geneCount = 0;
    int markerCount = 0;
    std::string prevGene;
    for (int i = 1; i < ttm.size(); i++, count++) {

      // Skip commented out rows as well as header rows.
      if (ttm.ignoreRow(i) || ( !((ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)] == "Y") || (ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)] == "N")))) {
        continue;
      }


      // Initialize

      std::string gene = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENE)];

      if (prevGene.empty() || (prevGene != gene)) {
        geneCount++;
        prevGene = gene;
        markerCount = 0;
        geneSummaryWritten[gene] = false;
      }

      ExperimentGeneResults *egr = NULL;

      if (itSER->second->m_geneResults.find(gene) != itSER->second->m_geneResults.end()) {
        egr = itSER->second->m_geneResults[gene];
      }


      TranslationRow cRowNew;
      TranslationRow * cRow = &cRowNew;


      std::string probeSet = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)];

      // Multiallelic Probe Set?
      bool isMultiAllelicRow = probeSetRow.count(probeSet) > 0;
      if (isMultiAllelicRow) {
        APT_ERR_ASSERT(probeSetRow[probeSet] < cRows.size(), "");
        cRow = &(cRows[probeSetRow[probeSet]]);
        APT_ERR_ASSERT(probeSet == cRow->m_probeSet, "");
      } else {
        markerCount++;
        probeSetRow[probeSet] = cRows.size();
        markerIndex[probeSet] = markerCount;
        cRow->m_isSummaryRow = false;
        cRow->m_isGeneSummary = false;
      }

      bool isProbeSetPresent = false;
      if (egr != NULL) {
        isProbeSetPresent
        = egr->m_probeSetToMarkerCallResults.count(probeSet) > 0;
      }

      int markerCallResultsIndex = 0;
      bool hasAlleleCalls = false;
      if (isProbeSetPresent && (egr != NULL)) {
        markerCallResultsIndex = egr->m_probeSetToMarkerCallResults[probeSet];
        hasAlleleCalls =  egr->m_callResults[markerCallResultsIndex].m_alleleCalls.size() > 0;
      }

      // Report dat

      // Index

      char buf[BUFSIZ];
      sprintf(buf, "%4.4d-%4.4d-%2.2d", m_experimentCount, geneCount, markerIndex[probeSet]);

      cRow->m_index       = buf;

      // CHP file name
      cRow->m_chpFileName = Fs::basename(itSER->second->m_experiment);

      // Gene
      cRow->m_gene       = gene;
      cRow->m_geneCopyNumberIndicator = ttm.m_geneCopyNumberIndicator[gene];

      // Haplotype Marker

      cRow->m_haplotypeMarker = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)];

      //bool isHaplotypeMarker = cRow->m_haplotypeMarker == "Y";

      // Known Calls, Unknown Calls
      std::map<std::string, GeneHaplotypeCallMapElement>::iterator itGHCM;
      itGHCM = ghcm.find(gene);

      int numKnownCalls = 0;
      int numUnknownCalls = 0;
      int geneCopyNumber = 2;

      if ( itGHCM != ghcm.end() ) {
        geneCopyNumber = itGHCM->second.m_copyNumber;
      }
      
      if (!isMultiAllelicRow && (itGHCM != ghcm.end())) {
        numKnownCalls = itGHCM->second.m_knownCall.size();
        for (size_t i = 0; i <  itGHCM->second.m_knownCall.size(); i++) {
          if (cRow->m_knownCall.empty()) {
            cRow->m_knownCall = itGHCM->second.m_knownCall[i];
          } else {
            cRow->m_knownCall = cRow->m_knownCall  + "," + itGHCM->second.m_knownCall[i];
          }
        }
        numUnknownCalls = itGHCM->second.m_unknownCall.size();
        for (size_t i = 0; i <  itGHCM->second.m_unknownCall.size(); i++) {
          if (cRow->m_unknownCall.empty()) {
            cRow->m_unknownCall = itGHCM->second.m_unknownCall[i];
          } else {
            cRow->m_unknownCall = cRow->m_unknownCall  + "," + itGHCM->second.m_unknownCall[i];
          }
        }
      }

      // Summary Flag
      if (isMultiAllelicRow) {
        if (cRow->m_summaryFlag.empty()) {
          cRow->m_summaryFlag = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)];
        } else  if (!ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)].empty() && (cRow->m_summaryFlag != ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)])) {
          cRow->m_summaryFlag = cRow->m_summaryFlag + "," + ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)];
        }
      } else {
        cRow->m_summaryFlag = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)];

      }

      // Relevant Alelles
      _generateRelevantAlleles(ttm, cRow, i, isMultiAllelicRow);

      // Marker Name (External Id)
      cRow->m_markerName = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_COMMON_NAME)];

      // Probe Set
      cRow->m_probeSet = probeSet;

      // Reference Base
      _generateReferenceBase( i, cRow, ttm);

      // Variant Base
      _generateVariantBase( i, cRow, ttm, isMultiAllelicRow);

      // Call
      cRow->m_call = "";
      if (isProbeSetPresent && (egr != NULL) && (egr->m_callResults[markerCallResultsIndex].m_alleleCalls.size() > 0)) {
        APT_ERR_ASSERT(markerCallResultsIndex < egr->m_callResults.size(), "");
        cRow->m_call = egr->m_callResults[markerCallResultsIndex].m_alleleCalls[0].m_allele;

      }

      // Base Call
      _generateBaseCall(egr, i, markerCallResultsIndex, cRow, ttm, isMultiAllelicRow, isProbeSetPresent);


      /* HACK: m_copyNumber is set in _generateBaseCall */
      if (( egr != NULL) && (cRow->m_copyNumber == 0) &&
          (egr->m_callResults[markerCallResultsIndex].m_alleleCalls.size() > 0)) {
        cRow->m_call = egr->m_callResults[markerCallResultsIndex].m_alleleCalls[0].m_allele;
      }
      else if (cRow->m_baseCall == "NotAvailable") {
        cRow->m_call = "";
      }

      // Interpretation Code

      if ((egr != NULL) && geneInterpretationCode.count(gene) == 0) {
        _generateInterpretationCode(*egr, ttm, cRow, i, numKnownCalls, numUnknownCalls, geneInterpretationCode, geneCopyNumber);
      }

      if ( geneInterpretationCode.count(gene ) > 0 ) {
        cRow->m_interpretationCode = geneInterpretationCode[gene];
      }
      else {
        cRow->m_interpretationCode = "NC/PRA/NA";
      }

      if (isMultiAllelicRow) {
        if (cRow->m_changeForVariant.empty()) {
          cRow->m_changeForVariant = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)];
        } else if (! ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)].empty() && (cRow->m_changeForVariant != ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)])) {

          cRow->m_changeForVariant = cRow->m_changeForVariant + "," + ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)];
        }

      } else {
        cRow->m_changeForVariant = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)];
      }

      // CDNA Change
      if (isMultiAllelicRow) {
        if (cRow->m_CDNAChange.empty()) {
          cRow->m_CDNAChange = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CDNA)];
        } else if (!ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CDNA)].empty()) {
          cRow->m_CDNAChange = cRow->m_CDNAChange + "," + ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CDNA)];
        }
      } else  {
        cRow->m_CDNAChange = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CDNA)];
      }

      // Genomic Position
      if (isMultiAllelicRow) {
        if (cRow->m_genomicPosition.empty()) {
          cRow->m_genomicPosition = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENOME_POSITION)];
        } else if (!ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENOME_POSITION)].empty()) {
          cRow->m_genomicPosition = cRow->m_genomicPosition + "," + ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENOME_POSITION)];
        }
      } else {
        cRow->m_genomicPosition = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENOME_POSITION)];
      }

      // dbSNP Id
      cRow->m_dbSNPId = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DBSNP)];

      // Override
      // Check egr
      if ((m_gotm != NULL) && (egr != NULL)) {
        std::vector< std::string > overrideInfo = m_gotm->getOverrideInfo(egr->m_experimentName, cRow->m_probeSet);

        if (overrideInfo.size() == 2){
          cRow->m_originalBasecall = _generateOriginalBasecall( cRow->m_probeSet, overrideInfo[0], ttm);
          cRow->m_overrideComment  = overrideInfo[1];
          std::string invalidOriginalBasecall =  m_gotm->getOverrideInfo(egr->m_experimentName, cRow->m_probeSet, true)[0];
          if ( !(cRow->m_overrideComment.find("ZeroCopyNumber") > 0) && ((cRow->m_baseCall == invalidOriginalBasecall) || (cRow->m_originalBasecall == cRow->m_baseCall )) ){
            cRow->m_originalBasecall.clear();
            cRow->m_overrideComment.clear();
          }
        }
      }
      // Fix up ZeroCopyNumber 
      if ( (egr != NULL) && (geneCopyNumber == 0 ) ) {
        if ( cRow->m_haplotypeMarker == "N" )  {
          if ( cRow->m_interpretationCode == INTERPRETATION_CODES[TranslationInterpretationCodeMap::NoHAP].m_reportString) {
            cRow->m_interpretationCode = INTERPRETATION_CODES[TranslationInterpretationCodeMap::UNIQ].m_reportString;
          }
          cRow->m_knownCall = egr->m_geneZeroCopyHaplotype + "/" + egr->m_geneZeroCopyHaplotype;
        }
        else if ( cRow->m_knownCall.empty() ) {
          cRow->m_knownCall = egr->m_geneZeroCopyHaplotype + "/" + egr->m_geneZeroCopyHaplotype;
        }
      }

      // Summary Report determination.
      if (!cRow->m_summaryFlag.empty() && !cRow->m_call.empty() &&
          (cRow->m_summaryFlag != "N") && 
          RE_VAR.PartialMatch(cRow->m_call)) {
        cRow->m_isSummaryRow = true;
      } else if (RE_WILDCARD.FullMatch(cRow->m_baseCall)) {
        cRow->m_isSummaryRow = true;
      }
      if (cRow->m_isSummaryRow) {
        geneSummaryWritten[gene] = true;
      }

      // Uncalled Report determination
      _setIsUncalledRow(rte, *cRow);

      if (!isMultiAllelicRow) {
        cRows.push_back(*cRow);
      }


    } //for each translation table model marker

    for (size_t j = 0; j < cRows.size(); j++) {

      if (!geneSummaryWritten[cRows[j].m_gene]) {
        cRows[j].m_isGeneSummary = true;
        geneSummaryWritten[cRows[j].m_gene] = true;
      }
      _report(rte,cRows[j], sampleInfo);
    }

  } // for each gene, experiment result

  if ( rte.m_adtOpts.m_summaryReportSort ) {
    _reportSummarySorted(sampleInfo);
  }

  return okReport;

}
// end TranslationExperimentReport::generate
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::getComprehensiveReportName
 * Synopsis:
 *
 * Create a Translation report file name from various constant and
 * runtime environment settings. The name is cached in a private
 * variable after the first call.
 *
 * Note that the temporary extension of ".tmp" is not added here but
 * rather a requirement for the callee. 
 *
 * @param  rte       - the single instance RunTimeEnvironment with name options in rte.m_adtOpts that were passed in from the console or command line. 
 * 
 * @return comprehensivefileName - the name of the comprehensive report file.
 *
 */
/*****************************************************************************/
std::string TranslationExperimentReport::getComprehensiveReportName(const RunTimeEnvironment & rte)
{


  if (!m_comprehensiveReportName.empty()) {
    return m_comprehensiveReportName;
  }

  std::string comprehensiveFileName;
  if (!rte.m_adtOpts.m_outputReportPrefix.empty()) {
    comprehensiveFileName = rte.m_adtOpts.m_outputReportPrefix;
  } else if (rte.m_adtOpts.m_inputExperimentFiles.size() > 0) {
    comprehensiveFileName = _getFirstCHPFileRoot(rte);
  } else if (rte.m_adtOpts.m_inputGenoFile != "") {
    comprehensiveFileName = Fs::basename(rte.m_adtOpts.m_inputGenoFile);
    comprehensiveFileName = comprehensiveFileName.substr(0, comprehensiveFileName.find('.'));
  } else {
    comprehensiveFileName = name();
  }

  comprehensiveFileName = Fs::join(rte.m_adtOpts.m_outputDir,
                                       comprehensiveFileName+COMPREHENSIVE_REPORT_FILE_EXT);

  m_comprehensiveReportName = comprehensiveFileName;

  return comprehensiveFileName;

}
// end TranslationExperimentReport::getComprehensiveReportName
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::getSummaryReportName
 * Synopsis:
 *
 * Create a Translation report file name from various constant and
 * runtime environment settings. The name is cached in a private
 * variable after the first call.
 *
 * Note that the temporary extension of ".tmp" is not added here but
 * rather a requirement for the callee. 
 *
 * @param  rte      - the single instance RunTimeEnvironment which contains name options in rte.m_adtOpts passed in from the console or command line.
 * 
 * @return summaryfileName - the name of the summary report file.
 *
 */
/*****************************************************************************/
std::string TranslationExperimentReport::getSummaryReportName(const RunTimeEnvironment & rte)
{

  if (!m_summaryReportName.empty()) {
    return m_summaryReportName;
  }

  std::string summaryFileName;
  if (! rte.m_adtOpts.m_outputReportPrefix.empty()) {
    summaryFileName = rte.m_adtOpts.m_outputReportPrefix;
  } else if (rte.m_adtOpts.m_inputExperimentFiles.size() > 0) {
    summaryFileName = _getFirstCHPFileRoot(rte);
  } else if (!rte.m_adtOpts.m_inputGenoFile.empty()) {
    summaryFileName = Fs::basename(rte.m_adtOpts.m_inputGenoFile);
    summaryFileName = Fs::noextname1(summaryFileName);
  } else {
    summaryFileName = name();
  }

  summaryFileName = Fs::join(rte.m_adtOpts.m_outputDir,
                                 summaryFileName+SUMMARY_REPORT_FILE_EXT);

  m_summaryReportName = summaryFileName;

  return summaryFileName;

}
// end TranslationExperimentReport::getSummaryReportName
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::getUncalledReportName
 * Syopsis:
 *
 * Create a Translation report file name from various constant and
 * runtime environment settings. The name is cached in a private
 * variable after the first call.
 *
 * Note that the temporary extension of ".tmp" is not added here but
 * rather a requirement for the callee. 
 *
 * @param  rte      - the single instance RunTimeEnvironment which contains name options in rte.m_adtOpts passed in from the console or command line.
 * 
 * @return uncalledfileName - the name of the uncalled report file.
 *
 */
/*****************************************************************************/
std::string TranslationExperimentReport::getUncalledReportName(const RunTimeEnvironment & rte)
{

  if (!m_uncalledReportName.empty()) {
    return m_uncalledReportName;
  }

  std::string uncalledFileName;
  if (! rte.m_adtOpts.m_outputReportPrefix.empty()) {
    uncalledFileName = rte.m_adtOpts.m_outputReportPrefix;
  } else if (rte.m_adtOpts.m_inputExperimentFiles.size() > 0) {
    uncalledFileName = _getFirstCHPFileRoot(rte);
  } else if (rte.m_adtOpts.m_inputGenoFile != "") {
    uncalledFileName = Fs::basename(rte.m_adtOpts.m_inputGenoFile);
    uncalledFileName = Fs::noextname1(uncalledFileName);
  } else {
    uncalledFileName = name();
  }

  uncalledFileName = Fs::join(rte.m_adtOpts.m_outputDir,
                                  uncalledFileName+UNCALLED_REPORT_FILE_EXT);

  m_uncalledReportName = uncalledFileName;

  return uncalledFileName;

}
// end TranslationExperimentReport::getUncalledReportName
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_closeComprehensiveFile
 * Synopsis:
 * 
 * Close a TSV file of comprehensive report results.
 * 
 * !!!WARNING!!! under no circumstances should messaging via Verbose
 * take place when abort is passed as true. This will crash the
 * Windows console due to a bug in the Windows COM code. 
 * 
 * @param abort - if this close is happening in the midst of an abort exception.
 *
 * @return true - if the temporary file is moved to the non-temporary file name.
 */
/*****************************************************************************/
bool TranslationExperimentReport::_closeComprehensiveFile(bool abort)
{


  bool okMove = true;

  m_comprehensiveTsv.close();
  m_comprehensiveTsv.clear();

  if (!abort) {
    std::string tmpFile = m_comprehensiveReportName + TEMP_REPORT_FILE_EXT;
    std::string completeFile = m_comprehensiveReportName + COMPLETE_REPORT_FILE_EXT ;
    okMove = Fs::fileRename(tmpFile,  completeFile, false);

    if (!okMove) {
      Verbose::warn(ADT_VERBOSE_NORMAL, "Can't rename file " + tmpFile + " as " + completeFile + ".");
    }
  }

  return okMove;
}
// end TranslationExperimentReport::_closeComprehensiveFile
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_closeSummaryFile
 * Synopsis:
 * 
 * Close a TSV file for outputing summary report results.
 *
 * !!!WARNING!!! under no circumstances should messaging via Verbose
 * take place when abort is passed as true. This will crash the
 * Windows console due to a bug in the Windows COM code. 
 *
 * @param abort - if this close is happening in the midst of an abort exception.
 *
 * @return true - if the temporary file is moved to the non-temporary file name.
 */
/*****************************************************************************/
bool TranslationExperimentReport::_closeSummaryFile(bool abort)
{

  bool okMove = true;

  m_summaryTsv.close();
  m_summaryTsv.clear();

  if (!abort) {
    std::string tmpFile = m_summaryReportName + TEMP_REPORT_FILE_EXT;
    std::string completeFile = getSummaryFileName();
    okMove = Fs::fileRename(tmpFile,  completeFile, false);

    if (!okMove) {
      Verbose::warn(ADT_VERBOSE_NORMAL, "Can't rename file " + tmpFile + " as " + completeFile + ".");
    }
  }


  return okMove;

}
// end TranslationExperimentReport::_closeTranslationFile

/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::getSummaryFileName
 * Synopsis:
 *
 * Returns the name of the file storing the summary results.
 *
 * @return the file name.
 *
 */
/*****************************************************************************/
std::string TranslationExperimentReport::getSummaryFileName()
{
    std::string completeFile = m_summaryReportName + COMPLETE_REPORT_FILE_EXT ;
	return completeFile;
}
// end TranslationExperimentReport::getSummaryFileName

/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::getComprehensiveFileName
 * Synopsis:
 *
 * Returns the name of the file storing the summary results.
 *
 * @return the file name.
 *
 */
/*****************************************************************************/
std::string TranslationExperimentReport::getComprehensiveFileName()
{
    std::string completeFile = m_comprehensiveReportName + COMPLETE_REPORT_FILE_EXT ;
	return completeFile;
}
// end TranslationExperimentReport::getComprehensiveFileName

/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_closeUncalledFile
 * Synopsis:
 *
 * Close a TSV file for outputing uncalled marker report results.
 *
 * !!!WARNING!!! under no circumstances should messaging via Verbose
 * take place when abort is passed as true. This will crash the
 * Windows console due to a bug in the Windows COM code. 
 *
 * @param abort - if this close is happening in the midst of an abort exception.
 *
 * @return true - if the temporary file is moved to the non-temporary file name.
 *
 */
/*****************************************************************************/
bool TranslationExperimentReport::_closeUncalledFile(bool abort)
{

  bool okMove = true;

  m_uncalledTsv.close();
  m_uncalledTsv.clear();

  if (!abort) {
    std::string tmpFile = m_uncalledReportName + TEMP_REPORT_FILE_EXT;
    std::string completeFile = m_uncalledReportName + COMPLETE_REPORT_FILE_EXT ;
    okMove = Fs::fileRename(tmpFile,  completeFile, false);

    if (!okMove) {
      Verbose::warn(ADT_VERBOSE_NORMAL, "Can't rename file " + tmpFile + " as " + completeFile + ".");
    }
  }

  return okMove;

}
// end TranslationExperimentReport::_closeUncalledFile
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_generateBaseCall
 * Synopsis:
 * 
 * Helper function to TranslationExperimentReport::generate, broken
 * out for readability only and not for reuse. 
 *
 * Set the m_baseCall field in the TranslationRow structure to something like
 * "*1/ *1". 
 *
 * 
 * @param egr - the ExperimentGeneResult for the cRow being processed. 
 * @param ttmRow - the translation table model row index. 
 * @param markerCallResultsIndex - egr->m_callResults index for the cRow being passed. 
 * @param cRow - the output row being reported on, prepopulated with gene and  probe set id
 * @param ttm - the single instance translation table model
 * @param isMultiAllelicRow - if this probe set id being processed was previously processed. 
 * @param isProbeSetPresent - if the genotype data contained this probe set. 
 *
 * @return - cRow.m_baseCall filled in.
 */
/*****************************************************************************/
void TranslationExperimentReport::_generateBaseCall(ExperimentGeneResults *egr,  int ttmRow,  int markerCallResultsIndex,  TranslationRow * cRow, TranslationTableModel & ttm, bool isMultiAllelicRow, bool isProbeSetPresent)
{


  if (!isMultiAllelicRow) {
    cRow->m_copyNumber = 2;
    cRow->m_baseCall = std::string("NotAvailable");
  }

  if (isProbeSetPresent && (egr != NULL)) {
    if ((egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.size() > 0) &&
        (egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.size() > 0) &&
        pcrecpp::RE("PRA|PossibleRareAllele").FullMatch(egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0])) {
      cRow->m_baseCall = "PossibleRareAllele";
    } else if (egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_copyNumber == 0) {
      if ((egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.size() == 0) ||
          (pcrecpp::RE("NA|NotAvailable").FullMatch(egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0]))) {
        cRow->m_baseCall = std::string("NotAvailable");
      } else {
        cRow->m_baseCall = std::string("ZeroCopyNumber");
        cRow->m_copyNumber = 0;
      }
    } else if (egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_copyNumber == 1) {

      cRow->m_copyNumber = 1;
      std::string definingValue = cRow->m_call.substr(cRow->m_call.find('/') + 1);

      if (definingValue  == "Ref") {
        cRow->m_baseCall = cRow->m_referenceReportBase;
      } else if (definingValue  == "Var") {
        cRow->m_baseCall = egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0];
      } else if (pcrecpp::RE("NC|NoCall").FullMatch(egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0])) {
        cRow->m_baseCall = "NoCall";
      }
    } else if (cRow->m_call == "Ref/Ref") {
      cRow->m_baseCall = cRow->m_referenceReportBase + "/" + cRow->m_referenceReportBase;
    } else if (cRow->m_call == "Ref/Var") {
      std::string variance =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0];
      if (variance ==  cRow->m_referenceBase) {
        variance =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid2.m_ceSet.begin()->second.m_bases[0];
      }
      if (variance == ttm.m_rows[ttmRow][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)]) {
        cRow->m_baseCall = cRow->m_referenceReportBase + "/" + ttm.getReportAllele( cRow->m_probeSet, variance);
      }
    } else if (cRow->m_call == "Var/Var") {
      std::string variance1 =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0];
      if (variance1 == ttm.m_rows[ttmRow][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)]) {
        std::string variance2 =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid2.m_ceSet.begin()->second.m_bases[0];
        cRow->m_baseCall = ttm.getReportAllele( cRow->m_probeSet, variance1) + "/" + ttm.getReportAllele(cRow->m_probeSet, variance2) ;
      }
    } else if (pcrecpp::RE("NC|NoCall").FullMatch(egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0])) {
      cRow->m_baseCall = "NoCall";
    }
  }

  return;
}
// end TranslationExperimentReport::_generateBaseCall
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_generateGeneHaplotypeCallMap:
 * Synopsis:
 *
 * Helper function to TranslationExperimentReport::generate, broken
 * out for readability only and not for reuse. To be called once per gene
 * and not once per row. 
 *
 *  Haplotype calls are per gene, not per marker. Therefore we need to compile
 *  a list of known and unkown calls per gene. These calls then get repeated
 *  per every assay id in the haplotype group for that gene.
 *
 * 
 * @param rte  - the single instance run time environment.
 * @param er   - the ExperimentResults container of all the ExperimentGeneResults
 * @param geneHaplotypeMap - the index of entire gene saved for all rows of geno data for this 
 *
 * @return - geneHaplotypeMap created for this gene
 * 
 */
/*****************************************************************************/
void TranslationExperimentReport::_generateGeneHaplotypeCallMap(RunTimeEnvironment & rte, std::map<std::string, ExperimentResults *>&er, std::map< std::string, GeneHaplotypeCallMapElement > & geneHaplotypeMap)
{

  pcrecpp::RE reUnknown("UNK");

  std::map<std::string, ExperimentResults*>::iterator itSER;

  for (itSER = er.begin(); (itSER != er.end()) ; itSER++) {

    if (itSER->second->m_markerCallCount == 0) {
      continue;
    }

    std::map<std::string, ExperimentGeneResults*>::iterator itSEGR;

    for (itSEGR = itSER->second->m_geneResults.begin();
         itSEGR != itSER->second->m_geneResults.end(); itSEGR++) {

      ExperimentGeneResults & egr = *(itSEGR->second) ;

      for (size_t i = 0; i < egr.size(); i++) {

        std::string gene        = egr.m_callResults[i].m_geneName;
        geneHaplotypeMap[gene].m_gene = gene;
        if ( egr.m_callResults[i].m_experimentChromatid1.m_copyNumber == 0 ) {
          geneHaplotypeMap[gene].m_copyNumber = 0;
        }
        if (egr.m_callResults[i].getCallType() != ADT_CALL_TYPE_HAPLOTYPE_GROUP) {
          continue;
        }

        for (int j = 0; j < egr.m_callResults[i].size() ; j++) {

          std::string call =  egr.m_callResults[i].m_alleleCalls[j].m_allele;

          geneHaplotypeMap[gene].m_gene = gene;
          if (reUnknown.PartialMatch(call)) {
            geneHaplotypeMap[gene].m_unknownCall.push_back(call);
          } else {
            geneHaplotypeMap[gene].m_knownCall.push_back(call);
          }

        } // Foreach  call, if they are more than one (wild cards )

      } // foreach experiment gene call list

    } // for each experiment gene

  } // foreach experiment

  return;

}
// end TranslationExperimentReport::_generateGeneHaplotypeCallMap
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_generateInterpretationCode:
 * Synopsis:
 *
 * Helper function to TranslationExperimentReport::generate, broken
 * out for readability only and not for reuse. 
 *
 *  Contains the busines logic for generating the interpretation code
 *  Report field
 *
 * Interpretation Code Description:
 * UNIQ: Unique haplotype pair ( 1 known call, 0 unknown calls )
 * UNIQ+UNK: Unique annotated haplotype pair other haplotype pairs involving
 * anannotated haplotypes also possible
 *       (1 known call, 1+ unkown calls )
 * MULT: Multiple haplotype pairs possible due to phase ambiguity
 *           ( 2+ known calls, 0 unkowns )
 * MULT+UNK: Multiple annotated haplotype calls possible other haplotype
 * pairs involving unannotated haplotypes also possible
 *           (2+ known calls, 1+ unkowns)
 * UNDH: Only haplotype pairs with undefined haplotypes possible
 *         (0 known calls, 1+ unkown calls)
 * NC/PRA: NoCall or PossibleRareAllele call for one or more markers
 * resulting in multiple haplotype pairs
 *         (1+ assay id is NC or PCRA)
 * NoHAP: No haplotypes defined at this gene
 *        ( translation table "Haplotype" column is all 'N' for assay ids)
 * NA: No genotyping data available
 *        (Markers ignored by override or filter operations).
 *
 * @param egr - all call results for the experiment gene.
 * @param ttm - the translation table
 * @param cRow - has the already created known and unkown haplotype calls
 * @param ttmRow - the row of the current assay id in question.
 * @param numKnownCalls - the number of calls in cRow->m_knownCall std::string.
 * @param numUnknownCalls - the number of calls in cRow->m_unknownCall std::string.
 * @param geneInterpretationCode - returned.
 _
 * @return - geneInterpretationCode is filled in.
 */
/*****************************************************************************/
void  TranslationExperimentReport::_generateInterpretationCode(
  const ExperimentGeneResults & egr,
  TranslationTableModel & ttm,
  TranslationRow * cRow,
  int ttmRow,
  int numKnownCalls,
  int numUnknownCalls,
  std::map< std::string, std::string> & geneInterpretationCode,
  int geneCopyNumber)
{

  int headerRow = ttm.getHeaderRow(ttmRow);
  std::map< std::string , std::vector<std::string> > haplotypeDefiningBases;
  std::map< std::string, std::vector<std::string> >  probeSetDefiningBases;

  // NoHAP, scan the translation table for 'Haplotype' = Y.
  bool geneHasHaplotype = false;
  bool hasNCOrPRA       = false;
  bool isCopyNumberZero = (geneCopyNumber == 0);
  bool zeroCopyCallHasNotAvailable  = false;
  
  for (int row = headerRow + 1;
       (row < ttm.size()) && (ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_GENE)] == cRow->m_gene); row++) {

    if ( ttm.ignoreRow(row) ) {
      continue;
    }
    std::string probeSet = ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)];
    std::map< std::string, int>::const_iterator itSI = egr.m_probeSetToMarkerCallResults.find(probeSet);
    
    if (ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)] == "Y") {

      geneHasHaplotype = true;

      if (itSI !=  egr.m_probeSetToMarkerCallResults.end()) {

        int callIndex = itSI->second;

        // NOTE: Only need to check one chromatid because if one is
        // a wild card it is the case they will both be wildcard.
        // Wildcard is defined as either No Call, Not Available or Possible Rare Allele.

        zeroCopyCallHasNotAvailable = zeroCopyCallHasNotAvailable ||
          egr.m_callResults[callIndex].m_experimentChromatid1.m_hasZeroCopyNumberNotAvailable;
        if (egr.m_callResults[callIndex].m_experimentChromatid1.m_ceSet.size() > 0) {
          if (zeroCopyCallHasNotAvailable || CallElement::isWildCardBase(egr.m_callResults[callIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0])) {
            hasNCOrPRA = true;
          }
        }
      }
    }
  }

  if (!geneHasHaplotype && !isCopyNumberZero){
    geneInterpretationCode[cRow->m_gene] =  INTERPRETATION_CODES[TranslationInterpretationCodeMap::NoHAP].m_reportString;
    return;
  }

  // NC/PCRA, overrides all other codes. If there is a No Call or
  // Possible Rare Allele then that wins.
  if (hasNCOrPRA ) {
    //    if ( !geneHasHaplotype && isCopyNumberZero ) {
    if ( isCopyNumberZero ) {
      geneInterpretationCode[cRow->m_gene] = INTERPRETATION_CODES[TranslationInterpretationCodeMap::UNIQ].m_reportString;
    }
    else {
      geneInterpretationCode[cRow->m_gene] = INTERPRETATION_CODES[TranslationInterpretationCodeMap::NC_PRA_NA].m_reportString;

    }
    return;
  }


  // Call count. The rest of the codes are dependent on the count of
  // known and unkown calls.

  // UNIQ

  if (((numKnownCalls == 1) && (numUnknownCalls == 0)) ||
      (!geneHasHaplotype && isCopyNumberZero)) {
    geneInterpretationCode[cRow->m_gene] = INTERPRETATION_CODES[TranslationInterpretationCodeMap::UNIQ].m_reportString;
    return;
  }

  // UNIQ+UNK

  if (numKnownCalls == 1)  {
    geneInterpretationCode[cRow->m_gene] = INTERPRETATION_CODES[TranslationInterpretationCodeMap::UNIQ_UNK].m_reportString;
   return;
  }

  //MULTI
  if ((numKnownCalls > 1) && (numUnknownCalls == 0)) {
    geneInterpretationCode[cRow->m_gene] = INTERPRETATION_CODES[TranslationInterpretationCodeMap::MULT].m_reportString;
    return;
  }

  //MULTI+UNK

  if ((numKnownCalls > 1)) {
    geneInterpretationCode[cRow->m_gene] = INTERPRETATION_CODES[TranslationInterpretationCodeMap::MULT_UNK].m_reportString;
    return;
  }

  // UNDH

  if (numUnknownCalls > 0) {
    geneInterpretationCode[cRow->m_gene] = INTERPRETATION_CODES[TranslationInterpretationCodeMap::UNDH].m_reportString;
    return;
  }


  return;

}
// end TranslationExperimentReport::_generateInterpretationCode
/*****************************************************************************/
/*****************************************************************************/
/***
 * TranslationExperimentReport::_generateReferenceBase:
 * Synopsis:
 *
 * Helper function to TranslationExperimentReport::generate, broken
 * out for readability only and not for reuse. 
 *
 * The reference base is simply the translation table base in the
 * "Reference" column except when a report allele doesn't agree.
 *
 *
 * @param ttmRow - the row in question being reported in.
 * @param ttm - the in memory copy of the Translation Table
 * @param cRow - to be returned output row data structure
 *
 * @return - sets cRow->m_referenceBase, cRow->m_referenceReportBase
 */
/*****************************************************************************/
void TranslationExperimentReport::_generateReferenceBase(
    int ttRow,
    TranslationRow * cRow,
    TranslationTableModel & ttm )
{

  cRow->m_referenceBase = ttm.m_rows[ttRow][ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE)];



  cRow->m_referenceReportBase = ttm.getReportAllele( cRow->m_probeSet, cRow->m_referenceBase);

  return;

}
// end TranslationExperimentReport::_generateReferenceBase
/*****************************************************************************/
/*****************************************************************************/
/***
 * TranslationExperimentReport::_getFirstCHPFileRoot
 * Synopsis:
 *
 * Helper function for the various "getXXXXReportName" methods.
 *
 * If the base name option "-b" is not passed in then the
 * first CHP file name root is used as the root for all report files.
 * 
 * This method is used to generate report base names when no report
 * base name is provided.
 *
 * @param rte - contains the option 'inputExperimentListFiles'.
 *
 *
 * @return - the root of the file, or the empty std::string if not found.
 */
/*****************************************************************************/
std::string TranslationExperimentReport::_getFirstCHPFileRoot(const RunTimeEnvironment & rte)
{

  std::string firstCHPFileRoot = "";

  if (rte.m_adtOpts.m_inputExperimentFiles.size() > 0) {

    pcrecpp::RE("^(.*)(?:\\.[Cc][Hh][Pp]$)").FullMatch(Fs::basename(rte.m_adtOpts.m_inputExperimentFiles[0]), &firstCHPFileRoot);
  }


  return firstCHPFileRoot;


}
//end TranslationExperimentReport::_getFirstCHPFileRoot
/*****************************************************************************/
/*****************************************************************************/
/***
 * TranslationExperimentReport::_generateVariantBase:
 * Synopsis:
 *
 * Helper function to TranslationExperimentReport::generate, broken
 * out for readability only and not for reuse.
 * 
 * The variant base starts as the "Variant" column in the translation
 * table. For multi-allelic probe sets then the variant base is
 * a comma separated list. In addition report alleles override
 * the translation table. 
 *
 *
 * 
 * 
 * @param ttmRow - the row in question being reported in.
 * @param ttm - the in memory copy of the Translation Table
 * @param cRow - to be returned output row data structure
 * @param isMultiAllelicRow - if this is not the first time this probe set has been seen. 
 *
 * @return - sets cRow->m_variantBase, cRow->m_variantReportBase
 */
/*****************************************************************************/
void TranslationExperimentReport::_generateVariantBase(
  int ttRow,
  TranslationRow * cRow,
  TranslationTableModel & ttm,
  bool isMultiAllelicRow )
{


  std::string variantBase = ttm.m_rows[ttRow][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)];

  if (isMultiAllelicRow) {
    cRow->m_variantBase = cRow->m_variantBase + "," + variantBase;
  } else {
    APT_ERR_ASSERT(cRow->m_variantBase.empty(), "");
    cRow->m_variantBase = variantBase;
  }

  cRow->m_variantUncalledBase = variantBase;

  variantBase = ttm.getReportAllele( cRow->m_probeSet, variantBase );

  if (isMultiAllelicRow) {
    cRow->m_variantReportBase = cRow->m_variantReportBase + "," + variantBase;
  } else {
    APT_ERR_ASSERT(cRow->m_variantReportBase.empty(), "");
    cRow->m_variantReportBase = variantBase;
  }

  return;

}
// end TranslationExperimentReport::_generateVariantBase
/*****************************************************************************/
/***
 * TranslationExperimentReport::_generateRelevantAlleles:
 * Synopsis:
 *
 * Set the "m_relevantAlleles" attribute of the TranslationRow object.
 * Broken out of "generate" for readability purposes only.
 *
 * The relevant alleles is the list of all alleles that this probe set
 * does not participate in. 
 *
 * @param ttmRow - the row in question being reported in.
 * @param ttm - the in memory copy of the Translation Table
 * @param cRow - to be returned output row data structure
 * @param isMultiAllelicRow - if this is not the first time this probe set has been seen. 
 *
 * @return - cRow.m_relevantAlleles filled in.
 */
/*****************************************************************************/
void TranslationExperimentReport::_generateRelevantAlleles(
  TranslationTableModel & ttm,
  TranslationRow * cRow,
  int ttmRow,
  bool isMultiAllelicRow)
{
  // RELEVANT ALLELES
  int headerRow = ttm.getHeaderRow(ttmRow);

  std::string relevantAlleles;
  for (int j = ttm.getColumnIndex(ADT_DMET3_TT_ALLELE_START) + 1; j < ttm.m_rows[ttmRow].size(); j++) {
    if (!ttm.m_rows[ttmRow][j].empty() && ((cRow->m_haplotypeMarker == "N") || (ttm.m_rows[headerRow][j].c_str())[0] != '#')) {
      relevantAlleles = relevantAlleles + ttm.m_rows[headerRow][j] + ",";
    }
  }
  // Chop off the trailing comma.
  if (!relevantAlleles.empty()) {
    relevantAlleles.resize(relevantAlleles.size() - 1);
  }
  if (isMultiAllelicRow) {
    if (! relevantAlleles.empty()) {
      cRow->m_relevantAlleles = cRow->m_relevantAlleles  + "," + relevantAlleles;
    }
  } else {
    cRow->m_relevantAlleles = relevantAlleles;
  }

  if (!cRow->m_relevantAlleles.empty() &&
      cRow->m_relevantAlleles[ cRow->m_relevantAlleles.size() - 1 ] == ',') {
    cRow->m_relevantAlleles.resize(cRow->m_relevantAlleles.size() - 1);
  }
  return;

}
// end TranslationExperimentReport::_generateRelevantAlleles
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_getNumColumns:
 * Synopsis:
 *
 *  Helper function to various other helper functions.
 *
 *  Returns either the number of fixed columns for the comprehensive or
 *  summary reports only. Not to be used for the uncalled report.
 *  Pass a boolean as true to  include the number of sample columns.
 *
 * 
 * @param includeSampleColumns - boolean indicates whether to include the sample info columsn. 
 * 
 * @return - the number of columns
 */
/*****************************************************************************/
int TranslationExperimentReport::_getNumColumns(bool includeSampleColumns)
{

  if (includeSampleColumns) {
    APT_ERR_ASSERT(m_totalColumnsWithSampleInfo > 0, "");
    return m_totalColumnsWithSampleInfo;
  }

  return m_numFixedColumns;
}
// end TranslationExperimentReport::_getNumColumns
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_open:
 * Synopsis:
 *
 * Helper function to TranslationExperimentReport::generate that wraps
 * the open call to summary, comprehensive and the uncalled reports.
 *
 * The private variable "m_openCalled" protects against calling this
 * multiple times. 
 * 
 *
 * @param rte - the single instance RunTimeEnvironment. 
 *
 * @return - opens the report files or throws exception.
 */
/*****************************************************************************/
void TranslationExperimentReport::_open(const RunTimeEnvironment & rte)
{

  if (m_openCalled) {
    return;
  }

  std::vector< std::string > sampleHeaders;

  if (m_sitm != NULL) {
    sampleHeaders = m_sitm->getHeaderAsVector() ;
  }

  m_guid = affxutil::Guid::GenerateNewGuid();

  _openComprehensiveFile(rte, sampleHeaders);

  _openSummaryFile(rte, sampleHeaders);

  _openUncalledFile(rte);

  m_openCalled = true;

  m_experimentFileCount++;

  return;
}
// end TranslationExperimentReport::_open
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_openComprehensiveFile
 * Synopsis:
 * 
 * Helper function to TranslationExperimentReport::generate broken out
 * for readability.
 *
 * Open a TSV file for outputing comprehensive report results.
 * Opens the TSV file, writes file headers and column headers.
 *
 *
 * @param rte - the single instance RunTimeEnvironment.
 * @param sampleInfoHeaders - the dynamic sample headers that are not parsed
 *
 */
/*****************************************************************************/
void TranslationExperimentReport::_openComprehensiveFile(const RunTimeEnvironment & rte, std::vector< std::string > & sampleInfoHeaders)
{

  std::string completeName = getComprehensiveReportName(rte) + COMPLETE_REPORT_FILE_EXT;
  if (Fs::fileExists(completeName))
    Fs::rm(completeName, false);

  // describe our columns.
  for (int i = 0; i < _getNumColumns(); i++) {
    m_comprehensiveTsv.defineColumn(0, i, COMPREHENSIVE_TSV_COLUMNS[i]);
  }


  m_totalColumnsWithSampleInfo = _getNumColumns() + sampleInfoHeaders.size();

  std::vector< std::string >::iterator itS  = sampleInfoHeaders.begin();
  if (itS != sampleInfoHeaders.end()) {
    itS++;
  }

  for (int i = _getNumColumns();
       itS != sampleInfoHeaders.end(); i++, itS++) {
    m_comprehensiveTsv.defineColumn(0, i, *itS);
  }

  m_comprehensiveTsv.addHeaderComment(" For research use only. Not for diagnostic purposes.");
  m_comprehensiveTsv.addHeader("dmet3-report-guid", m_guid);
  m_comprehensiveTsv.addHeader("Program", Fs::basename(rte.m_programName));
  m_comprehensiveTsv.addHeader("Version", ADT_VERSION);
  m_comprehensiveTsv.addHeader("Date", Util::getTimeStamp());
  m_comprehensiveTsv.addHeader("TranslationFile", rte.m_adtOpts.m_inputTTableFile);
  m_comprehensiveTsv.addHeader("MarkerList", rte.m_adtOpts.m_inputMarkerListFile);
  m_comprehensiveTsv.addHeader("HaplotypeReportOption", rte.m_adtOpts.m_useFirstDupAlleleDefHeaderText);

  m_comprehensiveTsv.addHeader("GenotypeOverrideFile", rte.m_adtOpts.m_inputGenotypeOverrideFile);
  m_comprehensiveTsv.addHeaderComment("Interpretation Code Description:");
  m_comprehensiveTsv.addHeaderComment("UNIQ: Unique haplotype pair");
  m_comprehensiveTsv.addHeaderComment("MULT: Multiple haplotype pairs possible due to phase ambiguity");
  m_comprehensiveTsv.addHeaderComment("UNIQ+UNK: Unique annotated haplotype pair, with other haplotype pairs requiring unannotated haplotypes also possible");
  m_comprehensiveTsv.addHeaderComment("MULT+UNK: Multiple annotated haplotype calls possible, with other haplotype pairs requiring unannotated haplotypes also possible");
  m_comprehensiveTsv.addHeaderComment("UNDH: Only haplotype pairs with undefined haplotypes possible");
  //  m_comprehensiveTsv.addHeaderComment("6_CN: Copy Number state indicates deletion or duplication alleles");

  m_comprehensiveTsv.addHeaderComment("NC/PRA/NA: NoCall, PossibleRareAllele, or NotAvailable call for one or more markers resulting in multiple haplotype pairs");
  m_comprehensiveTsv.addHeaderComment("NoHAP: No haplotypes defined at this gene");
  m_comprehensiveTsv.writeTsv_v1(getComprehensiveReportName(rte) + TEMP_REPORT_FILE_EXT);


  return;

}
// end TranslationExperimentReport::_openComprehensiveFile
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_openSummaryFile
 * Synopsis:
 *
 * Helper function to TranslationExperimentReport::generate broken out
 * for readability.
 *
 * Open a TSV file for outputing summary report results.
 * Opens the TSV file, writes file headers and column headers.
 *
 *
 * @param rte - the single instance RunTimeEnvironment.
 * @param sampleInfoHeaders - the dynamic sample headers that are not parsed.
 */
/*****************************************************************************/
void TranslationExperimentReport::_openSummaryFile(const RunTimeEnvironment & rte, std::vector< std::string > & sampleInfoHeaders)
{

  std::string completeName = getSummaryReportName(rte) + COMPLETE_REPORT_FILE_EXT;
  if (Fs::fileExists(completeName))
    Fs::rm(completeName, false);

  // describe our columns.
  for (int i = 0; i < _getNumColumns(); i++) {
    m_summaryTsv.defineColumn(0, i, COMPREHENSIVE_TSV_COLUMNS[i]);
  }


  m_totalColumnsWithSampleInfo = _getNumColumns() + sampleInfoHeaders.size();

  std::vector< std::string >::iterator itS  = sampleInfoHeaders.begin();
  if (itS != sampleInfoHeaders.end()) {
    itS++;
  }

  for (int i = _getNumColumns();
       itS != sampleInfoHeaders.end(); i++, itS++) {
    m_summaryTsv.defineColumn(0, i, *itS);
  }

  m_summaryTsv.addHeaderComment(" For research use only. Not for diagnostic purposes.");
  m_summaryTsv.addHeader("dmet3-report-guid", m_guid);
  m_summaryTsv.addHeader("Program", Fs::basename(rte.m_programName));
  m_summaryTsv.addHeader("Version", ADT_VERSION);
  m_summaryTsv.addHeader("Date", Util::getTimeStamp());
  m_summaryTsv.addHeader("TranslationFile", rte.m_adtOpts.m_inputTTableFile);
  m_summaryTsv.addHeader("MarkerList", rte.m_adtOpts.m_inputMarkerListFile);
  m_summaryTsv.addHeader("HaplotypeReportOption", rte.m_adtOpts.m_useFirstDupAlleleDefHeaderText);
  m_summaryTsv.addHeader("GenotypeOverrideFile", rte.m_adtOpts.m_inputGenotypeOverrideFile);
  m_summaryTsv.addHeaderComment("Interpretation Code Description:");
  m_summaryTsv.addHeaderComment("UNIQ: Unique haplotype pair");
  m_summaryTsv.addHeaderComment("MULT: Multiple haplotype pairs possible due to phase ambiguity");
  m_summaryTsv.addHeaderComment("UNIQ+UNK: Unique annotated haplotype pair, with other haplotype pairs requiring unannotated haplotypes also possible");
  m_summaryTsv.addHeaderComment("MULT+UNK: Multiple annotated haplotype calls possible, with other haplotype pairs requiring unannotated haplotypes also possible");
  m_summaryTsv.addHeaderComment("UNDH: Only haplotype pairs with undefined haplotypes possible");
  //  m_summaryTsv.addHeaderComment("6_CN: Copy Number state indicates deletion or duplication alleles");

  m_summaryTsv.addHeaderComment("NC/PRA/NotAvailable: NoCall, PossibleRareAllele, or NotAvailable call for one or more markers resulting in multiple haplotype pairs");
  m_summaryTsv.addHeaderComment("NoHAP: No haplotypes defined at this gene");
  m_summaryTsv.writeTsv_v1(getSummaryReportName(rte) + TEMP_REPORT_FILE_EXT);


  return;

}
// end RegressionExperimentReport::_openSummaryFile
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_openUncalledFile
 *
 * Helper function to TranslationExperimentReport::generate broken out
 * for readability.
 *
 * Open a TSV file for outputing uncalled report results.
 * Opens the TSV file, writes file headers and column headers.
 *
 *
 * @param rte - the single instance RunTimeEnvironment.
 * @param sampleInfoHeaders - the dynamic sample headers that are not parsed.
 *
 */
/*****************************************************************************/
void TranslationExperimentReport::_openUncalledFile(const RunTimeEnvironment & rte)
{


  std::string completeName = getUncalledReportName(rte) + COMPLETE_REPORT_FILE_EXT;
  if (Fs::fileExists(completeName))
    Fs::rm(completeName, false);

  // describe our columns.
  int numUncalledColumns = sizeof(UNCALLED_TSV_COLUMNS) / sizeof(UNCALLED_TSV_COLUMNS[0]);

  for (int i = 0; i < numUncalledColumns; i++) {
    m_uncalledTsv.defineColumn(0, i, UNCALLED_TSV_COLUMNS[i]);
  }

  m_uncalledTsv.addHeaderComment(" For research use only. Not for diagnostic purposes.");
  m_uncalledTsv.addHeader("dmet3-report-guid", m_guid);
  m_uncalledTsv.addHeader("Program", Fs::basename(rte.m_programName));
  m_uncalledTsv.addHeader("Version", ADT_VERSION);
  m_uncalledTsv.addHeader("Date", Util::getTimeStamp());
  m_uncalledTsv.addHeader("TranslationFile", rte.m_adtOpts.m_inputTTableFile);
  m_uncalledTsv.addHeader("MarkerList", rte.m_adtOpts.m_inputMarkerListFile);
  m_uncalledTsv.addHeader("HaplotypeReportOption", rte.m_adtOpts.m_useFirstDupAlleleDefHeaderText);
  m_uncalledTsv.addHeader("GenotypeOverrideFile", rte.m_adtOpts.m_inputGenotypeOverrideFile);

  m_uncalledTsv.writeTsv_v1(getUncalledReportName(rte) + TEMP_REPORT_FILE_EXT);

  return;

}
// end RegressionExperimentReport::_openUncalledFile
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_report
 * Synopsis:
 *
 * Helper function to TranslationExperimentReport::generate broken out
 * for readability.
 *
 * Write or accumulate a tab delimited, output TSV file record
 * of marker data. The comprehensive report and uncalled report
 * are output in translation table marker order and therefore
 * are not further sorted.
 *
 * The summary report is sorted and therefore the records are accumulated
 * for a separate _report call after all the records have been accumulated
 * for sort. 
 *
 * @param rte - the single instance run time environment.
 * @param cRow - the current row to emit, complete will all data populated. 
 * @param sampleInfo -the same for all rows, appended to the row. 
 *
 * @return - writes the record to disk or accumulates the summary row. 
 */
/*****************************************************************************/
void TranslationExperimentReport::_report(const RunTimeEnvironment & rte, const TranslationRow & cRow, std::vector< std::string > & sampleInfo)
{

  // here we expect things to be in the order defined above.
  // otherwise we could say: m_tsv.set(0,"experiment",m_experiment);
  m_comprehensiveTsv.set(0, 0,  cRow.m_index);
  m_comprehensiveTsv.set(0, 1,  cRow.m_chpFileName);
  m_comprehensiveTsv.set(0, 2,  cRow.m_gene);
  m_comprehensiveTsv.set(0, 3,  cRow.m_knownCall);
  m_comprehensiveTsv.set(0, 4,  cRow.m_unknownCall);
  m_comprehensiveTsv.set(0, 5,  cRow.m_interpretationCode);
  m_comprehensiveTsv.set(0, 6,  cRow.m_summaryFlag);
  m_comprehensiveTsv.set(0, 7,  cRow.m_relevantAlleles);
  m_comprehensiveTsv.set(0, 8,  cRow.m_markerName);
  m_comprehensiveTsv.set(0, 9,  cRow.m_probeSet);
  m_comprehensiveTsv.set(0, 10, cRow.m_baseCall);
  m_comprehensiveTsv.set(0, 11, cRow.m_referenceReportBase);
  m_comprehensiveTsv.set(0, 12, cRow.m_variantReportBase);
  m_comprehensiveTsv.set(0, 13, cRow.m_call);
  m_comprehensiveTsv.set(0, 14, cRow.m_haplotypeMarker);
  m_comprehensiveTsv.set(0, 15, cRow.m_changeForVariant);
  m_comprehensiveTsv.set(0, 16, cRow.m_CDNAChange);
  m_comprehensiveTsv.set(0, 17, cRow.m_genomicPosition);
  m_comprehensiveTsv.set(0, 18, cRow.m_dbSNPId);
  m_comprehensiveTsv.set(0, 19, cRow.m_originalBasecall);
  m_comprehensiveTsv.set(0, 20, cRow.m_overrideComment);

  if (cRow.m_isUncalledRow) {
    m_uncalledTsv.set(0, 0, cRow.m_chpFileName);
    m_uncalledTsv.set(0, 1, cRow.m_gene);
    m_uncalledTsv.set(0, 2, cRow.m_markerName);
    m_uncalledTsv.set(0, 3, cRow.m_probeSet);
    m_uncalledTsv.set(0, 4, cRow.m_baseCall);
    m_uncalledTsv.set(0, 5, ""); // Override
    m_uncalledTsv.set(0, 6, cRow.m_referenceBase);
    m_uncalledTsv.set(0, 7, cRow.m_variantBase);

  }

  if (cRow.m_isSummaryRow || cRow.m_isGeneSummary) {
    // Save the summary rows so they can be sorted before output.
    if ( rte.m_adtOpts.m_summaryReportSort ) {
      SummaryRowIndex sri;
      sri.m_gene             = cRow.m_gene;
      sri.m_probeSet         = cRow.m_probeSet;
      sri.m_knownCall        = cRow.m_knownCall;
      sri.m_summaryFlag      = cRow.m_summaryFlag;
      sri.m_markerName       = cRow.m_markerName;
      sri.m_isSummaryRow     = cRow.m_isSummaryRow;
      sri.m_isGeneSummary    = cRow.m_isGeneSummary;
      sri.m_copyNumber       = cRow.m_copyNumber;
      m_summaryRows[sri]     = cRow;
    }
    else {
      _reportSummary( cRow, sampleInfo );
    }
  }

  for (size_t i = 0; i < sampleInfo.size() ; i++) {
    m_comprehensiveTsv.set(0, 21 + i, sampleInfo[i]);
  }

  // Comprehensive file write.
  m_comprehensiveTsv.writeLevel(0);


  // Uncalled file write.
  if (cRow.m_isUncalledRow) {
    m_uncalledTsv.writeLevel(0);
  }

  return;
}
// end TranslationExperimentReport::_report
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_reportSummary
 *
 * Synopsis:
 *
 * Helper function to TranslationExperimentReport::_reportSummarySorted
 * broken out for readability.
 
 * 
 * @param rte - the single instance run time environment.
 * @param cRow - the current row to emit, complete will all data populated. 
 * @param sampleInfo -the same for all rows, appended to the row. 
 *
 * @return - writes the summary record to disk in a TSV file. 
 *
 */
/*****************************************************************************/
void TranslationExperimentReport::_reportSummary(const TranslationRow & cRow, std::vector< std::string > & sampleInfo) {

  
  APT_ERR_ASSERT(cRow.m_isSummaryRow || cRow.m_isGeneSummary , "_reportSummary got a non-summary record.");

  if (cRow.m_isSummaryRow ) {
    
    m_summaryTsv.set(0, 0,  cRow.m_index);
    m_summaryTsv.set(0, 1,  cRow.m_chpFileName);
    m_summaryTsv.set(0, 2,  cRow.m_gene);
    m_summaryTsv.set(0, 3,  cRow.m_knownCall);
    m_summaryTsv.set(0, 4,  cRow.m_unknownCall);
    m_summaryTsv.set(0, 5,  cRow.m_interpretationCode);
    m_summaryTsv.set(0, 6,  cRow.m_summaryFlag);
    m_summaryTsv.set(0, 7,  cRow.m_relevantAlleles);
    m_summaryTsv.set(0, 8,  cRow.m_markerName);
    m_summaryTsv.set(0, 9,  cRow.m_probeSet);
    m_summaryTsv.set(0, 10, cRow.m_baseCall);
    m_summaryTsv.set(0, 11, cRow.m_referenceReportBase);
    m_summaryTsv.set(0, 12, cRow.m_variantReportBase);
    m_summaryTsv.set(0, 13, cRow.m_call);
    m_summaryTsv.set(0, 14, cRow.m_haplotypeMarker);
    m_summaryTsv.set(0, 15, cRow.m_changeForVariant);
    m_summaryTsv.set(0, 16, cRow.m_CDNAChange);
    m_summaryTsv.set(0, 17, cRow.m_genomicPosition);
    m_summaryTsv.set(0, 18, cRow.m_dbSNPId);
    m_summaryTsv.set(0, 19, cRow.m_originalBasecall);
    m_summaryTsv.set(0, 20, cRow.m_overrideComment);
  } else  {
    m_summaryTsv.set(0, 0,  cRow.m_index);
    m_summaryTsv.set(0, 1,  cRow.m_chpFileName);
    m_summaryTsv.set(0, 2,  cRow.m_gene);
    m_summaryTsv.set(0, 3,  cRow.m_knownCall);
    m_summaryTsv.set(0, 4,  cRow.m_unknownCall);
    m_summaryTsv.set(0, 5,  cRow.m_interpretationCode);
    if ( cRow.m_copyNumber == 0 ) {
      m_summaryTsv.set(0, 6, "Deletion alleles have been detected for this gene, CN=0");
    }
    else {
      m_summaryTsv.set(0, 6, "All markers responsible for functional changes are Ref/Ref");
    }
    for (int i = 7; i < 21; i++) {
      m_summaryTsv.set(0, i, "");
    }

  }

  for (size_t i = 0; i < sampleInfo.size() ; i++) {
    m_summaryTsv.set(0, 21 + i, sampleInfo[i]);
  }

  m_summaryTsv.writeLevel(0);

}
// end TranslationExperimentReport::_reportSummary
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_reportSummarySorted   
 *
 * Synopsis:
 * 
 * Helper function to TranslationExperimentReport::generate broken out
 * for readability.
 *
 * Separate from _report() because the sumary records need to be accumlated
 * and then sorted for output. This method outputs all the summary records
 * previously accumulated. 
 * 
 * @param sampleInfo - the sample info for all rows.
 *
 *
 */
/*****************************************************************************/
void TranslationExperimentReport::_reportSummarySorted(std::vector< std::string > & sampleInfo) {

  std::map< SummaryRowIndex, TranslationRow>::iterator itSRI;

  for ( itSRI = m_summaryRows.begin();
        itSRI != m_summaryRows.end();
        itSRI ++ ) {
    _reportSummary( itSRI->second, sampleInfo);
    
  }

  m_summaryRows.clear();
  

}
// end TranslationExperimentReport::_reportSummarySorted
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationExperimentReport::_setIsUncalledRow:
 * Synopsis:
 *
 * Helper function to TranslationExperimentReport::generate broken out
 * for readability.
 * 
 * Determines if a row should be reported on in the uncalled report.
 * From the SRS:
 *
 *  Conditions:
 *  basecall = NoCall, NotAvailaable or PRA 
 *
 *
 * @param cRow -  With baseCall already filled in data structure of row data.
 *
 * @return - void, cRow.m_isUncalledRow is set.
 */
/*****************************************************************************/
void TranslationExperimentReport::_setIsUncalledRow(const RunTimeEnvironment & rte, TranslationRow & cRow)
{

  cRow.m_isUncalledRow = false;

  // No Call or PRA
  if ( rte.m_adtOpts.m_uncalledReportAllMarkers || pcrecpp::RE("NoCall|NC|PRA|PossibleRareAllele|NotAvailable").FullMatch(cRow.m_baseCall)) {
    cRow.m_isUncalledRow = true;

  }

  return;

}
//end TranslationExperimentReport::_setIsUncalledRow
/*****************************************************************************/
