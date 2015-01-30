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
 * @file   MarkerExperimentReport.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  DEPRECATED DMET2: Class for the marker report as report files.
 */

//
#include "translation/MarkerExperimentReport.h"
//
#include "translation/ExperimentGeneResults.h"
#include "translation/ExperimentResults.h"
#include "translation/RunTimeEnvironment.h"
#include "translation/TranslationTableModel.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include "pcrecpp.h"

using namespace std;
using namespace affx;


/*****************************************************************************/
/**
 * MarkerRowComp::operato()
 *
 * Used in the STL std::sort for sorting MarkerRow objects.
 * Sort order is experiment -> gene -> probeSet.
 *
 *
 * @param  a,b - the two MarkerRow objects to compare.
 * @return true - if ( a < b);
 *
 */
/*****************************************************************************/
bool MarkerRowComp::operator()(const MarkerRow &a,
                               const MarkerRow &b) const
{

  if (a.m_experiment != b.m_experiment) {
    return(a.m_experiment.compare(b.m_experiment) <= 0);
  }

  if (a.m_gene != b.m_gene) {
    return(a.m_gene.compare(b.m_gene) <= 0);
  }

  return (a.m_probeSet.compare(b.m_probeSet) <= 0);

}
// end MarkerRowComp::operator
/*****************************************************************************/
/*****************************************************************************/
/**
 * MarkerExperimentReport::close:
 * Synopsis:
 *
 *   Close the Marker report file.
 *
 * @return - void, closes report files.
 */
/*****************************************************************************/
bool MarkerExperimentReport::close(bool abort)
{

  if (! abort) {
    Verbose::out(ADT_VERBOSE_NORMAL, "MarkerExperimentReport::generate called: ", false);
    Verbose::out(ADT_VERBOSE_NORMAL,  m_reportName);
  }

  if (m_open_called) {
    m_marker_tsv.close();
    m_open_called = false;
  }

  return true;

}
// end MarkerExperimentReport::close
/*****************************************************************************/
/*****************************************************************************/
/**
 * _report_marker
 * Synopsis:
 *
 * Helper function for one time use by MarkerExperimentReport::generate.
 *
 * Write a tab delimited record to the output TSV file for marker data.
 *
 * @param marker_tsv - the TsvFile reference
 * @param mRow - the MarkerRow class defined just for this report.
 *
 */
/*****************************************************************************/
static void _report_marker(affx::TsvFile & marker_tsv, const MarkerRow & mRow)
{

  // here we expect things to be in the order defined above.
  // otherwise we could say: m_tsv.set(0,"experiment",m_experiment);
  marker_tsv.set(0, 0,  mRow.m_experiment);
  marker_tsv.set(0, 1,  mRow.m_gene);
  marker_tsv.set(0, 2,  mRow.m_sample);
  marker_tsv.set(0, 3,  mRow.m_functionalChange);
  marker_tsv.set(0, 4,  mRow.m_externalId);
  marker_tsv.set(0, 5,  mRow.m_probeSet);
  marker_tsv.set(0, 6,  mRow.m_baseCall);
  marker_tsv.set(0, 7,  mRow.m_referenceBase);
  marker_tsv.set(0, 8,  mRow.m_variantBase);
  marker_tsv.set(0, 9,  mRow.m_call);
  marker_tsv.set(0, 10, mRow.m_haplotypeMarker);
  marker_tsv.set(0, 11, mRow.m_changeForVariant);
  marker_tsv.set(0, 12, mRow.m_variantCDNAChange);
  marker_tsv.set(0, 13, mRow.m_variantDNAChange);
  marker_tsv.set(0, 14, mRow.m_dbSNPId);
  marker_tsv.set(0, 15, mRow.m_validated);
  marker_tsv.set(0, 16, mRow.m_alleleDefiningMarker);
  marker_tsv.set(0, 17, mRow.m_relevantAlleles);

  marker_tsv.writeLevel(0);

  return;
}
// end _report_marker
/*****************************************************************************/
/*****************************************************************************/
/**
 * MarkerExperimentReport::generate
 *
 * Create the report this object was made for.
 * As with all other report objects the data is generated one experiment
 * at a time.
 *
 *
 * @param rte - the single instnace RunTimeEnvironment which contains options
 * @param ttm - the single instance translation table model
 * @param er  - the container for all experiment results.
 *
 * @return true - when report was generated as expected.
 *
 */
/*****************************************************************************/
bool MarkerExperimentReport::generate(class RunTimeEnvironment & rte,
                                      class TranslationTableModel & ttm,
                                      std::map<std::string, ExperimentResults *>&er)
{

  bool okReport = true;


  if (! m_open_called) {
    open(rte);
  }


  //std::sort(m_results.begin(), m_results.end());
  int count = 0;


  std::map<std::string, ExperimentResults*>::const_iterator itSER;

  // For each experiment, recapitulate the translation table.

  for (itSER = er.begin(); (itSER != er.end()) && okReport ;  itSER++) {

    // IGNORE Experiments with no marker calls for DMET2.
    if (rte.m_adtOpts.m_dmet2Calling && (itSER->second->m_markerCallCount == 0)) {
      continue;
    }

    std::vector<MarkerRow> mRows;
    std::map<std::string, int> probeSetRow;

    for (int i = 1; i < ttm.size(); i++, count++) {

      // Skip commented out rows as well as header rows.
      if (ttm.ignoreRow(i) || (!((ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)] == "Y") || (ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)] == "N")))) {
        continue;
      }


      std::string gene = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENE)];

      ExperimentGeneResults *egr = NULL;

      if (itSER->second->m_geneResults.find(gene) != itSER->second->m_geneResults.end()) {
        egr = itSER->second->m_geneResults[gene];
      }


      MarkerRow mRowNew;
      MarkerRow * mRow = &mRowNew;
      std::string probeSet = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)];

      bool isMultiAllelicRow = probeSetRow.count(probeSet) > 0;
      if (isMultiAllelicRow) {
        APT_ERR_ASSERT(probeSetRow[probeSet] < mRows.size(), "");
        mRow = &(mRows[probeSetRow[probeSet]]);
      } else {
        probeSetRow[probeSet] = mRows.size();
      }

      bool isProbeSetPresent = false;
      if (egr != NULL) {
        isProbeSetPresent
        = egr->m_probeSetToMarkerCallResults.count(probeSet) > 0;
      }

      int markerCallResultsIndex = 0;
      bool hasAlleleCalls = false;
      if (isProbeSetPresent) {
        markerCallResultsIndex = egr->m_probeSetToMarkerCallResults[probeSet];
        hasAlleleCalls =  egr->m_callResults[markerCallResultsIndex].m_alleleCalls.size() > 0;
      }

      // EXPERIMENT
      mRow->m_experiment  = itSER->first;

      // GENE
      mRow->m_gene       = gene;

      // SAMPLE
      mRow->m_sample     = itSER->second->m_sample;

      // FUNCTIONAL CHANGE
      mRow->m_functionalChange = "";

      // EXTERNAL ID
      mRow->m_externalId = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_COMMON_NAME)];

      // PROBE SET
      mRow->m_probeSet = probeSet;

      // REFERENCE BASE
      mRow->m_referenceBase = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE)];
      // VARIANT BASE
      if (isMultiAllelicRow) {
        mRow->m_variantBase = mRow->m_variantBase + "," + ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)];
      } else {
        mRow->m_variantBase = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)];
      }

      // CALL
      mRow->m_call = "";
      if (isProbeSetPresent && (egr->m_callResults[markerCallResultsIndex].m_alleleCalls.size() > 0)) {
        APT_ERR_ASSERT(markerCallResultsIndex < egr->m_callResults.size(), "");
        APT_ERR_ASSERT(egr->m_callResults[markerCallResultsIndex].m_alleleCalls.size() > 0, "");
        mRow->m_call = egr->m_callResults[markerCallResultsIndex].m_alleleCalls[0].m_allele;
      }

      // BASE CALL

      if (!isMultiAllelicRow) {
        mRow->m_baseCall = "NotAvailable";
      }

      if (isProbeSetPresent) {
        if ((egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.size() > 0) && (egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.size() > 0) && egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0] == "PRA") {
          mRow->m_baseCall = "PossibleRareAllele";
        }

        if (mRow->m_call == "Ref/Ref") {
          mRow->m_baseCall = mRow->m_referenceBase + "/" + mRow->m_referenceBase;
        } else if (mRow->m_call == "Ref/Var") {
          std::string variance =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0];
          if (variance ==  mRow->m_referenceBase) {
            variance =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid2.m_ceSet.begin()->second.m_bases[0];
          }
          if (variance == ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)]) {
            mRow->m_baseCall = mRow->m_referenceBase + "/" + variance;
          }
        } else if (mRow->m_call == "Var/Var") {
          std::string variance1 =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0];
          if (variance1 == ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)]) {
            std::string variance2 =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid2.m_ceSet.begin()->second.m_bases[0];
            mRow->m_baseCall = variance1 + "/" + variance2 ;
          }
        } else if (egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0] == "NC") {
          mRow->m_baseCall = "NoCall";
        }
      }

      // HAPLOTYPE MARKER
      mRow->m_haplotypeMarker = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)];

      // CHANGE FOR VARIANT
      if (isMultiAllelicRow) {
        if (mRow->m_changeForVariant.empty()) {
          mRow->m_changeForVariant = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)];
        } else if (! ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)].empty() && (mRow->m_changeForVariant != ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)])) {

          mRow->m_changeForVariant = mRow->m_changeForVariant + "," + ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)];
        }

      } else {
        mRow->m_changeForVariant = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CHANGE)];
      }

      // VARIANT CDNA CHANGE
      if (isMultiAllelicRow) {
        if (mRow->m_variantCDNAChange.empty()) {
          mRow->m_variantCDNAChange = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CDNA)];
        } else if (!ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CDNA)].empty()) {
          mRow->m_variantCDNAChange = mRow->m_variantCDNAChange + "," + ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CDNA)];
        }
      } else  {
        mRow->m_variantCDNAChange = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_CDNA)];
      }

      // VARIANT DNA CHANGE
      if (isMultiAllelicRow) {
        if (mRow->m_variantDNAChange.empty()) {
          mRow->m_variantDNAChange = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENOME_POSITION)];
        } else if (!ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENOME_POSITION)].empty()) {
          mRow->m_variantDNAChange = mRow->m_variantDNAChange + "," + ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENOME_POSITION)];
        }
      } else {
        mRow->m_variantDNAChange = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENOME_POSITION)];
      }

      // DBSNP ID
      mRow->m_dbSNPId = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DBSNP)];

      // VALIDATED
      mRow->m_validated = ttm.m_rows[i][ADT_DMET2_TT_VALIDATED];
      // ALLELE DEFINING MARKER
      if (isMultiAllelicRow) {
        if (mRow->m_alleleDefiningMarker.empty()) {
          mRow->m_alleleDefiningMarker = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)];
        } else  if (!ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)].empty() && (mRow->m_alleleDefiningMarker != ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)])) {
          mRow->m_alleleDefiningMarker = mRow->m_alleleDefiningMarker + "," + ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)];
        }
      } else {
        mRow->m_alleleDefiningMarker = ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_DEFINING)];
      }

      // RELEVANT ALLELES
      int headerRow = 0;

      for (headerRow = i - 1; headerRow >= 0; headerRow--) {
        if (headerRow == 0) break;
        if (ttm.m_rows[headerRow][0] != mRow->m_gene)  {
          headerRow++;
          break;
        }
      }

      std::string relevantAlleles;
      int relevantAllelesIndex = ttm.getColumnIndex(ADT_DMET3_TT_VARIANT);
      for (int j = ttm.getColumnIndex(ADT_DMET3_TT_ALLELE_START) + 1; j < ttm.m_rows[i].size(); j++) {
        if ((mRow->m_haplotypeMarker == "N") || (ttm.m_rows[headerRow][j].c_str())[0] != '#') {
          if (ttm.m_rows[i][j].empty()) {
            if (j == ttm.getColumnIndex(ADT_DMET3_TT_ALLELE_START))  {
              relevantAlleles = ttm.m_rows[headerRow][j] + ",";
            }
          } else   {
            relevantAllelesIndex = j;
            relevantAlleles = relevantAlleles + ttm.m_rows[headerRow][j] + ",";
          }
        }
      }
      if (!relevantAlleles.empty()) {
        relevantAlleles.resize(relevantAlleles.size() - 1);
      }
      if (isMultiAllelicRow) {
        if (! relevantAlleles.empty()) {
          mRow->m_relevantAlleles = mRow->m_relevantAlleles  + "," + relevantAlleles;
        }
      } else {
        mRow->m_relevantAlleles = relevantAlleles;
      }

      if (isProbeSetPresent && !relevantAlleles.empty() && (mRow->m_haplotypeMarker == "N")  && (ttm.m_rows[i][relevantAllelesIndex] ==  ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)])) {

        std::string variance1 =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid1.m_ceSet.begin()->second.m_bases[0];
        std::string variance2 =  egr->m_callResults[markerCallResultsIndex].m_experimentChromatid2.m_ceSet.begin()->second.m_bases[0];
        if ((mRow->m_call == "Var/Var")) {

          if (mRow->m_functionalChange.empty()) {
            if (mRow->m_relevantAlleles.find(",") != std::string::npos) {
              std::string firstAllele = mRow->m_relevantAlleles.substr(0, mRow->m_relevantAlleles.find(","));
              mRow->m_functionalChange = firstAllele + "/" + firstAllele;
            } else {
              mRow->m_functionalChange = mRow->m_relevantAlleles  + "/" + mRow->m_relevantAlleles;
            }
          } else {
            mRow->m_functionalChange = mRow->m_functionalChange.substr(0, mRow->m_functionalChange.find("/"))  + "/" + relevantAlleles;
          }
          if ((variance1 != variance2) &&
              (mRow->m_relevantAlleles.find(",") != std::string::npos)) {
            mRow->m_functionalChange = mRow->m_relevantAlleles;
          }
        } else if (mRow->m_call == "Ref/Var") {
          mRow->m_functionalChange = mRow->m_relevantAlleles;
          if (mRow->m_relevantAlleles.find(",") != std::string::npos) {
            mRow->m_functionalChange =  mRow->m_relevantAlleles.substr(0 , mRow->m_relevantAlleles.find(","));
          }

        }
      }

      if (!isMultiAllelicRow) {
        mRows.push_back(*mRow);
      }


    } //for each translation table model marker

    for (size_t j = 0; j < mRows.size(); j++) {
      _report_marker(m_marker_tsv, mRows[j]);
    }

  } // for each gene, experiment result

  return okReport;

}
// end MarkerExperimentReport::generate
/*****************************************************************************/
/*****************************************************************************/
/**
 * MarkerExperimentReport::getMarkerReportName
 * Synopsis:
 *
 * Create the Marker file name from various constant and runtime
 * environment settings.
 *
 *
 * @param rte - the single instance RunTimeEnvironment which contains options
 *
 * @return fileName - the name of the marker report file to use.
 *
 */
/*****************************************************************************/
std::string MarkerExperimentReport::getMarkerReportName(const class RunTimeEnvironment & rte)
{

  std::string markerFileName;
  if (rte.m_adtOpts.m_inputGenoFile != "") {
    markerFileName = Fs::basename(rte.m_adtOpts.m_inputGenoFile);
    markerFileName = markerFileName.substr(0, markerFileName.find('.'));
  } else {
    markerFileName = name();
  }


  markerFileName = Fs::join(rte.m_adtOpts.m_outputDir,
                            markerFileName+".dmet3_"+ADT_VERSION+"."+MARKER_FILE_EXT);

  m_reportName = markerFileName;

  return markerFileName;

}
// end MarkerExperimentReport::getMarkerReportName
/*****************************************************************************/
/*****************************************************************************/
/**
 * MarkerExperimentReport::open
 * Synopsis:
 *
 * Open the TSV file for output.
 *
 *
 * @param rte - the single instance run time environment.
 */
/*****************************************************************************/
void MarkerExperimentReport::open(const class RunTimeEnvironment & rte)
{

  m_guid = affxutil::Guid::GenerateNewGuid();

  // describe our columns.
  for (int i = 0; i < NUM_MARKER_TSV_COLUMNS; i++) {
    m_marker_tsv.defineColumn(0, i, MARKER_TSV_COLUMNS[i]);
  }

  m_marker_tsv.addHeaderComment(" For research use only. Not for diagnostic purposes.");
  m_marker_tsv.addHeader("Version", ADT_VERSION);
  m_marker_tsv.addHeader("Date", Util::getTimeStamp());
  m_marker_tsv.addHeader("CopyNumberFile", rte.m_adtOpts.m_inputCopyFile);
  m_marker_tsv.addHeader("TranslatTableFile", rte.m_adtOpts.m_inputTTableFile);
  m_marker_tsv.addHeader("DMET3File", rte.m_adtOpts.m_inputGenoFile);
  m_marker_tsv.addHeader("dmet3-report-guid", m_guid);

  std::string markerName = getMarkerReportName(rte);

  m_marker_tsv.writeTsv_v1(markerName);
  m_open_called = true;
  return;

}
// end MarkerExperimentReport::open
/*****************************************************************************/

