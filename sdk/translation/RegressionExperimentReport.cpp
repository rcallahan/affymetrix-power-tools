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
 * @file   RegressionExperimentReport.cpp
 * @author Mybrid Spalding
 * @date   Wed Apr 30 11:43:47 PDT 2008
 * @brief  Take advantage of the reporting mechanism to output regression results. Command line only, not used by the Console. 
 */

//
#include "translation/RegressionExperimentReport.h"
//
#include "translation/ExperimentGeneResults.h"
#include "translation/ExperimentResults.h"
#include "translation/TranslationTableModel.h"
//
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <algorithm>

using namespace std;
using namespace affx;

class RegressionHaplotypeRowComp
{
public:
  bool operator()(const RegressionHaplotypeRow &a, const RegressionHaplotypeRow  &b) const;

} rhrcObject;


/*****************************************************************************/
/**
 * RegressionMarkerRowComp::operator()
 * Synopsis:
 * 
 * Used in the STL std::sort for sorting RegressionMarkerRow objects.
 * We are guaranteed that results are derived from ExperimentGeneResults object
 * so the first discriminating field is that ProbeSet. Just compare that.
 *
 *
 * @param  a,b - the two RegressionMarkerRow objects to compare.
 *
 * @return true - if ( a < b);
 *
 */
/*****************************************************************************/
bool RegressionMarkerRowComp::operator()(const RegressionMarkerRow &a,
    const RegressionMarkerRow &b) const
{

  return (a.m_probeSet.compare(b.m_probeSet) <= 0);

}
// end RegressionMarkerRowComp::operator
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionHaplotypeRowComp::operator()
 * Synopsis:
 *
 * Used in the STL std::sort for sorting RegressionHaplotypeRow objects.
 * We are guaranteed that results are derived from ExperimentGeneResults object
 * so the first discriminating field is that Call field. Just compare that.
 * We sort on a number found in the first call OR the std::string comparison of
 * the two call std::strings if there is no number in the first call.
 *
 *
 * @param  a,b - the two RegressionHaplotypeRow objects to compare.
 *
 * @return true - if ( a < b);
 *
 */
/*****************************************************************************/
bool RegressionHaplotypeRowComp::operator()(const RegressionHaplotypeRow &a,
    const RegressionHaplotypeRow &b) const
{


  pcrecpp::StringPiece aConsume(a.m_call);
  std::string aDigits;
  int aInt = 0;

  if (pcrecpp::RE("(\\d+)").FindAndConsume(&aConsume, &aDigits)) {
    aInt = atoi(aDigits.c_str());
  }

  pcrecpp::StringPiece bConsume(b.m_call);
  std::string bDigits;
  int bInt = 0;

  if (pcrecpp::RE("(\\d+)").FindAndConsume(&bConsume, &bDigits)) {
    bInt = atoi(bDigits.c_str());
  }

  if (aInt && bInt && (aInt != bInt)) {
    return (aInt < bInt);
  }

  bool test = a.m_call.compare(b.m_call) <= 0;
  return (test);

}
// end RegressionHaplotypeRowComp::operator
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::close
 * Synopsis:
 *
 *   Close the report files.
 *
 * @param abort - true if close is happening due to an abort event. 
 *
 * @return - void, closes files on disk.
 */
/*****************************************************************************/
bool RegressionExperimentReport::close(bool abort)
{

  if (!abort) {
    Verbose::out(ADT_VERBOSE_NORMAL, "RegressionExperimentReport::generate called:" , false);
    Verbose::out(ADT_VERBOSE_NORMAL, m_markerReportName + ", ", false);
    Verbose::out(ADT_VERBOSE_NORMAL, m_haplotypeReportName);
  }

  if (m_open_called) {
    _close_marker_file();
    _close_haplotype_file();
    m_open_called = false;
  }

  return true;

}
// end RegressionExperimentReport::close
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::generate
 * Synopsis:
 *
 * Create the report this object was made for. This should be called
 * once per experiment. 
 * 
 * 
 * @param rte  - the single instance run time environment.
 * @param ttm  - the single instance of the translation table.
 * @param er   - the container of all ExperimentGeneResults per experiment.
 * 
 * @return true - when report was generated as expected.
 *
 */
/*****************************************************************************/
bool RegressionExperimentReport::generate(RunTimeEnvironment & rte,
    TranslationTableModel & ttm,
    std::map<std::string, ExperimentResults*> & er)
{

  bool okToContinue = true;

  if (! m_open_called) {
    _open_marker_file(rte);
    _open_haplotype_file(rte);
  }

  std::map<std::string, ExperimentResults*>::iterator itSER;

  for (itSER = er.begin(); (itSER != er.end()) && okToContinue ; itSER++) {

    if (itSER->second->m_markerCallCount > 0) {

      std::map<std::string, ExperimentGeneResults*>::iterator itSEGR;

      for (itSEGR = itSER->second->m_geneResults.begin();
           itSEGR != itSER->second->m_geneResults.end(); itSEGR++) {

        okToContinue = _generateMarkerCallResults(rte, *(itSEGR->second));
        okToContinue = okToContinue &&
                       _generateHaplotypeCallResults(rte, *(itSEGR->second));
      }
    }
  }


  return okToContinue;

}
// end RegressionExperimentReport::generate
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::getHaplotypeReportName
 * Synopsis:
 *
 * Create the Haplotype regression file name from various constant and runtime
 * environment settings. Unlike the Console reports created 
 * there is now ".tmp" version created. 
 *
 *
 * @param  rte - the single instance run time environment
 *
 * @return fileName - the name of the regression haplotype file
 *
 */
/*****************************************************************************/
std::string RegressionExperimentReport::getHaplotypeReportName(const RunTimeEnvironment & rte)
{

  std::string fileName =  _getBaseReportName(rte) +
                     getDmet3HaplotypeFileNameExt();

  m_haplotypeReportName = fileName;

  return fileName;

}
// end getMarkerReportName
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::getMarkerReportName
 * Synopsis:
 *
 * Create the Marker regression file name from various constant and runtime
 * environment settings. Unlike the Console reports created 
 * there is now ".tmp" version created. 
 *
 *
 * @param rte - the sngle instance run time environment.
 *
 * @return fileName - the name of the regression marker file.
 *
 */
/*****************************************************************************/
std::string RegressionExperimentReport::getMarkerReportName(const RunTimeEnvironment & rte)
{

  std::string fileName = _getBaseReportName(rte) + getDmet3MarkerFileNameExt();

  m_markerReportName = fileName;

  return fileName;

}
// end getMarkerReportName
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::_close_marker_file
 * Synopsis:
 *
 * Close a TSV file for outputing marker regression file results.
 *
 *
 */
/*****************************************************************************/
void RegressionExperimentReport::_close_marker_file()
{

  m_marker_tsv.close();

}
// end RegressionExperimentReport::_close_marker_file
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::_close_haplotype_file
 * Synopsis:
 *
 * Close a TSV file for outputing haplotype regression file results.
 *
 *
 */
/*****************************************************************************/
void RegressionExperimentReport::_close_haplotype_file()
{

  m_haplotype_tsv.close();


}
// end RegressionExperimentReport::_close_haplotype_file
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::_generateHaplotypeCallResults
 * Synopsis:
 *
 * Helper function broken out for readability.
 *
 * Output the HAPLOTYPE regression results for the CallResults set for
 * a particular ExperimentGene. Note that for multi-allelic probeSets, or
 * markers, this means merging the CallResults from all the various
 * probeSet calls in the Translation Table for that probeSet.
 *
 *
 * @param  rte - the single instance run time environment 
 * @param  egr - the ExperimentGeneResults with per gene CallResults array
 *
 * @return true - when report was generated as expected.
 *
 */
/*****************************************************************************/
bool RegressionExperimentReport::_generateHaplotypeCallResults(RunTimeEnvironment & rte,  ExperimentGeneResults & egr)
{


  pcrecpp::RE reUnknown("UNK");
  std::vector<RegressionHaplotypeRow> rhrRows;

  for (size_t i = 0; i < egr.size(); i++) {

    if (egr.m_callResults[i].getCallType() != ADT_CALL_TYPE_HAPLOTYPE_GROUP) {
      continue;
    }

    std::string calls;
    int unknownCount = 0;
    for (int j = 0; j < egr.m_callResults[i].size() ; j++) {


      RegressionHaplotypeRow rhr;

      rhr.m_experiment  = egr.m_callResults[i].m_experimentName;
      rhr.m_gene        = egr.m_callResults[i].m_geneName;

      rhr.m_call =  egr.m_callResults[i].m_alleleCalls[j].m_allele;

      if (reUnknown.PartialMatch(rhr.m_call)) {
        unknownCount++;
      }
      rhrRows.push_back(rhr);
    }

    // Generate Summary Statistics
    std::string callCount = ToStr(egr.m_callResults[i].size());
    std::string knownCount = ToStr(egr.m_callResults[i].size() - unknownCount);

    for (size_t j = 0; j < rhrRows.size() ; j++) {
      rhrRows[j].m_callCount  = callCount;
      rhrRows[j].m_knownCount = knownCount;
    }

  } // for each CallResults, report or aggregate.


  std::sort(rhrRows.begin(), rhrRows.end(), rhrcObject);

  for (size_t i = 0; i < rhrRows.size(); i++) {
    _report_haplotype(rhrRows[i]);
  }


  return true;

}
// end RegressionExperimentReport::_generateHaplotypeCallResults
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::_generateMarkerCallResults
 * Synopsis:
 *
 * Helper function broken out for clarity.
 *
 * Output the MARKER regression results for the CallResults set for
 * a particular ExperimentGene. Note that for multi-allelic probeSets, or
 * markers, this means merging the CallResults from all the various
 * probeSet calls in the Translation Table for that probeSet.
 *
 *
 * @param  rte - the single instance run time environment. 
 * @param  egr - the ExperimentGeneResults with per gene  CallResults array.
 *
 * @returns true - when report was generated as expected.
 *
 */
/*****************************************************************************/
bool RegressionExperimentReport::_generateMarkerCallResults(RunTimeEnvironment & rte,    ExperimentGeneResults & egr)
{

  set< std::string > seenProbeSets;
  std::vector<RegressionMarkerRow> rmrRows;

  for (size_t i = 0; i < egr.size(); i++) {

    if (egr.m_callResults[i].getCallType() != ADT_CALL_TYPE_MARKER) {
      continue;
    }

    RegressionMarkerRow rmr;

    // A marker, even with wildcards, should only ever have one call
    // per input experiment set of data. Assert this.

    if (egr.m_callResults[i].m_alleleCalls.size() > 1) {
      for (size_t i = 0; i < egr.m_callResults[i].m_alleleCalls.size() ; i++) {
        cerr << egr.m_callResults[i].m_alleleCalls[i].m_allele << endl;
      }
    }
    APT_ERR_ASSERT(egr.m_callResults[i].m_alleleCalls.size() <= 1, "");

    bool possibleRareAllele = false;
    bool emptyResults = (egr.m_callResults[i].m_alleleCalls.size()  == 0);

    if (emptyResults) {


      possibleRareAllele = (((*egr.m_callResults[i].m_experimentChromatid1.m_ceSet.begin()).second.m_bases[0] == "PRA") || ((*egr.m_callResults[i].m_experimentChromatid2.m_ceSet.begin()).second.m_bases[0] == "PRA"));


    }

    std::string probeSet = egr.m_callResults[i].m_experimentChromatid1.getMarkerProbeSet();


    set<std::string>::iterator iSADit = seenProbeSets.find(probeSet);

    if (iSADit != seenProbeSets.end()) {
      continue;
    }


    rmr.m_experiment = egr.m_callResults[i].m_experimentName;
    rmr.m_gene       = egr.m_callResults[i].m_geneName;
    rmr.m_probeSet    = probeSet;

    if (possibleRareAllele) {

      rmr.m_a1 = "PossibleRareAllele";

    } else if (emptyResults) {
      rmr.m_a1 = "NoCall";

    } else {
      rmr.m_a1 =  "A1|" + (*egr.m_callResults[i].m_experimentChromatid1.m_ceSet.begin()).second.m_bases[0];
      rmr.m_a2 = "A2|" + (*egr.m_callResults[i].m_experimentChromatid2.m_ceSet.begin()).second.m_bases[0];

    }

    APT_ERR_ASSERT(egr.m_probeSetReference[ probeSet ].size() == 1, "");
    APT_ERR_ASSERT(egr.m_probeSetVariants[ probeSet ].size() > 0, "");

    rmr.m_ref = "R|" + (*egr.m_probeSetReference[ probeSet ].begin());



    std::string variantList = egr.m_probeSetVariants[probeSet][0];

    if (egr.m_probeSetVariants[probeSet].size() > 1) {

      std::sort(egr.m_probeSetVariants[probeSet].begin(), egr.m_probeSetVariants[probeSet].end());

      variantList = egr.m_probeSetVariants[probeSet][0];
      for (size_t j = 1;  j < egr.m_probeSetVariants[probeSet].size();  j++) {
        variantList = variantList + "," + egr.m_probeSetVariants[probeSet][j];
      }
    }

    rmr.m_var =  "V|" + variantList;


    if (! possibleRareAllele && !emptyResults) {
      rmr.m_call =  egr.m_callResults[i].m_alleleCalls[0].m_allele;
    }

    rmrRows.push_back(rmr);

    seenProbeSets.insert(probeSet);

  } // for each CallResults, report or aggregate.

  std::sort(rmrRows.begin(), rmrRows.end(), RegressionMarkerRowComp());

  for (size_t i = 0; i < rmrRows.size(); i++) {
    _report_marker(rmrRows[i]);
  }


  return true;

}
// end RegressionExperimentReport::_generateMarkerCallResults
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::_getBaseReportName
 * Synopsis:
 *
 * Common function for generation of all report names.
 *
 * Create the common base name for all regression files.
 *
 * @param rte - the single instance run time enviroment
 *
 * @return baseName - the common base report name for regression. 
 *
 */
/*****************************************************************************/
std::string RegressionExperimentReport::_getBaseReportName(const RunTimeEnvironment & rte) const
{
  std::string baseName = Fs::join(rte.m_adtOpts.m_outputDir,
                                  Fs::basename(rte.m_adtOpts.m_inputGenoFile)+".");
  // FsPath p1(rte.m_adtOpts.m_inputGenoFile);
  // p1.dump();
  // std::string baseName = FsPath::join(rte.m_adtOpts.m_outputDir,p1.getFileNameWoExt()+".");
  return baseName;
}

// end RegressionExperimentReport::_getBaseReportName
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::_open_marker_file
 * Synopsis:
 *
 * Open a TSV file for outputing marker regression file results, includes
 * adding any necessary headers. 
 * 
 * @param rte - the single instance run time environment 
 */
/*****************************************************************************/
void RegressionExperimentReport::_open_marker_file(const RunTimeEnvironment & rte)
{

  // describe our columns.
  for (int i = 0; i < REGRESSION_NUM_MARKER_TSV_COLUMNS; i++) {
    m_marker_tsv.defineColumn(0, i, REGRESSION_MARKER_TSV_COLUMNS[i]);
  }

  // No headers, not used during regression comparing.

  m_marker_tsv.writeTsv_v1(getMarkerReportName(rte));
  m_open_called = true;
  return;

}
// end RegressionExperimentReport::_open_marker_file
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::_open_haplotype_file
 * Synopsis:
 *
 * Open a TSV file for outputing haplotype regression file results, includes
 * any required headers. 
 *
 * @param rte - the single instance run time environment. 
 *
 */
/*****************************************************************************/
void RegressionExperimentReport::_open_haplotype_file(const RunTimeEnvironment & rte)
{

  for (int i = 0; i < REGRESSION_NUM_HAPLOTYPE_TSV_COLUMNS; i++) {
    m_haplotype_tsv.defineColumn(0, i,  REGRESSION_HAPLOTYPE_TSV_COLUMNS[i]);
  }

  // No headers, not used during regression comparing.
  m_haplotype_tsv.writeTsv_v1(getHaplotypeReportName(rte));


  m_open_called = true;
  return;

}
// end RegressionExperimentReport::_open_haplotype_file
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::_report_marker
 * Synopsis:
 *
 * Write a tab delimited marker record to the output TSV file for marker data.
 *
 * @param rmr - the row as RegressionMarkerRow class defined just for this report.
 *
 */
/*****************************************************************************/
void RegressionExperimentReport::_report_marker(const RegressionMarkerRow & rmr)
{

  // here we expect things to be in the order defined above.
  // otherwise we could say: m_tsv.set(0,"experiment",m_experiment);
  m_marker_tsv.set(0, 0, rmr.m_experiment);
  m_marker_tsv.set(0, 1, rmr.m_gene);
  m_marker_tsv.set(0, 2, rmr.m_probeSet);
  m_marker_tsv.set(0, 3, rmr.m_a1);
  m_marker_tsv.set(0, 4, rmr.m_a2);
  m_marker_tsv.set(0, 5, rmr.m_ref);
  m_marker_tsv.set(0, 6, rmr.m_var);
  m_marker_tsv.set(0, 7, rmr.m_call);

  m_marker_tsv.writeLevel(0);

  return;
}
// end RegressionExperimentReport::_report_marker
/*****************************************************************************/
/*****************************************************************************/
/**
 * RegressionExperimentReport::_report_haplotype
 * Synopsis:
 *
 * Write a tab delimited haplotype record to the output TSV file for
 * haplotype data.
 *
 * @param rhr - the row as RegressionHaplotypeRow class defined just for this report.
 *
 */
/*****************************************************************************/
void RegressionExperimentReport::_report_haplotype(const RegressionHaplotypeRow & rhr)
{

  // here we expect things to be in the order defined above.
  // otherwise we could say: m_tsv.set(0,"experiment",m_experiment);
  m_haplotype_tsv.set(0, 0, rhr.m_experiment);
  m_haplotype_tsv.set(0, 1, rhr.m_gene);
  m_haplotype_tsv.set(0, 2, rhr.m_call);
  m_haplotype_tsv.set(0, 3, rhr.m_callCount);
  m_haplotype_tsv.set(0, 4, rhr.m_knownCount);

  m_haplotype_tsv.writeLevel(0);

  return;

}
// end RegressionExperimentReport::_report_haplotype
/*****************************************************************************/
