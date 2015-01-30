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
 * @file   TranslationTable.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  The primary business object class that manipulates the various input data models,
 *         creating other business objects appropriate for the translation algorithm and invokes translation. 
 */


#include "translation/TranslationTable.h"
//
#include "translation/CopyNumberTableModel.h"
#include "translation/ExperimentGeneResults.h"
#include "translation/GeneCall.h"
#include "translation/GenotypeTableModel.h"
#include "translation/TranslationTableModel.h"
//
#include "util/Err.h" // includes Verbose.h
//

using namespace std;

/*****************************************************************************/
/**
 * _errorMessageDuplicateProbeSets
 * Synopsis:
 * Helper error message convienance function to make the code easier to read.
 *
 * @param rte - RunTimeEnviornment
 * @param ttm - TranslationTableModel.
 * @param row - the row in question, 0 based.
 *
 * @return false - indicates it is not ok to perform analysis.
 */
/*****************************************************************************/
static bool _errorMessageDuplicateProbeSets(const RunTimeEnvironment & rte,
    TranslationTableModel & ttm,
    const int row)
{


  bool okToContinue = false;
  std::stringstream msgSStr;
  msgSStr << rte.m_adtOpts.m_inputGenoFile << ": invalid ProbeSet with Allele [A1-A29] or Reference column for row: " << endl;
  std::string rowAsString = ttm.getRowAsString(row);
  msgSStr << "gene row [" << row << "]: " << rowAsString << endl;
  Verbose::warn(ADT_VERBOSE_EXCEPTION, msgSStr.str());

  return okToContinue;

}
// end TranslationTable::_errorMessageDuplicateProbeSets
/*****************************************************************************/
/*****************************************************************************/
/**
 * _errorMessageGeneRowOrder
 * Synopsis:
 * Helper error message convienance function to make the code easier to read.
 *
 * @param rte - RunTimeEnviornment
 * @param ttm - TranslationTableModel.
 * @param row - the row in question, 0 based.
 *
 * @return false - indicates it is not ok to perform analysis.
 */
/*****************************************************************************/
static bool _errorMessageGeneRowOrder(const RunTimeEnvironment & rte,
    TranslationTableModel & ttm,
    const int row,
    const std::string & geneInProcess)
{

  bool okToContinue = false;
  std::stringstream msgSStr;
  msgSStr << rte.m_adtOpts.m_inputGenoFile << ": invalid format, gene row is not under the correct gene header: " << endl;
  std::string rowAsString = ttm.getRowAsString(row);
  msgSStr << "gene row [" << row << "]: " << rowAsString << endl;
  msgSStr << "gene header : "  << geneInProcess ;
  Verbose::out(ADT_VERBOSE_INPUT_FILES,  msgSStr.str());

  return okToContinue;

}
// end _errorMessageGeneRowOrder
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTable::TranslationTable
 * Synopsis:
 * Main constructor.
 * Truth be told, this object should really be refactored to
 * just "Translation" and not "TranslationTable". There was a time where
 * object was closely coupled with the translation table but no more.
 *
 * This object is primarily comprised of a hash of gene call business objects
 * derived from the translation table. 
 * The key is the gene name. The value is a std::vector of GenCalls corresponding
 * to the "Reference", "Variant", and "A1" ... "AN" columns of the
 * translation table. Non-haplotype markers will only have GeneCalls
 * that correspond to the "Reference" and "Variant" columns. 
 * There will be a GeneCall per each non-haplotype marker plus one for the
 * single possible Haplotype set of markers.
 *
 * During translation, each GeneCall will be attempted to match
 * the CHP experiment data and if a match is find then the underlying
 * CallSet contained within the GeneCall will be attached to
 * the CallResults object. 
 *
 * 1.) GeneCall with CallSets for making translation calls. 
 * 2.) Zero Copy Number CallSets per gene. 
 *
 * @param rte - the single instance run time environment
 * @param ttm - the single instance TranslationTableModel 
 * @param cntm - CopyNumberTableModel which carries the override alleles.
 */
/*****************************************************************************/
TranslationTable::TranslationTable(const RunTimeEnvironment & rte,
                                   TranslationTableModel & ttm,
                                   const CopyNumberTableModel & cntm)
{


  std::string    geneInProcess;
  std::string    priorGene;
  int       geneHeaderRow           = 0;
  GeneCall *nextMarkerGeneCall      = NULL;
  GeneCall *nextHaplotypeGeneCall   = NULL;
  bool      okToContinue            = true;
  bool      rowIsHaplotype          = false;
  int       haplotypeCount          = 0;

  //Profile profile1;
  int numRows = ttm.size();

  for (int row = 0; (row < numRows) && okToContinue ;  row++) {

    std::string rowAsString = ttm.getRowAsString(row);

    if (ttm.ignoreRow(row)) {
      Verbose::out(ADT_VERBOSE_INPUT_FILES, "IGNORING: " + ToStr(row + 2) + ":" + rowAsString);
      continue;
    } else {
      Verbose::out(ADT_VERBOSE_INPUT_FILES, "PROCESSING: " + ToStr(row + 2) + ":" + rowAsString);
    }

    std::string gene = ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_GENE)];
    std::string probeSet = ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)];


    if (gene != geneInProcess)  {
      // Starting a new GENE.

      if (ttm.validateRowIsHeader(row)) {

        m_geneGroupCompleteCallSet[gene] = CallSet(gene + ":CompleteSet");
        priorGene = geneInProcess;
        geneInProcess = gene;
        rowIsHaplotype = false;
        geneHeaderRow = row;
        haplotypeCount = 0;

        if (nextHaplotypeGeneCall != NULL) {
          nextHaplotypeGeneCall =
            _addGeneCallToGroup(rte, priorGene, nextHaplotypeGeneCall);

        }

        // Copy Number is hard code for DMET2
        // Copy Number is data driven for DMET3
        if (rte.m_adtOpts.m_dmet2Calling) {
          if (isGeneCopyNumberSensitive(gene)) {
            std::map<std::string, std::string>::const_iterator itSS;

            for (itSS =  cntm.m_experimentCopyNumberCall.begin();
                 itSS !=  cntm.m_experimentCopyNumberCall.end(); itSS++) {
              std::string geneExperiment = gene + ":" + itSS->first;
              m_geneExperimentCopyNumberCall[geneExperiment] = itSS->second;
            }
          }
        }
      } else {
        okToContinue
        = _errorMessageGeneRowOrder(rte, ttm, row, geneInProcess) ;
      }

      continue;

    } // if starting new gene and encounted the gene header row

    rowIsHaplotype = ttm.isRowHaplotype(row);

    // We currently only support one Haplotype per gene.
    // The haplotype rows may be out of order so just continue to append
    // the current haplotype gene call.
    if (rowIsHaplotype && ! nextHaplotypeGeneCall) {
      nextHaplotypeGeneCall = new GeneCall(rte, ttm, geneHeaderRow, rowIsHaplotype);
    }

    // For multi-allelic probeSets then the nextMarkerGeneCall is the
    // previous version.
    bool okToAddMarkerGeneCall = true;
    if ((nextMarkerGeneCall = _findMarkerGeneCall(geneInProcess,  probeSet))  != NULL) {

      okToAddMarkerGeneCall = false;

    } else {
      nextMarkerGeneCall = new GeneCall(rte, ttm, geneHeaderRow, false);
    }

    bool okCallElement = true;
    // Add a marker GeneCall for both Haplotype and Marker.
    int allMarkerAlleleIndex[] = { ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE), ttm.getColumnIndex(ADT_DMET3_TT_VARIANT) };

    // Make sure the reference and variant are not the same.
    std::string ref = ttm.m_rows[row][allMarkerAlleleIndex[0]];
    std::string var = ttm.m_rows[row][allMarkerAlleleIndex[1]];

    if (ref == var) {
      APT_ERR_ABORT(gene + " " + probeSet + ": translation file reference and variant base are identical.");
    }
    for (int i = 0; (i < 2) && okCallElement; i++) {

      okCallElement = okCallElement && nextMarkerGeneCall->
                      addCallElementToAlleleCallSet(ttm, geneHeaderRow, row, allMarkerAlleleIndex[i], true, false);

      m_geneGroupCompleteCallSet[gene].
      addCallElement(ttm, row, allMarkerAlleleIndex[i], true);

    }

    if (rowIsHaplotype) {
      // Add to Allele call sets
      std::vector<int> alleleColumns = ttm.getHaplotypeAlleleColumns(geneHeaderRow);
      for (int j = 0; okCallElement && (j < alleleColumns.size()); j++) {
        APT_ERR_ASSERT(nextHaplotypeGeneCall, "");

        if (ttm.m_rows[row][alleleColumns[j]].empty()) {
          continue;
        }

        okCallElement = okCallElement &&
                        nextHaplotypeGeneCall->
                        addCallElementToAlleleCallSet(ttm, geneHeaderRow, row, alleleColumns[j], true, true);
        m_geneGroupCompleteCallSet[gene].
        addCallElement(ttm, row, alleleColumns[j], true);
      }
    } else if ( ttm.getType() == ADT_TRANSLATION_TABLE_TYPE_DMET3 ){
      // For non-haplotype markers the copy number base needs to be added.
      int copyNumberColumn =  ttm.getCopyNumberColumn(row);
      if (copyNumberColumn > -1) {

        okCallElement = okCallElement && nextMarkerGeneCall->
                        addCallElementToAlleleCallSet(ttm, geneHeaderRow, row, copyNumberColumn, true, false);

        m_geneGroupCompleteCallSet[gene].
        addCallElement(ttm, row, copyNumberColumn, true);

      }
    }

    // Add a single marker call
    if (okToAddMarkerGeneCall && okCallElement) {
      nextMarkerGeneCall = _addGeneCallToGroup(rte, gene, nextMarkerGeneCall);
    } else if (!okCallElement) {
      okToContinue = _errorMessageDuplicateProbeSets(rte, ttm, row);
    }

  } // for each row

  if (okToContinue && (nextHaplotypeGeneCall != NULL)) {
    nextHaplotypeGeneCall =
      _addGeneCallToGroup(rte, geneInProcess, nextHaplotypeGeneCall);

  }
  std::map<std::string, CallSet>::iterator itSCS;
  itSCS = m_geneGroupCompleteCallSet.begin();

  for ( ; itSCS != m_geneGroupCompleteCallSet.end(); itSCS++ ) {
    std::map<std::string, CallElement >::iterator itCE;
    for (  itCE = itSCS->second.m_ceSet.begin();
           itCE != itSCS->second.m_ceSet.end();   itCE++ ) {
      itCE->second.setCompleteSetSize();
    }
  }

  if (ttm.getType() == ADT_TRANSLATION_TABLE_TYPE_DMET3) {
    okToContinue = okToContinue && _duplicateCallSetCheck(rte);
    okToContinue = okToContinue && _reconcileCopyNumberCallSet(rte, ttm);
  }

  if (!okToContinue) {
    std::stringstream msgSStr;
    msgSStr << rte.m_adtOpts.m_inputGenoFile ;
    msgSStr << ": business logic errors detected.";
    APT_ERR_ABORT(msgSStr.str());
  }


}
// end TranslationTable::TranslationTable
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTable::getGeneCallGroup
 * Synopsis:
 *
 * A selector for the m_geneGroup.
 * 
 * Typically a GeneCallGroup will be made up of SNPs (single
 * markers) and one possilbe Haplotype CallSet. In any case, the
 * entire set is returned with this method. 
 *
 * @param geneName - TranslationTableModel gene.
 * @returns std::vector<GeneCall> - the group of GeneCalls found. a
 */
/*****************************************************************************/
std::vector<GeneCall> TranslationTable::getGeneCallGroup(const std::string & geneName)
{
  return m_geneGroup[geneName];

}
// end TranslationTable::getGeneCallGroup
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTable::hasGeneCall
 * Synopsis:
 * 
 * Interrogate the std::vector of GeneCalls if a particular gene name
 * is available for making a call. If not then there is no experiment
 * analysis to be done.
 *
 * @return true - if the passed in gene name has any corresponding GeneCall.
 */
/*****************************************************************************/
bool TranslationTable::hasGeneCall(const std::string & geneName)
{

  return (m_geneGroup[geneName].size() != 0);

}
// end TranslationTable::hasGeneCall
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTable::translateExperimentGene
 * Snyopsis:
 * Perform translation and invoke the translation algorithm found in the
 * GeneCall business object.
 * Marshalls all the requisite data for making each translation
 * call and then again packages up the results in the CallResults
 * business object whenever a call is made.
 * 
 * To re-state this routine is mostly
 * just a wrapper around the GeneCall method of similar name,
 * translateExperimentCall. The difference being is that this
 * method iterates over all the GeneCall's in a group and summarizes
 * the results in the return ExperimentGeneResults.
 *
 * So the way  this works is that each GeneCall acts as a filter on the
 * input genotype experiment data. Genotype data is grouped into genes within
 * the GenotypeTableModel being passed in. Each GeneCall corresponds to
 * precisely one gene marker or probe set, or the gene's
 * single Haplotype marker
 * set. The input genotype data (CHP) represented by the
 * GenotypeTableModel is searched for the probe sets required for
 * GeneCall. The two chromatid CallSets are built from the
 * GenotypTableModel for the GeneCall in question and then the
 * translation algorithm is invoked to find a match.
 *
 * If the input genotype data (CHP) does not contain the probe sets
 * required for the GeneCall in question then the GeneCall is skipped.
 * There is an exception with respect to the Haplotype set of markers.
 * If the program option is given to enforce that entire Haplotype
 * sets be complete for translation and markers are found missing
 * within the input genotype data (CHP), then an exception will be thrown.
 *
 * In the common case, an incomplete haplotype group results in the
 * haplotype GeneCall being skipped, but realize that each individual
 * marker has an individual marker GeneCall and so will still have its
 * marker "Ref/Var" call made. 
 * 
 * 
 * 
 * @param rte - the single instance run time environment
 * @param gtm - the experiment Genotype Tabel model filtered for a single gene. 
 * @param ttm - the single instance translation table model. 
 * @param cntm - DMET2 ONLY: the copy number table model
 * @param egr - the experiment gene results for return
 *
 * @return - egr, experiment gene results is filled in with the calls for the entire experiment gene. 

 */
/*****************************************************************************/
void TranslationTable::translateExperimentGene(RunTimeEnvironment & rte,  GenotypeTableModel  & gtm, const TranslationTableModel &ttm, const CopyNumberTableModel & cntm, ExperimentGeneResults *egr)
{

  // Copy the static stuff from the model.

  APT_ERR_ASSERT(egr, "");

  if (! hasGeneCall(gtm.m_geneName)) {

    // This data should be passed through for reporting purposes.
    Verbose::out(ADT_VERBOSE_TMI, "Experiment Gene with no call, ignoring: " + gtm.m_experimentName + ", " + gtm.m_geneName);
    return;
  }

  // Each CallSet will have a pair of CallElements corresponding to a
  // ALLELE1 and ALLELE2 for the row.
  if (rte.m_adtOpts.m_profile) {
    rte.m_profiles["call_loop1"]->begin();
  }

  std::vector<CallSet> invalidGTMRows = _invalidateExperimentGeneProbeSetBaseValues(rte, gtm, ttm);

  if (invalidGTMRows.size() > 0) {
    Verbose::warn(ADT_VERBOSE_NORMAL, gtm.m_geneName + ": Invalid gene markers not found in translation file");
    for (int i = 0; i < invalidGTMRows.size(); i++) {
      invalidGTMRows[i].describeVerbose(ADT_VERBOSE_NORMAL);
    }
    APT_ERR_ABORT("Invalid experiment gene marker values detected.");
  }
  if (rte.m_adtOpts.m_profile) {
    rte.m_profiles["call_loop1"]->end();
  }

  std::string geneExperiment = gtm.getGeneName() + ":" + gtm.m_experimentName;
  std::string dmet2CopyNumberCall;

  if (rte.m_adtOpts.m_dmet2Calling) {
    std::map<std::string, std::string>::const_iterator itSS
    = m_geneExperimentCopyNumberCall.find(geneExperiment);

    if (itSS != m_geneExperimentCopyNumberCall.end()) {
      dmet2CopyNumberCall = itSS->second;
    }
  }

  // This is a std::vector because a GeneCall represents one Haplotype grouping
  // call and many marker calls.

  std::vector<GeneCall> geneCallToMatch = getGeneCallGroup(gtm.getGeneName());

  std::vector<GeneCall>::iterator matchCall;

  for (matchCall = geneCallToMatch.begin();
       matchCall != geneCallToMatch.end(); matchCall++) {

    Verbose::out(ADT_VERBOSE_CALL, "CALL: " + gtm.m_geneName);

    //matchCall->describeVerbose( rte, ADT_VERBOSE_CALL );

    // Ok, the gtm will contain all the
    //  markers needed for all GeneCalls
    // This constructer filters the gtm with only elements found in the
    // referenceSet
    if (rte.m_adtOpts.m_profile) {
      rte.m_profiles["call_loop2"]->begin();
    }
    CallSet chromatid1(rte, gtm, matchCall->getReferenceCallSet(), 1);

    CallSet chromatid2(rte, gtm, matchCall->getReferenceCallSet(), 2);

    if (rte.m_adtOpts.m_profile) {
      rte.m_profiles["call_loop2"]->end();
    }

    if (chromatid1.isEmpty()) {
      // This experiment didn't contain any data for this call,
      // this is expected, just continue.
      continue;
    }

    if (rte.m_adtOpts.m_profile) {
      rte.m_profiles["call_loop3"]->begin();
    }
    CallResults call =
      matchCall->
      translateExperimentCall(rte,
                              m_geneGroupCompleteCallSet[gtm.m_geneName],
                              gtm,
                              chromatid1, chromatid2,
                              dmet2CopyNumberCall,
                              m_geneCopyNumberZeroCallSet[gtm.m_geneName]);

    if (rte.m_adtOpts.m_profile) {
      rte.m_profiles["call_loop3"]->end();
    }

    // Set the TranslatedExperimentResults data
    egr->appendCallResults(call, *matchCall);

  } // for each experiment gene call


  return;
}
// end   TranslationTable::translateExperimentGene
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTable::_addGeneCallToGroup
 * Synopsis:
 * Helper function for the constructor, TranslationTable::TranslationTable
 * 
 * Add a completed GeneCall to a gene's call Group.
 *
 * @param rte - RunTimeEnvironment
 * @param gene - std::string name of the Gene from the translation table.
 * @param addGeneCall - pointer to the GeneCall to add
 * @param markerProbeSet - 
 *
 * @returns NULL, to reset a the pointer if needed.
 */
/*****************************************************************************/
GeneCall * TranslationTable::_addGeneCallToGroup(
  const RunTimeEnvironment & rte,
  const std::string & gene,
  GeneCall * addGeneCall)
{

  APT_ERR_ASSERT(addGeneCall, "");

  if (addGeneCall->getReferenceCallSet().getCallType() == ADT_CALL_TYPE_HAPLOTYPE_GROUP) {
    _completeAlleleSets(addGeneCall);
  } else {
    // Erase the allele name CallSets for MARKERS.
    // The first two CallSets are always Reference and Variant.
    // The copy number 0 set will be reconciled in another API.
    for (int i = addGeneCall->m_alleleSet.size() - 1; i >= 2; i--) {
      if (addGeneCall->m_alleleSet[i].m_copyNumber != 0) {
        addGeneCall->m_alleleSet.erase(addGeneCall->m_alleleSet.begin() + i);
      }
    }
  }

  m_geneGroup[gene].push_back(*addGeneCall);

  //m_geneGroup[gene][m_geneGroup[gene].size() -1].describeVerbose(rte, ADT_VERBOSE_NORMAL);

  delete addGeneCall;

  return  NULL;

}
// end TranslationTable::_addGeneCallToGroup
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTable::_completeAlleleSets
 * Synopsis:
 * Helper function to TranslationTable::_addGeneCallToGroup
 * 
 * The allele bases  in the input translation table represent
 * changes with respect to the reference set. Fill in the blanks where the
 * reference set ProbeSet bases are added when not found in the allele set.
 * For SNP markers this means filling in the reference value for the
 * first column which matches the reference.
 *
 * @param  gc - pointer to the GeneCall class in question.
 *
 */
/*****************************************************************************/
void TranslationTable::_completeAlleleSets(GeneCall * gc)
{

  APT_ERR_ASSERT(gc, "");

  CallSet & referenceSet = gc->getReferenceCallSet();

  std::map<std::string, CallElement>::iterator iCEit;


  for (int i = 0; i < gc->m_alleleSet.size(); i++) {

    CallSet & alleleSet = gc->getAlleleCallSet(i);

    alleleSet.m_type = ADT_CALL_TYPE_HAPLOTYPE_GROUP;
    alleleSet.m_isDescriptive = true;

    // For each CallElement in the reference set,
    // search for the ProbeSet in the allele set and if
    // it is not there, then add the reference CallElement to the allele set.
    for (iCEit = referenceSet.m_ceSet.begin();
         iCEit != referenceSet.m_ceSet.end();
         iCEit++) {

      std::map<std::string, CallElement>::iterator jCEit;
      bool okToAdd = true;

      for (jCEit = alleleSet.m_ceSet.begin();
           (jCEit != alleleSet.m_ceSet.end()) &&
           okToAdd;
           jCEit++) {

        if (jCEit->first == iCEit->first)
          okToAdd = false;


      } // for each CallElement in the alleleSet, find the corresponding
      // CallElement probeSet in the referenceSet.

      if (okToAdd) {
        alleleSet.m_ceSet[iCEit->first] = iCEit->second;
      }

    } // for each CallElement in the referenceSet


  } // for each allele call, complete the allele set


  return;

}
// end TranslationTable::_completeAlleleSets
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTable::_duplicateCallSetCheck
 * Synopsis:
 * Helper function to the constructor, TranslationTable::TranslationTable
 *
 * The input translation table is created by humans and is
 * not guaranteed to have a unique set of
 * defining markers per haplotype alllele. The algorithm 
 * behavior during translation is to match the first
 * allele set of defining markers and ignore the others.
 * This API is used to make the determinate if any duplicates
 * exist and the calllee can use this to abort the program or
 * not depending on the business requirements. 
 * 
 *
 * @param rte  - the single run time environment instance
 * 
 * @returns true - the duplicates do not exist 
 */
/*****************************************************************************/
bool TranslationTable::_duplicateCallSetCheck(const RunTimeEnvironment & rte)
{

  std::map<std::string, std::vector<GeneCall> >::iterator itMSVG;

  ADT_VERBOSE_ENUM level = ADT_VERBOSE_NORMAL;

  bool okDuplicateCheck = true;

  // Foreach gene
  for (itMSVG = m_geneGroup.begin(); itMSVG != m_geneGroup.end(); itMSVG++) {

    std::string gene = itMSVG->first;

    // Foreach GeneCall allele call set
    for (int i = 0; i < itMSVG->second.size(); i++) {

      // Foreach call set
      for (int j = 1; j < itMSVG->second[i].m_alleleSet.size(); j++) {

        // Check subsequent call sets for duplicates.

        for (int k = j + 1; k < itMSVG->second[i].m_alleleSet.size(); k++) {

          if (!rte.m_adtOpts.m_audit.empty() && itMSVG->second[i].m_alleleSet[j] == itMSVG->second[i].m_alleleSet[k]) {
            Verbose::warn(level, "Invalid translation table, duplicate allele set of variants detected.");
            itMSVG->second[i].m_alleleSet[j].describeVerbose(level);
            itMSVG->second[i].m_alleleSet[k].describeVerbose(level);
            //Verbose::warn(level, "Translation file, duplicate allele set of variants detected.");
            okDuplicateCheck = false;;
          }
        }
      }
    }
  }


  return okDuplicateCheck;
}
// end TranslationTable::_duplicateCallSetCheck
/*****************************************************************************/
/**
 * TranslationTable::_findMarkerGeneCall
 * Synopsis:
 * Helper function to the constructor, TranslationTable::TranslationTable
 *
 * When constructing marker GeneCalls, each marker gets a new GeneCall
 * except for multiallelic markers. These markers will already have GeneCalls
 * in the m_geneGroup when the marker occurs for more than once in the
 * input translation table. 
 * This method finds the appropriate GeneCall for the multiallelic  marker
 * that already has a GeneCall. 
 * 
 *
 * @param gene - the multiallelic probeSet
 * @param probeSet - the marker
 * 
 * @return - The corresponding GeneCall for the inputs.
 */
/*****************************************************************************/
GeneCall * TranslationTable::_findMarkerGeneCall(std::string gene, std::string probeSet)
{

  if (! m_geneGroupCompleteCallSet[gene].hasProbeSet(probeSet)) {
    return NULL;
  }

  for (int i = 0; i < m_geneGroup[gene].size(); i++) {

    for (int j = 0; j < m_geneGroup[gene][i].m_alleleSet.size(); j++) {
      if (m_geneGroup[gene][i].m_alleleSet[j].hasProbeSet(probeSet) &&
          m_geneGroup[gene][i].m_alleleSet[j].getCallType() == ADT_CALL_TYPE_MARKER) {
        return &(m_geneGroup[gene][i]);
      }
    }
  }
  // Programming error. The complete set cannot have an probeSet that
  // is not in a GeneCall.

  APT_ERR_ASSERT(false, "");

  return NULL;

}
// end TranslationTable::_findMarkerGeneCall
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTable::_invalidateExperimentGeneProbeSetBaseValues
 * Synopsis:
 * Helper function to TranslationTable::translateExperimentGene
 * 
 * This is just and optimiaztion fucntion where the *MAJORITY* 
 * of inbound GenotypeTableModel probe sets will have ProbeSet/Base pairs which
 * are not recognized for any GeneCall. This helper method modifies
 * the passed in GenotypeTableModel by marking rows in the
 * "m_rows" attribute that are unknown with a special experiment name.
 * The invalid rows are also returned
 * in the form of a CallSet where each ALLELE1 and ALLELE2
 * CallElement corresponds to one CallSet and one GTM row.
 *
 * @param rte - the single instance RunTimeEnvironment
 * @param gtm - a GenotypeTableModel for a particulare gene experiment.
 * @param ttm - the single instance translation table model
 *
 * @return invalidRows - a std::vector<CallSet> where each CallSet represents a
 *                       row where ALLELE1 or ALLELE2 had an invalid base.
 */
/*****************************************************************************/
std::vector<CallSet> TranslationTable::_invalidateExperimentGeneProbeSetBaseValues(
  const RunTimeEnvironment & rte, GenotypeTableModel & gtm, const TranslationTableModel & ttm)
{


  std::vector<CallSet> invalidGTMRows;

  // Unknown gene names should be handled before calling this API.

  APT_ERR_ASSERT(m_geneGroupCompleteCallSet.find(gtm.getGeneName()) != m_geneGroupCompleteCallSet.end(), "");


  const CallSet & testCS = m_geneGroupCompleteCallSet[gtm.getGeneName()];

  for (int i = 0; i < gtm.getGeneRowSize(gtm.getGeneName()); i++) {

    int row = gtm.getGeneRow(gtm.getGeneName(), i);

    std::string probeSet = gtm.m_rows[row][GT_PROBE_SET_INDEX];
    std::string geneName = gtm.m_rows[row][GT_GENES_INDEX];
    std::string base1    = gtm.m_rows[row][GT_ALLELE1_INDEX];
    std::string base2    = gtm.m_rows[row][GT_ALLELE2_INDEX];

    if (base2.empty()) {
      base2 = base1;
    }

    CallElement allele1(probeSet, base1);
    CallElement allele2(probeSet, base2);

    if (testCS.hasProbeSet(probeSet)) {

      // Non-haplotype markers are always valid.
      if (ttm.isHaplotypeMarker(probeSet)) {


        // Wildcard markers are always valid.
        if ((! allele1.hasWildCards() && !testCS.contains(allele1)) ||
            (! allele2.hasWildCards() && !testCS.contains(allele2))) {

          // Just ignore the row by setting the Experiment Name to something diffferent.
          gtm.m_rows[row][GT_EXPERIMENT_INDEX] =
            "Invalid ALLELE Base For ProbeSet";
          //Flag an exception by returning an invalid CallSEt.
          CallSet invalidCS(geneName);
          invalidCS.m_ceSet[probeSet] =  allele1 ;
          invalidGTMRows.push_back(invalidCS);
        }
      }
    } else {

      APT_ERR_ASSERT(ttm.m_probeSetRowIndex.count(probeSet) == 0, "");
      // Just ignore the row by setting the Experiment Name to something diffferent.
      gtm.m_rows[row][GT_EXPERIMENT_INDEX] =
        "Invalid ProbeSet for Gene";
    }

  }

  return invalidGTMRows;

}
// end TranslationTable::_invalidateExperimentGeneProbeSetBaseValues
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTable::_reconcileCopyNumberCallSet:
 * Synopsis:
 * Helper function to the constructor, TranslationTable::TranslationTable.
 *
 * The format of the translation table was hacked to use an allele column
 * (A1-AN) to carry the copy number zero designation. This routines
 * reconciles the hack and creates the appropriate business object
 * and then marks "allele" column as one to be ignored. 
 *
 * How to know if an allele name is copy number zero?
 * There is no naming conention, instead the first row after the allele
 * name column is checked for a "0" value. If the row immediately below
 * the allele name contains a "0" then the allele is the single
 * copy number zero designation. 
 *
 * 1.) It is an exception to have more than one copy number 0 call set.
 * 2.) Remove the copy number 0 allele call set from the GeneCall set.
 * 3.) Instantiate the copy number 0 member CallSet set for the GeneCall.
 *
 * @params rte - the RunTimeEnvironment
 @ @params ttm - the translation table model for updates.
 *
 * @return - false on error, the parent is responsible for throwing an exception.
 *           ttm is updated to comment out the copy number allele
 *           designations.
 */
/*****************************************************************************/
bool TranslationTable::_reconcileCopyNumberCallSet(const RunTimeEnvironment & rte, TranslationTableModel & ttm)
{

  std::map<std::string, std::vector<GeneCall> >::iterator itMSVG;

  // Foreach gene
  for (itMSVG = m_geneGroup.begin(); itMSVG != m_geneGroup.end(); itMSVG++) {

    std::string gene = itMSVG->first;

    // Foreach GeneCall allele call set
    for (int i = 0; i < itMSVG->second.size(); i++) {


      // Skip over reference and variant sets, start at 2.
      for (int j = 2; j < itMSVG->second[i].m_alleleSet.size(); j++) {

        if (itMSVG->second[i].m_alleleSet[j].m_copyNumber == 0) {

          // Multiple check
          if (m_geneCopyNumberZeroCallSet[gene].getCallType() != ADT_CALL_TYPE_NULL) {
            if (m_geneCopyNumberZeroCallSet[gene].m_name == itMSVG->second[i].m_alleleSet[j].m_name) {
              continue;
            }
            Verbose::warn(ADT_VERBOSE_EXCEPTION, "Fatal error,  translation file has more than one copy number 0 designation for gene: " + itMSVG->first + ": " +  m_geneCopyNumberZeroCallSet[gene].m_name + ", "  + itMSVG->second[i].m_alleleSet[j].m_name);
            return false;
          }
          ttm.m_geneCopyNumberIndicator[gene] = true;
          m_geneCopyNumberZeroCallSet[gene] = itMSVG->second[i].m_alleleSet[j];
          m_geneCopyNumberZeroCallSet[gene].m_type = ADT_CALL_TYPE_COPY_NUMBER;
          m_geneCopyNumberZeroCallSet[gene].size() > 1 ?
          m_geneCopyNumberZeroCallSet[gene].m_type = ADT_CALL_TYPE_HAPLOTYPE_GROUP :
              m_geneCopyNumberZeroCallSet[gene].m_type = ADT_CALL_TYPE_MARKER ;
          itMSVG->second[i].m_alleleSet.erase(itMSVG->second[i].m_alleleSet.begin() + j);
          //m_geneCopyNumberZeroCallSet[gene].describeVerbose( ADT_VERBOSE_NORMAL );
          // Delete the copy number 0 so that reporting doesn't
          // use it in relevant alleles.
          std::map<std::string, CallElement >::iterator itMSCE = m_geneCopyNumberZeroCallSet[gene].m_ceSet.begin();
          std::string probeSet = itMSCE->second.m_probeSet;
          int probeSetRow = ttm.getProbeSetRowIndex(probeSet);
          int headerRow = ttm.getHeaderRow(probeSetRow);

          for (int k = ttm.getColumnIndex(ADT_DMET3_TT_ALLELE_START); i < ttm.columnCount(); k++) {

            if (ttm.m_rows[headerRow][k] == m_geneCopyNumberZeroCallSet[gene].m_name) {
              //ttm.m_rows[probeSetRow][k] = "#" + m_geneCopyNumberZeroCallSet[gene].m_name;
              ttm.m_rows[probeSetRow][k].clear();
              ttm.m_rows[headerRow][k].clear();
              break;
            }
          }
        } else if (itMSVG->second[i].m_alleleSet[j].getCallType() == ADT_CALL_TYPE_MARKER) {
          APT_ERR_ASSERT(itMSVG->second[i].m_alleleSet.size() == 2, "");
        }
      }
    }

  }

  return true;
}
// TranslationTable::_reconcileCopyNumberCallSet
/*****************************************************************************/
