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
 * @file   TranslationAudit.cpp
 * @author Mybrid Spalding
 * @date   Thu Sep 25 13:56:24 PDT 2008
 * @brief  Business object that contains the annotation file and translation file audit requirements.
 */


#include "translation/TranslationAudit.h"
//
#include "translation/AnnotationTableModel.h"
#include "translation/RunTimeEnvironment.h"
#include "translation/TranslationTableModel.h"
//
#include "util/Util.h"
//

using namespace std;

// TranslationAuditType
// An audit type is similar to a regression type where
// new classes of audits are completely independent of other
// audit classes. To implement a new audit type create the
// the code and register the top level method here. 
const TranslationAuditType TRANSLATE_AUDIT_TYPES[] = {
  // command line, callback pointer
  { std::string("annotation"), &TranslationAudit::auditAnnotation },
};

#define TRANSLATE_AUDIT_TYPES_SIZE() ( (int) sizeof(TRANSLATE_AUDIT_TYPES[0]) / sizeof( TRANSLATE_AUDIT_TYPES) )

/*****************************************************************************/
/**
 * TranslationAudit::TranslationAudit
 * Synopsis:
 *
 * Main constructor. Does not start the audit logic. To start the audit
 * see TranslationAudit::audit
 *
 *
 * @param rte - the single instance RunTimeEnvironment 
 */
/*****************************************************************************/
TranslationAudit::TranslationAudit(const class RunTimeEnvironment &rte) :
    m_numInvalidGeneProbeSets(0),  m_numInvalidSwitchDesignStrands(0)
{

  APT_ERR_ASSERT(! rte.m_adtOpts.m_audit.empty(), "");

  bool optionFound = false;

  int size = TRANSLATE_AUDIT_TYPES_SIZE();

  for (int i = 0 ; (i < size) && !optionFound; i++) {
    if (rte.m_adtOpts.m_audit == TRANSLATE_AUDIT_TYPES[i].m_commandLine) {
      optionFound = true;
    }

  }

  if (! optionFound) {
    APT_ERR_ABORT("invalid option: " + rte.m_adtOpts.m_audit + " is not a recognized audit option value.");
  }

}
// end TranslationAudit::TranslationAudit
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationAudit::audit
 * Synopsis:
 *
 * The "main" program if you will. Wrapper for the various call backs
 * registered in the TranslationAuditType TRANSLATE_AUDIT_TYPES at the
 * top of this file.
 *
 * At the command line the audit option is passed an audit type, i.e.
 * "annotation". This type is then looked up in the type table and if
 * found the audit is run. Only one audit can be run at a time.
 * Invoke the command line multiple times with the various names
 * to run various audits. 
 *
 * An audit type is similar to a regression type or unit test scenario where
 * new classes of audits are completely independent of other.
 * To implement a new audit type create the
 * the code as methods here and then register the top level method in
 * TranslationAuditType TRANSLATE_AUDIT_TYPES. Pass in the type name
 * at the command line "-a mynewtype" to run the audit. 
 * 
 * @param rte - the single instance RunTimeEnvironment
 * @param TRANSLATE_AUDIT_TYPES - with the registered call backs to invoke. 
 *
 * @return - exit value 0 if ok, -1 otherwise
 */
/*****************************************************************************/
int TranslationAudit::audit(class RunTimeEnvironment & rte)
{

  APT_ERR_ASSERT(! rte.m_adtOpts.m_audit.empty(), "");

  int size = TRANSLATE_AUDIT_TYPES_SIZE();
  
  for (int i = 0; i < size; i++) {
    int (TranslationAudit::*callback)(RunTimeEnvironment & rte) =
      TRANSLATE_AUDIT_TYPES[i].m_auditCallback;
    if ( TRANSLATE_AUDIT_TYPES[i].m_commandLine == rte.m_adtOpts.m_audit ) {
      return (this->*callback)(rte);
    }
  }


  APT_ERR_ABORT("Programmer error: audit called with " + rte.m_adtOpts.m_audit + ", an unknown audit option type.");

  return -1;

}
// TranslationAudit::audit
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationAudit::auditAnnotation
 * Synopsis:
 *
 * Audit the translation file. This is a work in progress as of
 * 11/04/2008 and needs completion.
 *
 *
 * @param rte - the single instance RunTimeEnvironment 
 *
 * @return - exit value, 0 when all is ok, -1 otherwise.
 */
/*****************************************************************************/
int TranslationAudit::auditAnnotation(class RunTimeEnvironment & rte)
{

  int exitValue = 0;

  Verbose::out(ADT_VERBOSE_NORMAL, "AUDIT: annotation file.");

  TranslationTableModel auditTT(rte, rte.m_adtOpts.m_inputTTableFile, rte.m_adtOpts.m_inputTTableType);

  m_msgSStr <<  "Input translation file records read:      " << auditTT.size() ;

  Verbose::out(rte.m_currentVerbosity, m_msgSStr.str());
  m_msgSStr.str("");

  AnnotationTableModel auditAT(rte, rte.m_adtOpts.m_inputAnnotationFile);

  m_msgSStr <<  "Input annotation table records read:       " << auditAT.size() ;

  Verbose::out(rte.m_currentVerbosity, m_msgSStr.str());
  m_msgSStr.str("");


  std::string gene;

  for (int row = 0; row < auditTT.size(); row++) {

    bool isHeaderRow = auditTT.validateRowIsHeader(row, false);

    if (isHeaderRow) {
      gene = auditTT.m_rows[row][auditTT.getColumnIndex(ADT_DMET3_TT_GENE)];
      continue;
    }

    // PROBE SET
    std::string probeSet = auditTT.m_rows[row][auditTT.getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)];
    if (auditAT.m_probeSetIndex.count(probeSet) == 0) {
      m_missingProbeSets.push_back(probeSet);
      continue;
    }

    for (int i = 0; i < auditAT.m_probeSetIndex[probeSet].size(); i++) {

      int annotationRow = auditAT.m_probeSetIndex[probeSet][i];

      // GENE
      std::string annotationGene = auditAT.m_rows[annotationRow][AnnotationTableModel::ASSOCIATED_GENE];

      if (annotationGene != gene) {
        exitValue = -1;
        m_numInvalidGeneProbeSets++;
        m_invalidGeneProbeSetsSStr << "translation file row " << row << ": " << probeSet << ", " <<  gene << " != annotation file row " << annotationRow << ": " << probeSet << ", " << annotationGene << endl;

      }
      // SWITCH STRAND
      std::string annotationSwitchStrand = auditAT.m_rows[annotationRow][AnnotationTableModel::SWITCH_DESIGN_STRAND];
      std::string translationSwitchStrand = auditTT.m_rows[row][auditTT.getColumnIndex(ADT_DMET3_TT_SWITCH_STRAND)];

      if ((annotationSwitchStrand == "1") &&
          !((translationSwitchStrand == "1") || (translationSwitchStrand == "Y"))) {
        m_numInvalidSwitchDesignStrands++;
        m_invalidSwitchDesignStrandSStr << "translation  file row " << row << ": " << probeSet << ", " <<  "switch design strand " << translationSwitchStrand << " != annotation table row " << annotationRow << ": " << probeSet << ", " << "switch design strand " << annotationSwitchStrand << endl;
      }

    } // for each annotation row per probe set

    // Reference and Variant values
    exitValue = _auditAnnotationReferenceVariant(rte) || exitValue;
  } // for each translation table row


  exitValue = _auditAnnotationReport(rte) || exitValue;

  return exitValue;

}
// end TranslationAudit::auditAnnotation
/*****************************************************************************/

/*****************************************************************************/
/**
 * TranslationAudit::_auditAnnotationReferenceVariant
 * Synopsis:
 *
 *  A helper function to auditAnnotation broken out strictly for
 *  readability purposes.
 *  Audits that the "Allele Reported Strand" values like "A // G" reflect
 *  the translation table where the first allele is the reference
 *
 *
 * @param rte - the RunTimeEnvironment
 *
 * @return - 0 if ok, -1 otherwise
 */
/*****************************************************************************/
int TranslationAudit::_auditAnnotationReferenceVariant(class RunTimeEnvironment & rte)
{


  return 0;

}
// end TranslationAudit::_auditAnnotationReferenceVariant
/*****************************************************************************/

/*****************************************************************************/
/**
 * TranslationAudit::_auditAnnotationReport
 * Synopsis:
 *
 *  A helper function to auditAnnotation broken out strictly for
 *  readability purposes.
 *
 *  Writes the audit report to log and console using Verbose.
 *
 * @param rte - the RunTimeEnvironment
 *
 * @return - 0 if ok, -1 otherwise.
 */
/*****************************************************************************/
int TranslationAudit::_auditAnnotationReport(class RunTimeEnvironment & rte)
{

  int exitValue = 0;

  // PROBE SET
  if (m_missingProbeSets.size() == 0) {
    Verbose::out(ADT_VERBOSE_NORMAL, "Probe Set Ids: AUDIT OK (all translation file probe sets accounted for in the annotation file).");
  } else {
    exitValue = -1;
    Verbose::warn(ADT_VERBOSE_NORMAL, "Probe Set Ids: AUDIT FAILED");
    m_msgSStr << "Total missing probe sets: " << m_missingProbeSets.size() << endl;
    m_msgSStr << "The following translation file probe sets were missing:";
    for (int i = 0; i < m_missingProbeSets.size(); i++) {
      m_msgSStr << m_missingProbeSets[i] << " ";
    }
    m_msgSStr << endl;
    Verbose::out(ADT_VERBOSE_NORMAL, m_msgSStr.str());
    m_msgSStr.str("");
  }

  // GENE

  if (m_numInvalidGeneProbeSets == 0) {
    Verbose::out(ADT_VERBOSE_NORMAL, "Genes: AUDIT OK (all probe sets have identical gene designations).");
  } else {
    Verbose::warn(ADT_VERBOSE_NORMAL, "Genes: AUDIT FAILED");
    Verbose::out(ADT_VERBOSE_NORMAL, m_invalidGeneProbeSetsSStr.str());
  }

  // DESIGN STRAND
  if (m_numInvalidSwitchDesignStrands == 0) {
    Verbose::out(ADT_VERBOSE_NORMAL, "Switch Design Strand: AUDIT OK (all probe sets have matching switch design strand indicators).");
  } else {
    Verbose::warn(ADT_VERBOSE_NORMAL, "Switch Design Strand: AUDIT FAILED");
    Verbose::out(ADT_VERBOSE_NORMAL, m_invalidSwitchDesignStrandSStr.str());
  }
  rte.m_programTime.end();

  Verbose::out(ADT_VERBOSE_NORMAL, "Elasped time: " +
               rte.m_programTime.getElapsedFormatedString() + ".");

  std::string timeStamp = Util::getTimeStamp();

  Verbose::out(ADT_VERBOSE_NORMAL, "[" + rte.m_programName  + " " + timeStamp + "]  END ");

  return exitValue;

}
// end TranslationAudit::_auditAnnotationReport
/*****************************************************************************/
