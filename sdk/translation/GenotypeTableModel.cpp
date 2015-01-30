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
 * @file   GenotypeTableModel.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 10:15:23 PDT 2008
 * @brief  Genotype Data (CHP) class
 */


#include "translation/GenotypeTableModel.h"
//
#include "translation/CallElement.h"
#include "translation/TranslationTableModel.h"
//
#include "util/Err.h"  // includes Verbose.h
//
#include <sstream>
//

using namespace affx;
using namespace std;

//Sample Name Experiment Name Gene External Id ProbeSet Id Allele 1 Allele 2
//001-000505 96dev_27_001-000505 ABCB1 MDR1_A-44G 367417 A G
const TittmColumnDefinition GTM_COLUMN_DEFINITIONS[] = { // column_name, empty_ok, index, valid_regex
  { std::string("Sample Name"),     0,  0, std::string(".*"), NULL },
  { std::string("Experiment Name"), 0,  1, std::string(".*"), NULL },
  { std::string("Gene"),            0,  2, std::string("^[A-Z\\d]+$"), NULL},
  { std::string("External Id"),     0,  3, std::string(".*"), NULL  },
  { std::string("Assay Id"),        0,  4, std::string("^(?:[_\\-\\w\\d]*\\d+)$"), NULL },
  { std::string("Allele 1"),        0,  5,
    std::string("^(?:[ACGT0]+)|(?:INS|Ins|DEL|Del|NC|PRA|NoCall|PossibleRareAllele|\\-)$"), NULL },
  { std::string("Allele 2"),        1,  6,
    std::string("^(?:[ACGT0]+)|(?:INS|Ins|DEL|Del|NC|PRA|NoCall|PossibleRareAllele|\\-)$"), NULL },
};

/*****************************************************************************/
/**
 * GenotypeTableModel::GenotypeTableModel
 * Synopsis:
 *
 * DMET2: TsvFile constructor.
 *
 * This is a single instance class where only one of these should ever
 * be in existence:
 *
 * genotype data file -> genotype table model.
 *
 * A quick glance at the column headers defined for this table shows
 * that the schema for this table looks nothing like what's in a CHP
 * file header.
 *
 * A decision was made early on that the only data to be used for DMET3
 * translation with the CHP file would be a subset of the DMET2 data.
 * Therefore the DMET2 table format was left in place.
 *
 * This table is the only object to stream data. This model then is the
 * closest to being a data access object. The transfer object is simply the
 * rows (m_rows) and the "readNextExperiment" accessor is the data
 * access method. This is unlike other "data access objects" using the
 * "TranslationInputTsvTableModel" where all data is read into memory
 * at construction.
 *
 * The business object that consumes this data access / transfer object
 * is the "TranslationTable" business object. The "TranslationTable" object
 * creates chromatid call sets from the genotype data transfered.
 *
 * When this object is constructed, no data is read in. To access
 * data use the "readNextExperiment" accessor and the next
 * experiment's worth of data will be read in, albiet from
 * a DMET2 genotype short report file or a DMET3 CHP file.
 *
 *
 *
 * @param rte - the single instance RunTimeEnvironment which contains options
 * @param genoFileName - the input file to read in.
 * @param gotm - the single instance of the genotype override file table model
 */
/*****************************************************************************/
GenotypeTableModel::GenotypeTableModel(const RunTimeEnvironment &rte,  const std::string &genoFileName, GenotypeOverrideTableModel *gotm) :
    TranslationInputStreamTableModel(rte, genoFileName, GTM_COLUMN_DEFINITIONS,
                                     (sizeof(GTM_COLUMN_DEFINITIONS) / sizeof(GTM_COLUMN_DEFINITIONS[0])), true, GT_EXPERIMENT_INDEX, GT_GENES_INDEX, GT_PROBE_SET_INDEX, gotm)
{


  return;

}
// end GenotypeTableModel::GenotypeTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * GenotypeTableModel::GenotypeTableModel
 * Synopsis
 *
 * DMET3 CHP constructor to coalesce the GenotypeTableModel with the
 * DMET2 files. See the TSV constructor for a complete discussion.
 *
 * @param rte - RunTimeEnvironment which contains options
 * @param gotm - single instance of the genotype override table model
 *
 */
/*****************************************************************************/
GenotypeTableModel::GenotypeTableModel(const RunTimeEnvironment &rte,  GenotypeOverrideTableModel *gotm) :
    TranslationInputStreamTableModel(rte, GTM_COLUMN_DEFINITIONS,
                                     (sizeof(GTM_COLUMN_DEFINITIONS) / sizeof(GTM_COLUMN_DEFINITIONS[0])), gotm)
{

  return;

}
// end GenotypeTableModel::GenotypeTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * GenotypeTableModel::describeVerboseExperimentGene
 * Synopsis:
 *
 * Run time debug routine that can be invoked with various verbose levels.
 * ADT_VERBOSE_INPUT_FILES is the default level.
 *
 * @param rte       - run time environment
 * @param overrideLevel  - ADT_VERBOSE_INPUT_FILES by default,
 *                    otherwise the ADT_VERBOSE_ENUM level to output
 * @returns - void
 */
/*****************************************************************************/
void GenotypeTableModel::describeVerboseExperimentGene(const RunTimeEnvironment & rte, ADT_VERBOSE_ENUM overrideLevel)
{

  if (!overrideLevel && (rte.m_currentVerbosity < ADT_VERBOSE_INPUT_FILES)) {
    return;
  }

  std::stringstream msgSStr;


  ADT_VERBOSE_ENUM outLevel = overrideLevel ? overrideLevel : rte.m_currentVerbosity;

  if (m_geneName.empty()) {

    msgSStr <<  "describeVerboseExperimentGene: no experiment gene data has been read yet." << endl;

  }

  else {

    msgSStr << "describeVerboseExperimentGene: " << m_experimentName;
    msgSStr << " " << m_geneName << endl;

    for (int row = 0; row < m_rows.size(); row++) {
      msgSStr << getRowAsString(row) << endl;

    }
  }

  Verbose::out(outLevel, msgSStr.str());

  return;
}

//end GenotypeTableModel::describeVerobseExperimentGene
/*****************************************************************************/
/*****************************************************************************/
/**
 * getGeneRowSize:
 * Synopsis:
 *
 *   Returns the number of rows in the experiment m_rows table for one gene.
 * Rows are guaranteed to be in gene order.
 *
 * @params gene - the gene std::string to test.
 *
 * @return - The number of gene related rows, 0 if the gene doesn't exist.
 */
/*****************************************************************************/
int GenotypeTableModel::getGeneRowSize(const std::string & gene)
{

  if (m_geneRowIndex.count(gene) == 0) {
    return 0;
  }

  return m_geneRowIndex[gene].size();

}
// end GenotypeTableModel::getGeneRowSize
/*****************************************************************************/
/*****************************************************************************/
/**
 * GenotypeTableModel::getGeneRow:
 * Synopsis:
 *
 *   Returns index std::map from gene row to experiment row index.
 *
 *
 * @return - the row integer, -1 if the gene or row do not exist.
 */
/*****************************************************************************/
int GenotypeTableModel::getGeneRow(const std::string & gene, int row)
{

  if ((m_geneRowIndex.count(gene) == 0) ||
      (m_geneRowIndex[gene].size() <= row)) {
    return -1;
  }

  return m_geneRowIndex[gene][row];

}
// end GenotypeTableModel::getGeneRow
/*****************************************************************************/
/*****************************************************************************/
/**
 * GenotypeTableModel::_validateGenotypes:
 * Synopsis:
 *
 * Helper function used one time by GenotypeTableModel::getNextExperiment.
 * Broken out for readability.
 *
 *   Validate that the experiment genotype corresponds to either the reference
 *  or one of the variant bases in the translation table.
 *
 *
 *
 *
 *
 * @param rte         - the single instance run time environment
 * @param ttm         - the single instance translation table model
 * @param probeSetId  - the probeSetId in question
 * @param base1       - allele
 * @param base2       - allele
 */
/*****************************************************************************/
static void _validateGenotypes(RunTimeEnvironment & rte,
                               const TranslationTableModel & ttm,
                               const std::string & probeSetId,
                               const std::string & base1,
                               const std::string & base2)
{

  std::map< std::string, std::vector< int >  >::const_iterator itSVI;

  itSVI = ttm.m_probeSetRowIndex.find(probeSetId);

  // It's given that experiments will contain probes not in the
  // translation table. Those are perfectly valid and dealt with
  // during translation.
  if (itSVI == ttm.m_probeSetRowIndex.end()) {
    return;
  }

  bool base1OK = CallElement::isWildCardBase(base1);
  bool base2OK = CallElement::isWildCardBase(base2);

  for (int i = 0; (i < itSVI->second.size()) && (!base1OK || !base2OK); i++) {

    int row = itSVI->second[i];

    if (!base1OK) {

      if ((ttm.getType() == ADT_TRANSLATION_TABLE_TYPE_DMET3) && (base1 == "0")) {
        base1OK = true;
      } else if (ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE)] == base1) {
        base1OK = true;
      } else if (ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)] == base1) {
        base1OK = true;
      }
    }
    if (!base2OK) {
      if ((ttm.getType() == ADT_TRANSLATION_TABLE_TYPE_DMET3) && (base2 == "")) {
        base2OK = true;
      }
      if (ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE)] == base2) {
        base2OK = true;
      } else if (ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)] == base2) {
        base2OK = true;
      }
    }
  }

  if (!base1OK) {
    std::stringstream msgSStr;
    msgSStr << probeSetId << " has unknown chromatid 1 base " << base1 << " not found in translation file.";
    if (rte.m_adtOpts.m_ignoreUnknownAlleles) {
      Verbose::warn(ADT_VERBOSE_EXCEPTION, msgSStr.str());
    } else {
      APT_ERR_ABORT(msgSStr.str());
    }
  } else if (!base2OK) {
    std::stringstream msgSStr;
    msgSStr << probeSetId << " has unknown chromatid 2 base " << base2 << " not found in translation file.";
    if (rte.m_adtOpts.m_ignoreUnknownAlleles) {
      Verbose::warn(ADT_VERBOSE_EXCEPTION, msgSStr.str());
    } else {
      APT_ERR_ABORT(msgSStr.str());
    }
  }

  return;
}
// end _validateGenotypes
/*****************************************************************************/
/*****************************************************************************/
/**
 * GenotypeTableModel::getNextExperiment:
 * Synopsis:
 *
 *  This API is simply a wrapper for
 *  TranslationInputStreamTableModel->readNextExperiment.
 *
 *  Records read in are then scanned again for validation and reformating.
 *
 *
 * @param rte - the single instance RunTimeEnvironment
 * @param ttm - the single instance translation table model
 * @param geneExperimentCopyNumberCall - the copy "allele" designation per gene
 *
 * @return - the number of records read.
 */
/*****************************************************************************/
int GenotypeTableModel::getNextExperiment(RunTimeEnvironment & rte,
    TranslationTableModel & ttm,
    const std::map<std::string, std::string> &geneExperimentCopyNumberCall)
{


  if (rte.m_adtOpts.m_profile) {
    rte.m_profiles["readNextExperiment"]->begin();
  }

  m_geneCopyNumber.clear();
  
  size_t numRecords = readNextExperiment(rte, ttm, geneExperimentCopyNumberCall, m_geneCopyNumber);

  if (rte.m_adtOpts.m_profile) {
    rte.m_profiles["readNextExperiment"]->end();
  }

  m_geneRowIndex.clear();
  m_experimentGenes.clear();

  if (numRecords  > 0) {

    for (int row = 0; row < m_rows.size(); row++) {

      if (row == 0) {
        m_experimentName = m_rows[row][GT_EXPERIMENT_INDEX];
        m_geneName       = m_rows[row][GT_GENES_INDEX];
      }
      m_geneRowIndex[m_rows[row][GT_GENES_INDEX]].push_back(row);

      if (m_experimentGenes.count(m_rows[row][GT_GENES_INDEX]) == 0) {
        m_experimentGenes.insert(m_rows[row][GT_GENES_INDEX]);
      }

      if (rte.m_adtOpts.m_dmet2Calling) {
        if (m_rows[row][GT_ALLELE1_INDEX] == "Ins") {
          m_rows[row][GT_ALLELE1_INDEX] = "INS";
        }
        if (m_rows[row][GT_ALLELE2_INDEX] == "Ins") {
          m_rows[row][GT_ALLELE2_INDEX] = "INS";
        }
        if (m_rows[row][GT_ALLELE1_INDEX] == "Del") {
          m_rows[row][GT_ALLELE1_INDEX] = "DEL";
        }
        if (m_rows[row][GT_ALLELE2_INDEX] == "Del") {
          m_rows[row][GT_ALLELE2_INDEX] = "DEL";
        }
      }

      _validateGenotypes(rte, ttm, m_rows[row][GT_PROBE_SET_INDEX], m_rows[row][GT_ALLELE1_INDEX], m_rows[row][GT_ALLELE2_INDEX]);


    }
  } else {
    m_geneName = "";
    m_experimentName = "";
  }

  return(numRecords);

}
// end GenotypeTableModel::getNextExperiment
/*****************************************************************************/
