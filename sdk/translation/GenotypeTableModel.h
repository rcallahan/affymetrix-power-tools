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
 * @file   GenotypeTableModel.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 10:15:23 PDT 2008
 * @brief  Genotype Data (CHP) class
 */

#ifndef TRANSLATION_GENOTYPE_TABLE_MODEL_H
#define TRANSLATION_GENOTYPE_TABLE_MODEL_H

#include "translation/GenotypeOverrideTableModel.h"
#include "translation/TranslationInputStreamTableModel.h"
//
#include <set>


extern const TittmColumnDefinition GTM_COLUMN_DEFINTIONS[];

class GenotypeTableModel : public TranslationInputStreamTableModel
{

public:

  // empty until the first experiment has been read in
  std::set< std::string >             m_experiments;
  std::set< std::string >             m_experimentGenes;
  std::string                         m_experimentName;
  std::map< std::string, int>         m_geneCopyNumber;
  std::string                         m_geneName;
  std::map< std::string , std::vector< int > >
  m_geneRowIndex;
  
  GenotypeTableModel(const RunTimeEnvironment & rte,
                     const std::string & genoFileName,
                     GenotypeOverrideTableModel *gotm = NULL);

  GenotypeTableModel(const RunTimeEnvironment & rte,
                     GenotypeOverrideTableModel *gotm = NULL);


  GenotypeTableModel() {};
  ~GenotypeTableModel() {};


  void                     describeVerboseExperimentGene(
    const RunTimeEnvironment & rte,
    ADT_VERBOSE_ENUM overrideLevel
    = ADT_VERBOSE_NULL);

  std::string              getGeneName() {
    return m_geneName;
  };

  int                      getGeneRow(const std::string & gene, int row) ;
  int                      getGeneRowSize(const std::string & gene) ;

  int                      getNextExperiment(RunTimeEnvironment & rte,
      TranslationTableModel & ttm,
      const std::map<std::string, std::string> &geneExperimentCopyNumberCall);


};

#endif /* TRANSLATION_GENOTYPE_TABLE_MODEL_H */
