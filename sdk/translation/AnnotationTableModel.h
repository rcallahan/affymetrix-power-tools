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
 * @file   AnnotationTableModel.h
 * @author Mybrid Spalding <
 * @date   Thu Oct 30 15:36:28 2008
 *
 * @brief  Audit input tranlation and annoatation files for errors.
 *
 *
 */


#ifndef TRANSLATION_ANNOTATION_TABLE_MODEL_H
#define TRANSLATION_ANNOTATION_TABLE_MODEL_H

#include "translation/TranslationInputTsvTableModel.h"

class AnnotationTableModel : public TranslationInputTsvTableModel
{

public:

  enum ANNOTATION_FIELDS {
    PROBE_SET_ID,
    ASSOCIATED_GENE,
    DBSNP,
    COMMON_NAME,
    CHROMOSOME,
    PHYSICAL_POSITION,
    GENOME_CONTEXT,
    DESIGN_STRAND,
    ALLELES_DESIGN_STRAND,
    ALLELE_CODE,
    GENE_STRAND,
    SWITCH_DESIGN_STRAND,
    ALLELES_REPORTED_STRAND,
    REPORTED_STRAND,
    TYPE,
    VALIDATED,
    DBSNP_ANNOT,
    PHARMGKB_CODE,
    ALLELES_ALIASES,
    CONTEXT_CODE,
    CONTEXT_SEQUENCE,
  };

  std::map< std::string, std::vector< int > > m_probeSetIndex;


  AnnotationTableModel(const RunTimeEnvironment &rte,
                       const std::string & annotationFileName);


  ~AnnotationTableModel() {};


};

#endif /* TRANSLATION_ANNOTATION_TABLE_MODEL_H */
