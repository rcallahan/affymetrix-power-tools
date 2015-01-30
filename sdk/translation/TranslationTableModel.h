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
 * @file   TranslationTableModel.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  Single instance data model class that wraps the translation file in the TranslationInputTsvModel.
 */


#ifndef TRANSLATION_TRANSLATION_TABLE_MODEL_H
#define TRANSLATION_TRANSLATION_TABLE_MODEL_H

#include "translation/MarkerListModel.h"
#include "translation/TranslationInputTsvTableModel.h"
//
#include "pcrecpp.h"
//
#include <list>
#include <map>
#include <set>
#include <vector>
//

extern const TittmColumnDefinition TTM_DMET2_COLUMN_DEFINTIONS[];
extern const TittmColumnDefinition TTM_DMET3_COLUMN_DEFINTIONS[];

enum ADT_DMET3_TT_COLUMNS_ENUM {
  ADT_DMET3_TT_GENE,
  ADT_DMET3_TT_REFERENCE_LINK,
  ADT_DMET3_TT_PROBE_SET_ID,
  ADT_DMET3_TT_SWITCH_STRAND,
  ADT_DMET3_TT_DBSNP,
  ADT_DMET3_TT_DEFINING,
  ADT_DMET3_TT_CDNA,
  ADT_DMET3_TT_GENOME_POSITION,
  ADT_DMET3_TT_CHANGE,
  ADT_DMET3_TT_COMMON_NAME,
  ADT_DMET3_TT_HAPLOTYPE,
  ADT_DMET3_TT_REFERENCE,
  ADT_DMET3_TT_VARIANT,
  ADT_DMET3_TT_ALLELE_START,
};

enum ADT_DMET2_TT_COLUMNS_ENUM {
  ADT_DMET2_TT_GENE,
  ADT_DMET2_TT_REFERENCE_LINK,
  ADT_DMET2_TT_ASSAYID,
  ADT_DMET2_TT_DBSNP,
  ADT_DMET2_TT_DEFINING,
  ADT_DMET2_TT_CDNA,
  ADT_DMET2_TT_GENOMIC,
  ADT_DMET2_TT_CHANGE,
  ADT_DMET2_TT_EXTERNAL_ID,
  ADT_DMET2_TT_VALIDATED,
  ADT_DMET2_TT_HAPLOTYPE,
  ADT_DMET2_TT_REFERENCE,
  ADT_DMET2_TT_VARIANT,
  ADT_DMET2_TT_ALLELE_START,
};

class TranslationTableModel : public TranslationInputTsvTableModel
{

public:

  // Cache indexes, indicators and derived data for performance
  std::map< int, int >                          m_rowCopyNumberColumn;
  std::map< std::string, bool >                 m_probeSetToIsHaplotypeMarker;
  std::map< std::string, std::vector< int >  >  m_probeSetRowIndex;
  std::map< std::string, int >                  m_probeSetHeaderRowIndex;
  std::map< int, int >                          m_rowHeaderRowIndex;
  std::map< std::string, bool>                  m_geneCopyNumberIndicator;
  std::map< std::string, std::string >          m_reportAlleles;
  std::map< std::string, std::vector< std::string > >
                                                m_geneMarkers;
  
  // Methods
  TranslationTableModel(const RunTimeEnvironment & rte,
                        const std::string & ttableFileName,
                        ADT_TRANSLATION_TABLE_TYPE_ENUM ttableFileType);


  //~TranslationTableModel() {};

  void             generateGeneMarkerList();
  int              getColumnIndex(ADT_DMET3_TT_COLUMNS_ENUM) const;
  pcrecpp::RE      getColumnRegex(int column);
  int              getCopyNumberColumn(int row);
  std::vector<int> getHaplotypeAlleleColumns(const int row, const int startColumn = -1);
  int              getHeaderRow(const int childRow);
  std::string      getHeaderRowColumnName(int headerRow, int column);
  int              getProbeSetRowIndex(const std::string & probeSet,
                                  int multiAllelicRow = 0) const;
  int              getProbeSetRowIndexSize(const std::string & probeSet) const;
  std::string      getReportAllele( const std::string probeSet, const std::string allele);
  static ADT_TRANSLATION_TABLE_TYPE_ENUM
                   getTranslationTableFileType(const std::string & ttableFileName);
  ADT_TRANSLATION_TABLE_TYPE_ENUM
                   getType() const { return m_ttableFileType; }
  bool             ignoreRow(const int row, bool assignIgnore = false);
  bool             isHaplotypeMarker(const std::string & probeSet) const ;
  bool             isRowHaplotype(const int row);
  bool             probeSetFilter(const RunTimeEnvironment & rte, MarkerListModel & mlm);
  void             setReportAlleles( class GenoCallCoder & gcc );
  bool             validateRowIsHeader(const int row, bool giveWarning  = true);

private:

  std::map<   ADT_DMET3_TT_COLUMNS_ENUM, ADT_DMET2_TT_COLUMNS_ENUM >
                                 m_dmet3ToDmet2ColumnIndex;
  ADT_TRANSLATION_TABLE_TYPE_ENUM m_ttableFileType;
  std::set< int >                m_validatedHeaderRows;

  bool             _reconcileDeletedHaplotypeRows(const RunTimeEnvironment & rte, std::list< int >  deletedRows);
  

};


#endif /* TRANSLATION_TRANSLATION_TABLE_MODEL_H */
