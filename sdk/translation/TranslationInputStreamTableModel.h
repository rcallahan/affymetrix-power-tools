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
 * @file   TranslationInputStreamTableModel.h
 * @author Mybrid Spalding
 * @date   Fri Jun 13 13:09:30 PDT 2008
 * @brief  Started as base class data model for either CHP or TSV files, but ended as an exclusive class to GenotypeTableModel and contains GenotypeTableModel data. 
 */


#ifndef TRANSLATION_TRANSLATION_INPUT_STREAM_TABLE_MODEL_H
#define TRANSLATION_TRANSLATION_INPUT_STREAM_TABLE_MODEL_H

#include "translation/GenotypeOverrideTableModel.h"
#include "translation/RunTimeEnvironment.h"
#include "translation/TranslationInputTsvTableModel.h"
//
#include "calvin_files/array/src/ArrayId.h"
#include "calvin_files/parsers/src/CHPMultiDataFileReader.h"
#include "calvin_files/utils/src/GenoCallCoder.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h" // includes "util/Verbose.h"
//
#include "affy-pcre.h"
//
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>
//


enum GenotypeTableIndex {
  GT_SAMPLE_INDEX,
  GT_EXPERIMENT_INDEX,
  GT_GENES_INDEX,
  GT_EXTERNAL_ID_INDEX,
  GT_PROBE_SET_INDEX,
  GT_ALLELE1_INDEX,
  GT_ALLELE2_INDEX,
  GT_ORIGINAL_ALLELE1_INDEX,
  GT_ORIGINAL_ALLELE2_INDEX,
  GT_END_FIELD_INDEX,
};

const affymetrix_calvin_io::MultiDataType ADT_CHP_DATA_TYPES[]
= {
    affymetrix_calvin_io::DmetCopyNumberMultiDataType,
    affymetrix_calvin_io::DmetMultiAllelicMultiDataType,
    affymetrix_calvin_io::DmetBiAllelicMultiDataType,
};


class TranslationInputStreamTableModel
{

public:

  enum { COPY_TYPE, MULTI_TYPE,  BI_TYPE, };


  std::vector< std::vector<std::string> >      m_rows;
  GenotypeOverrideTableModel                  *m_gotm;
  ADT_EXPERIMENT_STREAM_TYPE_ENUM              m_type;


  TranslationInputStreamTableModel() {}

  // TSV/DMET2 constructor
  TranslationInputStreamTableModel(const RunTimeEnvironment& rte,
                                   const std::string& inputTsvFileName,
                                   const TittmColumnDefinition columnDefs[],
                                   size_t tcdSize,
                                   bool strictTsvColumnOrder,
                                   int tvsExperimentColumn,
                                   int tsvGeneColumn,
                                   int tsvProbeSetColumn,
                                   GenotypeOverrideTableModel * gotm = NULL);

  // CHP/DMET3 constructor
  TranslationInputStreamTableModel(const RunTimeEnvironment & rte,
                                   const TittmColumnDefinition columnDefs[],
                                   size_t tcdSize,
                                   GenotypeOverrideTableModel *gotm = NULL);



  //
  ~TranslationInputStreamTableModel() {
    if (m_genoCallCoder != NULL) {
      delete m_genoCallCoder;
    }
    if (m_chpData != NULL) {
      delete m_chpData;
    }
  };

  // TSV or CHP specific methods have CHP or Tsv in the name. 
  void        clearData();
  int         columnCount() const;
  std::string getCHPGuid() { return m_chpGuid; }
  std::string getColumnName(int column);
  std::string getRowAsString(const int row) const;
  int         readNextExperiment(const RunTimeEnvironment & rte,
                                 class TranslationTableModel & ttm,
                                 const std::map<std::string, std::string> &gecnc,
                                 std::map<std::string, int> & geneCopyNumber);
  
  int         size()  const { return m_rows.size(); }
  void        sortTsvFileByExperiment(const RunTimeEnvironment & rte);

private:

#ifndef WIN32
	  int m_fd;
#endif
	  
  // Common data
  // AFter we've read in an experiment, check for duplicates
  std::string                                m_experimentName;
  std::map< std::string, std::vector<int> >  m_experimentRowIndex;
  std::map< int, int>                        m_fileRowIndex;

  // TSV data for Genotype Short Report input
  bool                          m_strictTsvColumnOrder;
  affx::TsvFile                 m_tsv;
  int                           m_tsvColumnCount;
  std::vector<TittmColumnDefinition> m_tsvColumnDefinition;
  std::string                   m_tsvFileName;
  int                           m_tsvExperimentColumn;
  int                           m_tsvGeneColumn;
  std::vector<int>              m_tsvIndexToColumnDefinitionIndex;
  std::vector<std::string>      m_tsvNextRow;
  int                           m_tsvProbeSetColumn;
  int                           m_tsvRow;

  // CHP type for CHP files input
  affymetrix_calvin_io::CHPMultiDataData  *m_chpData;
  std::string                        m_chpDataEncodingSize;
  std::map< affymetrix_calvin_io::MultiDataType, int32_t >
                                     m_chpDataTypeSize;
  std::string                        m_chpGuid;
  std::string                        m_chpExperimentFile;
  std::list< std::string >           m_chpExperimentFileQueue;
  std::string                        m_chpMaxAlleles;
  std::string                        m_chpVersion;
  GenoCallCoder *                    m_genoCallCoder;
  std::map< std::string, int >       m_probeSetRowIndex;


  // HELPER methods. 

  void   _closeTsvFile();
  void   _openTsvFile(const RunTimeEnvironment & rte);

  void   _initializeColumnDefinitions(const RunTimeEnvironment & rte,
                                      const TittmColumnDefinition tsvColumnDefinitions[],
                                      int tcdSize);


  bool   _getAllelesFromProbeSet(const RunTimeEnvironment& rte,
                                 TranslationTableModel & ttm,
                                 affymetrix_calvin_data::DmetBiAllelicData & dbad,
                                 std::string & allele1,
                                 std::string & allele2);

  std::string _getCHPGuid(const RunTimeEnvironment & rte, affymetrix_calvin_io::GenericData & gd);
  std::string _getGeneFromProbeSet(const std::string probeSet,
                                   const TranslationTableModel & ttm);

  std::string _getExternalIdFromProbeSet(const std::string probeSet,
                                         const TranslationTableModel & ttm);

  bool   _initializeNextExperimentCHPFile(const RunTimeEnvironment & rte,
                                          TranslationTableModel &ttm);

  bool   _initializeTsvColumns(const std::string & inputTsvFileName,
                               affx::TsvFile & tsv);

  bool   _matchTsvColumnsWithColumnDefinitions(const std::string & inputTsvFileName,
      affx::TsvFile & tsv);
  std::string _readCopyNumber(const RunTimeEnvironment & rte,
                              int row,
                              affymetrix_calvin_data::DmetBiAllelicData & dbad);

  int    _readNextExperimentCHP(const RunTimeEnvironment & rte,
                                TranslationTableModel & ttm,
                                const std::map<std::string, std::string> &gecnc,
                                std::map< std::string, int> & geneCopyNumber);

  int    _readNextExperimentTsv(const RunTimeEnvironment & rte,
                                const TranslationTableModel & ttm,
                                const std::map<std::string, std::string> &gecnc);

  void   _updateGeneCopyNumberState(std::vector< std::string> & newRow,
                                    const TranslationTableModel & ttm,
                                    const std::map<std::string, std::string> & gecnc);

};

#endif /* TRANSLATION_TRANSLATION_INPUT_STREAM_TABLE_MODEL_H */
