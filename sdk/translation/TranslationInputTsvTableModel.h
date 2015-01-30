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
 * @file   TranslationInputTsvTableModel.h
 * @author Mybrid Spalding
 * @date   Wed Apr 23 10:09:02 PDT 2008
 * @brief  Base class for the multitude of single instance models of input TsvFile or Console PgOptions abstracted as a 'vector of std::vector of std::strings' table in memory. 
 */


#ifndef TRANSLATION_TRANSLATION_INPUT_TSV_TABLE_MODEL_H
#define TRANSLATION_TRANSLATION_INPUT_TSV_TABLE_MODEL_H

#include "translation/RunTimeEnvironment.h"
//
#include "file/TsvFile/TsvFile.h"
//
#include "affy-pcre.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//


class TittmColumnDefinition
{
public:
  std::string  m_columnName;
  bool         m_emptyOk;
  int          m_index;
  pcrecpp::RE  m_validRE;
  pcrecpp::RE  *m_ignoreRE;
};


class TranslationInputTsvTableModel
{

public:

  //
  std::vector<TittmColumnDefinition>      m_columnDefinition;
  std::map< int, int>                     m_fileRowIndex;
  std::vector< std::vector<std::string> > m_rows;
  bool                                    m_strictTsvColumnOrder;
  std::vector<int>                        m_tsvIndexToColumnDefinitionIndex;
  bool                                    m_validateColumnRegex;

  TranslationInputTsvTableModel() {}

  TranslationInputTsvTableModel(const RunTimeEnvironment& rte,
                                const std::string& inputTsvFileName,
                                const TittmColumnDefinition columnDefs[],
                                size_t cdSize,
                                bool strictTsvColumnOrder,
                                bool (*createdynamicColumnDefinitions)(const RunTimeEnvironment & rte, std::vector<TittmColumnDefinition> & columnDefs, const std::string & inputTsvFile) = NULL, bool (*filterRow)(std::vector< std::string > & row) = NULL,
                                bool validateColumnRegex = true);

  //
  ~TranslationInputTsvTableModel() {};

  void        clearData() {m_rows.resize(0);  }
  int         columnCount() const;
  std::string getColumnName(int column);
  std::string  getHeaderAsString(const std::string delimiter = ",") const;
  std::vector< std::string > getHeaderAsVector() const;
  std::string  getRowAsString(const int row, const std::string delimiter = ",") const;
  void        readTsvFile(const RunTimeEnvironment & rte,
                          const std::string& inputTsvFileName,
                          bool (*filterRow)(std::vector<std::string> & row));

  int         size()  const { return m_rows.size(); }

protected:
  std::map<std::string, int>              m_columnNameIndex;

private:

  void _initializeColumnDefinitions(const RunTimeEnvironment & rte,
                                    const std::vector<TittmColumnDefinition> c);


  bool _initializeTsvColumns(const std::string & inputTsvFileName,
                             affx::TsvFile & tsv);

  bool _matchTsvColumnsWithColumnDefinitions(const std::string & inputTsvFileName,
      affx::TsvFile & tsv);


};

#endif /* TRANSLATION_TRANSLATIONINPUTTSVTABLEMODEL_H */
