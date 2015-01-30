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
 * @file   ComprehensiveReportTableModel.h
 * @author Mybrid Spalding
 * @date   Fri Nov 14 14:16:02 PST 2008
 * @brief  Comprehensive report wrapper for testing purposes.
 */

#ifndef TRANSLATION_COMPREHENSIVE_REPORT_TABLE_MODEL_H
#define TRANSLATION_COMPREHENSIVE_REPORT_TABLE_MODEL_H

#include "translation/TranslationInputTsvTableModel.h"

extern const TittmColumnDefinition CRTM_COLUMN_DEFINTIONS[];

class ComprehensiveReportTableModel : public TranslationInputTsvTableModel
{

public:




  ComprehensiveReportTableModel(const class RunTimeEnvironment & rte,
                                const std::string & reportFileName );


  ComprehensiveReportTableModel() {};
  ~ComprehensiveReportTableModel() {};


  int         getColumnIndex( std::string columnName );
  std::string getColumnName( int row);
  
private:
  std::map<  std::string, int> m_columnIndex;
  
};

#endif /* TRANSLATION_COMPREHENSIVE_REPORT_TABLE_MODEL_H */
