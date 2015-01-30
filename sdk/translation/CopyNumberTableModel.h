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
 * @file   CopyNumberTableModel.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  DEPRECATED DMET2 ONLY: Class for wrapping the DMET2 copy number report TSV file.
 */

#ifndef TRANSLATION_COPY_NUMBER_TABLE_MODEL_H
#define TRANSLATION_COPY_NUMBER_TABLE_MODEL_H

#include "translation/TranslationInputTsvTableModel.h"

bool isGeneCopyNumberSensitive(std::string gene);

class CopyNumberTableModel : public TranslationInputTsvTableModel
{

public:

  // DATA
  enum  {
    EXPERIMENT_INDEX,
    SAMPLE_INDEX,
    PREDICTION_INDEX,
    NOTES_INDEX,
    END_CNTM_FIELD_INDEX,
  };

  std::map<std::string, std::string> m_experimentCopyNumberCall;

  CopyNumberTableModel(const RunTimeEnvironment &rte,
                       const std::string & copyNumberReportFileName);

  CopyNumberTableModel() {}

  ~CopyNumberTableModel() {};

};

#endif /* TRANSLATION_COPY_NUMBER_TABLE_MODEL_H */
