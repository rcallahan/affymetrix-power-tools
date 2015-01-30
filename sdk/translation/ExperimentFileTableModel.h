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
 * @file   ExperimentFileTableModel.h
 * @author Mybrid Spalding
 * @date   Wed Sep 24 10:03:50 PDT 2008
 * @brief  Class for wrapping the CHP list of experiment files.
 */

#ifndef TRANSLATION_EXPERIMENT_FILE_TABLE_MODEL_H
#define TRANSLATION_EXPERIMENT_FILE_TABLE_MODEL_H

#include "translation/TranslationInputTsvTableModel.h"


class ExperimentFileTableModel : public TranslationInputTsvTableModel
{

public:

  // DATA
  enum  {
    EXPERIMENT_FILE_INDEX,
    END_INDEX,
  };

  std::vector< std::string > m_experimentFiles;

  ExperimentFileTableModel(const RunTimeEnvironment &rte,
                           const std::string & probesetFilterFileName);


  ~ExperimentFileTableModel() {};

  bool validateFiles();


};

#endif /* TRANSLATION_EXPERIMENT_FILE_TABLE_MODEL_H */
