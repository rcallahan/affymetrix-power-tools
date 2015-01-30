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
 * @file   SampleInfoTableModel.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  A single instance class that wraps sample information either passed in by the console as a std::vector or when read from a file. 
 */

#ifndef TRANSLATION_SAMPLE_INFO_TABLE_MODEL_H
#define TRANSLATION_SMAPLE_INFO_TABLE_MODEL_H

#include "translation/TranslationInputTsvTableModel.h"


class SampleInfoTableModel : public TranslationInputTsvTableModel
{

public:

  std::string m_sampleInfoFileName;

  //Command Line
  SampleInfoTableModel(const class RunTimeEnvironment &rte,
                       const std::string & sampleInfoFileName);

  //Console
  SampleInfoTableModel(const class RunTimeEnvironment &rte);

  // STL
  SampleInfoTableModel() {}

  ~SampleInfoTableModel() {};

  std::vector< std::string > getSampleInfo(const std::string & experimentId);

};

#endif /* TRANSLATION_SAMPLE_INFO_TABLE_MODEL_H */
