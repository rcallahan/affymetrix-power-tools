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
 * @file   ExperimentResults.h
 * @author Mybrid Spalding
 * @date   Tue Sep 16 15:45:30 PDT 2008
 * @brief  Wrapper class for ExperimentGeneResults, all gene experiment results. There should ever only ever be one of these.
 */

#ifndef TRANSLATION_EXPERIMENT_RESULTS_H
#define TRANSLATION_EXPERIMENT_RESULTS_H

#include "translation/ExperimentGeneResults.h"

class ExperimentResults
{
public:

  std::string m_experimentGuid;
  std::string m_experiment;
  std::map< std::string, ExperimentGeneResults*>
  m_geneResults;
  int         m_markerCallCount;
  std::string m_sample;

  ExperimentResults() {
    m_markerCallCount = 0;
  }

  void Clear();

};

#endif /* TRANSLATION_EXPERIMENT_RESULTS_H */
