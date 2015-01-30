////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
#ifndef _PSBOARD_H_
#define _PSBOARD_H_

class PsBoard;

#include "bboard/Bboard.h"
#include "chipstream/DataStore.h"
#include "util/Options.h"

class PmAdjuster;
class AnalysisInfo;
class QuantMethodReport;
class GenotypeInfo;

class PsBoard : public Bboard {

public:

  PsBoard() {
    m_ProbeAnnotations = NULL;
    m_Options = NULL;
    m_PmAdjuster = NULL;
    m_AnalysisInfo = NULL;
    m_GTypeReporters = NULL;
    m_ExpressionReporters = NULL;
    m_GenotypeInfo = NULL;
  }

  ~PsBoard();

  void setOptions(Options *opts) {
    m_Options = opts;
  }

  Options *getOptions() {
    return m_Options;
  }

  void setPmAdjuster(PmAdjuster *pmAdj) {
    m_PmAdjuster = pmAdj;
  }

  PmAdjuster *getPmAdjuster() {
    return m_PmAdjuster;
  }

  void setAnalysisInfo(AnalysisInfo *aInfo) {
    m_AnalysisInfo = aInfo;
  }

  AnalysisInfo *getAnalysisInfo() {
    return m_AnalysisInfo;
  }

  void setProbeInfo(DataStore *ds) {
    m_ProbeAnnotations = ds;
  }

  DataStore *getProbeInfo() {
    return m_ProbeAnnotations;
  }

  void setGenotypeInfo(GenotypeInfo *gInfo) {
    m_GenotypeInfo = gInfo;
  }

  GenotypeInfo *getGenotypeInfo() {
    return m_GenotypeInfo;
  }

  void setQuantGTypeReporters(std::vector<QuantMethodReport *> *reporters) {
    m_GTypeReporters = reporters;
  }

  std::vector<QuantMethodReport *> *getQuantGTypeReporters() {
    return m_GTypeReporters;
  }

  void setQuantExpressionReporters(std::vector<QuantMethodReport *> *reporters) {
    m_ExpressionReporters = reporters;
  }

  std::vector<QuantMethodReport *> *getQuantExpressionReporters() {
    return m_ExpressionReporters;
  }


private: 
  DataStore *m_ProbeAnnotations;
  Options *m_Options;
  PmAdjuster *m_PmAdjuster;
  AnalysisInfo *m_AnalysisInfo;
  std::vector<QuantMethodReport *> *m_GTypeReporters;
  std::vector<QuantMethodReport *> *m_ExpressionReporters;
  GenotypeInfo *m_GenotypeInfo;
};

#endif /* _PSBOARD_H_ */
