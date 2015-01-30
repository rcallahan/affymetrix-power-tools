////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#include "chipstream/QuantMethodTransform.h"
#include "chipstream/SpfReader.h"
#include "chipstream/QuantLabelZ.h"
#include "chipstream/QuantBRLMM.h"
#include "chipstream/QuantBirdseed.h"
#include "chipstream/QuantMethodGTypeReport.h"
#include "chipstream/QuantMethodExprReport.h"
#include "util/Util.h"
#include "util/Guid.h"

using namespace std;

QuantMethodTransform::QuantMethodTransform(QuantMethod *qMethod) {
  m_QuantMethod = qMethod;
}

QuantMethodTransform::~QuantMethodTransform() {
  Freez(m_QuantMethod);
}

bool QuantMethodTransform::doAnalysis(ProbeSetGroup &psGroup,
                                      const IntensityMart &iMart,
                                      PmAdjuster &pmAdj,
                                      std::vector<QuantMethodReport *> &reporters,
                                      bool doReport) {
  bool success = true;
  vector<ChipStream *> streams; // all the normalization should be done so we're using empty
  if(m_QuantMethod->setUp(psGroup, iMart, streams, pmAdj)) {
    m_QuantMethod->computeEstimate();
    if(doReport) {
      for(unsigned int i = 0; i < reporters.size(); i++) {
        reporters[i]->report(psGroup, *m_QuantMethod, iMart, streams, pmAdj);
      }
    }
  }
  else {
    Verbose::out(5, "Warning setup failed for name: " + ToStr(psGroup.name));
    success = false;
  }
  if(!success && doReport) {
    for(unsigned int i = 0; i < reporters.size(); i++) {
      reporters[i]->reportFailure(psGroup, *m_QuantMethod, iMart, streams, pmAdj);
    }
  }
  return success;
}


void QuantMethodTransform::writeHeaders(PsBoard &board,
                                        const IntensityMart &iMart,
                                        QuantMethod *qMethod,
                                        std::vector<QuantMethodReport *> &reporters,
                                        std::vector<QuantMethodReport *> &exprReporters,
                                        AnalysisInfo &info) {
  Options *o = board.getOptions();
  // @todo should the guid be set someplace else?
  string reportGuid = affxutil::Guid::GenerateNewGuid();
  string execGuid =  o->getOpt("exec-guid");
  string timeStart = o->getOpt("time-start");
  string cmdLine = o->getOpt("command-line");
  string version = o->getOpt("version-to-report");
  for(unsigned int i = 0; i < reporters.size(); i++) {
    QuantMethodReport *reporter = reporters[i];
    reporter->addStdHeaders(reporter,
                            execGuid,
                            reportGuid,
                            timeStart,
                            cmdLine, 
                            version,
                            info);
    reporter->prepare(*qMethod, iMart);
  }
  
  for(size_t i = 0; i < exprReporters.size(); i++) {
    QuantMethodReport *reporter = exprReporters[i];
    reporter->addStdHeaders(reporter,
                            execGuid,
                            reportGuid,
                            timeStart, 
                            cmdLine,
                            version,
                            info);
    //    reporter->prepare(*qMethod, iMart);
  }


  
}  

bool QuantMethodTransform::transformData(PsBoard &board, const DataStore &in, DataStore &out) {
  std::vector<QuantMethodReport *> *gtypeReporters = board.getQuantGTypeReporters();
  APT_ERR_ASSERT(gtypeReporters != NULL, "QuantMethodtransform::transformData() - Reporters can't be NULL.");
  /* @todo refactor - finalize the blackboard interface... */
  AnalysisInfo *info = board.getAnalysisInfo();
  PmAdjuster *pmAdj = board.getPmAdjuster();
  ChipLayout *layout = NULL;
  board.get("chiplayout", &layout);
  std::vector<QuantMethodReport *> *exprReporters = board.getQuantExpressionReporters();
  if (InstanceOf(m_QuantMethod, QuantGTypeMethod)) {
    QuantGTypeMethod *qMethod = static_cast<QuantGTypeMethod *>(m_QuantMethod);
    float min = qMethod->getMinThresh();
    info->addParam("min-thresh",ToStr(min));
    float max = qMethod->getMaxThresh();
    info->addParam("max-thresh",ToStr(max));
  }
  info->addParam("quantification-name", m_QuantMethod->getType());
  info->addParam("quantification-version", m_QuantMethod->getVersion());
  info->addParam("quantification-scale", QuantMethod::scaleToTxt(m_QuantMethod->getScale()));
  info->addParam("quantification-type", QuantMethod::quantTypeToTxt(m_QuantMethod->getQuantType()));
  info->m_AlgVersion = m_QuantMethod->getVersion();
  
  for (int i = 0; i < exprReporters->size(); i++) {
    // @todo - refactor make this vector be expr reporters only to avoid cast...
    QuantMethodExprReport *report = static_cast<QuantMethodExprReport *>((*exprReporters)[i]);
    if(InstanceOf(m_QuantMethod,QuantLabelZ)) {
      report->m_qMethod = static_cast<QuantLabelZ *>(m_QuantMethod)->getQuantExprMethod();
      static_cast<QuantLabelZ *>(m_QuantMethod)->addExprReporter(report);
    } else if (InstanceOf(m_QuantMethod,QuantBRLMM)) {
      report->m_qMethod = static_cast<QuantBRLMM *>(m_QuantMethod)->getQuantExprMethod();
      static_cast<QuantBRLMM *>(m_QuantMethod)->addExprReporter(report);
    } else if (InstanceOf(m_QuantMethod,QuantBirdseed)) {
      report->m_qMethod = static_cast<QuantBirdseed *>(m_QuantMethod)->getQuantExprMethod();
      static_cast<QuantBirdseed *>(m_QuantMethod)->addExprReporter(report);
    }
  }

  int numPsSets = layout->getProbeSetCount();
  unsigned int dotMod = max(numPsSets/20, 1);
  writeHeaders(board, in, m_QuantMethod, *gtypeReporters, *exprReporters, *info);
  m_QuantMethod->prepare(in);

  Verbose::progressBegin(1, "Processing Probesets" , 20, (int)dotMod, numPsSets);
  for (unsigned int i = 0; i < layout->getProbeSetCount(); i++) {
    Verbose::progressStep(1);
    ProbeListPacked pList = layout->getProbeListAtIndex(i);
    ProbeSet *ps = ProbeListFactory::asProbeSet(pList);
    ProbeSetGroup psGroup(ps); 
    doAnalysis(psGroup, in, *pmAdj, *gtypeReporters, true);
    // psGroup should delete the memory for ps...
    //    ps = spfReader.readNextProbeSet();
  }
  Verbose::progressEnd(1, ToStr("Done."));  
  m_QuantMethod->finish();
  for(unsigned int i = 0; i < gtypeReporters->size(); i++) {
    (*gtypeReporters)[i]->finish(*m_QuantMethod);
  }
  return false;
}
