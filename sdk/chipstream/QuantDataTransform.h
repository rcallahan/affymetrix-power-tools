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
#ifndef _QUANTDATATRANSFORM_H_
#define _QUANTDATATRANSFORM_H_

//
#include "chipstream/DataStore.h"
#include "chipstream/ChipStream.h"
#include "chipstream/DataTransform.h"
#include "chipstream/SketchQuantNormTran.h"
#include "util/Util.h"
//

class QuantDataTransform {
  
public:


  QuantDataTransform(QuantMethodFactory::QuantType type, PsBoard &board, 
                     std::vector<QuantMethodReport *> &reporters,
                     const std::string &pmSpec, 
                     const std::string &qSpec,
                     const std::string &spfFile) {
    PmAdjusterFactory pmFactory;
    QuantMethodFactory qmFactory(type);
    m_QuantMethod = qmFactory.quantMethodForString(qSpec, board, type);
    m_PmAdjuster = pmFactory.pmAdjusterForString(pmSpec, board);
    m_Reporters = reporters;
    m_SpfFile = spfFile;
  }

  ~QuantDataTransform() {
    Freez(m_PmAdjuster);
    Freez(m_QuantMethod);
  }

  void setReporters(std::vector<QuantMethodReport *> &reporters) {
    m_Reporters = reporters;
  }

  QuantMethod *getQuantMethod() {
    return m_QuantMethod;
  }

  bool doAnalysis(ProbeSetGroup &psGroup,
                  IntensityMart &iMart,
                  bool doReport) {
    bool success = true;
    vector<ChipStream *> streams; // all the normalization should be done so we're using empty
    if(m_QuantMethod->setUp(psGroup, iMart, streams, *m_PmAdjuster)) {
      m_QuantMethod->computeEstimate();
      if(doReport) {
        for(unsigned int i = 0; i < m_Reporters.size(); i++) {
          m_Reporters[i]->report(psGroup, *m_QuantMethod, iMart, streams, *m_PmAdjuster);
        }
      }
    }
    else {
      Verbose::out(5, "Warning setup failed for name: " + ToStr(psGroup.name));
      success = false;
    }
    if(!success && doReport) {
      for(unsigned int i = 0; i < m_Reporters.size(); i++) {
        m_Reporters[i]->reportFailure(psGroup, *m_QuantMethod, iMart, streams, *m_PmAdjuster);
      }
    }
    return success;
  }
  
  void transformData(PsBoard &board, DataStore &in, DataStore &out) {
      SpfReader spfReader;
      spfReader.openSpf(m_SpfFile);
      
      int numPsSets = 0;
      board.get("num-probesets", &numPsSets);
      unsigned int dotMod = max(numPsSets/20, 1);
      Verbose::progressBegin(1, ToStr("Processing Probesets"), 20, (int)dotMod, numPsSets);
      ProbeSet *ps = spfReader.readNextProbeSet();
      while(ps != NULL) {
        Verbose::progressStep(1);
        ProbeSetGroup psGroup(ps); 
        doAnalysis(psGroup, in, true);
        // psGroup should delete the memory for ps...
        ps = spfReader.readNextProbeSet();
      }
      Verbose::progressEnd(1, ToStr("Done."));
  }

private:
  std::string m_SpfFile;
  PmAdjuster *m_PmAdjuster;
  QuantMethod *m_QuantMethod;
  std::vector<QuantMethodReport *> m_Reporters;

};

#endif  /* _QUANTDATATRANSFORM_H_ */
