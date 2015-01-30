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

#include "chipstream/ChipStreamDataTransform.h"
#include "chipstream/ArtifactReduction.h"

bool ChipStreamDataTransform::transformData(PsBoard &board, const DataStore &in, DataStore &out) {
  bool needTwoPass = false;
  if((InstanceOf(m_ChipStream, SketchQuantNormTran) && board.getOptions()->getOpt("target-sketch").empty()) || 
      InstanceOf(m_ChipStream, ArtifactReduction)) {
    Verbose::out(1, "Doing two passes.");
    needTwoPass = true;
  }
  
  /* Send the chip data through the chipstream object. */
  int numDataSets = in.getCelDataSetCount();
  Verbose::progressBegin(1, ToStr("Processing Chips"), numDataSets, 1, numDataSets);
  std::vector<float> data;
  for(int chipIx = 0; chipIx < in.getCelDataSetCount(); chipIx++) {
    Verbose::progressStep(1);
    data.clear();
    //    Verbose::out(3, "Reading " + m_ChipStream->getDocName() + " for chip: " + ToStr(chipIx));
    in.fillInCelData(chipIx, data);
    m_ChipStream->newChip(data);
    if(!needTwoPass) {
      for(int probeIx = 0; probeIx < data.size(); probeIx++) {
        data[probeIx] = m_ChipStream->transform(probeIx, chipIx, data[probeIx], board);
      }
      out.writeColumn(chipIx, in.getCelName(chipIx), data);
    }
  }
  Verbose::progressEnd(1, ToStr("Done."));

  if(needTwoPass) {
    m_ChipStream->finishedChips();
    Verbose::progressBegin(1, ToStr("Processing Stage Two Chips"), numDataSets, 1, numDataSets);
    /* Do the chipstream transformation and store the new intensities. */
    for(int chipIx = 0; chipIx < in.getCelDataSetCount(); chipIx++) {
      Verbose::progressStep(1);
    //      Verbose::out(3, "Applying " + m_ChipStream->getDocName() + " to chip: " + ToStr(chipIx));
      std::vector<float> data;
      in.fillInCelData(chipIx, data);
        
      for(int probeIx = 0; probeIx < data.size(); probeIx++) {
        data[probeIx] = m_ChipStream->transform(probeIx, chipIx, data[probeIx]);
      }
      out.writeColumn(chipIx, in.getCelName(chipIx), data);
    }
    Verbose::progressEnd(1, ToStr("Done."));
  }
  return true;
}
