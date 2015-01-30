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

#include "chipstream/DTFactory.h"
#include "chipstream/ChipStreamDataTransform.h"
#include "chipstream/QuantMethodTransform.h"
#include "chipstream/PmAdjustTransform.h"

using namespace std;

DataTransform *DTFactory::dataTransformForString(const std::string &spec, 
                                                 PsBoard &board) {

  if (m_CSFactory.canBuild(spec)) {
    //  stage.spec is something like 'quant-norm.sketch=0'
    ChipStream *trans = m_CSFactory.chipStreamForStringBoard(board, spec, "");
    
    // Wrap the chipstream object as a data transform stage
    ChipStreamDataTransform *csStage = new ChipStreamDataTransform(trans);
    return csStage;
  }
  else if (m_QMFactory.canBuild(spec)) {
    QuantMethod *qMethod = m_QMFactory.quantMethodForString(spec, board);
    QuantMethodTransform *qTrans = new QuantMethodTransform(qMethod);
    return qTrans;
  }
  else if (m_PmAdjFactory.canBuild(spec)) {
    PmAdjuster *pmAdj = m_PmAdjFactory.pmAdjusterForString(spec, board);
    PmAdjusterTransform *pmAdjTrans = new PmAdjusterTransform(pmAdj);
    return pmAdjTrans;
  }
  else {
    Err::errAbort("DTFactory::dataTransformForString() - Don't know how to build specification: " + spec);
  }
  return NULL;
}
