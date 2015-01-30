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

/**
 * @file   QuantMethodTransform.h
 * @author Chuck Sugnet
 * @date   Mon Dec 14 14:39:27 2009
 * 
 * @brief  DataTransform wrapper for quantification methods
 * 
 * 
 */


#ifndef _QUANTMETHODTRANSFORM_H_
#define _QUANTMETHODTRANSFORM_H_

#include "chipstream/DataTransform.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodReport.h"
#include "chipstream/AnalysisInfo.h"

class QuantMethodTransform : public DataTransform {

public:

  QuantMethodTransform(QuantMethod *qMethod);

  ~QuantMethodTransform();

  virtual bool transformData(PsBoard &board, const DataStore &in, DataStore &out);

  /** 
   * Get a unique identifier for this data transform.
   * 
   * @return string that identifies the type of transform
   */
  virtual std::string getName() { return m_QuantMethod->getType(); }

private:

  bool doAnalysis(ProbeSetGroup &psGroup,
                  const IntensityMart &iMart,
                  PmAdjuster &pmAdj,
                  std::vector<QuantMethodReport *> &reporters,
                  bool doReport);

  void writeHeaders(PsBoard &board,
                    const IntensityMart &iMart,
                    QuantMethod *qMethod,
                    std::vector<QuantMethodReport *> &reporters,
                    std::vector<QuantMethodReport *> &exprReporters,
                    AnalysisInfo &info);

  QuantMethod *m_QuantMethod;

};

#endif /* _QUANTMETHODTRANSFORM_H_ */
