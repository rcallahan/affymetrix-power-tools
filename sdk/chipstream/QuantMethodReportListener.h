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
 * @file   QuantMethodReportListener.h
 * @author Chuck Sugnet
 * @date   Tue Mar 28 10:04:33 2006
 * 
 * @brief Listens for reporting signals from QuantMethod and then passes
 * them off to another coordinating object. This is needed as the analysis
 * stream "owns" the reporters and will destroy them when the stream is
 * destroyed. This is a problem if you want multiple analysis streams
 * reporting into the same thing.
 */

#ifndef _QUANTMETHODEXPRREPORTLISTENER_H_
#define _QUANTMETHODEXPRREPORTLISTENER_H_

//
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodReport.h"
//
#include "util/Util.h"
//
#include <cstdlib>
#include <iostream>
//

class QuantMethodReportListener : public QuantMethodReport {

public:
  
  QuantMethodReportListener();

  ~QuantMethodReportListener();

  virtual bool prepare(QuantMethod &qMethod, const IntensityMart &iMart);

  virtual bool report(ProbeSetGroup &psGroup, 
                      QuantMethod &qMethod,
                      const IntensityMart &iMart, 
                      std::vector<ChipStream *> &iTrans, 
                      PmAdjuster &pmAdjust);

  virtual bool reportFailure(ProbeSetGroup &psGroup, 
                             QuantMethod &qMethod,
                             const IntensityMart &iMart, 
                             std::vector<ChipStream *> &iTrans, 
                             PmAdjuster &pmAdjust);

  virtual bool finish(QuantMethod &qMethod);


  void registerListener(QuantMethodReport *listener);

private:
  /// Listeners which we pass stuff off to - mem owned elsewhere
  std::vector<QuantMethodReport *> m_Listeners;
};

#endif /* _QUANTMETHODEXPRREPORTLISTENER_H_ */
