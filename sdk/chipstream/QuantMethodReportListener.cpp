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

//
#include "chipstream/QuantMethodReportListener.h"

QuantMethodReportListener::QuantMethodReportListener()
{}

QuantMethodReportListener::~QuantMethodReportListener()
{}

bool QuantMethodReportListener::prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
  bool result = true;
  for(int i=0; i<m_Listeners.size(); i++)
    result = result && m_Listeners[i]->prepare(qMethod, iMart);
  return result;
}

bool QuantMethodReportListener::report(ProbeSetGroup &psGroup, 
                                       QuantMethod &qMethod,
                                       const IntensityMart &iMart, 
                                       std::vector<ChipStream *> &iTrans, 
                                       PmAdjuster &pmAdjust) { 
  bool result = true;
  for(int i=0; i<m_Listeners.size(); i++)
    result = result && m_Listeners[i]->report(psGroup, qMethod, iMart, iTrans, pmAdjust);
  return result;
}

bool QuantMethodReportListener::reportFailure(ProbeSetGroup &psGroup, 
                                              QuantMethod &qMethod,
                                              const IntensityMart &iMart, 
                                              std::vector<ChipStream *> &iTrans, 
                                              PmAdjuster &pmAdjust) { 
  bool result = true;
  for(int i=0; i<m_Listeners.size(); i++)
    result = result && m_Listeners[i]->reportFailure(psGroup, qMethod, iMart, iTrans, pmAdjust);
  return result;
}

bool QuantMethodReportListener::finish(QuantMethod &qMethod) {
  bool result = true;
  for(int i=0; i<m_Listeners.size(); i++)
    result = result && m_Listeners[i]->finish(qMethod);
  return result;
}

//  virtual void out(const std::string &s, bool comment=false) {
//      for(int i=0; i<m_Listeners.size(); i++)
//          m_Listeners[i]->out(s,comment);
//  }

//  virtual void commentDelim() {
//      for(int i=0; i<m_Listeners.size(); i++)
//          m_Listeners[i]->commentDelim();
//  }

//  virtual void printHeader(const std::vector<std::string> &s, const std::string& timeStr) {
//      for(int i=0; i<m_Listeners.size(); i++)
//          m_Listeners[i]->printHeader(s, timeStr);
//  }

//  virtual void finishTextHeader() {
//      for(int i=0; i<m_Listeners.size(); i++)
//          m_Listeners[i]->finishTextHeader();
//  }

//  virtual void writeHeader(QuantMethodReport *qReport,
//                           const std::string& execGuid, 
//                           const std::string& reportGuid,
//                           const std::string& timeStr,
//                           const std::string& commandLine,
//                           const std::vector<std::string> &colNames,
//                           const std::string& execVersion, 
//                           const AnalysisInfo &info) {
//    for(int i=0; i<m_Listeners.size(); i++)
//      m_Listeners[i]->writeHeader(qReport, execGuid, reportGuid, timeStr, commandLine, colNames, execVersion, info);
//  }

void QuantMethodReportListener::registerListener(QuantMethodReport *listener) {
  m_Listeners.push_back(listener);
}


