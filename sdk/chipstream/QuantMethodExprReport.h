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
 * @file   QuantMethodExprReport.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 12:09:42 2005
 *
 * @brief  Class for reporting results of quantification methods.
 */

#ifndef _QUANTMETHODEXPRREPORT_H_
#define _QUANTMETHODEXPRREPORT_H_

//
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantMethodReport.h"
#include "chipstream/TsvReport.h"
//
#include "util/Util.h"
//
#include <cstdlib>
#include <iostream>
//

/**
 *   Class for reporting results of quantification methods.
 */
class QuantMethodExprReport : public QuantMethodReport {

public:
  
  QuantMethod* m_qMethod;

  void setCompactFile5Format(int maxNameLength=TSVREPORT_PROBESET_STRLEN);
 
  bool m_DoSummary;
  affx::TsvReport m_summary_tsv;
  int m_summary_tsv_precision;
  affx::File5_Group* m_summary_a5_group;

  bool m_DoFeatureEffects;
  affx::TsvReport m_feffects_tsv;
  int m_feffects_tsv_precision;
  affx::File5_Group* m_feffects_a5_group;
  bool m_WriteOldStyleFeatureEffectsFile;


  bool m_DoResiduals;
  affx::TsvReport m_residuals_tsv;
  int m_residuals_tsv_precision;
  int m_NumCol;
  bool m_unmungePSName;  // strip off allele suffix added by apt-probeset-genotype.

  affx::File5_Group* m_residuals_a5_group;

  QuantMethodExprReport(int numCol);
  virtual ~QuantMethodExprReport();

  void setDoResiduals(bool b) {m_DoResiduals = b;}

  void setWriteOldStyleFeatureEffectsFile(bool inputValue) { m_WriteOldStyleFeatureEffectsFile = inputValue;}
  bool getWriteOldStyleFeatureEffectsFile(){ return m_WriteOldStyleFeatureEffectsFile;}

  bool prepare(QuantMethod &qMethod,
               const IntensityMart &iMart);

  bool report(ProbeSetGroup &psGroup,
              QuantMethod &qMethod,
              const IntensityMart &iMart,
              std::vector<ChipStream *> &iTrans,
              PmAdjuster &pmAdjust);

  bool finish(QuantMethod &qMethod);

  void registerProbeSetsToReport(std::set<const char *, Util::ltstr> *probeSetsToReport) {
      m_ProbeSetsToReport = probeSetsToReport;
  }

protected: 
  QuantMethodExprReport();

private:
  void init();
  std::set<const char *, Util::ltstr> *m_ProbeSetsToReport;

  affx::File5_File *m_File5;
  affx::File5_Group* m_F5Group;
  int m_MaxNameLength;
  /// Should we output confidences as a float rather than a double and calls as a char
  bool m_DoCompact;
};

#endif /* _QUANTMETHODEXPRREPORT_H_ */
