////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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

#ifndef FILE5_EQUIVALENT_CHECK
#define FILE5_EQUIVALENT_CHECK

#include "util/RegressionCheck.h"

#include <string>
#include <set>

class File5EquivalentCheck : public RegressionCheck {

public:

  File5EquivalentCheck(const std::string &file1, const std::string &file2,
                       const std::string &group, const std::string &tsvName,
                       const std::set<std::string> &setIgnore,
                       double epsilon,
                       double corrCutoff,
                       bool bAllowNegation=false);

  bool check(std::string &msg);

private:
  std::string m_File1;
  std::string m_File2;
  std::string m_Group;
  std::string m_TsvName;
  std::set<std::string> m_SetIgnore;
  double m_Epsilon;
  double m_CorrCutoff;
  bool m_bAllowNegation;
};

#endif /* FILE5_EQUIVALENT_CHECK */
