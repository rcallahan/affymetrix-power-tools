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

#include "copynumber/apt-copynumber-cyto/regression/File5EquivalentCheck.h"
#include "file5/File5.h"

File5EquivalentCheck::File5EquivalentCheck(const std::string &file1, const std::string &file2,
                                           const std::string &group, const std::string &tsvName,
                                           const std::set<std::string> &setIgnore,
                                           double epsilon,
                                           double corrCutoff,
                                           bool bAllowNegation) :
    m_File1(file1),
    m_File2(file2),
    m_Group(group),
    m_TsvName(tsvName),
    m_SetIgnore(setIgnore),
    m_Epsilon(epsilon),
    m_CorrCutoff(corrCutoff),
    m_bAllowNegation(bAllowNegation) {
    }

bool File5EquivalentCheck::check(std::string &msg) {
    return affx::File5_File::equivalent(m_File1, m_File2, m_Group, m_TsvName, m_SetIgnore, m_Epsilon, m_CorrCutoff, m_bAllowNegation);
}
