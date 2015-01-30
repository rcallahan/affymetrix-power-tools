////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////
/**
 * @file  CnchpCheck.cpp 
 * 
 * @brief Class for doing a check of two matrix files after regression run.
 */
//

#include "calvin_files/utils/src/Calvin.h"
#include "util/CnchpCheck.h"
#include "util/Fs.h"
#include "util/RegressionCheck.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <set>
#include <map>
//


bool CnchpCheck::check(std::string &msg) {
    msg = "";
    // Fix pathnames to work for this platform.
    m_NewDirectory = Fs::convertToUncPath(m_NewDirectory);
    m_GoldDirectory = Fs::convertToUncPath(m_GoldDirectory);

    bool bFilePass = true; 
    for (int iIndex=0; iIndex<m_vFileNames.size(); iIndex++ )
    {

        std::string goldFile = Fs::join(m_GoldDirectory , m_vFileNames[iIndex]) + "." + m_infix + ".CNCHP"; 
        std::string newFile = Fs::join(m_NewDirectory , m_vFileNames[iIndex]) + "." + m_infix + ".CNCHP"; 
        if (!Fs::isReadable(goldFile)) 
        {
            msg += "Can't open file: "+FS_QUOTE_PATH(goldFile)+" to read.";
            return false;
        }
        if (!Fs::isReadable(newFile)) 
        {
            msg += "Can't open file: "+FS_QUOTE_PATH(newFile)+" to read.";
            return false;
        }

        if(! Calvin::equivalent(
				goldFile,
                                newFile,
                                m_SetIgnore,
                                m_SetSetIgnore,
                                m_MapEpsilon,
                                m_dEpsilon,
                                m_dCorrelation,
                                m_bCheckHeader)) 
        {
            bFilePass = false;
        }
    }
    return bFilePass;
}
