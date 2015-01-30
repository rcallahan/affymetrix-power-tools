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
#ifndef CNCHPCHECK_H
#define CNCHPCHECK_H

//

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

/**
 * Class for testing that two cnchp files are the same +/- some
 * epsilon.
 */
class CnchpCheck : public RegressionCheck {

private:
    std::vector<std::string> m_vFileNames;
    std::string m_GoldDirectory;
    std::string m_NewDirectory;
    std::string m_infix;
    std::set<std::string> m_SetIgnore;
    std::set<std::string> m_SetSetIgnore;
    std::map<std::string, float> m_MapEpsilon;
    double  m_dEpsilon;
    double  m_dCorrelation;
    bool    m_bCheckHeader;   



public:

  CnchpCheck(   const std::vector<std::string> &vFileNames,
                const std::string &goldDirectory, 
                const std::string &newDirectory, 
                const std::string &infix,
                std::set<std::string> setIgnore, 
                std::set<std::string> setSetIgnore, 
                std::map<std::string, float> mapEpsilon,
                double epsilon,
                double correlation,
                bool bCheckHeader) :  
                    m_vFileNames(vFileNames), m_GoldDirectory(goldDirectory), m_NewDirectory(newDirectory), 
                    m_infix(infix), m_SetIgnore(setIgnore),
                    m_SetSetIgnore(setSetIgnore), m_MapEpsilon(mapEpsilon), m_dEpsilon(epsilon),
                    m_dCorrelation(correlation), m_bCheckHeader(bCheckHeader){};
    

   // @return - Returns true if files are close enough, false otherwise.
  bool check(std::string &msg);
};


#endif /* CNCHPCHECK_H */
