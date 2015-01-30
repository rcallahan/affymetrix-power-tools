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

#ifndef INBREDSTATUS_H
#define INBREDSTATUS_H

//
#include "chipstream/BioTypes.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @brief Class for determining inbreeding covariate from a cel set of cel files.
 */
class InbredStatus {

public:

  /** 
   * Get the inbreeding status for the cel files that have been seen.
   * @return - vector of inbreeding covariates, one for each cel file in order seen.
   */
  virtual std::vector<double> getInbredStatus() = 0;

  std::string getInbredName() { return m_InbredName; }
  std::string getInbredDescription() { return m_InbredDescription; }

  virtual ~InbredStatus(){}

protected:
  std::string m_InbredName;
  std::string m_InbredDescription;
  
};

/**
 * @brief Instance of InbredStatus which provides null covariates 
 */
class NoInbredStatus : public InbredStatus {

public:
    NoInbredStatus(int numberChips) {
        m_InbredStatus.resize(numberChips);
        for(size_t i=0; i<m_InbredStatus.size(); i++)
            m_InbredStatus[i] = 0;
        m_InbredName = "none";
        m_InbredDescription = "no inbreeding status, no modifications to likelihood";
    }
    std::vector<double> getInbredStatus() { return m_InbredStatus; }

private:
    std::vector<double> m_InbredStatus;
};

#endif /* GENDERCALLS_H */
