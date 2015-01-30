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

#ifndef GENDERCALLS_H
#define GENDERCALLS_H

//
#include "chipstream/BioTypes.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @brief Class for determining gender from a cel set of cel files.
 */
class GenderCalls {

public:

  /** 
   * Get the genders for the cel files that have been seen.
   * @return - vector of gender calls, one for each cel file in order seen.
   */
  virtual std::vector<affx::Gender> getGenders() = 0;

  std::string getGenderName() { return m_GenderName; }
  std::string getGenderDescription() { return m_GenderDescription; }

  virtual ~GenderCalls(){}

protected:
  std::string m_GenderName;
  std::string m_GenderDescription;
  
};

/**
 * @brief Instance of GenderCalls which provides unknown genders
 */
class NoGenderCalls : public GenderCalls {

public:
    NoGenderCalls(int numberChips) {
        m_GenderCalls.resize(numberChips);
        for(size_t i=0; i<m_GenderCalls.size(); i++)
            m_GenderCalls[i] = affx::UnknownGender;
        m_GenderName = "none";
        m_GenderDescription = "no gender calls - assume unknown";
    }
    std::vector<affx::Gender> getGenders() { return m_GenderCalls; }

private:
    std::vector<affx::Gender> m_GenderCalls;
};

#endif /* GENDERCALLS_H */
