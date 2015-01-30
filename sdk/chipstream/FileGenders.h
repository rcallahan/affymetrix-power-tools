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

#ifndef FILEGENDERS_H
#define FILEGENDERS_H

//
#include "chipstream/BioTypes.h"
#include "chipstream/GenderCalls.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @brief Read gender assignments from file
 */
class FileGenders : public GenderCalls {

public:

  /**
   * Constructor. Reads file with gender calls.
   * @param genderFile - the gender file to read. must have cel_files, and gender columns
   * @param celFiles   - the cel files names expected
   */
  FileGenders(const std::string &genderFile, const std::vector<std::string> &celFiles);

  /** 
   * Get the genders for the cel files that have been seen.
   * @return - vector of gender calls, one for each cel file in order seen.
   */
  std::vector<affx::Gender> getGenders() { return m_Genders; }

private:
    std::vector<affx::Gender> m_Genders;
  
};

#endif /* FILEGENDERS_H */
