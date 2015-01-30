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

#ifndef FILEINBRED_H
#define FILEINBRED_H

//
#include "chipstream/BioTypes.h"
#include "chipstream/InbredStatus.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @brief Read sample covariate relative inbreeding level assignments from file
 */
class FileInbred : public InbredStatus {

public:

  /**
   * Constructor. Reads file with inbreeding level.
   * @param inbredFile - the inbreeding file to read. must have cel_files, and inbred_het_penalty columns
   * @param celFiles   - the cel files names expected
   */
  FileInbred(const std::string &inbredFile, const std::vector<std::string> &celFiles);

  /** 
   * Get the  penalty for inbreeding levels the cel files that have been seen.
   * @return - vector of inbred covariates, one for each cel file in order seen.
   */
  std::vector<double> getInbredStatus() { return m_InbredStatus; }

private:
    std::vector<double> m_InbredStatus;
  
};

#endif /* FILEGENDERS_H */
