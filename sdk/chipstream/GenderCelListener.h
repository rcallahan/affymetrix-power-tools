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

#ifndef GENDERCELLISTENER_H
#define GENDERCELLISTENER_H

//
#include "chipstream/BioTypes.h"
#include "chipstream/CelListener.h"
#include "chipstream/ChipSummary.h"
#include "chipstream/GenderCalls.h"
//
#include <cstring>
#include <string>
#include <utility>
#include <vector>
//

/**
 * @brief Class for determining gender from a cel set of cel files.
 */
class GenderCelListener : public ChipSummary, public CelListener, public GenderCalls {

public:

  /** 
   * Get the names for the cel files that have been seen.
   * @return - vector of cel file names, one for each cel file in order seen.
   */
  virtual std::vector<std::string> getCelNames() = 0;

};

#endif /* GENDERCELLISTENER_H */
