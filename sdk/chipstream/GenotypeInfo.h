////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
#ifndef _GENOTYPEINFO_H_
#define _GENOTYPEINFO_H_

#include "chipstream/GenderCalls.h"
#include "chipstream/InbredStatus.h"
#include "chipstream/SpecialSnps.h"

#include <map>
#include <string>

class GenotypeInfo {

public: 
  GenderCalls *m_Genders;
  InbredStatus *m_Inbred;
  std::map<std::string,bool> *m_HaploidSnps;
  SpecialSnpMap *m_SpecialSnps;
  CopyNumberMap *m_CopyNumberMap;
    
};

#endif /* _GENOTYPEINFO_H_ */

