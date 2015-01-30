////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

/**
 * @file   QuantSea.cpp
 * @author Pete Klosterman
 * @date   Fri Jun 23 12:48:12 2006
 *
 * @brief  Class for doing probe set quantification using the SEA, or
 * Simplified Expression Analysis, method.
 */

#include "chipstream/QuantSea.h"

using namespace std;

/** Constructor. */
QuantSea::QuantSea(QuantPlierParams &plierParams)
  : QuantPlier(plierParams) {
  setupSelfDoc(*this);
  m_Type = getDocName();
}
