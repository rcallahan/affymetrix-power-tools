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
//
// ~/Affy-projects/apt2/trunk/PnodeFactory.h ---
//
// $Id: PnodeFactory.h,v 1.1 2009-10-28 18:25:58 harley Exp $
//

#ifndef _PNODEFACTORY_H_
#define _PNODEFACTORY_H_

//
#include "bboard/Bboard.h"
#include "bboard/pnode/Pnode.h"
//
#include <string>

class PnodeFactory {
public:
  static Pnode* fromString(const std::string& name);
};

#endif // _PNODEFACTORY_H_
