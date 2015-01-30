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
// ~/Affy-projects/apt2/trunk/PnodeFactory.cpp ---
//
// $Id: PnodeFactory.cpp,v 1.1 2009-10-28 18:25:58 harley Exp $
//

//
#include "bboard/pnode/PnodeFactory.h"
//
#include "bboard/pnode/PN_Null.h"
#include "bboard/pnode/PN_TestVecStat.h"
//
#include <assert.h>

// This translates names to newly created Pnodes.
// Not all Pnodes need to be listed here, just the ones we want to
// create on the fly.
Pnode* PnodeFactory::fromString(const std::string& name) {
  Pnode* pnode;
  //
  if (name=="PN_Null") {
    pnode=new PN_Null();
    return pnode;
  }
  //
  if (name=="PN_TestVecStat") {
    pnode=new PN_TestVecStat();
    return pnode;
  }
  //
  printf("PnodeFactory::fromString('%s')==NULL\n",name.c_str());
  assert(0);
  return NULL;
}
