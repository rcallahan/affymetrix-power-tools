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
// ~/Affy-projects/apt2/trunk/PN_TestVecStats.h ---
//
// $Id: PN_TestVecStat.h,v 1.1 2009-10-28 18:25:58 harley Exp $
//

#ifndef _PN_TESTVECSTATS_H_
#define _PN_TESTVECSTATS_H_

#include "Pnode.h"

/// This file defines three Pnodes: PN_TestVecStat,
/// PN_TestVecStatSize and PN_TestVecStatSum

// this class is registered with PnodeFactory as "PN_TestVecStat".
// It can be created by name from the factory.
class PN_TestVecStat : public Pnode {
public:
  PN_TestVecStat();
  //
  virtual AptErr_t doRunNodePre(Bboard* bb);
  virtual AptErr_t doRunNodePost(Bboard* bb);
};

// these arent registered with PnodeFactory.
// but you can still make them directly.
class PN_TestVecStatSize : public Pnode {
public:
  PN_TestVecStatSize();
  //
  virtual AptErr_t doRunNode(Bboard* bb);
};
class PN_TestVecStatSum : public Pnode {
public:
  PN_TestVecStatSum();
  //
  virtual AptErr_t doRunNode(Bboard* bb);
};

#endif // _PN_TESTVECSTATS_H_
