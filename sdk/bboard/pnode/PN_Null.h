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
// ~/Affy-projects/apt2/trunk/PN_Null.h ---
//
// $Id: PN_Null.h,v 1.1 2009-10-28 18:25:58 harley Exp $
//

#ifndef _PN_NULL_H_
#define _PN_NULL_H_

#include "Pnode.h"

/// This is a null Pnode. It doesnt do anything, but 
/// what is inherited from "Pnode".
class PN_Null : public Pnode {
public:
  // 
  PN_Null();
  virtual ~PN_Null();
};

#endif // _PN_NULL_H_
