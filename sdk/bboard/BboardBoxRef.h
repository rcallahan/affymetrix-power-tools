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
// ~/Affy/apt2/trunk/BboardBoxRef.h ---
// 
// $Id: BboardBoxRef.h,v 1.4 2009-10-31 23:23:02 harley Exp $
// 

#ifndef _BBOARDBOXREF_H_
#define _BBOARDBOXREF_H_

#include "BboardTypes.h"
#include "Bboard.h"
#include "BboardBox.h"

class BboardBoxRef {

private:

  BboardBox* m_boxptr;

public:

  BboardBoxRef();
  BboardBoxRef(const BboardBox* boxptr);
  BboardBoxRef(const BboardBoxRef& boxref);
  virtual ~BboardBoxRef();

  //
  void init();
  void clear();
  void dump();

  //
  bool isNull();
  bool isNullRef();

  //
  const BboardBoxRef& assign(const BboardBox* valptr);
  const BboardBoxRef& assign(const BboardBoxRef& valref);

  //
  const BboardBoxRef& operator=(const BboardBox* vptr);
  const BboardBoxRef& operator=(const BboardBoxRef& vref);
  // dereferencing it
  BboardBox* operator->();

  //
  void swapBoxPtrs(BboardBoxRef& boxref);
  BboardBox* boxPtr();
  BboardType_t getBbType();
};

#endif // _BBOARDBOXREF_H_
