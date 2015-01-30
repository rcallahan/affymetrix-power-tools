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
// ~/Affy-projects/apt2/trunk/PN_Null.cpp ---
//
// $Id: PN_Null.cpp,v 1.1 2009-10-28 18:25:58 harley Exp $
//

#include "PN_Null.h"

/// say hello when we are created.
PN_Null::PN_Null() {
  printf("new PN_Null()  %p\n",this);
  init();
}

/// be sure to call doDeleteChildren() in the destructor if you write one.
PN_Null::~PN_Null() {
  doDeleteChildren();
};
