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
// ~/CL/trunk/CL_private.h ---
// 
// $Id: CL_private.h,v 1.1 2009-10-27 16:53:53 harley Exp $
// 

#ifndef _CL_PRIVATE_H_
#define _CL_PRIVATE_H_

// stuff which doesnt need to be seen by users of CalvinLite.

// delete a vector of pointers.
template<typename T1>
void delete_ptr_vec(std::vector<T1>& vec) {
  for (int i=0;i<vec.size();i++) {
    delete vec[i];
    vec[i]=NULL;
  }
  vec.clear();
}

#endif // _CL_PRIVATE_H_
