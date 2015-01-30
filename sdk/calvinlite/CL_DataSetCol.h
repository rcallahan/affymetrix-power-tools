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
// ~/CalvinLite/CL_DataSetCol.h ---
//
// $Id: CL_DataSetCol.h,v 1.2 2009-10-29 22:28:59 harley Exp $
//

#ifndef _CL_DATASETCOL_H_
#define _CL_DATASETCOL_H_

//
#include "calvinlite/CL_string.h"
//
#include <string>

class CL_DataSetCol {
public:
  CL_string m_name;
  int m_type_code;
  int m_byte_len;
  //
  int m_src_cidx;
  //
  int m_byte_offset;

  //
  CL_DataSetCol();
  void clear();
  void dump(const std::string& prefix_in);
};

#endif // _CL_DATASETCOL_H_
