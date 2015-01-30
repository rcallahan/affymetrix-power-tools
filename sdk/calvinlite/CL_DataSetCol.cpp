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
// ~/CalvinLite/CL_DataSetCol.cpp ---
//
// $Id: CL_DataSetCol.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
//

#include "calvinlite/CL_DataSetCol.h"
#include <stdio.h>

CL_DataSetCol::CL_DataSetCol()
{
  clear();
}

void CL_DataSetCol::clear()
{
  m_name.clear();
  m_type_code=CL_TC_UNSET;
  m_byte_len=0;
  m_byte_offset=0;
  m_src_cidx=-1;
}

void CL_DataSetCol::dump(const std::string& prefix_in)
{
  std::string prefix=prefix_in;

  CL_DUMP_PREFIX("DCOL ");
  CL_DUMP_MEMB_STR(m_name);
  CL_DUMP_MEMB_INT(m_type_code);
  CL_DUMP_MEMB_INT(m_byte_len);
  CL_DUMP_MEMB_INT(m_src_cidx);
}
