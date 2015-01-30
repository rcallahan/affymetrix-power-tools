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
// ~/CL/trunk/CL_Param.h ---
// 
// $Id: CL_Param.h,v 1.2 2009-10-29 22:28:59 harley Exp $
// 

#ifndef _CL_PARAM_H_
#define _CL_PARAM_H_

//
#include "calvinlite/CalvinLite.h"
#include "calvinlite/CL_string.h"
#include "calvinlite/CL_Object.h"
//
#include <string>

class CL_Param : public CL_Object {
public:
  int m_start_fpos;
  //
  CL_string m_name;
  //
  int m_val_int;
  double m_val_double;
  CL_string m_val_string;
  // the code
  CL_TypeCode_t  m_type_code;
  // the *suggested* file length.
  int m_suggested_byte_len;
  
  //
  CL_Param();
  CL_Param(const std::string& key,const std::string& val);
  CL_Param(const std::string& key,int val);
  CL_Param(const std::string& key,double val);
  //
  CL_Param(const std::string& key,
           const std::string& type,
           const std::string& value,
           int byte_size);
  //
  void clear();
  //
  std::string valueAsString();
  CL_Err_t setValueFromString(const std::string& value);
  int valueByteLen();
  int suggestedByteLen();
  //
  void dump(const std::string& prefix_in);
};

#endif // _CL_PARAM_H_
