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
// ~/CL/trunk/CL_ObjectWithParams.h ---
// 
// $Id: CL_ObjectWithParams.h,v 1.2 2009-10-29 22:28:59 harley Exp $
// 

#ifndef _CL_OBJECTWITHPARAMS_H_
#define _CL_OBJECTWITHPARAMS_H_

//
#include "calvinlite/CalvinLite.h"
//
#include "calvinlite/CL_Object.h"
#include "calvinlite/CL_Param.h"
//
#include <vector>

class CL_ObjectWithParams: public CL_Object {
public:
  CL_ObjectWithParams();
  ~CL_ObjectWithParams();
  //
  std::vector<CL_Param> m_params;
  //
  void clear();
  void clearParams();
  //
  int getParamCount() const;
  CL_Param* getParam(int idx);

  //
  CL_Param* addParam(const CL_Param& param);
  CL_Param* addParam(const std::string& key,const std::string& val);
  CL_Param* addParam(const std::string& key,int val);
  CL_Param* addParam(const std::string& key,double val);

  //
  CL_Param* addParam(const std::string& key,
                     const std::string& type,
                     const std::string& value,
                     int byte_size);

  // delete all the params with key and return the num deleted.
  int delParam(const std::string& key);
};

#endif // _CL_OBJECTWITHPARAMS_H_
