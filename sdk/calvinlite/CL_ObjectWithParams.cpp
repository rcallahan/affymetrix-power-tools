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
// ~/CL/trunk/CL_ObjectWithParams.cpp ---
// 
// $Id: CL_ObjectWithParams.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
// 

#include "calvinlite/CL_ObjectWithParams.h"

CL_ObjectWithParams::CL_ObjectWithParams()
{
  clear();
}
CL_ObjectWithParams::~CL_ObjectWithParams()
{
  clear();
}
void CL_ObjectWithParams::clear()
{
  CL_Object::clear();
  clearParams();
}
void CL_ObjectWithParams::clearParams()
{
  m_params.clear();
}
int CL_ObjectWithParams::getParamCount() const
{
  return m_params.size();
}
CL_Param* CL_ObjectWithParams::getParam(int idx)
{
  if ((0<=idx)&&(idx<m_params.size())) {
    return &m_params[idx];
  }
  return NULL;
}
CL_Param* CL_ObjectWithParams::addParam(const CL_Param& param)
{
  m_params.push_back(param);
  return &m_params[m_params.size()-1];
}

CL_Param* CL_ObjectWithParams::addParam(const std::string& key,const std::string& val)
{
  return addParam(CL_Param(key,val));
}
CL_Param* CL_ObjectWithParams::addParam(const std::string& key,int val)
{
  return addParam(CL_Param(key,val));
}
CL_Param* CL_ObjectWithParams::addParam(const std::string& key,double val)
{
  return addParam(CL_Param(key,val));
}

CL_Param* CL_ObjectWithParams::addParam(const std::string& key,
                                        const std::string& type,
                                        const std::string& value,
                                        int byte_len)
{
  return addParam(CL_Param(key,type,value,byte_len));
}

int CL_ObjectWithParams::delParam(const std::string& key)
{
  int cnt=0;
  for (int i=0;i<m_params.size();i++) {
    if (key==m_params[i].m_name.Nstr()) {
      m_params.erase(m_params.begin()+i);
      i--; // step back
      cnt++;
    }
  }
  return cnt;
}
