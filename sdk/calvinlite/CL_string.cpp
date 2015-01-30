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
// CalvinLite/CL_string.cpp ---
//
// $Id: CL_string.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
//

//
#include "calvinlite/CL_string.h"
//
#include <assert.h>
#include <stdio.h>

CL_string::CL_string()
{
  clear();
}

void CL_string::clear()
{
  m_calvin_fmt=CL_STRING_NULL;
  //
  m_nstr_valid=true;
  m_nstr.clear();
  m_wstr_valid=true;
  m_wstr.clear();
}

void CL_string::ensureNstrValid()
{
  if (m_nstr_valid) {
    return;
  }
  if (m_wstr_valid) {
    WstrToNstr(m_wstr,m_nstr);
    m_nstr_valid=true;
    return;
  }
  clear();
}

void CL_string::ensureWstrValid()
{
  if (m_wstr_valid) {
    return;
  }
  if (m_nstr_valid) {
    NstrToWstr(m_nstr,m_wstr);
    m_wstr_valid=true;
    return;
  }
  clear();
}

std::string CL_string::Nstr()
{
  ensureNstrValid();
  return m_nstr;
}

std::string CL_string::Wstr()
{
  ensureWstrValid();
  return m_wstr;
}

const char* CL_string::c_str()
{
  ensureNstrValid();
  return m_nstr.c_str();
}

const char* CL_string::blob_str()
{
  ensureWstrValid();
  return m_wstr.c_str();
}

void CL_string::setNstr(const std::string& val)
{
  m_wstr.clear();
  m_wstr_valid=false;
  m_nstr=val;
  m_nstr_valid=true;
  m_calvin_fmt=CL_STRING_NSTR;
}

void CL_string::setWstr(const std::string& val)
{
  m_nstr.clear();
  m_nstr_valid=false;
  // val could be a real wide string or a fake wide string.
  // if val starts with a '0', then it is a real wide string.
  if ((m_wstr.size()>0) && (m_wstr[0]==0)) {
    m_wstr=val;
  }
  else {
    NstrToWstr(val,m_wstr);
  }
  //
  m_wstr_valid=true;
  m_calvin_fmt=CL_STRING_WSTR;
}

char* CL_string::resizeNstr(int size)
{
  m_wstr.clear();
  m_wstr_valid=false;
  m_nstr.resize(size,0);
  m_nstr_valid=true;
  m_calvin_fmt=CL_STRING_NSTR;
  return &m_nstr[0];
}

char* CL_string::resizeWstr(int size)
{
  m_nstr.clear();
  m_nstr_valid=false;
  m_wstr.resize(size);
  m_wstr_valid=true;
  m_calvin_fmt=CL_STRING_WSTR;
  return &m_wstr[0];
}

char* CL_string::resizeBlob(int size)
{
  m_nstr.clear();
  m_nstr_valid=false;
  m_wstr.resize(size);
  m_wstr_valid=true;
  //m_calvin_fmt
  return &m_wstr[0];
}

void CL_string::NstrToWstr(const std::string& nin,std::string& wout)
{
  wout.resize(nin.size()*2);
  for (int i=0;i<nin.size();i++) {
    int ii=i*2;
    wout[ii]=0;
    wout[ii+1]=nin[i];
  }
}

void CL_string::WstrToNstr(const std::string& win,std::string& nout)
{
  // we copy the all the non null chars to turn a wide
  // string into a narrow one.
  nout.clear();
  for (int i=0;i<win.size();i++) {
    char c=win[i];
    if (c!=0) {
      nout.push_back(c);
    }
  }
}

int CL_string::byteLen()
{
  if (m_calvin_fmt==CL_STRING_NULL) {
    return 0;
  }
  if (m_calvin_fmt==CL_STRING_NSTR) {
    ensureNstrValid();
    return m_nstr.size();
  }
  if (m_calvin_fmt==CL_STRING_WSTR) {
    ensureWstrValid();
    return m_wstr.size();
  }
  //
  assert(0);
  return -1;
}

void CL_string::dump() const
{
  printf("\n");
  printf("Nstr: %d : %3d : '%s'\n",int(m_nstr_valid),int(m_nstr.size()),m_nstr.c_str());
  printf("Wstr: %d : %3d : '",int(m_wstr_valid),int(m_wstr.size()));
  for (int i=0;i<m_wstr.size();i++) {
    printf("%02x ",m_wstr[i]);
  }
  printf("'\n");
}
