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
// CalvinLite/CL_Gdh.cpp ---
//
// $Id: CL_Gdh.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
//

//
#include "calvinlite/CL_Gdh.h"
//
#include "calvinlite/CL_private.h"
#include "calvinlite/CL_util.h"
//
#include "util/Convert.h"
//
#include <assert.h>
#include <stdio.h>

//////////

CL_Gdh::CL_Gdh()
{
  CL_ObjectWithParams::clear();
  defaultValues();
}
CL_Gdh::~CL_Gdh()
{
  clear();
}
void CL_Gdh::clear()
{
  delete_ptr_vec(m_gdh_leaves);
}

void CL_Gdh::defaultValues()
{
  // give this gdh sensible defaults.
  m_locale_str.setNstr("en-US");
  m_uuid_str.setNstr(CL_gen_uuid());
  m_datetime_str.setNstr(CL_gen_timestr());
}

void CL_Gdh::setDataType(const std::string& val) {
  m_datatype_str.setNstr(val);
}
void CL_Gdh::setUuid(const std::string& val) {
  m_uuid_str.setNstr(val);
}
void CL_Gdh::setDateTime(const std::string& val) {
  m_datetime_str.setWstr(val);
}
void CL_Gdh::setLocale(const std::string& val) {
  m_locale_str.setWstr(val);
}

//////////

CL_Gdh* CL_Gdh::newLeafGdh()
{
  CL_Gdh* leaf=new CL_Gdh();
  m_gdh_leaves.push_back(leaf);
  return leaf;
}

int CL_Gdh::getLeafCount()
{
  return m_gdh_leaves.size();
}

CL_Gdh* CL_Gdh::getLeafGdh(int i)
{
  assert((i<=0)&&(i<m_gdh_leaves.size()));
  return m_gdh_leaves[i];
}

CL_Gdh* CL_Gdh::getGdhByPath(std::vector<std::string>& path_vec) {
  // ourselves?
  if ((path_vec.size()==0)||(path_vec[0]=="")||(path_vec[0]==".")) {
    return this;
  }
  //
  int idx=Convert::toInt(path_vec[0]);

  // stick on leaves
  while (getLeafCount()<=idx) {
    newLeafGdh();
  }
  CL_Gdh* leaf=getLeafGdh(idx);
  
  path_vec.erase(path_vec.begin());
  return leaf->getGdhByPath(path_vec);
}

//////////

void CL_Gdh::dump()
{
  dump("");
}

void CL_Gdh::dump(const std::string& prefix_in)
{
  std::string prefix=prefix_in;

  CL_DUMP_PREFIX("GDH ");
  CL_DUMP_MEMB_INT(m_fpos_start);
  CL_DUMP_MEMB_INT(m_fpos_size);
  CL_DUMP_MEMB_STR(m_datatype_str);
  CL_DUMP_MEMB_STR(m_uuid_str);
  CL_DUMP_MEMB_STR(m_datetime_str);
  CL_DUMP_MEMB_STR(m_locale_str);

  if (1) {
    for (int i=0;i<m_params.size();i++) {
      m_params[i].dump(prefix);
    }
  }
  if (1) {
    for (int i=0;i<m_gdh_leaves.size();i++) {
      m_gdh_leaves[i]->dump(prefix);
    }
  }
}

//////////

static std::string prefix_init(const std::string& prefix_in) {
  if (prefix_in=="") {
    return "0";
  }
  return prefix_in;
}

static std::string prefix_idx(const std::string& prefix_in,int idx)
{
  char buf[20];
  sprintf(buf,"%d",idx);
  return prefix_in+":"+buf;
}

//////////

#define GDH_PRINT(_key,_val) { printf("%s | %s='%s'\n",prefix.c_str(),_key,_val.c_str()); }

void CL_Gdh::print(const std::string& prefix_in)
{
  std::string prefix=prefix_in;
  std::string prefix_out;

  prefix=prefix_init(prefix_in);

  GDH_PRINT("datatype",m_datatype_str);
  GDH_PRINT("uuid"    ,m_uuid_str);
  GDH_PRINT("date"    ,m_datetime_str);
  GDH_PRINT("locale"  ,m_locale_str);

  for (int i=0;i<m_params.size();i++) {
    CL_Param* param=getParam(i);
    GDH_PRINT(param->m_name.c_str(),param->valueAsString());
  }

  if (1) {
    for (int i=0;i<m_gdh_leaves.size();i++) {
      prefix_out=prefix_idx(prefix,i);
      m_gdh_leaves[i]->print(prefix_out);
    }
  }
}

//////////

void CL_Gdh::printMatchingKeys(const std::vector<std::string>& key_vec,
                               const std::string& prefix_in)
{
  std::set<std::string> key_set;
  key_set.insert(key_vec.begin(),key_vec.end());

  printMatchingKeys(key_set,prefix_in);
}

#define GDH_PRINT_IF_MATCH(_key_set,_key,_val) { \
  if (key_set.find(_key)!=key_set.end()) { \
    GDH_PRINT(_key,_val); \
  } \
  }

void CL_Gdh::printMatchingKeys(const std::set<std::string>& key_set,
                               const std::string& prefix_in)
{
  std::string prefix_out;
  std::string prefix=prefix_init(prefix_in);

  GDH_PRINT_IF_MATCH(key_set,"datatype",m_datatype_str);
  GDH_PRINT_IF_MATCH(key_set,"uuid"    ,m_uuid_str);
  GDH_PRINT_IF_MATCH(key_set,"date"    ,m_datetime_str);
  GDH_PRINT_IF_MATCH(key_set,"locale"  ,m_locale_str);

  for (int i=0;i<m_params.size();i++) {
    CL_Param* param=getParam(i);
    if (key_set.find(param->m_name.c_str())!=key_set.end()) {
      GDH_PRINT(param->m_name.c_str(),param->valueAsString());
    }
  }

  if (1) {
    for (int i=0;i<m_gdh_leaves.size();i++) {
      prefix_out=prefix_idx(prefix,i);
      m_gdh_leaves[i]->printMatchingKeys(key_set,prefix_out);
    }
  }
}
