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
// ~/CL/trunk/CL_DataGroup.cpp ---
//
// $Id: CL_DataGroup.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
//

#include "calvinlite/CL_DataGroup.h"
#include "calvinlite/CL_private.h"
#include <stdio.h>

CL_DataGroup::CL_DataGroup()
{
  init();
  clear();
}
CL_DataGroup::~CL_DataGroup()
{
  clear();
}

void CL_DataGroup::init()
{
  CL_Object::init();
}

int CL_DataGroup::close()
{
  return CL_OK;
}

void CL_DataGroup::clear()
{
  delete_ptr_vec(m_datasets);
  //
  m_next_datagroup=NULL;
  m_next_datagroup_fpos=0;
  m_first_dataset_fpos=0;
  //
  m_name.clear();
  //
  CL_Object::clear();
}

CL_DataSet* CL_DataGroup::newDataSet(const std::string& ds_name)
{
  CL_DataSet* ds=new CL_DataSet();
  // this DataGroup remembers the dataset...
  m_datasets.push_back(ds);
  // ... and the dataset remember its parents.
  ds->m_parent_file=m_parent_file;
  ds->m_parent_group=this;
  // and its name
  ds->setName(ds_name);
  //
  return ds;
}

std::string CL_DataGroup::getName()
{
  return m_name.Nstr();
}
void CL_DataGroup::setName(const std::string& name)
{
  m_name.setNstr(name);
}

void CL_DataGroup::chainDataSets()
{
  int ds_cnt=m_datasets.size();
  for (int dsi=0;dsi<ds_cnt;dsi++) {
    m_datasets[dsi]->m_dataset_next=NULL;
    if (0<dsi) {
      m_datasets[dsi-1]->m_dataset_next=m_datasets[dsi];
    }
  }
}

CL_DataGroup* CL_DataGroup::nextDataGroup()
{
  return m_next_datagroup;
}

int CL_DataGroup::nextDataGroupFposStart()
{
  if (m_next_datagroup==NULL) {
    return 0;
  }
  return m_next_datagroup->getFposStart();
}

CL_DataSet* CL_DataGroup::firstDataSet()
{
  if (m_datasets.size()==0) {
    return NULL;
  }
  return m_datasets[0];
}

int CL_DataGroup::firstDataSetFposStart()
{
  if (m_datasets.size()==0) {
    return 0;
  }
  return m_datasets[0]->getFposStart();
}

void CL_DataGroup::dump(const std::string& prefix_in)
{
  std::string prefix=prefix_in;

  CL_DUMP_PREFIX("DG  ");
  // m_hdr_gdh.dump(prefix);
  //CL_DUMP_MEMB_INT(m_next_datagroup);
  //CL_DUMP_MEMB_INT(m_first_dataset_fpos);
  CL_DUMP_MEMB_STR(m_name);

  for (int i=0;i<m_datasets.size();i++) {
    m_datasets[i]->dump(prefix);
  }
}

int CL_DataGroup::getDataSetCount() const
{
  return m_datasets.size();
}

CL_DataSet* CL_DataGroup::getDataSet(int ds_idx) const
{
  if ((ds_idx<0)||(ds_idx>=m_datasets.size())) {
    return NULL;
  }
  return m_datasets[ds_idx];
}

CL_DataSet* CL_DataGroup::getDataSet(const std::string& ds_name) const
{
  for (int ds_idx=0;ds_idx<m_datasets.size();ds_idx++) {
    if (m_datasets[ds_idx]->m_name.Nstr()==ds_name) {
      // found the first matching name.
      return m_datasets[ds_idx];
    }
  }
  return NULL;
}
