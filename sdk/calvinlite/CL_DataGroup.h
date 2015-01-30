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
// CalvinLite/CL_DataGroup.h ---
// 
// $Id: CL_DataGroup.h,v 1.2 2009-10-29 22:28:59 harley Exp $
// 

#ifndef _CALVINLITE_DATAGROUP_H_
#define _CALVINLITE_DATAGROUP_H_

//
#include "calvinlite/CalvinLite.h"
//
#include "calvinlite/CL_DataSet.h"
#include "calvinlite/CL_Gdh.h"
#include "calvinlite/CL_Object.h"
#include "calvinlite/CL_string.h"
//
#include <vector>

class CL_DataGroup : public CL_Object {
public:
  // in file order
  CL_DataGroup* m_next_datagroup;
  int m_next_datagroup_fpos;
  int m_first_dataset_fpos;
  CL_string m_name;
  //
  std::vector<CL_DataSet*> m_datasets;
  
  //
  CL_DataGroup();
  ~CL_DataGroup();
  //
  void init();
  void clear();
  int close();
  void dump(const std::string& prefix_in);

  void setName(const std::string& name);
  std::string getName();

  CL_DataSet* newDataSet(const std::string& ds_name);

  //
  void chainDataSets();
  CL_DataGroup* nextDataGroup();
  int nextDataGroupFposStart();
  CL_DataSet* firstDataSet();
  int firstDataSetFposStart();

  //
  int getDataSetCount() const;
  //
  CL_DataSet* getDataSet(int ds_idx) const;
  CL_DataSet* getDataSet(const std::string& ds_name) const;
};

#endif
