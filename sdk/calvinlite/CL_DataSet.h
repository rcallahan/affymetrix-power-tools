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
// CalvinLite/CL_DataSet.h ---
// 
// $Id: CL_DataSet.h,v 1.2 2009-10-29 22:28:59 harley Exp $
// 

#ifndef _CALVINLITE_DATASET_H_
#define _CALVINLITE_DATASET_H_

//
#include "calvinlite/CalvinLite.h"
//
#include "calvinlite/CL_DataGroup.h"
#include "calvinlite/CL_DataSetCol.h"
#include "calvinlite/CL_Gdh.h"
#include "calvinlite/CL_Param.h"
#include "calvinlite/CL_ObjectWithParams.h"
#include "calvinlite/CL_string.h"
//
#include <string>
#include <vector>

class CL_DataSet : public CL_ObjectWithParams {
public:
  // The DataGroup which contains this DataSet
  CL_DataGroup* m_parent_group;
  //
  int m_fpos_data_start;
  int m_fpos_data_end;
  // this might be used later to force padding.
  // int m_fpos_pad_bytesize;
  //
  int m_next_dataset_fpos;
  CL_string m_name;
  //
  std::vector<CL_DataSetCol> m_col_vec;
  //
  int m_row_cnt;
  mutable int m_row_bytesize;
  //
  CL_DataSet* m_pushed_dataset;
  //
  CL_DataSet* m_dataset_next;
  //
  int m_buffer_size;
  void* m_buffer_ptr;

  //
  CL_DataSet();
  ~CL_DataSet();
  //
  void init();
  void clear();
  int  close();
  void dump(const std::string& prefix_in);

  std::string getName();
  void setName(const std::string& name);

  //
  void setParentGroup(CL_DataGroup* cl_group);
  CL_DataGroup* getParentGroup();

  //
  void pushDataInto(CL_DataSet* dstset);
  void maybePushData();
  void copyDataFrom(CL_DataSet* src);
  static void copyData(CL_DataSet* dst,CL_DataSet* src);
  void maybePullData();

  //
  int colCount() const;
  int rowCount() const;
  int setRowCount(int rcnt);
  int rowBytesize();
  int computeRowBytesize();
  //
  CL_DataSetCol* getColPtr(int idx);
  void resizeCols(int cnt);
  int getColName(int cidx,std::string* colName);

  // @todo rename to "defineColumn"
  int newColumn(const std::string& cname,CL_TypeCode_t ctype);
  int newColumn(const std::string& cname,CL_TypeCode_t ctype,int bytelen);

  //
  void* getBufferPtr();
  int getBufferSize();
  void* allocBuffer();
  void freeBuffer();

  //
  unsigned char* getValuePtr(int row,int col);

  //
  int get(int row,int col,int* val);
  int get(int row,int col,float* val);
  int get(int row,int col,std::string* val);
  int getAsString(int ridx,int cidx,std::string* val);

  //
  CL_Err_t set(int ridx,int cidx,int   val);
  CL_Err_t setFloat(int ridx,int cidx,float val);
  CL_Err_t set(int ridx,int cidx,const std::string& val);

  //
  CL_Err_t setFromArray(int ridx_start,int cidx,char* arrayData,int arrayLen);
  CL_Err_t setFromArray(int ridx_start,int cidx,int* arrayData,int arrayLen);
  CL_Err_t setFromArray(int ridx_start,int cidx,float* arrayData,int arrayLen);

  //
  int getDataFposStart() const;
  int setDataFposStart(int fpos);
  int getDataFposEnd() const;
  int setDataFposEnd(int fpos);
};

#endif
