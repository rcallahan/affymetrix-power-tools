////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

//
// affy/sdk/file5/File5_Matrix.h ---
//
// $Id: File5_Matrix.h,v 1.11 2009-09-25 17:49:03 mspald Exp $
//


#ifndef _FILE5_MATRIX_H_
#define _FILE5_MATRIX_H_

//
#include "file5/File5_File.h"
#include "file5/File5_types.h"
#include "file5/File5_util.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

//////////

class affx::File5_Matrix : public affx::File5_Object {
public:
  int m_rank;
  std::vector<hsize_t> m_dims;
  //
  File5_Matrix();
  File5_Matrix(affx::File5_File* file5_ptr);
  virtual ~File5_Matrix();
 
  //
  int init();
  int close();
  //int dump();
  int flush();
  virtual char file5_kind_char();

  int open(const std::string& name,
           const affx::File5_dtype_t& dtype,
           const std::vector<int> dims,
           int flags);

  //
  int resize(const std::vector<int>& dims);

  //
  int set_void(const std::vector<int>&dims,void* val);
  //
  int set(const std::vector<int>&dims,int val);
  int set(const std::vector<int>&dims,float val);
  int set(const std::vector<int>&dims,double val);

  //
  int get_void(const std::vector<int>&dims,void* val);
  //
  int get(const std::vector<int>&dims,int* val);
  int get(const std::vector<int>&dims,float* val);
  int get(const std::vector<int>&dims,double* val);
};

#endif // _FILE5_MATRIX_H_
