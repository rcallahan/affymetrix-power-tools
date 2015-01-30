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
// affy/sdk/file5/File5_util.h ---
// 
// $Id: File5_util.h,v 1.13 2009-09-18 03:37:27 mspald Exp $
// 

#ifndef _FILE5_UTIL_H_
#define _FILE5_UTIL_H_

//
#include "file5/File5_types.h"
//
#include "file/TsvFile/TsvFile.h"
//
#include <vector>

// We are dead.
#define FILE5_ABORT(xmsg) \
  { \
    Err::errAbort("File5: " __FILE__ ":"+ToStr(__LINE__)+" "+ToStr(xmsg)); \
  }
// xcond must be true or we are dead.
#define FILE5_ASSERT(xcond) \
  { \
    if (!(xcond)) { \
      FILE5_ABORT(#xcond); \
    } \
  }
// the return value must be "0" or we are dead.
#define FILE5_CHECKRV(xrv,xmsg) \
  { \
  int tmp_xrv=xrv; \
  if (tmp_xrv!=0) { \
    FILE5_ABORT("rv="+ToStr(tmp_xrv)+": "+xmsg); \
  } \
}
// the ID must be positive.
#define FILE5_CHECKID(xid,xmsg) \
{ \
  int tmp_xid=xid; \
  if (tmp_xid<0) { \
    FILE5_ABORT("id="+ToStr(tmp_xid)+": "+xmsg);  \
  } \
}

//
#define FILE5_DIM1(NAME,D0) \
  hsize_t NAME[1]; \
  NAME[0]=D0;

#define FILE5_DIM2(NAME,D0,D1) \
  hsize_t NAME[2]; \
  NAME[0]=D0; \
  NAME[1]=D1;

#define FILE5_DIM3(NAME,D0,D1,D2) \
  hsize_t NAME[3]; \
  NAME[0]=D0; \
  NAME[1]=D1; \
  NAME[2]=D2;

//
namespace affx {
  //
  typedef union File5_type_pun {
    int i;
    float f;
    double d;
  } File5_type_pun_t;

  //
  std::vector<int> file5_dims(int d0);
  std::vector<int> file5_dims(int d0,int d1);
  std::vector<int> file5_dims(int d0,int d1,int d2);
  std::vector<int> file5_dims(int d0,int d1,int d2,int d3);
  //
  void file5_dims_hsize_copy(std::vector<hsize_t>& dst,const std::vector<int>& src);
  //
  affx::File5_dtype_t as_file5_dtype(affx::tsv_type_t tsv_dtype);
  affx::File5_dtype_t as_file5_dtype(hid_t hdf5_dtype);
  //
  std::string file5_dtype_to_string(affx::File5_dtype_t dtype);
  affx::File5_dtype_t file5_string_to_dtype(const std::string& str);
  //
  hid_t as_hdf5_dtype(affx::tsv_type_t tsv_dtype);
  hid_t as_hdf5_dtype(affx::File5_dtype_t file5_dtype);
  hid_t as_h5t_native_type(hid_t dtype);
  
  //
  void dump_h5_dset(hid_t dset_id,std::ostream& ostm=std::cout);
  void dump_h5_dspace(hid_t dspace_id,std::ostream& ostm=std::cout);
  void dump_h5_dtype(hid_t dtype_id,std::ostream& ostm=std::cout);
};

#endif // _FILE5_UTIL_H_
