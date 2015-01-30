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
// affy/sdk/file5/File5Item.h ---
// 
// $Id: File5_Object.h,v 1.21 2009-09-25 17:49:03 mspald Exp $
// 

#ifndef _FILE5_OBJECT_H_
#define _FILE5_OBJECT_H_

// #define FILE5_DEBUG_USAGE 1

//
#include "file5/File5_types.h"
//
#include <cstring>
#include <string>
#include <vector>
//

namespace affx {
  // These functions start and stop the HDF5 library.
  void File5_open();
  void File5_close();

#ifdef FILE5_DEBUG_USAGE
  void file5_usecnt_inc(affx::File5_Object* ptr);
  void file5_usecnt_dec(affx::File5_Object* ptr);
  void file5_usecnt_print();
#endif
}

class affx::File5_Object {
private:
  // which group do we belong to?
  affx::File5_Group* m_parent;
  affx::File5_File*  m_file;

public:
  //
  std::string   m_name;
  int           m_errno;
  int           m_flags;
  int           m_dirty;
  int           m_readonly;
  //
  File5_kind_t  m_kind;
  File5_state_t m_state;

  // a char to mark the type with.
  char          m_kind_char;

  hid_t m_h5_obj;
  hid_t m_h5_dspace;
  herr_t m_h5_status;
  //
  affx::File5_dtype_t m_dtype;
  hid_t m_h5_dtype;
  hid_t m_h5_dtype_tofree;
  size_t m_h5_dtype_size;

  //
  int           m_refcnt;
  std::vector<affx::File5_Object*> m_refcnt_vec;

  // options
  int m_opt_compress;
  int m_opt_chunksize;
  int m_opt_crc;

  //
  File5_Object();
  virtual ~File5_Object();

  //
  affx::File5_return_t open(const std::string& filename,int flags);
  affx::File5_return_t init();
  affx::File5_return_t dump();
  affx::File5_return_t dump1();
  affx::File5_return_t flush();
  affx::File5_return_t close();
  //
  bool is_open();
  //
  bool ReadOnly();
  bool setReadOnly(bool val);

  //
  void setOptChunkSize(int size);
  void setOptCompress(int level);
  void setOptCrc(int state);

  //
  void setParent(affx::File5_Group* parent);
  affx::File5_Group* getParent() { return m_parent; };
  void setFile(affx::File5_File* file) { m_file=file; };
  affx::File5_File*  getFile() { return m_file; };

  //
  std::string file5_kind_string();
  int set_file5_kind_string(const std::string& kind_str);
  virtual char file5_kind_char();

  //
  int refcnt();
  int refcnt_inc(affx::File5_Object* referrer);
  int refcnt_dec(affx::File5_Object* referrer);
  int refcnt_has_referrer(affx::File5_Object* referrer);
  int refcnt_print();
  int refcnt_print(int indent);

  //
  hid_t h5_parent_id();
  int parent_refcnt_inc();
  int parent_refcnt_dec();
  //int set_parent(affx::File5_Object* parent);

  //
  int find_attrib_idx(const std::string& attrib_name);
  
  //
  int attrib_set_void(const std::string& attrib_name,hid_t h5_dtype,void* val);
  //
  int attrib_set(const std::string& attrib_name,int val);
  int attrib_set(const std::string& attrib_name,size_t val);
  int attrib_set(const std::string& attrib_name,double val);
  int attrib_set(const std::string& attrib_name,const std::string& val);
  //
  int attrib_get_void(const std::string& attrib_name,hid_t h5_dtype,void* val);
  //
  int attrib_get(const std::string& attrib_name,int* val);
  int attrib_get(const std::string& attrib_name,size_t* val);
  int attrib_get(const std::string& attrib_name,double* val);
  int attrib_get(const std::string& attrib_name,std::string* val);
  //
  int attrib_delete(const std::string& attrib_name);

  //
  int get_objinfo(H5G_stat_t* statbuf);
  int get_name(std::string* name);

  //
  static void close_oid(hid_t* oid);
};

#endif // _FILE5_OBJECT_H_
