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
// affy/sdk/file5/File5_Matrix.cpp ---
// 
// $Id: File5_Matrix.cpp,v 1.20 2009-08-12 21:54:08 harley Exp $
// 

//
#include "file5/File5_Matrix.h"
//
#include "util/Err.h"

//
affx::File5_Matrix::File5_Matrix()
{
  m_kind_char='M';
  init();
}

//
affx::File5_Matrix::~File5_Matrix()
{
  close();
}

int affx::File5_Matrix::init() 
{
  File5_Object::init();
  //
  m_kind=FILE5_KIND_MATRIX;
  m_kind_char='M';
  //
  m_rank=0;
  //
  return 0;
}

int affx::File5_Matrix::flush() 
{
  if (m_state!=affx::FILE5_STATE_OPEN) {
    return 0;
  }
    
  return 0;
}

char affx::File5_Matrix::file5_kind_char() {
  return 'M';
}

/////

int
affx::File5_Matrix::open(const std::string& name,
                         const affx::File5_dtype_t& dtype,
                         const std::vector<int> dims,
                         int flags)
{
  close();
  //parent_refcnt_inc();
  //
  m_name=name;
  m_flags=flags;
  m_dtype=dtype;
  m_h5_dtype=as_hdf5_dtype(dtype);
  // copy from int to hsize_t
  affx::file5_dims_hsize_copy(m_dims,dims);
  m_rank=m_dims.size();
  FILE5_ASSERT(m_rank>0);

  if (getParent()->name_exists(name)==1) {
    FILE5_ABORT(": name exists: '"+name+"'");
  }

  //
  // affx::File5_Object::open();

  //
  hid_t h5_cparms;
  std::vector<hsize_t> chunk_dims;
  chunk_dims.resize(m_rank,100);
  h5_cparms=H5Pcreate(H5P_DATASET_CREATE);
  // chunk the data
  m_h5_status=H5Pset_chunk(h5_cparms,m_rank,&chunk_dims[0]);
  // turn on crc checks.
  if (m_opt_crc==1) {
    m_h5_status=H5Pset_fletcher32(h5_cparms);
  }

  std::vector<hsize_t> dims_max;
  dims_max.resize(m_rank,H5S_UNLIMITED);

  m_h5_dspace=H5Screate_simple(m_rank,&m_dims[0],&dims_max[0]);
  m_h5_obj=H5Dcreate(h5_parent_id(),
                      m_name.c_str(),
                      m_h5_dtype,
                      m_h5_dspace,
                      h5_cparms);
  FILE5_CHECKID(m_h5_obj,"H5Dcreate");

  //
  set_file5_kind_string("file5-matrix");
  
  //
  m_h5_status=H5Pclose(h5_cparms);
  m_dirty=1;
  m_state=affx::FILE5_STATE_OPEN;
  return 0;
}

int affx::File5_Matrix::close() 
{
  flush();
  //
  affx::File5_Object::close();
  //
  return 0;
}

int
affx::File5_Matrix::resize(const std::vector<int>& dims)
{
  // convert to the hsize_t
  std::vector<hsize_t> dims_cpy;
  affx::file5_dims_hsize_copy(dims_cpy,dims);
  //
  m_h5_status=H5Dextend(m_h5_obj,&dims_cpy[0]);
  FILE5_CHECKRV(m_h5_status,"H5Dextend");

  return 0;
}

//////////

int
affx::File5_Matrix::set_void(const std::vector<int>& dims,void* val_ptr)
{
  // the memory dspace
  FILE5_DIM1(m_dims,   1);
  FILE5_DIM1(m_offset, 0);
  FILE5_DIM1(m_count,  1);

  hid_t m_dspace=H5Screate_simple(1,m_dims,NULL);
  H5Sselect_hyperslab(m_dspace,H5S_SELECT_SET,m_offset,NULL,m_count,NULL);

  // the file dspace
  std::vector<hsize_t> f_offset;
  affx::file5_dims_hsize_copy(f_offset,dims);
  std::vector<hsize_t> f_count;
  f_count.resize(m_rank,1);

  hid_t f_dspace=H5Dget_space(m_h5_obj);
  H5Sselect_hyperslab(f_dspace,H5S_SELECT_SET,&f_offset[0],NULL,&f_count[0],NULL);
  //
  m_h5_status=H5Dwrite(m_h5_obj,m_h5_dtype,m_dspace,f_dspace,H5P_DEFAULT,val_ptr);
  FILE5_CHECKRV(m_h5_status,"H5Dwrite");
  //
  H5Sclose(m_dspace);
  H5Sclose(f_dspace);
  //
  return 0;
}

int
affx::File5_Matrix::set(const std::vector<int>& dims,int val)
{
  File5_type_pun_t pun;
  pun.i=val;
  return set_void(dims,(void*)&pun);
}
int
affx::File5_Matrix::set(const std::vector<int>& dims,float val)
{
  File5_type_pun_t pun;
  pun.f=val;
  return set_void(dims,(void*)&pun);
}
int
affx::File5_Matrix::set(const std::vector<int>& dims,double val)
{
  File5_type_pun_t pun;
  pun.d=val;
  return set_void(dims,(void*)&pun);
}

//////////

int
affx::File5_Matrix::get_void(const std::vector<int>& dims,void* val_ptr)
{
  // the memory dspace
  FILE5_DIM1(m_dims,   1);
  FILE5_DIM1(m_offset, 0);
  FILE5_DIM1(m_count,  1);

  hid_t m_dspace=H5Screate_simple(1,m_dims,NULL);
  H5Sselect_hyperslab(m_dspace,H5S_SELECT_SET,m_offset,NULL,m_count,NULL);

  // the file dspace
  std::vector<hsize_t> f_offset;
  affx::file5_dims_hsize_copy(f_offset,dims);
  std::vector<hsize_t> f_count;
  f_count.resize(m_rank,1);

  hid_t f_dspace=H5Dget_space(m_h5_obj);
  H5Sselect_hyperslab(f_dspace,H5S_SELECT_SET,&f_offset[0],NULL,&f_count[0],NULL);

  //
  m_h5_status=H5Dread(m_h5_obj,m_h5_dtype,m_dspace,f_dspace,H5P_DEFAULT,val_ptr);
  FILE5_CHECKRV(m_h5_status,"H5Dread");
  
  //
  return 0;
}

int
affx::File5_Matrix::get(const std::vector<int>& dims,int* val)
{
  return get_void(dims,(void*)val);
}
int
affx::File5_Matrix::get(const std::vector<int>& dims,float* val)
{
  return get_void(dims,(void*)val);
}
int
affx::File5_Matrix::get(const std::vector<int>& dims,double* val)
{
  return get_void(dims,(void*)val);
}
