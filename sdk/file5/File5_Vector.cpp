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
// affy/sdk/file5/File5_Vector.cpp ---
//
// $Id: File5_Vector.cpp,v 1.48 2009-11-04 19:07:38 harley Exp $
//

//
#include "file5/File5_Vector.h"
//
#include "file5/File5.h"
//
#include "util/Convert.h"
//
#include <cstring>
#include <string.h>
//

#define VECTOR_FILL_INDEX "fill-index"

//

affx::File5_Vector::File5_Vector()
{
  m_kind_char='V';
  init();
}

affx::File5_Vector::~File5_Vector()
{
  close();
}

int affx::File5_Vector::init() {
  affx::File5_Object::init();
  //
  m_kind=affx::FILE5_KIND_VECTOR;
  m_kind_char='V';
  //
  m_vec_end_idx=0;
  m_vec_fill_idx=0;
  //
  m_buf_ptr=NULL;
  //
  m_opt_autoresize=0;
  m_opt_buffer_size=-1;
  m_opt_string_size=-1; // default to unlimited strings
  m_opt_resize_on_flush=0;
  //
  m_h5_dtype=-1;
  m_h5_dtype_size=0;
  m_dtype=affx::FILE5_DTYPE_ERR;
  m_dtype_size=0;
  //
  buffer_free();
  //
  return 0;
}

//
affx::File5_Vector* affx::File5_Vector::allocate(affx::File5_dtype_t dtype)
{
  affx::File5_Vector* vec=new affx::File5_Vector();
  vec->m_dtype=dtype;
  return vec;
}

char affx::File5_Vector::file5_kind_char() {
  return 'V';
}

/////

void affx::File5_Vector::setOptAutoResize(int flag)
{
  m_opt_autoresize=flag;
}
void affx::File5_Vector::setOptStringSize(int size)
{
  m_opt_string_size=size;
}

int affx::File5_Vector::getOptStringSize()
{
  return m_opt_string_size;
}

/////

int affx::File5_Vector::close()
{
  flush();
  // truncate the file to the vector length.
  resizeFile(m_vec_end_idx);
  //
  buffer_free();
  //
  affx::File5_Object::close();
  //
  return 0;
}

#define F5V_dump(_memb) { printf("  %30s = %d\n",#_memb,(int)_memb); }
int affx::File5_Vector::dump()
{
  printf("File5_Vector(%p)==\n",this);
  F5V_dump(m_dtype);
  F5V_dump(m_vec_fill_idx);
  F5V_dump(m_file_end_idx);
  F5V_dump(m_buf_start_idx);
  F5V_dump(m_buf_end_idx);
  F5V_dump(m_buf_max_cnt);
  //
  return FILE5_OK;
}

void affx::File5_Vector::check_sanity() {
  // check the sanity of our vector.
  FILE5_ASSERT(m_buf_start_idx<=m_buf_end_idx);
  FILE5_ASSERT(0<=m_vec_end_idx);
  FILE5_ASSERT(0<=m_file_end_idx);
  FILE5_ASSERT(m_vec_fill_idx<=m_vec_end_idx);
  FILE5_ASSERT(0<=m_buf_max_cnt);
  //FILE5_ASSERT(m_buf_start_idx+m_buf_cnt<=m_vec_end_idx);
}

int affx::File5_Vector::flush()
{
  if ((m_state!=affx::FILE5_STATE_OPEN) || (m_dirty==0)) {
    return 0;
  }

  check_sanity();

  // affx::dump_h5_dset(m_h5_obj); //dbg

  //  printf("flush() = m_vec_fill_idx=%lu m_vec_end_idx=%lu m_file_end_idx=%lu\n",
  //         m_vec_fill_idx,
  //         m_vec_end_idx,
  //         m_file_end_idx);

  // make sure the file is at least m_vec_end_idx long.
  if (m_opt_resize_on_flush==1) {
    resizeFile(m_vec_end_idx);
  } else {
    reserveFile(m_vec_end_idx);
  }

  // start should always be less than end.
  FILE5_ASSERT(m_buf_start_idx<=m_buf_end_idx);

  // write the buffer contents
  size_t cnt=m_buf_end_idx-m_buf_start_idx;

  // nothing to do.
  if (cnt==0) {
    return 0;
  }

  // put the data out.
  write_array_io(m_buf_start_idx,cnt,m_buf_ptr);
  // keep the meta-data up to date.
  attrib_set(VECTOR_FILL_INDEX,m_vec_fill_idx);
  m_dirty=0;

  //
  getFile()->flush();

  return 0;
}

//////////

size_t affx::File5_Vector::resizeFile(size_t size)
{
  if (m_state!=affx::FILE5_STATE_OPEN) {
    return 0;
  }
  // printf("resizeFile(%lu)\n",size); // dbg
  //
  FILE5_DIM1(dims,size);
  m_h5_status=H5Dextend(m_h5_obj,dims);
  FILE5_CHECKRV(m_h5_status,"H5Dextend");
  //
  m_file_end_idx=size;
  //if (m_file_end_idx<m_vec_end_idx) {
  //resize(m_file_end_idx);
  //}
  //
  return size;
}

size_t affx::File5_Vector::reserveFile(size_t size)
{
  if (size>m_file_end_idx) {
    resizeFile(size);
  }
  return size;
}

size_t affx::File5_Vector::resize(size_t size)
{
  // put the data out.
  flush();
  // adjust the size of the vector.
  m_vec_end_idx=size;
  // clamp the fill pointer to the end.
  if (m_vec_fill_idx>m_vec_end_idx) {
    m_vec_fill_idx=m_vec_end_idx;
  }
  // mark us as dirty as, our size is different.
  m_dirty=1;
  //
  return size;
}

size_t affx::File5_Vector::setFillIndex(size_t idx)
{
  flush();
  //
  m_dirty=1;
  //
  m_vec_fill_idx=idx;
  // check for in bounds.
  if (m_vec_fill_idx>m_vec_end_idx) {
    FILE5_ABORT("Fill index is past vector end.");
  }
  return idx;
}

//////////

void affx::File5_Vector::setBufferSize(int size)
{
  FILE5_ASSERT(size>0);
  buffer_alloc(size);
}

void affx::File5_Vector::buffer_alloc(int cnt)
{
  buffer_free();
  //
  FILE5_ASSERT(m_dtype_size!=0);
  //
  m_buf_max_cnt=cnt;
  m_buf_ptr=(char*)malloc(cnt*m_dtype_size);
  FILE5_ASSERT(m_buf_ptr!=NULL);
  // zero it out, we know it to be empty.
  buffer_zero();
}

void affx::File5_Vector::buffer_free()
{
  flush();
  //
  if (m_buf_ptr!=NULL) {
    free(m_buf_ptr);
    m_buf_ptr=NULL;
  }
  //
  //m_buf_cnt=0;
  m_buf_max_cnt=0;
  m_buf_start_idx=0;
  m_buf_end_idx=0;
}

// dont call this: use buffer_clear() instead.
// Pointers need to be freed first.
void affx::File5_Vector::buffer_zero()
{
  if (m_buf_ptr!=NULL) {
    memset(m_buf_ptr,0,m_buf_max_cnt*m_dtype_size);
  }
}

void affx::File5_Vector::buffer_clear() {
  // was this varlen data?
  if ((m_dtype==affx::FILE5_DTYPE_STRING)&&(m_dtype_size==-1)) {
    // yes, reclaim it.
    FILE5_DIM1(m_dims, (m_buf_end_idx-m_buf_start_idx));
    hid_t m_dspace=H5Screate_simple(1,m_dims,NULL);
    H5Sselect_all(m_dspace);
    H5Dvlen_reclaim(m_h5_dtype,m_dspace,H5P_DEFAULT,&m_buf_ptr[m_buf_start_idx]);
  }
  //
  buffer_zero();
}

void affx::File5_Vector::buffer_print()
{
  char* ptr_c;

  for (int i=0;i<m_buf_max_cnt;i++) {
    printf("%3d: ",i);
    //
    switch (m_dtype) {
    case affx::FILE5_DTYPE_CHAR:
      printf("%d",(int)((unsigned char*)m_buf_ptr)[i]);
      break;
    case affx::FILE5_DTYPE_SHORT:
      printf("%d",(int)((short*)m_buf_ptr)[i]);
      break;
    case affx::FILE5_DTYPE_INT:
      printf("%d",((int*)m_buf_ptr)[i]);
      break;
    case affx::FILE5_DTYPE_FLOAT:
      printf("%f",((float*)m_buf_ptr)[i]);
      break;
    case affx::FILE5_DTYPE_DOUBLE:
      printf("%f",((float*)m_buf_ptr)[i]);
      break;
    case affx::FILE5_DTYPE_STRING:
      if (m_opt_string_size==-1) {
        ptr_c=(((char**)m_buf_ptr)[i]);
        printf("%p => '%s'",ptr_c,ptr_c);
      }
      else {
        ptr_c=&(((char*)m_buf_ptr)[i*m_dtype_size]);
        printf("'%s'",ptr_c);
      }
      break;
    default:
      printf("?");
      break;
    }
    //
    printf(" %s %s %s\n",
           (((i+m_buf_start_idx)==m_vec_fill_idx)?"<=v-fill":""),
           (((i+m_buf_start_idx)==m_vec_end_idx)?"<=v-end":""),
           (((i+m_buf_start_idx)==m_file_end_idx)?"<=f-end":"")
           );
    //
    if (i>=11) {
      break;
    }
  } // for
  printf("v-fill=%d  v-end=%d   f-end=%d\n",
         int(m_vec_fill_idx),
         int(m_vec_end_idx),
         int(m_file_end_idx));
}

void affx::File5_Vector::buffer_seek(size_t idx)
{
  // printf("buffer_seek(%lu)\n",idx); // dbg
  //
  if (m_vec_end_idx<idx) {
    FILE5_ABORT("buffer_seek(idx) is past end.");
  }
  //
  FILE5_ASSERT(0<m_buf_max_cnt);
  FILE5_ASSERT(m_buf_ptr!=NULL);
  check_sanity();
  //
  // affx::dump_h5_dset(m_h5_obj); // dbg
  //
  flush();
  //
  buffer_clear();
  //
  if (idx>m_vec_end_idx) {
    FILE5_ABORT("beyond vector end.");
  }

  // set the start
  m_buf_start_idx=idx;
  // now set where we think the end should be.
  m_buf_end_idx=m_buf_start_idx+m_buf_max_cnt;
  // trim the end of the buffer to the vector...
  if (m_buf_end_idx>m_vec_end_idx) {
    m_buf_end_idx=m_vec_end_idx;
  }
  // ...and to the file.
  //if (m_buf_end_idx>m_file_end_idx) {
  //m_buf_end_idx=m_file_end_idx;
  //}
  // set the count.
  int cnt=m_buf_end_idx-m_buf_start_idx;
  if (idx+cnt>m_file_end_idx) {
    cnt=m_file_end_idx-idx;
  }
  //
  if (cnt>0) {
    read_array_io(idx,cnt,m_buf_ptr);
  }
  //
  // buffer_print(); // dbg
}

char* affx::File5_Vector::buffer_idx2ptr(size_t idx)
{
  // printf("buffer_idx2ptr(%lu)\n",idx); // dbg
  // should
  FILE5_ASSERT(m_buf_ptr!=NULL);
  FILE5_ASSERT(m_buf_max_cnt!=0);

  // in bounds
  if ((m_buf_start_idx<=idx)&&(idx<m_buf_end_idx)) {
    return m_buf_ptr+((idx-m_buf_start_idx)*m_dtype_size);
  }
  // out of bounds... so move the buffer.
  buffer_seek(idx);
  // try again...
  if ((m_buf_start_idx<=idx)&&(idx<m_buf_end_idx)) {
    return m_buf_ptr+((idx-m_buf_start_idx)*m_dtype_size);
  }
  // failed!
  FILE5_ABORT("buffer_idx2ptr: internal error.");
  //
  return NULL;
}

char* affx::File5_Vector::buffer_pushback2ptr()
{
  // printf("buffer_pushback2ptr() >>>\n"); // dbg
  // buffer_print(); // dbg

  // Is there space at the end of the buffer?
  if ((m_vec_fill_idx==m_vec_end_idx) &&
      (m_vec_fill_idx==m_buf_end_idx) &&
      ((m_buf_end_idx-m_buf_start_idx)<m_buf_max_cnt)) {
    // bump up indexes.
    m_vec_fill_idx++;
    m_vec_end_idx++;
    m_buf_end_idx++;
    //m_buf_cnt++;
    //
    char* rv=buffer_idx2ptr(m_vec_fill_idx-1);
    // printf("buffer_pushback2ptr() <<<\n"); // dbg

    //
    m_dirty=1;
    return rv;
  }

  //
  m_vec_fill_idx++;
  if (m_vec_fill_idx>m_vec_end_idx) {
    m_vec_end_idx++;
  }

  char* rv=buffer_idx2ptr(m_vec_fill_idx-1);
  // printf("buffer_pushback2ptr() <<<\n"); // dbg
  return rv;
}

//////////

size_t affx::File5_Vector::size()
{
  return m_vec_fill_idx;
}

size_t affx::File5_Vector::reserved()
{
  return m_vec_end_idx;
}

bool  affx::File5_Vector::empty() {
  return (m_vec_fill_idx==0);
}

affx::File5_dtype_t affx::File5_Vector::file5_dtype() {
  return m_dtype;
}

//////////

int
affx::File5_Vector::open(const std::string& name)
{
  return open(name,affx::FILE5_DTYPE_ANY,affx::FILE5_OPEN);
}

int
affx::File5_Vector::open(const std::string& name,affx::File5_dtype_t dtype,int flags)
{
  hid_t h5_cparms;
  herr_t rv;
  //
  close();

  //
  m_vec_fill_idx=0;
  m_vec_end_idx=0;
  //
  m_dtype=affx::FILE5_DTYPE_UNKNOWN;
  m_dtype_size=0;

  //
  m_name=name;
  m_flags=flags;

  // if given the dtype of ANY, then we can only be an open.
  if ((dtype==affx::FILE5_DTYPE_ANY)&&(flags!=affx::FILE5_OPEN)) {
    FILE5_ABORT("May only use 'FILE5_DTYPE_ANY' with 'FILE5_OPEN'");
  }

  // remove the orginal
  if ((flags&affx::FILE5_REPLACE)==affx::FILE5_REPLACE) {
    //printf("file5: vector: %p: replace\n",this);
    H5E_BEGIN_TRY {
      H5Gunlink(h5_parent_id(),m_name.c_str());
    } H5E_END_TRY;
    // replace implies create...
    flags|=FILE5_CREATE;
  }

  //
  if (((flags&affx::FILE5_OPEN)==affx::FILE5_OPEN) &&
      (getParent()->name_exists(name))) {
    //printf("### File5_Vector::open(%s,%d,%d) => open\n",name.c_str(),dtype,flags);
    m_h5_obj=H5Dopen(h5_parent_id(),m_name.c_str());
    FILE5_CHECKID(m_h5_obj,"H5Dopen failed. Out of memory?");
  }
  //
  else if ((flags&affx::FILE5_CREATE)==affx::FILE5_CREATE) {
    // the options which will be set.
    //printf("file5: vector: %p: create\n",this);

    h5_cparms=H5Pcreate(H5P_DATASET_CREATE);
    FILE5_CHECKID(h5_cparms,"H5Pcreate failed? Not enough memory?");

    // attempt to create the vector.
    if ((flags&FILE5_CREATE)==FILE5_CREATE) {
      // set the type since we are creating
      m_dtype=dtype;
      // for strings we set compression & create our datatype...
      if (m_dtype==FILE5_DTYPE_STRING) {
        //
        if (m_opt_chunksize==-1) {
          m_opt_chunksize=FILE5_CHUNKSIZE_STRING;
        }
        if (m_opt_compress==-1) {
          m_opt_compress=FILE5_DEFAULT_COMPRESS;
        }
        if (m_opt_crc==-1) {
          m_opt_crc=1;
        }
        //
        m_h5_dtype=H5Tcopy(H5T_C_S1);
        m_h5_dtype_tofree=m_h5_dtype;
        // treat "0" the same as "-1"
        if (m_opt_string_size==0) {
          m_opt_string_size=-1;
        }
        if (m_opt_string_size==-1) {
          // size=-1 => Make a variable length string
          H5Tset_size(m_h5_dtype,H5T_VARIABLE);
        }
        else {
          H5Tset_size(m_h5_dtype,m_opt_string_size);
        }
      }
      else {
        // non string
        m_h5_dtype=as_hdf5_dtype(m_dtype);
        m_h5_dtype_tofree=-1;
        //
        if (m_opt_chunksize==-1) {
          m_opt_chunksize=FILE5_CHUNKSIZE;
        }
      }

      // chunk it?
      if (m_opt_chunksize>0) {
        FILE5_DIM1(chunk_dims,m_opt_chunksize);
        m_h5_status=H5Pset_chunk(h5_cparms,1,chunk_dims);
        FILE5_CHECKRV(m_h5_status,"H5Pset_chunk")
      }
      // compress it?
      if (m_opt_compress>=0) {
        m_h5_status=H5Pset_deflate(h5_cparms,m_opt_compress);
        FILE5_CHECKRV(m_h5_status,"H5Pset_deflate");
      }
      // turn on crc checks?
      if (m_opt_crc==1) {
        m_h5_status=H5Pset_fletcher32(h5_cparms);
        FILE5_CHECKRV(m_h5_status,"H5Pset_fletcher32");
      }

      //
      FILE5_DIM1(dims,0);
      FILE5_DIM1(dims_max,H5S_UNLIMITED);
      m_h5_dspace=H5Screate_simple(1,dims,dims_max);
      FILE5_CHECKID(m_h5_dspace,"H5Screate_simple");

      //printf("File5_Vector: create: string_size=%d chunksize=%d compress=%d\n",
      //m_opt_stringsize,m_opt_chunksize,m_opt_compressionlevel);

      m_h5_obj=H5Dcreate(h5_parent_id(),
                         m_name.c_str(),
                         m_h5_dtype,
                         m_h5_dspace,
                         h5_cparms);
      FILE5_CHECKID(m_h5_obj,"H5Dcreate");

      //
      m_h5_status=H5Pclose(h5_cparms);
      // cant destroy basic types/
      if (m_h5_dtype_tofree>0) {
        H5Tclose(m_h5_dtype_tofree);
      }
      m_h5_dtype_tofree=-1;
      m_h5_dtype=-1;
      // setup the meta-data.
      m_vec_fill_idx=0;
      attrib_set(VECTOR_FILL_INDEX,m_vec_fill_idx);
      //m_dtype=dtype;
      set_file5_kind_string("file5-vector");
      // meta-data needs update
      m_dirty=1;
    }
  }
  //  // open it
  //  else {
  //    //printf("file5: vector: %p: open\n",this);
  //    m_h5_obj=H5Dopen(h5_parent_id(),m_name.c_str());
  //    rv=attrib_get("fill-index",&m_vec_fill_idx);
  //    //
  //    m_h5_dtype=H5Dget_type(m_h5_obj);
  //    m_h5_dtype_tofree=m_h5_dtype;
  //    m_dtype=as_file5_dtype(m_h5_dtype);
  //    m_opt_chunksize=FILE5_CHUNKSIZE;
  //    //
  //    if (H5Tget_class(m_h5_dtype)==H5T_STRING) {
  //      m_opt_string_size=H5Tget_size(m_h5_dtype);
  //      // @todo: not quite right
  //      m_opt_chunksize=FILE5_CHUNKSIZE_STRING;
  //      //
  //      if (m_opt_string_size==H5T_VARIABLE) {
  //        m_opt_string_size=-1;
  //      }
  //    }
  //  }

  // how big are the chunks?
  h5_cparms=H5Dget_create_plist(m_h5_obj);
  FILE5_CHECKID(h5_cparms,"open");
  // allocates an array to be used
  FILE5_DIM1(chunk_dims,0);
  H5E_BEGIN_TRY {
    // ask for the chunk-size
    rv=H5Pget_chunk(h5_cparms,1,chunk_dims);
  } H5E_END_TRY;
  if (rv<0) {
    // an error: there wasnt a chunksize on the dataset, 
    // so lets just use this value for our buffering.
    m_opt_chunksize=5000;
  }
  else {
    // we will do our buffering in chunk-size chunks.
    m_opt_chunksize=chunk_dims[0];
  }
  H5Pclose(h5_cparms);

  // find the ends of this vector
  hid_t obj_space=H5Dget_space(m_h5_obj);
  FILE5_CHECKID(obj_space,"open");
  FILE5_DIM1(obj_dims,0);
  FILE5_DIM1(obj_max_dims,0);
  H5Sget_simple_extent_dims(obj_space,obj_dims,obj_max_dims);
  H5Sclose(obj_space);
  // set the file,vector and fill indexes.
  m_file_end_idx=obj_dims[0];
  m_vec_end_idx=m_file_end_idx;
  attrib_get("fill-index",&m_vec_fill_idx);

  // figure out the data type size
  hid_t tmp_dtype_file=H5Dget_type(m_h5_obj);
  hid_t tmp_dtype_native=affx::as_h5t_native_type(tmp_dtype_file);
  // we could convert
  if (tmp_dtype_native!=-1) {
    m_h5_dtype=tmp_dtype_native;
    m_h5_dtype_tofree=-1;
    H5Tclose(tmp_dtype_file);
  }
  else if (H5Tget_class(tmp_dtype_file)==H5T_STRING) {
    m_h5_dtype=tmp_dtype_file;
    m_h5_dtype_tofree=m_h5_dtype;
    m_h5_dtype_size=H5Tget_size(m_h5_dtype);
    m_dtype_size=m_h5_dtype_size;
    m_opt_string_size=m_h5_dtype_size;
    // if it is variable length, use the size of a pointer.
    if (m_dtype_size==H5T_VARIABLE) {
      m_dtype_size=sizeof(char*);
    }
  }
  else {
    FILE5_ABORT("Unhandled datatype.");
  }
  m_dtype=affx::as_file5_dtype(m_h5_dtype);
  m_h5_dtype_size=H5Tget_size(m_h5_dtype);
  m_dtype_size=m_h5_dtype_size;

  // if we havent been told what buffer size to use...
  if (m_opt_buffer_size==-1) {
    // make it the same as the chunksize.
    m_opt_buffer_size=m_opt_chunksize;
  }
  buffer_alloc(m_opt_buffer_size);

  // @todo: what happens if the data types dont match?
  if ((dtype!=affx::FILE5_DTYPE_ANY)&&(m_dtype!=dtype)) {
    FILE5_ABORT("Data type of the vector does not match the requested datatype to open.");
  }

  //
  //parent_refcnt_inc();
  m_state=affx::FILE5_STATE_OPEN;
  // write out meta-data if needed.
  flush();
  //printf("file5: vector: %p: fill-index=%d\n",this,m_fill_idx);
  return 0;
}

//////////

/// @brief     The internal read function. Read data from HDF5 without using our cache.
/// @param     idx       start
/// @param     cnt       count
/// @param     ptr       where to put the data
/// @return    number of items read
int affx::File5_Vector::read_array_io(size_t idx,size_t cnt,void* ptr)
{
  //printf("read_array(%lu,%lu,%p)\n",idx,cnt,ptr); // dbg

  // zero the mem before filling the buffer.
  // memset(ptr,0,cnt*m_h5_dtype_size);

  // clamp the end of the read to the end of the vector
  // by adjusting the length.
  if (idx+cnt>m_file_end_idx) {
    cnt=m_file_end_idx-idx;
  }
  // nothing to read?
  if (cnt<=0) {
    return 0;
  }

  // be sure!
  FILE5_ASSERT(idx+cnt<=m_file_end_idx);
  FILE5_ASSERT(idx+cnt<=m_vec_end_idx);

  // the memory dspace
  FILE5_DIM1(m_dims,  cnt);
  hid_t m_dspace=H5Screate_simple(1,m_dims,NULL);
  FILE5_CHECKID(m_dspace,"read_array_io");
  H5Sselect_all(m_dspace);

  // the file dspace
  hid_t f_dspace=H5Dget_space(m_h5_obj);
  FILE5_CHECKID(f_dspace,"read_array_io");
  FILE5_DIM1(f_offset,idx);
  FILE5_DIM1(f_count, cnt);
  H5Sselect_hyperslab(f_dspace,H5S_SELECT_SET,f_offset,NULL,f_count,NULL);

  //
  // printf("io: M: ");
  // affx::dump_h5_dspace(m_dspace);
  // printf("\n     F: ");
  // affx::dump_h5_dspace(f_dspace);
  // printf("\n");

  // affx::dump_h5_dset(m_h5_obj); // dbg
  // Do the actual IO
  m_h5_status=H5Dread(m_h5_obj,m_h5_dtype,m_dspace,f_dspace,H5P_DEFAULT,ptr);
  if (m_h5_status!=0) {
    printf("read_array_io: failed (idx=%d,cnt=%d,vec_end=%d,file_end=%d)",
           (int)idx,(int)cnt,(int)m_vec_end_idx,(int)m_file_end_idx);
  }
  FILE5_CHECKRV(m_h5_status,"HD5read");

  //
  H5Sclose(m_dspace);
  H5Sclose(f_dspace);
  //
  return cnt;
}

/// @brief     The internal write function. Write data to HDF5 without using our cache.
/// @param     idx       start
/// @param     cnt       count
/// @param     ptr       where to get the data
/// @return    number of items written
int
affx::File5_Vector::write_array_io(size_t idx,size_t cnt,const void* ptr)
{
  if (cnt==0) {
    return 0;
  }
  // need to flush, as this bypasses the cache.
  // flush();

  // buffer_print(); // dbg
  // make sure it is at least this big.
  reserveFile(idx+cnt);

  // the memory dspace
  FILE5_DIM1(m_dims, cnt);
  hid_t m_dspace=H5Screate_simple(1,m_dims,NULL);
  FILE5_CHECKID(m_dspace,"write_array_io");
  H5Sselect_all(m_dspace);

  // the file dspace
  FILE5_DIM1(f_start, idx);
  FILE5_DIM1(f_count, cnt);

  hid_t f_dspace=H5Dget_space(m_h5_obj);
  FILE5_CHECKID(f_dspace,"write_array_io");
  m_h5_status=H5Sselect_hyperslab(f_dspace,H5S_SELECT_SET,f_start,NULL,f_count,NULL);
  FILE5_CHECKRV(m_h5_status,"H5Sselect_hyperslab");
  //
  m_h5_status=H5Dwrite(m_h5_obj,m_h5_dtype,m_dspace,f_dspace,H5P_DEFAULT,ptr);
  FILE5_CHECKRV(m_h5_status,"H5Dwrite");
  //
  H5Sclose(m_dspace);
  H5Sclose(f_dspace);

  //
  return cnt;
}

/////

/// @brief     Read data into an array provided by the user.
/// @param     idx       idx to start the read from File5
/// @param     cnt       number of items to read
/// @param     ptr       where to put the data
/// @return    number of items read (==cnt)
int affx::File5_Vector::read_array(size_t idx,size_t cnt,void* ptr)
{
  FILE5_ASSERT(ptr!=NULL);
  // We flush our file5 cache to HDF5 first,
  // as we should see its data during our read.
  flush();
  buffer_clear();
  // now do our IO
  return read_array_io(idx,cnt,ptr);
}


/// @brief     Write data into an array provided by the user
/// @param     idx       idx to start the read from
/// @param     cnt       number of items to read
/// @param     ptr       where to put the data
/// @return    number of items written (==cnt)
int
affx::File5_Vector::write_array(size_t idx,size_t cnt,const void* ptr)
{
  FILE5_ASSERT(ptr!=NULL);
  // as above, we want our cache to be put into HDF5 space, before this IO.
  flush();
  // now do our IO
  int rv=write_array_io(idx,cnt,ptr);
  FILE5_ASSERT(rv==cnt);
  // move the fill_idx if needed.
  int write_end_idx=idx+cnt;
  if (write_end_idx>m_vec_end_idx) {
    m_vec_end_idx=write_end_idx;
  }
  if (write_end_idx>m_vec_fill_idx) {
    m_vec_fill_idx=write_end_idx;
  }
  //
  return rv;
}

//////////

int affx::File5_Vector::read_scalar(size_t idx,void* ptr)
{
  return read_array(idx,1,ptr);
}

int
affx::File5_Vector::write_scalar(size_t idx,void* ptr)
{
  return write_array(idx,1,ptr);
}

/////

void affx::File5_Vector::assert_dtype(affx::File5_dtype_t dtype,const std::string& msg)
{
  if (m_dtype!=dtype) {
    FILE5_ABORT(msg+": data type mismatch.");
  }
}

/////

// int
int affx::File5_Vector::read_vector(size_t idx,std::vector<char>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_CHAR,"File5_Vector::read_vector<char>");
  int cnt=read_array(idx,vec->size(),&(*vec)[0]);
  vec->resize(cnt);
  return cnt;
}
int affx::File5_Vector::write_vector(size_t idx,const std::vector<char>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_CHAR,"File5_Vector::write_vector<char>");
  return write_array(idx,vec->size(),&vec[0]);
}
// short
int affx::File5_Vector::read_vector(size_t idx,std::vector<short>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_SHORT,"File5_Vector::read_vector<short>");
  int cnt=read_array(idx,vec->size(),&(*vec)[0]);
  vec->resize(cnt);
  return cnt;
}
int affx::File5_Vector::write_vector(size_t idx,const std::vector<short>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_SHORT,"File5_Vector::write_vector<short>");
  return write_array(idx,vec->size(),&vec[0]);
}

// int
int affx::File5_Vector::read_vector(size_t idx,std::vector<int>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_INT,"File5_Vector::read_vector<int>");
  int cnt=read_array(idx,vec->size(),&(*vec)[0]);
  vec->resize(cnt);
  return cnt;
}
int affx::File5_Vector::write_vector(size_t idx,const std::vector<int>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_INT,"File5_Vector::write_vector<int>");
  return write_array(idx,vec->size(),&(*vec)[0]);
}
// float
int affx::File5_Vector::read_vector(size_t idx,std::vector<float>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_FLOAT,"File5_Vector::read_vector<float>");
  int cnt=read_array(idx,vec->size(),&(*vec)[0]);
  vec->resize(cnt);
  return cnt;
}
int affx::File5_Vector::write_vector(size_t idx,const std::vector<float>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_FLOAT,"File5_Vector::write_vector<float>");
  return write_array(idx,vec->size(),&(*vec)[0]);
}
// double
int affx::File5_Vector::read_vector(size_t idx,std::vector<double>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_DOUBLE,"File5_Vector::read_vector<double>");
  int cnt=read_array(idx,vec->size(),&(*vec)[0]);
  vec->resize(cnt);
  return cnt;
}
int affx::File5_Vector::write_vector(size_t idx,const std::vector<double>* vec)
{
  assert_dtype(affx::FILE5_DTYPE_INT,"File5_Vector::write_vector<double>");
  return write_array(idx,vec->size(),&(*vec)[0]);
}

//////////

// if (0) { printf("get(%d)=(%p)\n",(int)idx,ptr); }

int affx::File5_Vector::push_back(int val)
{
  switch (m_dtype) {
  case affx::FILE5_DTYPE_CHAR:
    return push_back_c(val);
    break;
  case affx::FILE5_DTYPE_SHORT:
    return push_back_s(val);
    break;
  case affx::FILE5_DTYPE_INT:
    return push_back_i(val);
    break;
  case affx::FILE5_DTYPE_FLOAT:
    return push_back_f(val);
    break;
  case affx::FILE5_DTYPE_DOUBLE:
    return push_back_d(val);
    break;
  default:
    FILE5_ABORT("File5_Vector::push_back(): wrong type for vector.");
    break;
  }
  FILE5_ABORT("File5_Vector::push_back(): internal error");
  return 0;
}

int affx::File5_Vector::get(size_t idx,int* val)
{
  void* ptr=buffer_idx2ptr(idx);

  *val=0;

  switch (m_dtype) {
  case affx::FILE5_DTYPE_CHAR:
    *val=*(char*)ptr;
    break;
  case affx::FILE5_DTYPE_SHORT:
    *val=*(short*)ptr;
    break;
  case affx::FILE5_DTYPE_INT:
    *val=*(int*)ptr;
    break;
  default:
    return affx::FILE5_ERR;
    break;
  }
  return affx::FILE5_OK;
}

int affx::File5_Vector::get(size_t idx,char* val)
{
  void* ptr=buffer_idx2ptr(idx);

  *val=0;

  switch (m_dtype) {
  case affx::FILE5_DTYPE_CHAR:
    *val=*(char*)ptr;
    break;
  default:
    return affx::FILE5_ERR;
    break;
  }
  return affx::FILE5_OK;
}

int affx::File5_Vector::get(size_t idx,float* val)
{
  void* ptr=buffer_idx2ptr(idx);

  *val=0.0;

  switch (m_dtype) {
  case affx::FILE5_DTYPE_FLOAT:
    *val=*(float*)ptr;
    break;
  case affx::FILE5_DTYPE_DOUBLE:
    *val=*(double*)ptr;
    break;
  default:
    return affx::FILE5_ERR;
    break;
  }
  return affx::FILE5_OK;
}

int affx::File5_Vector::get(size_t idx,double* val)
{
  void* ptr=buffer_idx2ptr(idx);

  *val=0.0;

  switch (m_dtype) {
  case affx::FILE5_DTYPE_FLOAT:
    *val=*(float*)ptr;
    break;
  case affx::FILE5_DTYPE_DOUBLE:
    *val=*(double*)ptr;
    break;
  default:
    return affx::FILE5_ERR;
    break;
  }
  return affx::FILE5_OK;
}

// @todo "pop_back()" isnt defined.
#define FILE5_ACCESSORS(TYPE,SUFFIX)                            \
  int affx::File5_Vector::get_##SUFFIX(size_t idx,TYPE* val)    \
  {                                                             \
    TYPE* ptr=(TYPE*)buffer_idx2ptr(idx);                       \
    *val=*ptr;                                                  \
    return affx::FILE5_OK;                                      \
  }                                                             \
  int affx::File5_Vector::assign_##SUFFIX(size_t idx,TYPE val)  \
  {                                                             \
    TYPE* ptr=(TYPE*)buffer_idx2ptr(idx);                       \
    FILE5_ASSERT(ptr!=NULL);                                    \
    *ptr=val;                                                   \
    return affx::FILE5_OK;                                      \
  }                                                             \
  int affx::File5_Vector::push_back_##SUFFIX(TYPE val)          \
  {                                                             \
    TYPE* ptr=(TYPE*)buffer_pushback2ptr();                     \
    FILE5_ASSERT(ptr!=NULL);                                    \
    *ptr=val;                                                   \
    return affx::FILE5_OK;                                      \
  }

// expand for the types.
FILE5_ACCESSORS(char,c);
FILE5_ACCESSORS(short,s);
FILE5_ACCESSORS(int,i);
FILE5_ACCESSORS(float,f);
FILE5_ACCESSORS(double,d);

// The string functions have two cases.
// * fixed size strings and
// * var length strings.

int affx::File5_Vector::get(size_t idx,std::string* val)
{
  // has to be the correct type.
  if (m_dtype!=affx::FILE5_DTYPE_STRING) {
    return FILE5_ERR;
  }

  if (m_opt_string_size==-1) {
    char** ptr=(char**)buffer_idx2ptr(idx);
    val->assign(*ptr);
  }
  else {
    char* ptr=(char*)buffer_idx2ptr(idx);
    //val->assign(ptr,m_opt_string_size);
    val->assign(ptr,m_opt_string_size);
    size_t pos=val->find_first_of(char(0));
    if (pos!=std::string::npos) {
      val->erase(pos);
    }
  }
  return affx::FILE5_OK;
}

//std::string& at_string(size_t idx);
int affx::File5_Vector::assign_string(size_t idx,const std::string& val)
{
  m_dirty=1;
  //
  if (m_opt_string_size==-1) {
    char** ptr=(char**)buffer_idx2ptr(idx);
    //
    if (*ptr!=NULL) {
      free(*ptr);
    }
    //
    *ptr=strdup(val.c_str());
  }
  else {
    char* ptr=(char*)buffer_idx2ptr(idx);
    strncpy(ptr,val.c_str(),m_opt_string_size);
  }
  return affx::FILE5_OK;
}

int affx::File5_Vector::push_back_string(const std::string& val)
{
  m_dirty=1;
  //
  if (m_opt_string_size==-1) {
    char** ptr=(char**)buffer_pushback2ptr();
    //
    if (*ptr!=NULL) {
      free(*ptr);
    }
    //
    *ptr=strdup(val.c_str());
  }
  else {
    ///@todo currently throwing errors in regression tests due to command line strings for apt-copynumber-workflow and apt-copynumber-cyto. If there really is a max string limit (eg 1024) then we probably want this check, but then we need to modify setting of meta info to truncate fields.
    /*
    if(m_opt_string_size >= 0 && val.length() >= m_opt_string_size) {
      FILE5_ABORT("String too long. Expecting " + ToStr(m_opt_string_size) + " max characters but got '" + val + "' with length: " + ToStr(val.size()));
    }
    */
    char* ptr=(char*)buffer_pushback2ptr();
    strncpy(ptr,val.c_str(),m_opt_string_size);
  }
  return affx::FILE5_OK;
}

////////// 

// @todo for now we do the IO as a series of "get" and "assign_string" calls.
//       we could do something better (block cache aware) in the future.
int affx::File5_Vector::read_vector(size_t file5_idx,std::vector<std::string>* vec)
{
  int rv;
  std::string* str_ptr;
  for (unsigned int i=0;i<vec->size();i++) {
    str_ptr=&((*vec)[file5_idx]);
    rv=get(file5_idx,str_ptr);
    if (rv!=FILE5_OK) {
      return rv;
    }
    file5_idx++;
  }
  return FILE5_OK;
}

int affx::File5_Vector::write_vector(size_t file5_idx,const std::vector<std::string>* vec)
{
  int rv;
  for (unsigned int i=0;i<vec->size();i++) {
    rv=assign_string(file5_idx,(*vec)[file5_idx]);
    if (rv!=FILE5_OK) {
      return rv;
    }
    file5_idx++;
  }
  return FILE5_OK;
}

/////

int affx::File5_Vector::getAsStr(size_t idx,std::string* str)
{
  int rv;

  switch (m_dtype) {

  case affx::FILE5_DTYPE_STRING:
    rv=get(idx,str);
    return rv;
    break;

  case affx::FILE5_DTYPE_CHAR:
    char tmp_char;
    rv=get_c(idx,&tmp_char);
    if (rv==affx::FILE5_OK) {
      *str=ToStr((int)tmp_char);
    }
    return rv;
    break;

  case affx::FILE5_DTYPE_SHORT:
    short tmp_short;
    rv=get_s(idx,&tmp_short);
    if (rv==affx::FILE5_OK) {
      *str=ToStr((int)tmp_short);
    }
    return rv;
    break;

  case affx::FILE5_DTYPE_INT:
    int tmp_int;
    rv=get_i(idx,&tmp_int);
    if (rv==affx::FILE5_OK) {
      *str=ToStr(tmp_int);
    }
    return rv;
    break;

  case affx::FILE5_DTYPE_DOUBLE:
    double tmp_double;
    rv=get_d(idx,&tmp_double);
    if (rv==affx::FILE5_OK) {
      *str=ToStr(tmp_double);
    }
    return rv;
    break;

  case affx::FILE5_DTYPE_FLOAT:
    float tmp_float;
    rv=get_f(idx,&tmp_float);
    if (rv==affx::FILE5_OK) {
      *str=ToStr(tmp_float);
    }
    return rv;
    break;

  default:
    FILE5_ABORT("Unable to convert this type to string.");
    break;
  }

  return affx::FILE5_ERR;
}

/////

int affx::File5_Vector::getAsStrVector(std::vector<std::string>* strvec)
{
  return getAsStrVector(strvec,0,m_vec_fill_idx);
}

int affx::File5_Vector::getAsStrVector(std::vector<std::string>* strvec,
                                       size_t read_start,
                                       size_t read_end)
{
  int strvec_len=read_end-read_start;
  strvec->resize(strvec_len);

  // @todo we could test for the type here and avoid the "select" in getAsStr()
  for (size_t strvec_idx=0;strvec_idx<strvec_len;strvec_idx++) {
    if (getAsStr(read_start+strvec_idx,&((*strvec)[strvec_idx]))!=affx::FILE5_OK) {
      FILE5_ABORT("Unable to convert string.");
    }
  }

  return affx::FILE5_OK;
}
