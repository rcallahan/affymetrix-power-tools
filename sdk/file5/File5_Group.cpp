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
// affy/sdk/file5/File5_Group.cpp ---
// 
// $Id: File5_Group.cpp,v 1.33 2009-09-25 17:49:03 mspald Exp $
// 


//
#include "file5/File5_Group.h"
//
#include "file5/File5.h"
#include "file5/File5_Matrix.h"
#include "file5/File5_Tsv.h"
#include "file5/File5_Vector.h"
//
#include "util/Convert.h"
#include "util/Err.h"
//
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
//

#define USER_META_NAME "_user_meta"

// #define FILE5_DEBUG_PRINT 1

//////////

affx::File5_Group::File5_Group()
{
  m_kind_char='G';
  init();
}

affx::File5_Group::~File5_Group()
{
  close();
}

//////////

affx::File5_return_t
affx::File5_Group::init()
{
  affx::File5_Object::init();
  //
  m_kind=FILE5_KIND_GROUP;
  m_kind_char='G';
  //
  m_h5_user_meta_id=-1;
  m_h5_user_meta_idx=-1;
  //
  return affx::FILE5_OK;
}

char affx::File5_Group::file5_kind_char() {
  return 'G';
}

affx::File5_return_t
affx::File5_Group::dump()
{
  printf("== File5_Group (%p)\n",this);
  dump1();
  return affx::FILE5_OK;
}

affx::File5_return_t
affx::File5_Group::dump1()
{
  //printf("   m_h5_group_id  = %d\n",m_h5_group_id);
  return affx::File5_Object::dump1();
}

hid_t
affx::File5_Group::usermeta_getid()
{
  if (m_h5_user_meta_id!=-1) {
    return m_h5_user_meta_id;
  }

  //
  if (name_exists(USER_META_NAME)) { // use the existing one
    m_h5_user_meta_id=H5Dopen(m_h5_obj,USER_META_NAME);
    FILE5_CHECKID(m_h5_user_meta_id,"usermeta_getid");
  }
  else if (m_readonly==1) {
    m_h5_user_meta_id=-1;
  }
  else { 
    // varlen string
    hid_t dtype=H5Tcopy(H5T_C_S1);
    FILE5_CHECKID(dtype,"H5Tcopy failed.");
    H5Tset_size(dtype,H5T_VARIABLE);
    // create a (2,inf) space
    FILE5_DIM2(dims,0,2);
    FILE5_DIM2(dims_max,H5S_UNLIMITED,2);
    hid_t h5_dspace=H5Screate_simple(2,dims,dims_max);
    FILE5_CHECKID(h5_dspace,"H5Screate_simple");
    //
    hid_t h5_cparms=H5Pcreate(H5P_DATASET_CREATE);
    FILE5_CHECKID(h5_cparms,"H5Pcreate failed.");
    FILE5_DIM2(chunk_dims,100,2);
    // say how big the chunk is
    m_h5_status=H5Pset_chunk(h5_cparms,2,chunk_dims);
    FILE5_CHECKRV(m_h5_status,"H5Pset_chunk");
    // turn on crc checks
    if (m_opt_crc==1) {
      m_h5_status=H5Pset_fletcher32(h5_cparms);
      FILE5_CHECKRV(m_h5_status,"H5Pset_fletcher32");
    }
    //
    m_h5_user_meta_id=H5Dcreate(m_h5_obj,
                                USER_META_NAME,
                                dtype,
                                h5_dspace,
                                h5_cparms);
    FILE5_CHECKID(m_h5_user_meta_id,"H5Dcreate");

    //
    H5Sclose(h5_dspace);
    H5Pclose(h5_cparms);
    H5Tclose(dtype);
  }
  //
  return m_h5_user_meta_id;
}

affx::File5_return_t
affx::File5_Group::usermeta_close()
{
  if (m_h5_user_meta_id!=-1) {
    H5Dclose(m_h5_user_meta_id);
  }
  m_h5_user_meta_id=-1;
  //
  return affx::FILE5_OK;
}

int
affx::File5_Group::usermeta_size()
{
  if (usermeta_getid()==-1) {
    return 0;
  }
  // resize
  hid_t dspace=H5Dget_space(m_h5_user_meta_id);
  FILE5_CHECKID(dspace,"H5Dget_space failed.");
  FILE5_DIM2(size_now,0,0);
  FILE5_DIM2(size_max,0,0);
  H5Sget_simple_extent_dims(dspace,size_now,size_max);
  //
  H5Sclose(dspace);
  //
  return size_now[0];
}

affx::File5_return_t
affx::File5_Group::usermeta_resize(int size)
{
  if (usermeta_getid()==-1) {
    return affx::FILE5_ERR_BADID;
  }
  // resize
  hid_t dspace=H5Dget_space(m_h5_user_meta_id);
  FILE5_CHECKID(dspace,"H5Dget_space failed.");
  FILE5_DIM2(size_now,0,0);
  FILE5_DIM2(size_max,0,0);
  H5Sget_simple_extent_dims(dspace,size_now,size_max);
  size_now[0]=size;
  m_h5_status=H5Dextend(m_h5_user_meta_id,size_now);
  FILE5_CHECKRV(m_h5_status,"H5Dextend failed.");
  
  //
  H5Sclose(dspace);
  //
  return affx::FILE5_OK;
}

// @todo: change to fixed or varlen & zero out buffer before write. -jhg
affx::File5_return_t
affx::File5_Group::usermeta_writeidx(int idx,const std::string& key,const std::string& val)
{
  if (usermeta_getid()==-1) {
    return affx::FILE5_ERR_BADID;
  }

  // varlen string
  hid_t dtype=H5Tcopy(H5T_C_S1);
  FILE5_CHECKID(dtype,"H5Tcopy failed.");
  H5Tset_size(dtype,H5T_VARIABLE);
  
  //
  const char* buf[2];
  buf[0]=key.c_str();
  buf[1]=val.c_str();
  //
  FILE5_DIM2(m_dims,1,2);
  hid_t m_dspace=H5Screate_simple(2,m_dims,NULL);
  FILE5_CHECKID(m_dspace,"H5Screate_simple failed.");
  H5Sselect_all(m_dspace);
  //
  FILE5_DIM2(f_start,idx,0);
  FILE5_DIM2(f_count,1,2);
  hid_t f_dspace=H5Dget_space(m_h5_user_meta_id);
  FILE5_CHECKID(f_dspace,"H5Dget_space failed.");
  m_h5_status=H5Sselect_hyperslab(f_dspace,H5S_SELECT_SET,f_start,NULL,f_count,NULL);
  FILE5_CHECKRV(m_h5_status,"H5Sselect_hyperslab");
  //
  m_h5_status=H5Dwrite(m_h5_user_meta_id,dtype,m_dspace,f_dspace,H5P_DEFAULT,&buf);
  FILE5_CHECKRV(m_h5_status,"H5Dwrite");
  //
  H5Sclose(m_dspace);
  H5Sclose(f_dspace);
  H5Tclose(dtype);
  //
  return affx::FILE5_OK;
}

affx::File5_return_t
affx::File5_Group::addHeader(const std::string& key,const std::string& val)
{
  // for now just punt.
#ifdef FILE5_DEBUG_PRINT
  printf("### File5_Group::addHeader('%s','%s')\n",key.c_str(),val.c_str());
#endif
  //
  if (usermeta_getid()==-1) {
    return affx::FILE5_ERR_BADID;
  }

  // resize
  int write_idx=usermeta_size();
  usermeta_resize(usermeta_size()+1);
  usermeta_writeidx(write_idx,key,val);
  //
  return affx::FILE5_OK;
}

affx::File5_return_t
affx::File5_Group::getHeaderByIdx(int idx,std::string* key,std::string* val)
{
  if (usermeta_getid()==-1) {
    return affx::FILE5_ERR_BADID;
  }
  
  // varlen string
  hid_t dtype=H5Tcopy(H5T_C_S1);
  H5Tset_size(dtype,H5T_VARIABLE);
  
  // spot for the varlen data.
  char* ptr_buf[2];

  //
  FILE5_DIM2(m_dims,1,2);
  hid_t m_dspace=H5Screate_simple(2,m_dims,NULL);
  H5Sselect_all(m_dspace);
  //
  FILE5_DIM2(f_start,idx,0);
  FILE5_DIM2(f_count,1,2);
  hid_t f_dspace=H5Dget_space(m_h5_user_meta_id);
  H5Sselect_hyperslab(f_dspace,H5S_SELECT_SET,f_start,NULL,f_count,NULL);

  //
  m_h5_status=H5Dread(m_h5_user_meta_id,dtype,m_dspace,f_dspace,H5P_DEFAULT,ptr_buf);
  FILE5_CHECKRV(m_h5_status,"H5Dread");

  *key=ptr_buf[0];
  *val=ptr_buf[1];
  //
  H5Dvlen_reclaim(dtype,m_dspace,H5P_DEFAULT,ptr_buf);
  H5Sclose(m_dspace);
  H5Sclose(f_dspace);
  H5Tclose(dtype);
  //
  return affx::FILE5_OK;
}

//
void affx::File5_Group::headersBegin()
{
  // go to the beginning.
  m_h5_user_meta_idx=0;
}

affx::File5_return_t
affx::File5_Group::headersNext(std::string* key,std::string* val)
{
  int h_cnt=usermeta_size();
  if (m_h5_user_meta_idx<h_cnt) {
    getHeaderByIdx(m_h5_user_meta_idx,key,val);
    m_h5_user_meta_idx++;
    return affx::FILE5_OK;
  }
  //
  return affx::FILE5_ERR;
}

affx::File5_return_t
affx::File5_Group::getHeader(const std::string& key,std::string* val)
{
  std::string h_key;
  std::string h_val;
  int h_cnt=usermeta_size();

  for (int h_idx=0;h_idx<h_cnt;h_idx++) {
    getHeaderByIdx(h_idx,&h_key,&h_val);
    if (h_key==key) {
      // found the first one, we are done.
      *val=h_val;
      return affx::FILE5_OK;
    }
  }
  return affx::FILE5_ERR;
}

affx::File5_return_t
affx::File5_Group::getHeader(const std::string& key,std::vector<std::string>* val)
{
  std::string h_key;
  std::string h_val;
  int h_cnt=usermeta_size();

  // start with none found
  val->clear();

  for (int h_idx=0;h_idx<h_cnt;h_idx++) {
    getHeaderByIdx(h_idx,&h_key,&h_val);
    if (h_key==key) {
      val->push_back(h_val);
    }
  }
  return affx::FILE5_OK;
}

//
affx::File5_return_t
affx::File5_Group::open(const std::string& name,int flags)
{
  close();
  //
  m_name=name;
  // Create parent groups if needed
  if((name.rfind("/") != std::string::npos) && (name.rfind("/") != 0)){
      affx::File5_Group::open(name.substr(0,name.rfind("/")),affx::FILE5_OPEN|affx::FILE5_CREATE);
  }
  //
  if ((flags&FILE5_REPLACE)==FILE5_REPLACE) {
    H5E_BEGIN_TRY {
      H5Gunlink(getParent()->m_h5_obj,name.c_str());
    } H5E_END_TRY;
    // replace implies create
    flags|=affx::FILE5_CREATE;
  }
  //
  if ((flags&FILE5_OPEN)==FILE5_OPEN) {
    hid_t tmp_id=-1;
    H5E_BEGIN_TRY {
      tmp_id=H5Gopen(h5_parent_id(),name.c_str());
    } H5E_END_TRY;
    if (tmp_id>0) {
      m_h5_obj=tmp_id;
      return affx::FILE5_OK;
    }
    if ((flags&FILE5_CREATE)!=FILE5_CREATE) {
      FILE5_ABORT("File5_Group::open('"+name+"'): Unable to open group!");
    }
  }
  //
  if ((flags&FILE5_CREATE)==FILE5_CREATE) {
    hid_t tmp_id;
    H5E_BEGIN_TRY {
      tmp_id=H5Gcreate(h5_parent_id(),name.c_str(),0);
    } H5E_END_TRY;
    if (tmp_id>=0) {
      m_h5_obj=tmp_id;
      set_file5_kind_string("file5-group");
      return affx::FILE5_OK;
    }
    FILE5_ABORT("File5_Group::open('"+name+"'): Unable to create group!");
  }
  //
  FILE5_ABORT("File5_Group::open('"+name+"'): Bad flags!");
  return affx::FILE5_OK;
}

affx::File5_return_t
affx::File5_Group::close()
{
  // abort before deallocating resources.
  if (refcnt()!=0) {
    // FILE5_ABORT("file5: Group: refcnt!=0");
	  return affx::FILE5_ERR_REFCNT;
  }
  //
  usermeta_close();
  //
  if (m_h5_obj!=-1) {
    H5Gclose(m_h5_obj);
  }
  m_h5_obj=-1;
  //
  affx::File5_Object::close();
  //
  return affx::FILE5_OK;
}

//////////

bool
affx::File5_Group::name_exists(const std::string& name)
{
  H5G_stat_t h5_stat_buf;
  herr_t rv;

  //printf("File5_Group_Mixin::name_exists('%s')\n",name.c_str());
  //printf("m_h5_group_id=%d\n",m_h5_group_id);

  H5E_BEGIN_TRY {
    rv=H5Gget_objinfo(m_h5_obj,name.c_str(),1,&h5_stat_buf);
  } H5E_END_TRY;

  // printf("rv=%d\n",rv);

  // translate the return code
  if (rv>=0) {
    return true;
  }
  else {
    return false;
  }
}

// the low level interface.
affx::File5_return_t
affx::File5_Group::unlink(const std::string& name)
{
  herr_t rv;

  FILE5_CHECKID(m_h5_obj,"unlink");

  H5E_BEGIN_TRY {
    rv=m_h5_status=H5Gunlink(m_h5_obj,name.c_str());
    if (rv!=0) {
      printf("File5_Group::unlink('%s'): rv=%d\n",name.c_str(),rv);
    }
  } H5E_END_TRY;

  return (rv==0)?affx::FILE5_OK:affx::FILE5_ERR;
}

// The higher level interface.
affx::File5_return_t
affx::File5_Group::deleteItem(const std::string& name)
{
  return File5_Group::unlink(name);
}

/// Delete the items in this group.
affx::File5_return_t
affx::File5_Group::deleteItemsIn(const std::string& name)
{
  if (!name_exists(name)) {
    //printf("File5_Group::deleteItemsIn('%s'): not found!\n",name.c_str());
    return affx::FILE5_ERR_NOTFOUND;
  }

  affx::File5_Group* group=this->openGroup(name,FILE5_OPEN);
  if (group==NULL) {
    return affx::FILE5_ERR;
  }

  int item_cnt=group->get_num_objs();
  std::string item_name;
  
  for (int i=0;i<item_cnt;i++) {
    item_name=group->get_name_by_idx(i);
    if (item_name!="") {
      group->deleteItem(item_name);
    }
  }
  group->close();
  delete group;
  return affx::FILE5_OK;
}

//
affx::File5_return_t
affx::File5_Group::deleteTsv(const std::string& name)
{
  if (!name_exists(name)) {
    //printf("File5_Group::deleteTsv('%s'): not found!\n",name.c_str());
    return affx::FILE5_ERR_NOTFOUND;
  }
  this->deleteItemsIn(name);
  this->deleteItem(name);
  //
  return affx::FILE5_OK;
}

//////////

int affx::File5_Group::get_num_objs()
{
  hsize_t num_objs;
  //herr_t rv=
  H5Gget_num_objs(m_h5_obj,&num_objs);
  return num_objs;
}

std::string affx::File5_Group::get_name_by_idx(int idx)
{
  std::string name="";  
  int name_len=0;
  int rv;

  name_len=H5Gget_objname_by_idx(m_h5_obj,idx,NULL,0);
  name_len++; // pad for null
  char *name_buf=(char*)malloc(name_len);
  FILE5_ASSERT(name_buf!=NULL);
  rv=H5Gget_objname_by_idx(m_h5_obj,idx,name_buf,name_len);
  FILE5_ASSERT((rv+1)==name_len);
  name.assign(name_buf,name_len);
  free(name_buf);
  //
  // printf("get_name_by_idx(%d)=='%s'\n",idx,name.c_str()); // debug
  //
  return name;
}

int affx::File5_Group::get_objtype_by_idx(int idx)
{
  return H5Gget_objtype_by_idx(m_h5_obj,idx);
}

affx::File5_Object* affx::File5_Group::get_object_by_idx(int idx)
{
  std::string obj_name=get_name_by_idx(idx);
  if (obj_name=="") {
    return NULL;
  }
  int h5_objtype=get_objtype_by_idx(idx);

  switch (h5_objtype) {
  case H5G_LINK: // 
    return NULL;
    break;
  case H5G_GROUP:
    return openGroup(obj_name,affx::FILE5_OPEN);
    break;
  case H5G_DATASET:
    return openDataset(obj_name,affx::FILE5_OPEN);
    break;
  case H5G_TYPE:
    return NULL;
    break;
  default:
    return NULL;
    break;
  }
  return NULL;
}

//////////

std::vector<std::string>
affx::File5_Group::listNames()
{
  std::vector<std::string> name_vec;
  int num_objs=get_num_objs();
  
  for (int i=0;i<num_objs;i++) {
    name_vec.push_back(get_name_by_idx(i));
  }

  return name_vec;
}
    
//////////

affx::File5_Group* 
affx::File5_Group::openGroup(const std::string& name,int flags)
{
  affx::File5_Group* group=new affx::File5_Group();
  group->setParent(this);
  group->open(name,flags);
  return group;
}

/////

affx::File5_Matrix* 
affx::File5_Group::openMatrix(const std::string& name,
                              affx::File5_dtype_t dtype,
                              const std::vector<int>& dims,
                              int flags)
{
  affx::File5_Matrix* matrix=new affx::File5_Matrix();
  matrix->setParent(this);
  matrix->open(name,dtype,dims,flags);
  return matrix;
}

affx::File5_Object*
affx::File5_Group::openDataset(const std::string& name,int flags)
{
  affx::File5_Object* dset=new affx::File5_Object();
  dset->setParent(this);
  dset->open(name,flags);
  return dset;
}

/////

affx::File5_Tsv*
affx::File5_Group::openTsv(const std::string& name)
{
  return openTsv(name,affx::FILE5_OPEN);
}

affx::File5_Tsv*
affx::File5_Group::openTsv(const std::string& name,int flags)
{
  affx::File5_Tsv* tsv=new affx::File5_Tsv();
  tsv->setParent(this);
  tsv->open(name,flags);
  return tsv;
}

/////

affx::File5_Vector*
affx::File5_Group::createVector()
{
  //affx::File5_Vector* vec=affx::File5_Vector::allocate(dtype);
  affx::File5_Vector* vec=new affx::File5_Vector();
  vec->setParent(this);
  return vec;
}

affx::File5_Vector*
affx::File5_Group::openVector(const std::string& name,
                              const affx::File5_dtype_t dtype,
                              int flags)
{
  //affx::File5_Vector* vec=affx::File5_Vector::allocate(dtype);
  affx::File5_Vector* vec=new affx::File5_Vector();
  vec->setParent(this);
  vec->open(name,dtype,flags);
  return vec;
}

affx::File5_Vector*
affx::File5_Group::openVector(const std::string& name)
{
  affx::File5_Vector* vec=new affx::File5_Vector();
  vec->setParent(this);
  vec->open(name,affx::FILE5_DTYPE_ANY,affx::FILE5_OPEN);
  return vec;
}

affx::File5_Vector*
affx::File5_Group::openVector_String(const std::string& name,int max_size,int flags)
{
  affx::File5_Vector* vec=new affx::File5_Vector();
  vec->setParent(this);
  vec->setOptStringSize(max_size);
  vec->open(name,affx::FILE5_DTYPE_STRING,flags);
  return vec;
}
