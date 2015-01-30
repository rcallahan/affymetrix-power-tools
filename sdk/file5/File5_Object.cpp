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
// affy/sdk/file5/File5_Object.cpp ---
//
// $Id: File5_Object.cpp,v 1.35 2009-09-25 17:49:03 mspald Exp $
//

#include "file5/File5_Object.h"
//
#include "file5/File5.h"
#include "file5/File5_File.h"
//
#include <algorithm>
#include <cstring>
#include <string.h>
#include <string>
#include <utility>
#include <vector>
//

//////////

void affx::File5_open()
{
  // printf("### File5_start():\n");
  H5open();
}

void affx::File5_close()
{
  // printf("### File5_close():\n");
#ifdef FILE5_DEBUG_USAGE
  // print the list if any objects havent been deleted.
  affx::file5_usecnt_print();
#endif
  //
  H5close();
}

//////////

affx::File5_Object::File5_Object()
{
  init();
#ifdef FILE5_DEBUG_USAGE
  file5_usecnt_inc(this);
#endif
}

affx::File5_Object::~File5_Object()
{
#ifdef FILE5_DEBUG_USAGE
  file5_usecnt_dec(this);
#endif
  //
  close();
  setParent(NULL);
}

affx::File5_return_t
affx::File5_Object::init()
{
  m_parent=NULL;
  m_file=NULL;
  //
  m_name="";
  m_errno=0;
  m_flags=0;
  m_dirty=0;
  m_refcnt=0;
  m_state=affx::FILE5_STATE_CLOSED;
  //
  m_readonly=true;
  //
  m_opt_chunksize=-1; // 
  m_opt_compress=-1; // our default
  m_opt_crc=1; // on
  //
  m_kind=FILE5_KIND_OBJECT;
  m_kind_char='O';
  //
  m_dtype=FILE5_DTYPE_UNKNOWN;
  m_h5_dspace=-1;
  m_h5_dtype=-1;
  m_h5_dtype_size=0;
  m_h5_dtype_tofree=-1;
  m_h5_obj=-1;
  m_h5_status=0;
  //
  return affx::FILE5_OK;
}

bool affx::File5_Object::ReadOnly()
{
  return m_readonly;
}
bool affx::File5_Object::setReadOnly(bool val)
{
  m_readonly=val;
  return val;
}

std::string affx::File5_Object::file5_kind_string() {
  std::string kind_str="";
  if (is_open()) {
    attrib_get("file5-kind-string",&kind_str);
  }
  return kind_str;
}

int affx::File5_Object::set_file5_kind_string(const std::string& kind_str) {
  if (is_open()) {
    return attrib_set("file5-kind-string",kind_str);
  }
  return -1;
}

char affx::File5_Object::file5_kind_char() {
  return m_kind_char;
}

bool
affx::File5_Object::is_open()
{
  return (m_state==affx::FILE5_STATE_OPEN);
}

affx::File5_return_t
affx::File5_Object::open(const std::string& name,int flags)
{
  close();
  //
  if (flags!=affx::FILE5_OPEN) {
    FILE5_ABORT("Only 'affx::FILE5_OPEN' is allowed.");
  }
  //
  m_name=name;
  hid_t tmp_id=H5Dopen(h5_parent_id(),name.c_str());
  FILE5_CHECKID(tmp_id,"H5Dopen failed");
  m_h5_obj=tmp_id;
  //
  //parent_refcnt_inc();
  m_state=affx::FILE5_STATE_OPEN;
  return affx::FILE5_OK;
}

void affx::File5_Object::close_oid(hid_t* oid)
{
  herr_t rv = 0;

  if (*oid==-1) {
    return;
  }
  //
  hid_t h5_type=H5Iget_type(*oid);
  switch (h5_type) {
  case H5I_ATTR:
    rv=H5Aclose(*oid);
    break;
  case H5I_DATASET:
    rv=H5Dclose(*oid);
    break;
  case H5I_DATASPACE:
    rv=H5Sclose(*oid);
    break;
  case H5I_DATATYPE:
    rv=H5Tclose(*oid);
    break;
  case H5I_FILE:
    rv=H5Fclose(*oid);
    break;
  case H5I_GROUP:
    rv=H5Gclose(*oid);
    break;
    //
  case H5I_BADID:
    FILE5_ABORT("H5I_BADID is unhandled.");
    break;
  default:
    FILE5_ABORT("unhandled case: oid="+ToStr(oid)+" type="+ToStr(h5_type));
    break;
  }
  //
  FILE5_CHECKRV(rv,"close_oid");
  *oid=-1;
}

affx::File5_return_t
affx::File5_Object::close()
{
  //
  m_h5_dtype=-1;
  close_oid(&m_h5_dtype_tofree);
  //
  close_oid(&m_h5_dspace);
  close_oid(&m_h5_obj);
  //
  m_state=affx::FILE5_STATE_CLOSED;
  //
  return affx::FILE5_OK;
}

/////

affx::File5_return_t
affx::File5_Object::flush()
{
  return affx::FILE5_OK;
}

//////////

affx::File5_return_t
affx::File5_Object::dump()
{
  printf("== File5_Object: %p\n",this);
  return dump1();
}

affx::File5_return_t
affx::File5_Object::dump1()
{
  printf("   m_name      : '%s'\n",m_name.c_str());
  printf("   m_parent    : '%p'\n",m_parent);
  printf("   m_file      : '%p'\n",m_file);
  printf("   m_errno     : '%d'\n",m_errno);
  printf("   m_flags     : '%d'\n",m_flags);
  printf("   m_dirty     : '%d'\n",m_dirty);
  printf("   m_h5_dspace : '%d'\n",m_h5_dspace);
  printf("   m_h5_obj    : '%d'\n",m_h5_obj);
  printf("   m_h5_dtype  : '%d'\n",m_h5_dtype);
  printf("   m_h5_dtype_tofree : '%d'\n",m_h5_dtype_tofree);
  //
  return affx::FILE5_OK;
}

//////////

void affx::File5_Object::setParent(affx::File5_Group* parent)
{
  // @todo: update refcnt.
  if (m_parent!=NULL) {
    m_parent->refcnt_dec(this);
  }
  m_parent=parent;
  setFile(NULL);
  if (m_parent!=NULL) {
    m_parent->refcnt_inc(this);
    m_readonly=m_parent->m_readonly;
    setFile(m_parent->getFile());
  }
}

int
affx::File5_Object::refcnt() {
  return m_refcnt;
}

int
affx::File5_Object::refcnt_inc(affx::File5_Object* referrer) {
  // dont add ourselves (circular) or null
  if ((referrer==NULL)||(referrer==this)) {
    return m_refcnt;
  }
  //
#ifdef FILE5_DEBUG_PRINT
  printf("== %c:%p->refcnt_inc(%c:%p)\n",
         file5_kind_char(),this,
         referrer->file5_kind_char(),referrer);
#endif
  // is it in the list already?
  if (refcnt_has_referrer(referrer)==1) {
    printf("ERROR: %c:%p: referrer %p already in the list!\n",file5_kind_char(),this,referrer);
    return m_refcnt;
  }
  //
  m_refcnt++;
  m_refcnt_vec.push_back(referrer);
  //
  //refcnt_print();
  //
  return m_refcnt;
}

int
affx::File5_Object::refcnt_has_referrer(affx::File5_Object* referrer)
{
  // std::vector<void*>::iterator it;
  // unsorted.
  // it=find(m_refcnt_vec.begin(),m_refcnt_vec.end(),referrer);
  for (int i=0;i<m_refcnt_vec.size();i++) {
    if (m_refcnt_vec[i]==referrer) {
      return 1;
    }
  }
  return 0;
}

int
affx::File5_Object::refcnt_dec(affx::File5_Object* referrer) {
  // 
  if ((referrer==NULL)||(referrer==this)) {
    return m_refcnt;
  }
  //
#ifdef FILE5_DEBUG_PRINT
  printf("== %c:%p->refcnt_dec(%c:%p)\n",
         file5_kind_char(),this,
         referrer->file5_kind_char(),referrer);
#endif
//  if ((referrer==NULL)||(referrer==this)) {
//    printf("ERROR: %c:%p: Dec of %c:%p\n",
//           file5_kind_char(),this,
//           referrer->file5_kind_char(),referrer);
//    return m_refcnt;
//  }
  //
  std::vector<affx::File5_Object*>::iterator it;
  int erase_cnt=0;
  //
  it=m_refcnt_vec.begin();
  for (int i=0;i<m_refcnt_vec.size();i++) {
    if (m_refcnt_vec[i]==referrer) {
      m_refcnt_vec.erase(it + i);
      i--;
      erase_cnt++;
      m_refcnt--;
    }
  }
  //
  if (erase_cnt==0) {
    printf("File5_Object::refcnt_dec(%p): Attempted to erase but not found.\n",referrer);
    refcnt_print();
    FILE5_ABORT("File5_Object::refcnt_dec()");
  }
  //
  return m_refcnt;
}


int
affx::File5_Object::refcnt_print() {
  return refcnt_print(0);
}

int
affx::File5_Object::refcnt_print(int indent) {
  std::string indent_str;
  indent_str.resize(indent,' ');

  printf("== %sReferences to %c:%p:  (cnt=%d) (size=%d)\n",
         indent_str.c_str(),
         file5_kind_char(),this,
         m_refcnt,int(m_refcnt_vec.size()));
  for (int i=0;i<m_refcnt_vec.size();i++) {
    affx::File5_Object* o=(affx::File5_Object*)m_refcnt_vec[i];
    printf("==   %s%2d: %c:%p '%s'\n",
           indent_str.c_str(),
           i,o->file5_kind_char(),o,o->m_name.c_str());
    if (o->refcnt()!=0) {
      o->refcnt_print(indent+3);
    }
  }
  //
  fflush(NULL);
  //
  return 0;
}

//////////

hid_t affx::File5_Object::h5_parent_id() {
  if (m_parent!=NULL) {
    return m_parent->m_h5_obj;
  }
  return m_h5_obj;
}

int affx::File5_Object::parent_refcnt_inc() {
  if (m_parent!=NULL) {
    return m_parent->refcnt_inc(this);
  }
  return 0;
}
int affx::File5_Object::parent_refcnt_dec() {
  if (m_parent!=NULL) {
    return m_parent->refcnt_dec(this);
  }
  return 0;
}

//////////

// FIXME: use H5E_BEGIN_TRY to just probe for it

int
affx::File5_Object::find_attrib_idx(const std::string& attrib_name)
{
  int num_attrs;
  hid_t attrib_id;
  char attrib_name_buf[200];

  memset(attrib_name_buf,0,sizeof(attrib_name));

  num_attrs=H5Aget_num_attrs(m_h5_obj);
  for (int idx=0;idx<num_attrs;idx++) {
    attrib_id=H5Aopen_idx(m_h5_obj,idx);
    int attrib_name_len=H5Aget_name(attrib_id,sizeof(attrib_name_buf)-1,attrib_name_buf);
    H5Aclose(attrib_id);
    if ((attrib_name_len==attrib_name.size()) &&
        (strcmp(attrib_name_buf,attrib_name.c_str())==0)) {
      return idx;
    }
  }
  return -1;
}

//////////

int
affx::File5_Object::attrib_set_void(const std::string& attrib_name,hid_t h5_dtype,void* val_ptr)
{
  hid_t attrib_space;
  hid_t attrib_id;

  //
  //printf("attrib_set_void('%s',%p)\n",attrib_name.c_str(),val_ptr);

  //
  FILE5_CHECKID(m_h5_obj,"attrib_set_void");

  // to overwrite, delete then create.
  attrib_delete(attrib_name);

  //
  attrib_space=H5Screate(H5S_SCALAR);
  attrib_id=H5Acreate(m_h5_obj,attrib_name.c_str(),h5_dtype,attrib_space,H5P_DEFAULT);
  m_h5_status=H5Awrite(attrib_id,h5_dtype,val_ptr);
  FILE5_CHECKRV(m_h5_status,"H5Awrite");
  //
  H5Sclose(attrib_space);
  H5Aclose(attrib_id);
  //
  return 0;
}

int
affx::File5_Object::attrib_set(const std::string& attrib_name,int val)
{
  affx::File5_type_pun_t pun;
  pun.i=val;
  return attrib_set_void(attrib_name,H5T_NATIVE_INT,&pun);
}

int
affx::File5_Object::attrib_set(const std::string& attrib_name,size_t val)
{
  return attrib_set(attrib_name,int(val));
}

int
affx::File5_Object::attrib_set(const std::string& attrib_name,double val)
{
  affx::File5_type_pun_t pun;
  pun.d=val;
  return attrib_set_void(attrib_name,H5T_NATIVE_DOUBLE,&pun);
}

// @todo fixme this is a array of variable length not a single item with a variable length.
int
affx::File5_Object::attrib_set(const std::string& attrib_name,const std::string& val)
{
  hid_t attrib_space;
  hid_t attrib_id;
  hid_t attrib_dtype;

  //
  FILE5_CHECKID(m_h5_obj,"attrib_set");

  //
  attrib_delete(attrib_name);

  //
  //printf("attrib_set_str('%s','%s')\n",attrib_name.c_str(),val.c_str());

  // make it just the right size.
  attrib_dtype=H5Tcopy(H5T_C_S1);
  H5Tset_size(attrib_dtype,val.size());
  // tmp ptr
  const char* val_ptr=val.c_str();
  //
  attrib_space=H5Screate(H5S_SCALAR);
  attrib_id=H5Acreate(m_h5_obj,attrib_name.c_str(),attrib_dtype,attrib_space,H5P_DEFAULT);
  m_h5_status=H5Awrite(attrib_id,attrib_dtype,val_ptr);
  FILE5_CHECKRV(m_h5_status,"H5Awrite");
  //
  H5Sclose(attrib_space);
  H5Aclose(attrib_id);
  H5Tclose(attrib_dtype);
  //
  return 0;
}

//////////

int
affx::File5_Object::attrib_get_void(const std::string& attrib_name,hid_t h5_dtype,void* val_ptr)
{
  hid_t h5_attrib_space;
  hid_t h5_attrib_id;

  //
  // printf("attrib_get_void('%s',%p)\n",attrib_name.c_str(),val_ptr); // dbg

  //
  FILE5_CHECKID(m_h5_obj,"Invalid hdf5 object.");

  //
  h5_attrib_space=H5Screate(H5S_SCALAR);
  h5_attrib_id=H5Aopen_name(m_h5_obj,attrib_name.c_str());
  if (h5_attrib_id<0) {
    return -1;
  }
  m_h5_status=H5Aread(h5_attrib_id,h5_dtype,val_ptr);
  FILE5_CHECKRV(m_h5_status,"H5Aread");
  //
  H5Sclose(h5_attrib_space);
  H5Aclose(h5_attrib_id);
  //
  return 0;
}

int
affx::File5_Object::attrib_get(const std::string& attrib_name,int* val)
{
  int rv;
  rv=attrib_get_void(attrib_name,H5T_NATIVE_INT,val);
  // printf("attrib_get('%s')=%d\n",attrib_name.c_str(),*val); // dbg
  return rv;
}

int
affx::File5_Object::attrib_get(const std::string& attrib_name,size_t* val)
{
  int val_int;
  int rv;
  rv=attrib_get_void(attrib_name,H5T_NATIVE_INT,&val_int);
  // printf("attrib_get('%s')=%d\n",attrib_name.c_str(),(int)*val); // dbg
  *val=val_int;
  return rv;
}

int
affx::File5_Object::attrib_get(const std::string& attrib_name,double* val)
{
  return attrib_get_void(attrib_name,H5T_NATIVE_DOUBLE,val);
}

int
affx::File5_Object::attrib_get(const std::string& attrib_name,std::string* val)
{
  hid_t attrib_id;
  hid_t attrib_dtype;
  int attrib_dtype_size;
  int attrib_idx;
  char* tmp_buf=NULL;

  attrib_idx=find_attrib_idx(attrib_name);

  if (attrib_idx<0) {
    *val="";
    return 1;
  }

  //
  attrib_id=H5Aopen_idx(m_h5_obj,attrib_idx);
  attrib_dtype=H5Aget_type(attrib_id);
  attrib_dtype_size=H5Tget_size(attrib_dtype);
  if (attrib_dtype_size==H5T_VARIABLE) {
    // variable size == HDF5 will allocate tmp_buf for us with malloc
    H5Aread(attrib_id,attrib_dtype,&tmp_buf);
  } 
  else {
    // known size == we allocate and HDF5 fills in.
    tmp_buf=(char*)malloc(attrib_dtype_size+1);
    FILE5_ASSERT(tmp_buf!=NULL);
    H5Aread(attrib_id,attrib_dtype,tmp_buf);
  }
  //
  val->assign(tmp_buf,attrib_dtype_size);
  //
  free(tmp_buf);
  H5Tclose(attrib_dtype);
  H5Aclose(attrib_id);

  return 0;
}

//////////

int
affx::File5_Object::attrib_delete(const std::string& attrib_name)
{
  // doing a delete without it being present seems to raise an error.

  // look for it...
  int attrib_idx;
  attrib_idx=find_attrib_idx(attrib_name);

  // found it?
  if (attrib_idx>=0) {
    // delete it by name.
    H5Adelete(m_h5_obj,attrib_name.c_str());
  }
  //
  return 0;
}

////

int affx::File5_Object::get_objinfo(H5G_stat_t* statbuf)
{
  m_h5_status=H5Gget_objinfo(m_h5_obj,".",0,statbuf);
  return m_h5_status;
}

int affx::File5_Object::get_name(std::string* str)
{
  size_t str_size;

  // once to get the size
  str_size=H5Iget_name(m_h5_obj,NULL,0)+1; // null
  if (str_size==0) {
    *str="";
  }
  else {
    char* buf=(char*)malloc(str_size);
    FILE5_ASSERT(buf!=NULL);
    str_size=H5Iget_name(m_h5_obj,buf,str_size);
    *str=buf;
    free(buf);
  }
  return str_size;
}

//////////

void affx::File5_Object::setOptCompress(int level)
{
  m_opt_compress=level;
}

void affx::File5_Object::setOptChunkSize(int size)
{
  m_opt_chunksize=size;
}

void affx::File5_Object::setOptCrc(int state)
{
  m_opt_crc=state;
}


#ifdef FILE5_DEBUG_USAGE

namespace affx {
  typedef std::map<affx::File5_Object*,int> file5_usecnt_t;
}

affx::file5_usecnt_t file5_usecnt;

void affx::file5_usecnt_inc(affx::File5_Object* ptr)
{
  affx::file5_usecnt_t::iterator i;
  int cnt=-1;

  i=file5_usecnt.find(ptr);
  if (i==file5_usecnt.end()) {
    file5_usecnt.insert(std::pair<affx::File5_Object*,int>(ptr,1));
    cnt=1;
  }
  else {
    i->second++;
    cnt=i->second++;
  }
  //
  printf("### file5_usecnt_inc: %c:%p=%3d '%s'\n",ptr->m_kind_char,ptr,cnt,ptr->m_name.c_str());
  fflush(NULL);
}

void affx::file5_usecnt_dec(affx::File5_Object* ptr)
{
  affx::file5_usecnt_t::iterator i;
  int cnt=-1;

  i=file5_usecnt.find(ptr);
  if (i==file5_usecnt.end()) {
    printf("### file5_usecnt_dec: %p: not found!\n",ptr);
  }
  else {
    i->second--;
    cnt=i->second;
    if (i->second==0) {
      file5_usecnt.erase(i);
    }
  }
  
  printf("### file5_usecnt_dec: %c:%p=%3d '%s'\n",ptr->m_kind_char,ptr,cnt,ptr->m_name.c_str());  
  fflush(NULL);
}

void affx::file5_usecnt_print()
{
  printf("========== File5 Usage ==========\n");
  int cnt=0;
  for (affx::file5_usecnt_t::iterator i=file5_usecnt.begin();i!=file5_usecnt.end();i++) {
    affx::File5_Object* ptr=i->first;
    printf("%3d: %c:%p=%3d '%s'\n",cnt++,ptr->m_kind_char,ptr,i->second,ptr->m_name.c_str());
  }
  printf("====================\n");
  fflush(NULL);
}

#endif
