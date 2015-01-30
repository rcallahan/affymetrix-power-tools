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
// ~/Affy/apt2/trunk/Bboard.cpp ---
//
// $Id: Bboard.cpp,v 1.9 2009-11-05 20:41:39 harley Exp $
//

//
#include "bboard/Bboard.h"
//
#include "bboard/dao/Dao.h"
//
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

Bboard::Bboard()
{
  init();
}
Bboard::Bboard(const std::string& name)
{
  init();
  m_name=name;
}

Bboard::Bboard(Bboard* bb_orig)
{
  m_name=bb_orig->m_name+"_copy";
  copyFrom(bb_orig);
}

Bboard::~Bboard()
{
  releaseAll();
  clear();
}

void Bboard::init()
{
  m_name="";
  m_parent=NULL;
}

void Bboard::clear()
{
  clearMap();
  clearArray();
}
void Bboard::clearArray()
{
  arrayResize(0);
}
void Bboard::clearMap()
{
  for (NameRefMap_t::iterator i=m_name_ref_map.begin();
       i!=m_name_ref_map.end();
       ++i) {
    // dont unbind() each one as we are going to clear this tree.
    // just NULL out the entries and then clear it.
    i->second.assign(NULL);
  }

  m_name_ref_map.clear();
}

void Bboard::dump()
{
  int cnt;
  printf("\n== Bboard::dump(): name='%s' this=%p\n",m_name.c_str(),this);
  printf("== Map:\n");
  if (m_name_ref_map.size()==0) {
    printf("  None.\n");
  }
  else {
    cnt=0;
    for (NameRefMap_t::iterator i=m_name_ref_map.begin();
         i!=m_name_ref_map.end();
         ++i) {
      printf("  %2d : %-20s ",cnt++,i->first.c_str());
      i->second->dump();
    }
  }

  //
  printf("== Array Values:\n");
  for (int i=0;i<m_ref_array.size();i++) {
    printf("  %2d : ",i);
    m_ref_array[i]->dump();
  }

  //
  printf("== Extra Values:\n");
  if (m_resource_refs.size()==0) {
    printf("  None.\n");
  }
  else {
    cnt=0;
    for (int i=0;i<m_resource_refs.size();i++) {
      printf(" %02d : ",cnt);
      m_resource_refs[i]->dump();
    }
  }
  //
  printf("\n");
}

//////////

int Bboard::arraySize() {
  return m_ref_array.size();
}

void Bboard::arrayClearAfterIdx(int idx) {
  for (int i=idx;i<m_ref_array.size();i++) {
    m_ref_array[i].assign(NULL);
  }
}
void Bboard::arrayResize(int newsize) {
  arrayClearAfterIdx(newsize);
  m_ref_array.resize(newsize);
}

//////////

BboardBox* Bboard::allocBboardBox()
{
  BboardBox* valptr=new BboardBox();
  m_resource_refs.push_back(BboardBoxRef());
  m_resource_refs[m_resource_refs.size()-1].assign(valptr);

  return valptr;
}

/////

void* Bboard::Malloc(int size)
{
  void* ptr=malloc(size);
  if (ptr!=NULL) {
    m_resource_mallocs.push_back(ptr);
  }
  //printf("Bboard::Malloc(%d)==%p\n",size,ptr);
  return ptr;
}

AptErr_t Bboard::Free(void* ptr) {
  // i=m_resource_mallocs.find(ptr)
  // if (i!=m_resource_mallocs.end()) {
  //   m_resource_mallocs.erase(i)
  // }
  for (ResMemVec_t::iterator i=m_resource_mallocs.begin();
       i!=m_resource_mallocs.end();
       i++) {
    if (*i==ptr) {
      m_resource_mallocs.erase(i);
      free(ptr);
      return APT_OK;
    }
  }
  // error: not found
  return APT_ERR;
}

//////////


void Bboard::releaseValRefs() {
  // zero our our value references.
  for (int i=0;i<m_resource_refs.size();i++) {
    m_resource_refs[i].assign(NULL);
  }
  m_resource_refs.clear();
}

void Bboard::releaseMalloc() {
  // free memory
  for (int i=0;i<m_resource_mallocs.size();i++) {
    if (m_resource_mallocs[i]!=NULL) {
      //printf("Bboard::releaseResources: free(%p)\n",m_resource_mallocs[i]);
      free(m_resource_mallocs[i]);
      m_resource_mallocs[i]=NULL;
    }
  }
  m_resource_mallocs.clear();
}

void Bboard::releaseAll()
{
  releaseValRefs();
  releaseMalloc();
}

//////////

Bboard* Bboard::getParentBboard() {
  return m_parent;
}
Bboard* Bboard::setParentBboard(Bboard* parent) {
  m_parent=parent;
  return m_parent;
}

Bboard* Bboard::makeChildBboard(const std::string& name) {
  Bboard* child=new Bboard(name);
  assert(child!=NULL);
  child->setParentBboard(this);
  // keep track of resources which might not be named,
  // but we want to track.
  BboardBox* val=allocBboardBox();
  val->setPtrAndType((void*)child,BBT_BBOARD);
  //
  return child;
}

//////////

AptErr_t Bboard::copyFrom(Bboard* bb_from)
{
  if (bb_from==this) {
    printf("copyFrom: wont copy from myself.");
    assert(0);
  }

  for (NameRefMap_t::iterator i=bb_from->m_name_ref_map.begin();
       i!=bb_from->m_name_ref_map.end();
       ++i) {
    BboardBox* bbval=i->second.boxPtr();
    bind(i->first,bbval);
  }
  return APT_OK;
}

//////////

AptErr_t Bboard::bind(const std::string& name,BboardBox* valptr)
{
  NameRefMap_t::iterator i;
  i=m_name_ref_map.find(name);
  if (i==m_name_ref_map.end()) {
    m_name_ref_map.insert(std::pair<std::string,BboardBoxRef>(name,valptr));
  }
  else {
    i->second.assign(valptr);
  }
  //
  return APT_OK;
}

AptErr_t Bboard::bind(const std::string& name,BboardBoxRef& boxref)
{
  return bind(name,boxref.boxPtr());
}

AptErr_t Bboard::unbind(const std::string& name)
{
  NameRefMap_t::iterator i;
  i=m_name_ref_map.find(name);
  if (i==m_name_ref_map.end()) {
    return APT_ERR;
  }
  // this does the dec
  i->second.assign(NULL);
  m_name_ref_map.erase(i);
  //
  return APT_OK;
}

//
AptErr_t Bboard::rename(const std::string& f_name,const std::string& t_name)
{
  BboardBox* valptr=findBoxPtr(f_name,0,0);
  if (valptr==NULL) {
    return APT_ERR;
  }
  //
  bind(t_name,valptr);
  unbind(f_name);
  //
  return APT_OK;
}

//////////

// @todo mkdir

BboardBox* Bboard::findBoxPtr(const std::string& name) {
  return findBoxPtr(name,0,1);
}

// @todo handle "one/two/three" names.
BboardBox* Bboard::findBoxPtr(const std::string& name,int doCreate,int doParents)
{
  BboardBox* vptr=NULL;

  NameRefMap_t::iterator i;
  i=m_name_ref_map.find(name);
  if (i!=m_name_ref_map.end()) {
    return i->second.boxPtr();
  }
  // search the parents if asked
  if (doParents==1) {
    if (m_parent!=NULL) {
      return m_parent->findBoxPtr(name,doCreate,doParents);
    }
  }
  // create always creates in this Bboard, no the parents.
  if (doCreate==1) {
    vptr=new BboardBox(name);
    bind(name,vptr);
    return vptr;
  }
  //
  return NULL;
}


AptErr_t Bboard::getPtrIfType(const std::string& name,void*& dataptr,BboardType_t tcode)
{
  BboardBox* vptr=findBoxPtr(name,0,1);
  dataptr=NULL;
  if (vptr==NULL) {
    return APT_ERR_NOTFOUND;
  }
  if (vptr->getBbType()!=tcode) {
    return APT_ERR_WRONGTYPE;
  }
  //
  dataptr=vptr->getPtr();
  //
  return APT_OK;
}

AptErr_t Bboard::setPtrAndType(const std::string& name,void* dataptr,BboardType_t tcode)
{
  BboardBox* vptr=findBoxPtr(name,1,0);
  assert(vptr!=NULL);
  return vptr->setPtrAndType(dataptr,tcode);
}

/////

#define BB_GETSET_SIMPLE(_bbt,_type,_unsetvalue)                        \
  AptErr_t Bboard::get(const std::string& name,_type* val,int abortOnErr) \
  {                                                                     \
    AptErr_t rv;                                                        \
    BboardBox* boxptr=findBoxPtr(name,0,1);                             \
    if (boxptr==NULL) {                                                 \
      *val=_unsetvalue;                                                 \
      rv=APT_ERR_NOTFOUND;                                              \
      APT_ABORT_ON_ERR(rv,"value for: '" + name + "'not found");        \
      return rv;                                                        \
    }                                                                   \
    rv=boxptr->getData(val);                                            \
    APT_ABORT_ON_ERR(rv,"Couldn't get data for: '" + name + "'");       \
    return rv;                                                          \
  }                                                                     \
  AptErr_t Bboard::set(const std::string& name,_type val,int abortOnErr) \
  {                                                                     \
    AptErr_t rv;                                                        \
    BboardBox* boxptr=findBoxPtr(name,1,0);                             \
    rv=boxptr->setData(val);                                            \
    APT_ABORT_ON_ERR(rv,"getData");                                     \
    return rv;                                                          \
  }                                                                     \
  AptErr_t Bboard::get(int idx,_type* val,int abortOnErr)               \
  {                                                                     \
    AptErr_t rv;                                                        \
    if (!(idx>=0)&&(idx<m_ref_array.size())) {                          \
      *val=_unsetvalue;                                                 \
      rv=APT_ERR_OUTOFBOUNDS;                                           \
      APT_ABORT_ON_ERR(rv,"out of bounds.");                            \
      return rv;                                                        \
    }                                                                   \
    if (m_ref_array[idx].isNullRef()) {                                 \
      *val=_unsetvalue;                                                 \
      rv=APT_ERR_ISNULL;                                                \
      APT_ABORT_ON_ERR(rv,"is null");                                   \
      return rv;                                                        \
    }                                                                   \
    rv=m_ref_array[idx]->getData(val);                                  \
    APT_ABORT_ON_ERR(rv,"getData");                                     \
    return rv;                                                          \
  }                                                                     \
  AptErr_t Bboard::set(int idx,_type val,int abortOnErr)                \
  {                                                                     \
    AptErr_t rv;                                                        \
    if (!(idx>=0)&&(idx<m_ref_array.size())) {                          \
      rv=APT_ERR_OUTOFBOUNDS;                                           \
      APT_ABORT_ON_ERR(rv,"out of bounds.");                            \
      return rv;                                                        \
    }                                                                   \
    rv=m_ref_array[idx]->setData(val);                                  \
    APT_ABORT_ON_ERR(rv,"setData");                                     \
    return rv;                                                          \
  }

//
BB_GETSET_SIMPLE(BBT_INT   ,int, 0  );
BB_GETSET_SIMPLE(BBT_DOUBLE,double ,0.0 );
BB_GETSET_SIMPLE(BBT_STRING,std::string,"");

// This will be redefined with extra checks later.
#define EXTRA_SET_CHECK() {}

// Like the above, but with pointers,
// which means we should check to see if there is a box
// for this pointer already.
#define BB_GETSET_PTR(_bbt,_type)                                       \
  AptErr_t Bboard::get(const std::string& name,_type* val,int abortOnErr) \
  {                                                                     \
    AptErr_t rv;                                                        \
    BboardBox* boxptr=findBoxPtr(name,0,1);                             \
    if (boxptr==NULL) {                                                 \
      *val=NULL;                                                        \
      rv=APT_ERR_NOTFOUND;                                              \
      APT_ABORT_ON_ERR(rv,"not found");                                 \
      return rv;                                                        \
    }                                                                   \
    rv=boxptr->getData(val);                                            \
    APT_ABORT_ON_ERR(rv,"getData");                                     \
    return rv;                                                          \
  }                                                                     \
  AptErr_t Bboard::set(const std::string& name,_type val,int abortOnErr) \
  {                                                                     \
    AptErr_t rv;                                                        \
    EXTRA_SET_CHECK();                                                  \
    BboardBox* boxptr=NULL;                                             \
    boxptr=BboardBox::findBoxForPtr((void*)val);                        \
    if (boxptr==NULL) {                                                 \
      boxptr=findBoxPtr(name,1,0);                                      \
    }                                                                   \
    rv=boxptr->setData(val);                                            \
    APT_ABORT_ON_ERR(rv,"setData");                                     \
    return rv;                                                          \
  }                                                                     \
  AptErr_t Bboard::get(int idx,_type* val,int abortOnErr)               \
  {                                                                     \
    AptErr_t rv;                                                        \
    if (!(idx>=0)&&(idx<m_ref_array.size())) {                          \
      *val=NULL;                                                        \
      rv=APT_ERR_OUTOFBOUNDS;                                           \
      APT_ABORT_ON_ERR(rv,"out of bounds");                             \
      return rv;                                                        \
    }                                                                   \
    if (m_ref_array[idx].isNullRef()) {                                 \
      *val=NULL;                                                        \
      rv=APT_ERR_ISNULL;                                                \
      APT_ABORT_ON_ERR(rv,"is null");                                   \
      return rv;                                                        \
    }                                                                   \
    rv=m_ref_array[idx]->getData(val);                                  \
    APT_ABORT_ON_ERR(rv,"getData");                                     \
    return rv;                                                          \
  }                                                                     \
  AptErr_t Bboard::set(int idx,_type val,int abortOnErr)                \
  {                                                                     \
    AptErr_t rv;                                                        \
    EXTRA_SET_CHECK();                                                  \
    if (!(idx>=0)&&(idx<m_ref_array.size())) {                          \
      rv=APT_ERR_OUTOFBOUNDS;                                           \
      APT_ABORT_ON_ERR(rv,"out of bounds");                             \
      return rv;                                                        \
    }                                                                   \
    BboardBox* boxptr=NULL;                                             \
    boxptr=BboardBox::findBoxForPtr((void*)val);                        \
    if (boxptr==NULL) {                                                 \
      m_ref_array[idx]=new BboardBox();                                 \
    }                                                                   \
    rv=m_ref_array[idx]->setData(val);                                  \
    APT_ABORT_ON_ERR(rv,"setData");                                     \
    return rv;                                                          \
  }

//
//BB_GETSET_PTR(BBT_CHAR_PTR,char*);
//
BB_GETSET_PTR(BBT_VEC_CHAR  ,std::vector<char>*);
BB_GETSET_PTR(BBT_VEC_INT   ,std::vector<int>*);
BB_GETSET_PTR(BBT_VEC_FLOAT ,std::vector<float>*);
BB_GETSET_PTR(BBT_VEC_DOUBLE,std::vector<double>*);

// This is special as we want to add a simple extra check.
// we are looking for simple loops in bboards.
#undef EXTRA_SET_CHECK
#define EXTRA_SET_CHECK() { if (val==this) { printf("you just made a simple BB loop!\n"); assert(0); } }
BB_GETSET_PTR(BBT_BBOARD,Bboard*);
#undef EXTRA_SET_CHECK
#define EXTRA_SET_CHECK() {}
BB_GETSET_PTR(BBT_BBOARDOBJ,BboardObj*);

//
BB_GETSET_PTR(BBT_BIOSPECIES,BioSpecies*);
BB_GETSET_PTR(BBT_CHIPLAYOUT,ChipLayout*);
BB_GETSET_PTR(BBT_PGOPTIONS,PgOptions*);

//
BB_GETSET_PTR(BBT_FSPATH,FsPath*);
BB_GETSET_PTR(BBT_DAO_FILE,Dao_File*);
BB_GETSET_PTR(BBT_DAO_GROUP,Dao_Group*);
BB_GETSET_PTR(BBT_DAO_TABLE,Dao_Table*);

// BBT__NEWTYPE
