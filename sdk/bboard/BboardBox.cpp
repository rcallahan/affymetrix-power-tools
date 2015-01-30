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
// ~/Affy/apt2/trunk/BboardBox.cpp ---
//
// $Id: BboardBox.cpp,v 1.6 2009-11-05 20:41:39 harley Exp $
//

//
#include "bboard/BboardBox.h"
//
#include "bboard/Apt2Types.h"
#include "bboard/Bboard.h"
//
// list of types which will get used.
#include "bboard/Bboard_ExampleTestObj.h"
//
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <utility>
#include <vector>


//////////

// allocate the static member.
BboardBox::boxForPtr_t BboardBox::m_boxForPtr_map;

//////////

void BboardBox::addBoxForPtr(void* ptr,BboardBox* box)
{
  // Yawn... trivial case.
  if (ptr==NULL) {
    return;
  }
  // Take the address of something on our stack.
  // if the ptr is bigger, then it is ABOVE our stack location, 
  // thus we think it is on our stack and NOT heap allocated.
  /*
    // This was commented out in r13234 by Chuck.
    // The comment was about getting this to work on windows.
    // Not sure why it wouldnt.
   int an_int_on_our_stack=0;
   if (ptr>&an_int_on_our_stack) {
    // which is bad bad bad! Whine about it...
    printf("BboardBox::addBoxForPtr(%p,%p): This pointer might be on the stack!\n",ptr,box);
    // ...and die.
    assert(0);
  }*/

  // We do not want to add a second box for the same pointer.
  // Check to see if there is a box with this pointer already.
  std::pair<boxForPtr_t::iterator,bool> rv_pair;
  rv_pair=m_boxForPtr_map.insert(std::pair<void*,BboardBox*>(ptr,box));
  // what we could do is reuse the box.
  if (rv_pair.second==false) {
    printf("BboardBox::addBoxForPtr(%p,%p): value already has a box. (%p)\n",ptr,box,rv_pair.first->second);
    assert(0);
  }
}

void BboardBox::delBoxForPtr(void* ptr,BboardBox* box)
{
  // trivial case.
  if (ptr==NULL) {
    return;
  }
  //
  boxForPtr_t::iterator i;
  i=m_boxForPtr_map.find(ptr);
  if (i==m_boxForPtr_map.end()) {
    printf("BboardBox::rmBoxForPtr(%p,%p): value does not exist.\n",ptr,box);
    assert(0);
  }
  if (i->second!=box) {
    printf("BboardBox::rmBoxForPtr(%p,%p):\n",ptr,box);
    assert(0);
  }
  m_boxForPtr_map.erase(i);
}

BboardBox* BboardBox::findBoxForPtr(void* ptr)
{
  boxForPtr_t::iterator i=m_boxForPtr_map.find(ptr);
  if (i==m_boxForPtr_map.end()) {
    return NULL;
  }
  return i->second;
}

//////////

BboardBox::BboardBox()
{
  init();
}
BboardBox::BboardBox(const std::string& name)
{
  init();
  m_name=name;
  m_data_owned=true;
}

BboardBox::BboardBox(const std::string& name,void* dataptr,BboardType_t tcode)
{
  init();
  m_name=name;
  setPtrAndType(dataptr,tcode);
}

BboardBox::~BboardBox()
{
  //printf("BboardBox::~BboardBox(): delete %p\n",this);
  clearData();
}

void BboardBox::init()
{
  m_name="";
  m_refcnt=0;
  m_data_ptr=NULL;
  m_data_tcode=BBT_UNSET;
}
  
// For simple cases, use this macro,
// it will delete the data right away.
#define CLEARDATA_CASE_FREE(_bbt)               \
  case _bbt:                                    \
  if (m_data_owned) {                           \
     free(m_data_ptr);                          \
  }                                             \
  break;

// For simple cases, use this macro,
// it will delete the data right away.
#define CLEARDATA_CASE_DEL(_bbt,_type) \
  case _bbt:                           \
  if (m_data_owned) {                  \
    delete (_type*)m_data_ptr;         \
  }                                    \
  break;

// Dont delete DAO object, just tell them they are no longer
// in a BboardBox.
#define CLEARDATA_CASE_DAO(_bbt,_type)     \
  case _bbt:                               \
  ((_type*)m_data_ptr)->delBoxRef(this); \
  break;

void BboardBox::clearData()
{
  if (m_data_ptr==NULL) {
    return;
  }
  switch (m_data_tcode) {
  case BBT_UNSET:
  case BBT_ERR:
      // Nothing to delete
    assert(m_data_ptr==NULL);
    break;
      //
    CLEARDATA_CASE_FREE(BBT_INT);
    CLEARDATA_CASE_FREE(BBT_DOUBLE);
    CLEARDATA_CASE_DEL(BBT_STRING    ,std::string);
    //
    CLEARDATA_CASE_DEL(BBT_VEC_CHAR  ,std::vector<char>);
    CLEARDATA_CASE_DEL(BBT_VEC_INT   ,std::vector<int>);
    CLEARDATA_CASE_DEL(BBT_VEC_FLOAT ,std::vector<float>);
    CLEARDATA_CASE_DEL(BBT_VEC_DOUBLE,std::vector<double>);
    //
    CLEARDATA_CASE_DEL(BBT_BIOSPECIES,BioSpecies);
    CLEARDATA_CASE_DEL(BBT_CHIPLAYOUT,ChipLayout);
    CLEARDATA_CASE_DEL(BBT_PGOPTIONS,PgOptions);
    //
    CLEARDATA_CASE_DEL(BBT_BBOARD     ,Bboard);
    //
    CLEARDATA_CASE_DEL(BBT_BBOARDOBJ  ,BboardObj);
    //
    CLEARDATA_CASE_DEL(BBT_FSPATH   ,FsPath);
    CLEARDATA_CASE_DAO(BBT_DAO_FILE   ,Dao_File);
    CLEARDATA_CASE_DAO(BBT_DAO_GROUP  ,Dao_Group);
    CLEARDATA_CASE_DAO(BBT_DAO_TABLE  ,Dao_Table);
    //
    //CLEARDATA_CASE_DEL(BBT_EXAMPLETESTOBJ,ExampleTestObj);
    
    // How should you get rid of the new type?
    // BBT__NEWTYPE
    
  default:
    printf("BboardBox::clearData(): unhandled type: '%d'\n",m_data_tcode);
    assert(0);
    break;
  }
  // do the moral equiv of:
  // setPtrAndType(NULL,BBT_UNSET);
  // but without the call
  delBoxForPtr(m_data_ptr,this);
  m_data_ptr=NULL;
  m_data_tcode=BBT_UNSET;
}

void BboardBox::dump()
{
  printf("Value: %p m_data_ptr=%p m_data_tcode=%d refcnt=%d\n",
         this,m_data_ptr,m_data_tcode,m_refcnt);
}

std::string BboardBox::asString()
{
  char buf[100];

  switch (m_data_tcode) {
    //
  case BBT_CHAR:
    return std::string((char*)m_data_ptr);
    break;
  case BBT_INT:
    sprintf(buf,"%d",*(int*)m_data_ptr);
    return buf;
    break;
  case BBT_FLOAT:
    sprintf(buf,"%f",*(float*)m_data_ptr);
    return buf;
    break;
  case BBT_DOUBLE:
    sprintf(buf,"%f",*(double*)m_data_ptr);
    return buf;
    break;
  case BBT_CHAR_PTR:
    return (char*)m_data_ptr;
    break;
  case BBT_STRING:
    return *(std::string*)m_data_ptr;
    break;
    //////////
  case BBT_PGOPTIONS:
    return "@todo filename-for-PgOptions";
    break;
    //////////
  case BBT_FSPATH:
    return ((FsPath*)m_data_ptr)->asString();
  case BBT_DAO_FILE:
    return std::string("Dao_File: ")+((Dao_File*)m_data_ptr)->m_dao_path.asString();
  case BBT_DAO_GROUP:
    return std::string("Dao_Group: ")+((Dao_Group*)m_data_ptr)->getName();
  case BBT_DAO_TABLE:
    return std::string("Dao_Table: ")+((Dao_Table*)m_data_ptr)->getName();
    //////////
  case BBT_EXAMPLETESTOBJ:
    return std::string("ExampleTestObj: ");
    // how should you convert this new type to a string?
    // BBT__NEWTYPE
  default:
    return "";
    break;
  }
  //
  return "";
}

int BboardBox::refcntStep(int cnt)
{
  m_refcnt+=cnt;
  if (m_refcnt<0) {
    assert(0);
  }
  return m_refcnt;
}

// @todo All sets of m_data_ptr should be funneled though here.
AptErr_t BboardBox::setPtrAndType(void* ptr,BboardType_t tcode)
{
  clearData();
  delBoxForPtr(m_data_ptr,this);
  m_data_ptr=ptr;
  m_data_tcode=tcode;
  // we own the data by default.
  m_data_owned=true;
  //
  addBoxForPtr(m_data_ptr,this);

  // should something special happen when we have a new reference? generally not.
  // however some collections of objects like DAOs, want to know when they have
  // a new reference.  So we tell them.
  // If we could derive all object from our own "Object" we would just call 
  // some virtual functions.
  // @todo Actually, we could do that now.
  
  switch (m_data_tcode) {
    // Tell the DAO we have a reference to it.
  case BBT_DAO_FILE:
    ((Dao_File*)m_data_ptr)->addBoxRef(this);
    break;
  case BBT_DAO_GROUP:
    ((Dao_Group*)m_data_ptr)->addBoxRef(this);
    break;
  case BBT_DAO_TABLE:
    ((Dao_Table*)m_data_ptr)->addBoxRef(this);
    break;
    // add new types here
    // BBT__NEWTYPE
  default:
    break;
  }
   
  //
  return APT_OK;
}

AptErr_t BboardBox::getPtrIfType(void*& dataptr,BboardType_t tcode)
{
  dataptr=NULL;
  if (m_data_tcode!=tcode) {
    return APT_ERR_WRONGTYPE;
  }
  dataptr=m_data_ptr;
  return APT_OK;
}

/////

AptErr_t BboardBox::getData(std::string* val)
{
  if (m_data_tcode!=BBT_STRING) {
    return APT_ERR_WRONGTYPE;
  }
  *val=*(std::string*)m_data_ptr;
  return APT_OK;
}
AptErr_t BboardBox::setData(const std::string& val)
{
  if (m_data_tcode!=BBT_STRING) {
    std::string* tmp_ptr=new std::string();
    setPtrAndType(tmp_ptr,BBT_STRING);
  }
  ((std::string*)m_data_ptr)->assign(val);
  return APT_OK;
}

#define BBBOX_GETSET_SIMPLE(_bbt,_type)             \
  AptErr_t BboardBox::getData(_type* val)           \
  {                                                 \
    if (m_data_tcode!=_bbt) {                       \
      return APT_ERR_WRONGTYPE;                     \
    }                                               \
    *val=*(_type*)m_data_ptr;                       \
    return APT_OK;                                  \
  }                                                 \
  AptErr_t BboardBox::setData(_type val)            \
  {                                                 \
    if (m_data_tcode!=_bbt) {                       \
      clearData();                                  \
      void* ptr=malloc(sizeof(_type));              \
      setPtrAndType(ptr,_bbt);                      \
    }                                               \
    *(_type*)m_data_ptr=val;                        \
    return APT_OK;                                  \
  }

BBBOX_GETSET_SIMPLE(BBT_INT,int);
BBBOX_GETSET_SIMPLE(BBT_DOUBLE,double);

////

#define BBBOX_GETSET_PTR(_bbt,_type)             \
  AptErr_t BboardBox::getData(_type* val)        \
  {                                              \
    *val=NULL;                                   \
    if (m_data_tcode!=_bbt) {                    \
      return APT_ERR_WRONGTYPE;                  \
    }                                            \
    *val=(_type)m_data_ptr;                      \
    return APT_OK;                               \
  }                                              \
  AptErr_t BboardBox::setData(_type val)         \
  {                                              \
    setPtrAndType(val,_bbt);                     \
    return APT_OK;                               \
  }

//BBBOX_GETSET_PTR(BBT_CHAR_PTR,char*);
BBBOX_GETSET_PTR(BBT_VEC_CHAR,std::vector<char>*);
BBBOX_GETSET_PTR(BBT_VEC_INT,std::vector<int>*);
BBBOX_GETSET_PTR(BBT_VEC_FLOAT,std::vector<float>*);
BBBOX_GETSET_PTR(BBT_VEC_DOUBLE,std::vector<double>*);
//
BBBOX_GETSET_PTR(BBT_BIOSPECIES,BioSpecies*);
BBBOX_GETSET_PTR(BBT_CHIPLAYOUT,ChipLayout*);
BBBOX_GETSET_PTR(BBT_PGOPTIONS,PgOptions*);
//
BBBOX_GETSET_PTR(BBT_BBOARD,Bboard*);
BBBOX_GETSET_PTR(BBT_BBOARDOBJ,BboardObj*);
//
BBBOX_GETSET_PTR(BBT_FSPATH,FsPath*);
BBBOX_GETSET_PTR(BBT_DAO_FILE,Dao_File*);
BBBOX_GETSET_PTR(BBT_DAO_GROUP,Dao_Group*);
BBBOX_GETSET_PTR(BBT_DAO_TABLE,Dao_Table*);

/// BBT__NEWTYPE
