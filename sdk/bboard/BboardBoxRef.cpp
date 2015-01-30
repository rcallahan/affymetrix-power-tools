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
//  ~/Affy/apt2/trunk/BboardBoxRef.cpp ---
//
//  $Id: BboardBoxRef.cpp,v 1.4 2009-10-31 23:23:02 harley Exp $
//

#include "BboardBoxRef.h"

BboardBoxRef::BboardBoxRef() {
  init();
}

BboardBoxRef::BboardBoxRef(const BboardBox* valptr) {
  init();
  assign(valptr);
}
BboardBoxRef::BboardBoxRef(const BboardBoxRef& boxref) {
  init();
  assign(boxref.m_boxptr);
}
BboardBoxRef::~BboardBoxRef() {
  assign(NULL);
};
//
void BboardBoxRef::init() {
  m_boxptr=NULL;
}
void BboardBoxRef::clear() {
  assign(NULL);
}
void BboardBoxRef::dump() {
  printf("BboardBoxRef(%p)==%p",this,m_boxptr);
}


bool BboardBoxRef::isNull() {
  if (m_boxptr==NULL) {
    return true;
  }
  if (m_boxptr->getBbType()<=BBT__START_OF_TYPES) {
    return true;
  }
  //
  return false;
}
bool BboardBoxRef::isNullRef() {
  if (m_boxptr==NULL) {
    return true;
  }
  return false;
}


// @todo
// all this should be renamed to make it clear:
// Ref -> Var -> Val
// Var* shouldnt be leaked to non-Refs.
const BboardBoxRef& BboardBoxRef::assign(const BboardBox* valptr) {
  if (0) {
    printf("%p: BboardBoxRef::assign(%p) refcnt=%d\n",
           this,valptr,
           ((valptr==NULL)?0:valptr->getRefcnt()));
  }
  //
  if (m_boxptr!=NULL) {
    // last one out, turns out the lights.
    if (m_boxptr->refcntDec()==0) {
      delete m_boxptr;
      m_boxptr=NULL;
    }
  }
  // this should be safe. Would mutable work?
  m_boxptr=const_cast<BboardBox*>(valptr);
  if (m_boxptr!=NULL) {
    m_boxptr->refcntInc();
  }
  //return m_boxptr;
  return *this;
}

const BboardBoxRef& BboardBoxRef::assign(const BboardBoxRef& boxref) {
  return this->assign(boxref.m_boxptr);
}
const BboardBoxRef& BboardBoxRef::operator=(const BboardBox* valptr) {
  return this->assign(valptr);
}
const BboardBoxRef& BboardBoxRef::operator=(const BboardBoxRef& boxref) {
  return this->assign(boxref.m_boxptr);
}

// dereferencing it
BboardBox* BboardBoxRef::operator->() {
  if (m_boxptr==NULL) {
    BboardBox* newboxptr=new BboardBox();
    assign(newboxptr);
    // now newboxptr->refcnt()==1
  }
  return m_boxptr;
}

void BboardBoxRef::swapBoxPtrs(BboardBoxRef& boxref) {
  // when doing this we dont need to inc/dec the refcnts.
  BboardBox* tmp_valptr1=m_boxptr;
  BboardBox* tmp_valptr2=boxref.m_boxptr;
  m_boxptr=tmp_valptr2;
  boxref.m_boxptr=tmp_valptr1;
}

BboardBox* BboardBoxRef::boxPtr() {
  return m_boxptr;
}

BboardType_t BboardBoxRef::getBbType() {
  if (m_boxptr==NULL) {
    return BBT_NULLREF;
  }
  return m_boxptr->getBbType();
}
