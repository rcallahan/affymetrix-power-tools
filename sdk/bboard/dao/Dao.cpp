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
// ~/Affy-projects/apt2/trunk/Dao.cpp ---
//
// $Id: Dao.cpp,v 1.5 2009-11-04 20:37:32 harley Exp $
//

//
#include "bboard/dao/Dao.h"
//
#include "bboard/dao/Dao_TsvFile.h"
//
#include <stdio.h>

//////////

#define DAO_REQIRE_DRIVER() { assert(m_driver!=NULL); }

//////////

Dao_Driver::Dao_Driver() {
}

Dao_Driver::~Dao_Driver() {
}

//////////

Dao_File::Dao_File() {
  init();
}

Dao_File::~Dao_File() {
  init();
}

void Dao_File::init() {
  m_driver=NULL;
  m_filefmt=FsPath::FILEFMT_NONE;
  m_box_refcnt=0;
  m_dao_path.clear();
  m_flags=0;
}

int Dao_File::getFormat() {
  return m_filefmt;
}
AptErr_t Dao_File::setFormat(int fmt)
{
  if (fmt==m_filefmt) {
    return APT_OK;
  }
  if (m_driver!=NULL) {
    this->close();
  }

  switch (fmt) {
  case FsPath::FILEFMT_TSVFILE:
    m_driver=new Dao_TsvFile_Driver();
    m_driver->attach(this);
    break;
  case FsPath::FILEFMT_FILE5:
    //m_driver=new Dao_Driver_File5();
    // m_driver->attach(this);
    assert(0);
    break;
  case FsPath::FILEFMT_NONE:
    m_driver=NULL;
    break;
  default:
    break;
  }
  //
  m_filefmt=fmt;
  return APT_OK;
}

AptErr_t Dao_File::create(const std::string& name,int flags) {
  return create(FsPath(name),flags);
}
AptErr_t Dao_File::create(const FsPath& dao_p,int flags) {
  this->close();
  //
  m_flags=flags;
  m_dao_path=dao_p;
  setFormat(m_dao_path.getFileFmt());
  //
  DAO_REQIRE_DRIVER();
  m_driver->create();
  //
  return APT_OK;
}
AptErr_t Dao_File::open(const std::string& path,int flags) {
  return open(FsPath(path),flags);
}
AptErr_t Dao_File::open(const FsPath& dao_p,int flags) {
  this->close();
  //
  m_flags=flags;
  m_dao_path=dao_p;
  setFormat(m_dao_path.getFileFmt());
  //
  DAO_REQIRE_DRIVER();
  m_driver->open();
  //
  return APT_OK;
}

AptErr_t Dao_File::close()
{
  if (m_driver!=NULL) {
    m_driver->close();
    delete m_driver;
    m_driver=NULL;
  }
  return APT_OK;
}

Dao_Group* Dao_File::createGroup(const std::string& name,int flags)
{
  DAO_REQIRE_DRIVER();
  return m_driver->openGroup(name,flags);
}

Dao_Group* Dao_File::openGroup(const std::string& name,int flags)
{
  DAO_REQIRE_DRIVER();
  return m_driver->openGroup(name,flags);
}

void Dao_File::addBoxRef(const BboardBox* box) {
  this->addBoxRef(box,this);
}
void Dao_File::delBoxRef(const BboardBox* box) {
  this->delBoxRef(box,this);
}
void Dao_File::addBoxRef(const BboardBox* box,void* dao_obj_ptr) {
  /// @todo
}
void Dao_File::delBoxRef(const BboardBox* box,void* dao_obj_ptr) {
  /// @todo
}
