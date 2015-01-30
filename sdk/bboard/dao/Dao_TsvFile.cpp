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
// ~/Affy-projects/apt2/trunk/Dao_TsvFile.cpp ---
//
// $Id: Dao_TsvFile.cpp,v 1.5 2009-11-05 20:41:39 harley Exp $
//

//
#include "bboard/Apt2Types.h"
#include "bboard/dao/Dao_TsvFile.h"
//
#include "assert.h"
#include <map>
#include <stdio.h>

enum TsvFile_State_t {
  TSV_STATE_UNSET=0,
  TSV_STATE_OPEN=1,
  TSV_STATE_CLOSED=2,
};

/// convert the TsvFile RV to a Dao RV.
#define CVT_RV(_rv) ((_rv==affx::TSV_OK)?APT_OK:APT_ERR)

#define CHECK_TSV_CLOSED() { }
#define ENSURE_TSV_OPEN() { }

/*
// #define DELETE_PTR_MAP(_iter,_mapname) {        \
//     for (_iter i=_mapname.begin();              \
//          i!=_mapname.end();                     \
//          ++i) {                                 \
//       delete i->second;                         \
//       i->second=NULL;                           \
//     }                                           \
//     _mapname.clear();                           \
//   }

// template<typename T1>
// AptErr_t delete_ptr_map(T1 map) {
//   T1 tmp_map;
//   //tmp_map=map;
//   T1::reverse_iterator i;
//   for (i=tmp_map.rbegin(); i!=tmp_map.rend(); ++i) {
//     delete i->second;
//   }
//   //
//   return APT_OK;
// };
// 
// template<typename T1>
// AptErr_t unbind_ptr_map(T1 map,T1::mapped_type key) {
//   for (T1::iterator i=map.rbegin(); i!=map.rend(); ++i) {
//     map.erase(i);
//   }
//   map.clear();
//   //
//   return APT_OK;
// }
*/

//////////

affx::tsv_type_t TsvDtype_From_DaoDtype(DaoDataType_t dao_dtype)
{
  switch (dao_dtype) {
  case DAO_CHAR:
  case DAO_STRING:
    return affx::TSV_TYPE_STRING;
    break;
  case DAO_SHORT:
  case DAO_INT:
    return affx::TSV_TYPE_INT;
    break;
  case DAO_FLOAT:
    return affx::TSV_TYPE_FLOAT;
    break;
  case DAO_DOUBLE:
    return affx::TSV_TYPE_DOUBLE;
    break;
  default:
    assert(0);
  }
  assert(0);
  return affx::TSV_TYPE_INT;
};

//////////

Dao_TsvFile_Driver::Dao_TsvFile_Driver() {
  init();
}
Dao_TsvFile_Driver::~Dao_TsvFile_Driver() {
  this->close();
}

void Dao_TsvFile_Driver::init() {
  //Dao_Driver::init();
  m_parent_file=NULL;
}

//
AptErr_t Dao_TsvFile_Driver::attach(Dao_File* dao_f) {
  m_parent_file=dao_f;
  return APT_OK;
}
AptErr_t Dao_TsvFile_Driver::create() {
  return APT_OK;
}
AptErr_t Dao_TsvFile_Driver::open() {
  return APT_OK;
}
AptErr_t Dao_TsvFile_Driver::close() {
  //return delete_ptr_map(name_group_map_iter_t,m_name_group_map);
  return APT_OK;
}
AptErr_t Dao_TsvFile_Driver::flush() {
  return APT_OK;
}

//
Dao_Group* Dao_TsvFile_Driver::createGroup(const std::string& name,int flags) {
  return openGroup(name,0);
}
Dao_Group* Dao_TsvFile_Driver::openGroup(const std::string& name,int flags) {
//  // have it already?
//  name_group_map_iter_t i=m_name_group_map.find(name);
//  if (i!=m_name_group_map.end()) {
//    return i->second;
//  }
  //
  Dao_TsvFile_Group* dg=new Dao_TsvFile_Group();
  dg->m_parent_file=m_parent_file;
  dg->m_parent_driver=this;
  dg->m_name=name;
  dg->m_flags=flags;
  //
  //m_name_group_map[name]=dg;
  //
  return dg;
}

//////////

Dao_TsvFile_Group::Dao_TsvFile_Group() {
  init();
}
Dao_TsvFile_Group::~Dao_TsvFile_Group() {
  this->close();
  //m_parent_file->unbindGroup(this);
}

void Dao_TsvFile_Group::init() {
  m_flags=0;
}

AptErr_t Dao_TsvFile_Group::close() {
  //return delete_ptr_map(name_table_map_iter_t,m_name_table_map);
  return APT_OK;
}
AptErr_t Dao_TsvFile_Group::flush() {
  return APT_OK;
}

void Dao_TsvFile_Group::delBoxRef(const BboardBox* box) {
  m_parent_file->delBoxRef(box,(void*)this);
}
void Dao_TsvFile_Group::addBoxRef(const BboardBox* box) {
  m_parent_file->delBoxRef(box,(void*)this);
}

Dao_Table* Dao_TsvFile_Group::createTable(const std::string& name,int flags) {
  Dao_TsvFile_Table* dt=new Dao_TsvFile_Table();
  dt->m_parent_file=this->m_parent_file;
  dt->m_parent_driver=this->m_parent_driver;
  dt->m_parent_group=this;
  dt->m_name=name;
  dt->m_flags=flags;
  //
  //m_name_table_map[name]=dt;
  //
  return dt;
}
Dao_Table* Dao_TsvFile_Group::openTable(const std::string& name,int flags) {
  //
  Dao_TsvFile_Table* dt=new Dao_TsvFile_Table();
  dt->m_parent_file=this->m_parent_file;
  dt->m_parent_driver=this->m_parent_driver;
  dt->m_parent_group=this;
  dt->m_name=name;
  dt->m_flags=flags;
  //
  //m_name_table_map[name]=dt;
  //
  dt->m_tsv_file.open(dt->generateFilename());
  //
  return dt;
}

//AptErr_t Dao_TsvFile_Group::bindTable(const std::string& name,Dao_TsvFile_Table* dt)
//{
//  //m_name_table_map[name]=dt;
//}

//AptErr_t Dao_TsvFile_Group::unbindTable(Dao_TsvFile_Table* dt)
//{
//  return unbind_ptr_map(m_name_table_map,dt);
//}



//////////

Dao_TsvFile_Table::Dao_TsvFile_Table() {
  init();
}

Dao_TsvFile_Table::~Dao_TsvFile_Table() {
  //m_parent_group->unbindTable(this);
  m_tsv_file.close();
}

void Dao_TsvFile_Table::init() {
  m_flags=0;
  m_tsv_state=0;
}

AptErr_t Dao_TsvFile_Table::close() {
  m_tsv_file.close();
  return APT_OK;
}
AptErr_t Dao_TsvFile_Table::flush() {
  m_tsv_file.flush();
  return APT_OK;
}

AptErr_t Dao_TsvFile_Table::clearSchema() {
  /// @todo
  return APT_OK;
}
AptErr_t Dao_TsvFile_Table::clearData() {
  /// @todo
  return APT_OK;
}


void Dao_TsvFile_Table::delBoxRef(const BboardBox* box) {
  m_parent_file->delBoxRef(box,(void*)this);
}
void Dao_TsvFile_Table::addBoxRef(const BboardBox* box) {
  m_parent_file->delBoxRef(box,(void*)this);
}

AptErr_t Dao_TsvFile_Table::addHeader(const std::string& key,const std::string& val)
{
  CHECK_TSV_CLOSED();
  // check that we are writing headers.
  m_tsv_file.addHeader(key,val);
  return APT_OK;
}
AptErr_t Dao_TsvFile_Table::addHeader(const std::string& key,int val)
{
  CHECK_TSV_CLOSED();
  // check that we are writing headers.
  m_tsv_file.addHeader(key,val);
  return APT_OK;
}

AptErr_t Dao_TsvFile_Table::setTsvFilename(const std::string& filename) {
  //m_tsv_filename=filename;
  return APT_OK;
}

AptErr_t Dao_TsvFile_Table::getTsvFilename(std::string* filename) {
  //*filename=m_tsv_filename;
  return APT_OK;
}

std::string Dao_TsvFile_Table::generateFilename() {
  
  FsPath dao_p;
  dao_p.copyFrom(m_parent_driver->m_parent_file->m_dao_path);
  dao_p.setInternalGroupName(m_parent_group->m_name);
  dao_p.setInternalTableName(m_name);
  
  std::string path=dao_p.asUnixPath();

  printf("generateFilename()=='%s'\n",path.c_str());

  return path;
}

AptErr_t Dao_TsvFile_Table::endHeaders() {
  m_tsv_file.writeTsv(generateFilename());
  return APT_OK;
}

AptErr_t Dao_TsvFile_Table::defineColumn(int clvl,int cidx,const std::string& col_name,DaoDataType_t col_type)
{
  defineColumn(clvl,cidx,col_name,col_type,-1);
  return APT_OK;
}

AptErr_t Dao_TsvFile_Table::defineColumn(int clvl,int cidx,const std::string& col_name,DaoDataType_t col_type,int col_size)
{
  CHECK_TSV_CLOSED();
  // check we arent writing data
  m_tsv_file.defineColumn(clvl,cidx,col_name,TsvDtype_From_DaoDtype(col_type));
  return APT_OK;
}

AptErr_t Dao_TsvFile_Table::rewind() {
  m_tsv_file.rewind();
  return APT_OK;
}
//
AptErr_t Dao_TsvFile_Table::nextLevel(int clvl)
{
  return CVT_RV(m_tsv_file.nextLevel(clvl));
}
AptErr_t Dao_TsvFile_Table::writeLevel(int clvl)
{
  return CVT_RV(m_tsv_file.writeLevel(clvl));
}

AptErr_t Dao_TsvFile_Table::nextRow()
{
  return CVT_RV(m_tsv_file.nextLevel(0));
}
AptErr_t Dao_TsvFile_Table::writeRow()
{
  return CVT_RV(m_tsv_file.writeLevel(0));
}

//////////

#define TSVFILE_GETSET(_type) \
  AptErr_t Dao_TsvFile_Table::get(int clvl,const std::string& cidx,_type* val) { \
    ENSURE_TSV_OPEN();                                                   \
    return CVT_RV(m_tsv_file.get(clvl,cidx,*val));                       \
  }                                                                     \
  AptErr_t Dao_TsvFile_Table::get(int clvl,int cidx,_type* val) {      \
    ENSURE_TSV_OPEN();                                                   \
    return CVT_RV(m_tsv_file.get(clvl,cidx,*val));  \
  }                                                                     \
AptErr_t Dao_TsvFile_Table::set(int clvl,const std::string& cidx,const _type& val) { \
  ENSURE_TSV_OPEN();                                                   \
  return CVT_RV(m_tsv_file.set(clvl,cidx,val));       \
} \
AptErr_t Dao_TsvFile_Table::set(int clvl,int cidx,const _type& val) { \
  ENSURE_TSV_OPEN();                                                   \
  return CVT_RV(m_tsv_file.set(clvl,cidx,val));      \
}

//
TSVFILE_GETSET(int);
TSVFILE_GETSET(std::string);
