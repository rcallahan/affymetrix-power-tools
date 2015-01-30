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
// ~/APT2/trunk/Dao.h ---
//
// $Id: Dao.h,v 1.6 2009-11-05 20:41:39 harley Exp $
//
#ifndef _DAO_H_
#define _DAO_H_

/// Notes:

// everything is 2d
// Calvin support mulitple 2d stores.
// TsvFile is only 1 2d store, but we can fake it with mutiple files.
// HDF5 is anything. (Nd, many groups.)

// The Dao is a cache in front of these three formats.
// * Does caching
// * Buffer management (wether to buffer or stream)
// * Reordering of data to "processing order".

// * generate a stable reference.

// * when creating a Dao, the column lifetimes should be similar.
//   if there is a short lived column, consider creating a separate DAO
//   and discarding that DAO sooner.
//   If a shortlived column is made in a long lived Table it will be around for a while.

// The Dao keeps a bunch of state about the current format
// of the data: the order, columns, containership, transient

// The Dao may keep several representations of the data about:
// orginal, reordered, selected, reordered and selected.

// The Dao uses drivers to handle the actual IO.
// Virtual functions to dispatch the IO calls to the backend.
// As virutal functions are expensive, we will want to do bulk IO
// when possible.

// Perhaps we will want to have the caching in the Dao itself
// rather than the drivers.  (or the buffers shared between them?)

// Possible backends:
//   Dao_Driver - base virtual class.
//   Dao_Driver_TsvFile
//   Dao_Driver_File5
//   Dao_Driver_Calvin
//   Dao_Driver_Memory

// Change to three levels in the Dao:
// * File
// * Group
// * Table

//
#include "bboard/dao/DaoTypes.h"
#include "util/FsPath.h"
#include "bboard/dao/DaoUtil.h"
//
#include "bboard/Apt2Types.h"
#include "bboard/BboardTypes.h"
#include "bboard/Bboard.h"
//
#include <string>
#include <vector>

// No virtual functions as we do our own dispatch.
// (Cause we might need to reorder the data into a new datafile.)

//////////

/// @todo should do the File5 inheritance thing

class Dao_File {
public:
  Dao_Driver* m_driver;
  //
  FsPath m_dao_path;
  //
  int m_filefmt;
  int m_flags;
  //
  std::vector<Dao_Group*> m_groups;

  // The number of boxes which refer to some part of this DAO structure.
  // we dont actually go about deleting our parts and deallocation ourselves
  // until this is zero.
  int m_box_refcnt;

  //
  Dao_File();
  virtual ~Dao_File();
  void init();
  //
  std::string getName() { return m_dao_path.getFileName(); }
  int getFormat();
  AptErr_t setFormat(int fmt);

  //
  AptErr_t create(const std::string& name,int flags);
  AptErr_t create(const FsPath& dp,int flags);
  AptErr_t open(const std::string& name,int flags);
  AptErr_t open(const FsPath& dp,int flags);
  //
  AptErr_t flush();
  AptErr_t close();
  //
  Dao_Group* createGroup(const std::string& name,int flags);
  Dao_Group* openGroup(const std::string& name,int flags);
  //
  // called when a ref goes away.
  void addBoxRef(const BboardBox* box);
  void delBoxRef(const BboardBox* box);
  void addBoxRef(const BboardBox* box,void* dao_obj_ptr);
  void delBoxRef(const BboardBox* box,void* dao_obj_ptr);

};

//////////

class Dao_Driver {
public:
  Dao_Driver();
  virtual ~Dao_Driver();
  virtual AptErr_t attach(Dao_File* df) = 0;
  virtual AptErr_t create() = 0;
  virtual AptErr_t open() = 0;
  virtual AptErr_t close() = 0;
  virtual AptErr_t flush() = 0;
  virtual void init()=0;
  //
  virtual std::string getName() = 0;
  //
  virtual Dao_Group* createGroup(const std::string& name,int flags) = 0;
  virtual Dao_Group* openGroup(const std::string& name,int flags) = 0;
};

/// Thse are the abstract classes which are returned.

class Dao_Group {
public:
  virtual ~Dao_Group() { };
  virtual AptErr_t close() = 0;
  virtual AptErr_t flush() = 0;
  virtual std::string getName() = 0;
  virtual Dao_File* getFile() = 0;
  //
  virtual void delBoxRef(const BboardBox* box) = 0;
  virtual void addBoxRef(const BboardBox* box) = 0;
  //
  virtual Dao_Table* openTable(const std::string& name,int flags) = 0;
  virtual Dao_Table* createTable(const std::string& name,int flags) = 0;
  //
};

//////////

#define DAO_TABLE_DEF_GETSET(_type) \
  virtual AptErr_t get(int clvl,const std::string& cidx,_type* val) = 0;  \
  virtual AptErr_t get(int clvl,int cidx,_type* val) = 0;                 \
  virtual AptErr_t set(int clvl,const std::string& cidx,const _type& val)  = 0; \
  virtual AptErr_t set(int clvl,int cidx,const _type& val) = 0;

class Dao_Table {
public:
  virtual ~Dao_Table() { } ;
  virtual AptErr_t close() = 0;
  virtual AptErr_t flush() = 0;
  virtual std::string getName() = 0;
  virtual Dao_File* getFile() = 0;
  virtual Dao_Group* getGroup() = 0;

  //
  virtual AptErr_t clearSchema() = 0;
  virtual AptErr_t clearData() = 0;

  //
  virtual void addBoxRef(const BboardBox* box) = 0;
  virtual void delBoxRef(const BboardBox* box) = 0;

  // sets the tsv filename for tsvfiles.
  virtual AptErr_t setTsvFilename(const std::string& filename) = 0;
  virtual AptErr_t getTsvFilename(std::string* filename) = 0;

  //
  virtual AptErr_t addHeader(const std::string& key,const std::string& val) = 0;
  virtual AptErr_t addHeader(const std::string& key,int val) = 0;
  //
  virtual AptErr_t defineColumn(int clvl,int cidx,const std::string& col_name,DaoDataType_t col_type) = 0;
  virtual AptErr_t defineColumn(int clvl,int cidx,const std::string& col_name,DaoDataType_t col_type,int col_size) = 0;

  //
  virtual AptErr_t endHeaders() = 0;

  //
  virtual AptErr_t rewind() = 0;

  //
  virtual AptErr_t nextLevel(int clvl=0) = 0;
  virtual AptErr_t writeLevel(int clvl=0) = 0;
  //
  virtual AptErr_t nextRow() = 0;
  virtual AptErr_t writeRow() = 0;

  //
  DAO_TABLE_DEF_GETSET(int);
  DAO_TABLE_DEF_GETSET(std::string);
};

#endif // _DAO_H_
