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
// ~/Affy-projects/apt2/trunk/Dao_TsvFile.h ---
//
// $Id: Dao_TsvFile.h,v 1.5 2009-11-05 20:41:39 harley Exp $
//

#ifndef _DAO_DRIVER_TSVFILE_H_
#define _DAO_DRIVER_TSVFILE_H_

#include "bboard/Apt2Types.h"
#include "bboard/dao/Dao.h"
#include "file/TsvFile/TsvFile.h"
//
#include <map>

//////////

class Dao_TsvFile_Driver : public Dao_Driver {
public:
  Dao_File* m_parent_file;
  //
  //typedef std::map<std::string,Dao_TsvFile_Group*>::iterator name_group_map_iter_t;
  //std::map<std::string,Dao_TsvFile_Group*> m_name_group_map;
  //
  Dao_TsvFile_Driver();
  virtual ~Dao_TsvFile_Driver();
  void init();
  //
  AptErr_t attach(Dao_File* df);
  AptErr_t create();
  AptErr_t open();
  AptErr_t close();
  AptErr_t flush();
  //
  Dao_Group* createGroup(const std::string& name,int flags);
  Dao_Group* openGroup(const std::string& name,int flags);
  //
  std::string getName() { return "Tsvfile_Driver"; };
  //
  //AptErr_t bindGroup(const std::string& name,Dao_TsvFile_Group* dt);
  //AptErr_t unbindGroup(Dao_TsvFile_Group* dt);
};

//////////

class Dao_TsvFile_Group : public Dao_Group {
public:
  Dao_File*           m_parent_file;
  Dao_TsvFile_Driver* m_parent_driver;
  std::string m_name;
  int m_flags;
  //
  //typedef std::map<std::string,Dao_TsvFile_Table*>::iterator name_table_map_iter_t;
  //std::map<std::string,Dao_TsvFile_Table*> m_name_table_map;
  //
  Dao_TsvFile_Group();
  virtual ~Dao_TsvFile_Group();
  void init();
  AptErr_t close();
  AptErr_t flush();
  //
  void delBoxRef(const BboardBox* box);
  void addBoxRef(const BboardBox* box);
  //
  Dao_Table* createTable(const std::string& name,int flags);
  Dao_Table* openTable(const std::string& name,int flags);
  //
  //AptErr_t bindTable(const std::string& name,Dao_TsvFile_Table* dt);
  //AptErr_t unbindTable(Dao_TsvFile_Table* dt);
  Dao_File* getFile() { return m_parent_file; };
  std::string getName() { return m_name; };
};

//////////

#define TSVFILE_DEF_GETSET(_type) \
  AptErr_t get(int clvl,const std::string& cidx,_type* val);       \
  AptErr_t get(int clvl,int cidx,_type* val);                      \
  AptErr_t set(int clvl,const std::string& cidx,const _type& val); \
  AptErr_t set(int clvl,int cidx,const _type& val);

class Dao_TsvFile_Table : public Dao_Table {
public:
  Dao_File*           m_parent_file;
  Dao_TsvFile_Driver* m_parent_driver;
  Dao_TsvFile_Group*  m_parent_group;
  //
  std::string m_name;
  int m_flags;
  // a filename to override with.
  std::string m_tsv_filename;
  //
  affx::TsvFile m_tsv_file;
  int m_tsv_state;
  //
  Dao_TsvFile_Table();
  virtual ~Dao_TsvFile_Table();
  void init();
  AptErr_t close();
  AptErr_t flush();
  //
  AptErr_t clearSchema();
  AptErr_t clearData();

  void delBoxRef(const BboardBox* box);
  void addBoxRef(const BboardBox* box);

  //
  Dao_File* getFile() { return m_parent_file; };
  Dao_Group* getGroup() { return m_parent_group; };
  std::string getName() { return m_name; };

  //
  std::string generateFilename();

  // sets the tsv filename for tsvfiles.
  AptErr_t setTsvFilename(const std::string& filename);
  AptErr_t getTsvFilename(std::string* filename);

  //
  AptErr_t addHeader(const std::string& key,const std::string& val);
  AptErr_t addHeader(const std::string& key,int val);
  //
  AptErr_t defineColumn(int clvl,int cidx,const std::string& col_name,DaoDataType_t col_type);
  AptErr_t defineColumn(int clvl,int cidx,const std::string& col_name,DaoDataType_t col_type,int col_size);
  
  //
  AptErr_t endHeaders();

  AptErr_t rewind();

  //
  AptErr_t nextLevel(int clvl=0);
  AptErr_t writeLevel(int clvl=0);
  AptErr_t nextRow();
  AptErr_t writeRow();

  //
  TSVFILE_DEF_GETSET(int);
  TSVFILE_DEF_GETSET(std::string);
};

#endif // _DAO_DRIVER_TSVFILE_H_
