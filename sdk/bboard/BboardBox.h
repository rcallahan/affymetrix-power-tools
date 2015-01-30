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
// ~/Affy/apt2/trunk/BboardBox.h ---
// 
// $Id: BboardBox.h,v 1.7 2009-11-05 20:41:39 harley Exp $
// 

#ifndef _BBOARDBOX_H_
#define _BBOARDBOX_H_

//
#include "bboard/Apt2Types.h"
#include "bboard/BboardTypes.h"
#include "bboard/BboardObj.h"
//
#include "bboard/dao/Dao.h"
#include "bboard/dao/DaoTypes.h"
//
// BBT__NEWTYPE
// add the ".h" file here for the class.
#include "bboard/extra/BioSpecies.h"
#include "chipstream/ChipLayout.h"
#include "util/PgOptions.h"

//
#include <map>
#include <string>
#include <vector>

// @todo keep the m_data_ptr values in a static class member
// so we can detect duplicates for debugging.

#define BBBOX_DEF_GETSET(_type) \
  AptErr_t getData(_type* val); \
  AptErr_t setData(_type  val);

class BboardBox {
public:
  // this is advisory only -- the name binding is given in the Bboard map.
  // (only used by the dump method.)
  std::string m_name;
private:
  /// The reference count.
  int m_refcnt;
  /// The type of the data stored n m_data_ptr
  BboardType_t m_data_tcode;
  /// the data we own.
  void* m_data_ptr;
  // do we own this memory?
  // unused, cause at the moment a BboardBox always owns its memory.
  bool m_data_owned;

  // to make sure there is only one Box per pointer.
  typedef std::map<void*,BboardBox*> boxForPtr_t;
  static boxForPtr_t m_boxForPtr_map;
  static void addBoxForPtr(void* ptr,BboardBox* box);
  static void delBoxForPtr(void* ptr,BboardBox* box);

public:
  static BboardBox* findBoxForPtr(void* ptr);

public:
  //
  BboardBox();
  BboardBox(const std::string& name);
  BboardBox(const std::string& name,void* dataptr,BboardType_t tcode);
  //
  virtual ~BboardBox();

  //
  void init();
  void clearData();
  void dump();
  std::string asString();

  /// how many refs referr to us?
  int  getRefcnt() const { return m_refcnt; };
  //
  int refcntStep(int cnt);
  int refcntInc() { return refcntStep(+1); };
  int refcntDec() { return refcntStep(-1); };

  //
  bool getDataOwned(bool val) { /* (val);  unused */ return m_data_owned; };
  bool setDataOwned(bool val) { m_data_owned=val; return val; };

  // 
  void* getPtr() { return m_data_ptr; } ;
  BboardType_t getBbType() { return m_data_tcode; };
  //
  AptErr_t setPtrAndType(void* dataptr,BboardType_t tcode);
  AptErr_t getPtrAndType(void*& dataptr,BboardType_t& tcode);
  //
  AptErr_t getPtrIfType(void*& dataptr,BboardType_t tcode);

  //
  BBBOX_DEF_GETSET(int);
  BBBOX_DEF_GETSET(double);
  // string is a special case.
  AptErr_t getData(std::string* val);
  AptErr_t setData(const std::string& val);
  //
  //BBBOX_DEF_GETSET(char*);
  //
  BBBOX_DEF_GETSET(std::vector<char>*);
  BBBOX_DEF_GETSET(std::vector<int>*);
  BBBOX_DEF_GETSET(std::vector<float>*);
  BBBOX_DEF_GETSET(std::vector<double>*);
  //
  BBBOX_DEF_GETSET(BioSpecies*);
  BBBOX_DEF_GETSET(ChipLayout*);
  BBBOX_DEF_GETSET(PgOptions*);
  //
  BBBOX_DEF_GETSET(Bboard*);
  BBBOX_DEF_GETSET(BboardObj*);
  //
  BBBOX_DEF_GETSET(FsPath*);
  BBBOX_DEF_GETSET(Dao_File*);
  BBBOX_DEF_GETSET(Dao_Group*);
  BBBOX_DEF_GETSET(Dao_Table*);
  //
  // Grep for this when adding a new type.
  // BBT__NEWTYPE
};

#endif // _BBOARDBOX_H_
