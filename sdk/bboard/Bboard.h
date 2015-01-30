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
// ~/Affy/apt2/trunk/Bboard.h ---
//
// $Id: Bboard.h,v 1.9 2009-11-05 20:41:39 harley Exp $
//

// Todo: (Meeting 20091102)
//
// * "dryrun" mode to check that the inputs and outputs
//   How to get the "dryrun" code to match up with the real code.
//   A flag to "doRun" and the method returns early?
//
// * resume a stage after a crash. "--resume"
//

// Features:
// * a means for pooling and releasing resources.
// * passing context through nodes which might not care about them.
// * BBoards can be use to implement perls "local" scoping.
// * Myscoping

// Could have a truly global BB with:
// "pushGlobalBboard()"
// "popGlobalBboard()"
// "getGlobalBboard()"

#ifndef _BBOARD_H_
#define _BBOARD_H_

//
#include "bboard/Apt2Types.h"
#include "bboard/BboardTypes.h"
#include "bboard/dao/DaoTypes.h"
//
#include "bboard/BboardBox.h"
#include "bboard/BboardBoxRef.h"
#include "bboard/BboardObj.h"
//
#include "bboard/dao/Dao.h"

// BBT__NEWTYPE: Add headers here for types.
// @todo there should be a common ".h" file for both ".h" files.
#include "chipstream/ChipLayout.h"
#include "bboard/extra/BioSpecies.h"
#include "util/PgOptions.h"


//
#include <map>
#include <string>
#include <vector>

/// There are three kinds of connections for child blackboards.
///
/// * If it is a child, it must be in m_name_ref_map.
///   It will be freed when the parent is freed.
/// * The child may or may not have m_parent set.
///   If set, it will look up vars which are in the parent.
/// * The child may or may not be bound to a name in the parent.
///   If bound, the child may be looked up from the parent.

/// Simple accessor methods for simple types.
#define BB_DEF_GETSET(_type)                                          \
  AptErr_t get(const std::string& name,_type* val,int abortOnErr=1);  \
  AptErr_t set(const std::string& name,_type val,int abortOnErr=1);   \
  AptErr_t get(int idx,_type* val,int abortOnErr=1);                  \
  AptErr_t set(int idx,_type val,int abortOnErr=1);

class Bboard {
public:
  // only used for debugging. (the orginal name.)
  std::string m_name;
  // when searching for names, where to search up the tree.
  Bboard* m_parent;
  // the (name->val) mapping
  typedef std::map<std::string,BboardBoxRef> NameRefMap_t;
  std::map<std::string,BboardBoxRef> m_name_ref_map;
  std::vector<BboardBoxRef> m_ref_array;

  // Unnamed resources.
  // memory with Malloc/Free
  typedef std::vector<void*> ResMemVec_t;
  ResMemVec_t m_resource_mallocs;
  // ValueRefs for things we own, but dont want to have looked up.
  // (normally BBT_BBOARD)
  typedef std::vector<BboardBoxRef> ResRefVec_t;
  ResRefVec_t m_resource_refs;

  //
  Bboard();
  Bboard(const std::string& name);
  Bboard(Bboard* bb_orig);
  ~Bboard();

  //
  void init();
  void clear();
  void clearMap();
  void clearArray();
  void dump();

  //
  int arraySize();
  void arrayResize(int newsize);
  void arrayClearAfterIdx(int idx);

  /// Memory pool functions.
  void* Malloc(int size);
  AptErr_t Free(void* ptr);

  ///
  BboardBox* allocBboardBox();

  ///
  void releaseValRefs();
  void releaseMalloc();
  void releaseAll();

  //
  Bboard* getParentBboard();
  Bboard* setParentBboard(Bboard* parent);
  Bboard* makeChildBboard(const std::string& name);

  //
  Bboard* mkdir(const std::string& name);

  /// Copy all the bindings from one BB to this BB
  AptErr_t copyFrom(Bboard* bb_from);

  /// bind (ie associate) a key/name to a value.
  AptErr_t bind(const std::string& name,BboardBox* vptr);
  AptErr_t bind(const std::string& name,BboardBoxRef& vptr);
  /// forget a binding.
  /// If this is the last reference, it is deleted.
  AptErr_t unbind(const std::string& name);
  /// changs the name of the binding.
  AptErr_t rename(const std::string& f_name,const std::string& t_name);

  /// Methods for getting a BboardBox.
  /// doParents = search the parents for the key.
  /// doCreate  = if not found, create an empty BboardBox and bind it.
  BboardBox* findBoxPtr(const std::string& name,int doCreate,int doParents);
  /// The same as getValuePtr(name, doParents=1, doCreate=0)
  /// returns NULL if not found.
  BboardBox* findBoxPtr(const std::string& name);

  //
  //AptErr_t getDataPtrNoCheck(const std::string& name,void*& dptr,BboardType_t& tcode);
  AptErr_t getPtrAndType(const std::string& name,void*& dptr,BboardType_t& tcode);
  AptErr_t setPtrAndType(const std::string& name,void* dptr,BboardType_t tcode);
  //
  AptErr_t getPtrIfType(const std::string& name,void*& dptr,BboardType_t tcode);

  /// Simple accessor methods for simple types.
  // API for each of the getter setters below looks like:
  // AptErr_t set(const string &s, TYPE t);
  // AptErr_t get(const string &s, TYPE *t);

  BB_DEF_GETSET(int);
  BB_DEF_GETSET(double);
  BB_DEF_GETSET(std::string);
  //BB_DEF_GETSET(char*);

  // vectors of values, user allocated, the bboard takes ownership when registered.
  // so dont delete them yourself.
  BB_DEF_GETSET(std::vector<char>*);
  BB_DEF_GETSET(std::vector<int>*);
  BB_DEF_GETSET(std::vector<float>*);
  BB_DEF_GETSET(std::vector<double>*);
  //
  BB_DEF_GETSET(BioSpecies*);
  BB_DEF_GETSET(ChipLayout*);
  BB_DEF_GETSET(PgOptions*);
  //
  BB_DEF_GETSET(Bboard*);
  BB_DEF_GETSET(BboardObj*);
  //
  BB_DEF_GETSET(FsPath*);
  // Drivers are not allowed to be put on the board.
  BB_DEF_GETSET(Dao_File*);
  BB_DEF_GETSET(Dao_Group*);
  BB_DEF_GETSET(Dao_Table*);
  //
  // Grep for this when adding a new type
  // BBT__NEWTYPE


  // These IO functions are just parked here until the APT1 merge.
  AptErr_t writePgOptionsToDaoTable(Dao_Table* dao_t,PgOptions* pgopts);

};

#endif // _BBOARD_H_
