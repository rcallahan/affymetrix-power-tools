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
// affy/sdk/file5/File5_Group.h ---
//
// $Id: File5_Group.h,v 1.18 2009-09-18 03:37:27 mspald Exp $
//

#ifndef _FILE5_GROUP_H_
#define _FILE5_GROUP_H_

//
#include "file5/File5_Object.h"
#include "file5/File5_types.h"
//


//
class affx::File5_Group : public affx::File5_Object {
public:
  //
  hid_t m_h5_user_meta_id;
  int   m_h5_user_meta_idx;

  //
  File5_Group();
  //File5_Group(affx::File5_File* file5);
  virtual ~File5_Group();

  //
  affx::File5_return_t open(const std::string& filename,int flags);

  //
  affx::File5_return_t init();
  affx::File5_return_t dump();
  affx::File5_return_t dump1();
  //int flush();
  affx::File5_return_t close();
  virtual char file5_kind_char();

  // @todo: move all this to File5_UserMeta
  hid_t usermeta_getid();
  affx::File5_return_t usermeta_close();
  int usermeta_size();
  affx::File5_return_t usermeta_resize(int newsize);
  affx::File5_return_t usermeta_writeidx(int idx,const std::string& key,const std::string& val);
  //
  affx::File5_return_t addHeader(const std::string& key,const std::string& val);
  // pretty name.
  int getHeaderCount() { return usermeta_size(); };
  affx::File5_return_t getHeaderByIdx(int idx,std::string* key,std::string* val);
  //
  void headersBegin();
  affx::File5_return_t headersNext(std::string* key,std::string* val);
  /// find the FIRST value for this header.
  affx::File5_return_t getHeader(const std::string& key,std::string* val);
  /// find ALL the values for this header.
  affx::File5_return_t getHeader(const std::string& key,std::vector<std::string>* val);

  //
  int get_num_objs();
  std::string get_name_by_idx(int idx);
  int get_objtype_by_idx(int idx);
  affx::File5_Object* get_object_by_idx(int idx);

  // the list of names in this group.
  std::vector<std::string> listNames();

  /// @brief     Is there an object (dataset,group,link,etc...) with this name?
  /// @param     name      name to check
  /// @return    true if there is an object with this name.
  bool name_exists(const std::string& name);
  /// @brief     Remove this name.  If no more names exist, the object is deleted.
  /// @param     name      name to remove.
  /// @return    ok if removed.
  affx::File5_return_t unlink(const std::string& name);

  //
  affx::File5_return_t deleteItem(const std::string& name);
  affx::File5_return_t deleteItemsIn(const std::string& name);
  affx::File5_return_t deleteTsv(const std::string& name);

  //
  affx::File5_Object* openDataset(const std::string& name,int flags);

  //
  affx::File5_Group* openGroup(const std::string& name,int flags);

  //
  affx::File5_Matrix* openMatrix(const std::string& name,
                                 affx::File5_dtype_t dtype,
                                 const std::vector<int>& dims,
                                 int flags);
  //
  affx::File5_Tsv* openTsv(const std::string& name,int flags);
  affx::File5_Tsv* openTsv(const std::string& name);

  // generic open vector...
  affx::File5_Vector* createVector();
  affx::File5_Vector* openVector(const std::string& name,
                                 const affx::File5_dtype_t dtype,
                                 const int flags);
  affx::File5_Vector* openVector(const std::string& name);

  // specialized opens.
  affx::File5_Vector* openVector_String(const std::string& name,
                                        int size,
                                        const int flags);
};

#endif // _FILE5_GROUP_H_
