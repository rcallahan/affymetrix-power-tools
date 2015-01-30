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
// /nfs/ss11/harley/Apt/work/apt-1/affy/sdk/bboard/BboardObj.h ---
//
// $Id$
//

#ifndef _BBOARDOBJ_H_
#define _BBOARDOBJ_H_

//
#include "bboard/Apt2Types.h"
#include "bboard/dao/Dao.h"
//
#include "assert.h"

class BboardObj {
public:
  virtual ~BboardObj();
  // returns "name-of-obj"
  virtual std::string getBboardTypeStr() = 0;
  //
  // The Bboard writer cand write a value as a string
  // or as a FsPath name to a table.
  virtual bool canWriteAsString() = 0;
  virtual std::string valueAsString() = 0;
  virtual AptErr_t writeToTable(Dao_Table& dao_t) = 0;
};


/// an example of a
class BboardObjExample: public BboardObj {
  // constructors are never virtual.
  BboardObjExample() {
    // no-op
  };

  // we always need a virtual destructor.
  virtual ~BboardObjExample() {
    // no-op
  };

  // this is what we are known as.
  // the factory will use this to make us.
  virtual std::string getBboardTypeStr() {
    return "bboard-obj-example";
  };

  //
  virtual bool canWriteAsString() {
    return true;
  };

  //
  virtual std::string valueAsString() {
    return "a bogus value";
  };

  virtual AptErr_t writeToTable(Dao_Table& dao_table);
};

#endif // _BBOARDOBJ_H_
