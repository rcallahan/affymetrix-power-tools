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
// ~/Affy-projects/apt2/trunk/DaoTypes.h ---
//
// $Id: DaoTypes.h,v 1.4 2009-11-04 20:37:32 harley Exp $
//

#ifndef _DAOTYPES_H_
#define _DAOTYPES_H_

//
class FsPath;
//
class Dao_File;
class Dao_Driver;
class Dao_Group;
class Dao_Table;
//
class Dao_TsvFile_Driver;
class Dao_TsvFile_Group;
class Dao_TsvFile_Table;

///
enum DaoFlags_t {
  DAO_NONE        = 0x000,
  //
  DAO_OPEN        = 0x001,
  DAO_CREATE      = 0x002,
  DAO_REPLACE     = 0x004,
  //
  DAO_RO          = 0x010,
  DAO_RW          = 0x020,
};

///
enum DaoDataType_t {
  //
  DAO_CHAR=1,
  DAO_SHORT,
  DAO_INT,
  DAO_FLOAT,
  DAO_DOUBLE,
  DAO_STRING,
};

#endif // _DAOTYPES_H_
