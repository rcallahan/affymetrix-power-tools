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
// ~/Affy-projects/apt2/trunk/Apt2Types.h ---
//
// $Id: Apt2Types.h,v 1.2 2009-11-01 01:02:07 harley Exp $
//

#ifndef _APT2TYPES_H_
#define _APT2TYPES_H_

//
#include "util/AptErrno.h"
#include "util/Convert.h"

/// Barf we see an error.
/// This assumes that "abortOnErr" is in scope;
#define APT_ABORT_ON_ERR(_rv,_msg)                            \
  {                                                           \
    if ((abortOnErr==1)&&(_rv!=APT_OK)) {                     \
      std::string file_line=__FILE__":"+ToStr(__LINE__)+":";  \
      Err::errAbort(file_line+_msg);                          \
    }                                                         \
  }

#endif // _APTTYPES_H_
