////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

/**
 * @file   AptTypes.h
 * @author harley
 * @date   Fri Mar  9 17:28:04 PST 2007
 * 
 * @brief  These types are common in the APT code.
 *         Collect them here for use by the rest of the code.
 *         No classes or functions here, just types.
 */

#ifndef _APTTYPES_H_
#define _APTTYPES_H_

//
#include "portability/affy-base-types.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

// These are ids of the probes/chips and such.
// As "ids" we know they are unsigned.
// @todo make them signed as we do have "-1" as an invalid value.
typedef uint32_t atomid_t;
typedef uint32_t chipid_t;
typedef uint32_t probeid_t;
typedef uint32_t probesetid_t;
// typedef int32_t atomid_t;
// typedef int32_t chipid_t;
// typedef int32_t probeid_t;
// typedef int32_t probesetid_t;

// An index is an an offset into a data structure.
// It *isnt* the same as the id and sometimes we use
// -1 to mark an invalid value.(Hence the signedness
typedef int32_t atomidx_t;
typedef int32_t chipidx_t;
typedef int32_t probeidx_t;
typedef int32_t probesetidx_t;

// Commonly use datastructures derived from above.
typedef std::map<probeid_t, std::vector<std::string> > probeidmap_t;

/// For row and column indexes.
typedef int colrow_t;

#endif /* __APTTYPES_H_ */
