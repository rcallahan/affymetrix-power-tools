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

/// @file   MidasCreateDirectory.h
/// @brief  Headers for midas directory creation function.

#ifndef MIDASCREATEDIRECTORY_H
#define MIDASCREATEDIRECTORY_H

/** midasCreateDirectory
 * @brief Create an output directory.
 * If an error occurs, return pointer to error message, else 0
 * @param outDirectory requested directory name
 * @return error message if any
 */
#include <cstring>
#include <string>
//
std::string midasCreateDirectory (const std::string& outDirectory);

#endif /* MIDASCREATEDIRECTORY_H */
