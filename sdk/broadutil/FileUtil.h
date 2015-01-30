////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

/// @file FileUtil.h  Header for FileUtil.cpp

#ifndef _FILEUTIL_H_
#define _FILEUTIL_H_

#include <cstdio>
#include <cstring>
#include <string>
//

namespace birdseed 
{
    /*
     * Read the input file until a line is found starting with desiredHeaderStartsWith.
     * Return that line.
     */
    std::string readHeader(FILE *dataStrm, const char *desiredHeaderStartsWith);

    std::string readWord(FILE *fp);
    bool readDouble(FILE *fp, double *val);
    bool readUnsignedInt(FILE *fp, unsigned int *val);
};

#endif /* _FILEUTIL_H_ */

/******************************************************************/
/**************************[END OF FileUtil.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
