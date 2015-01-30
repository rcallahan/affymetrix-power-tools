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

/// @file FileUtil.cpp  Birdseed functions to read and write strings and doubles.

#include "broadutil/FileUtil.h"
//
#include "broadutil/BroadException.h"
#include "broadutil/BroadUtil.h"
//
#include <cerrno>
#include <cstring>
#include <string.h>
//

using namespace std;

// Only one of these, so be a little careful.
static char buf[5000000];

/*
 * Read the input file until a line is found starting with desiredHeaderStartsWith.
 * Return that line.
 */
string birdseed::readHeader(FILE *dataStrm, const char *desiredHeaderStartsWith)
{
    while (fgets(buf, sizeof(buf), dataStrm) != NULL) {
        if (strlen(buf) == sizeof(buf) - 1) {
            throw BroadException("Line too long reading data file header", __FILE__, __LINE__);
        }
        if (startswith(buf, desiredHeaderStartsWith)) {
            return buf;
        }
    }
    throw BroadException("Problem encountered looking for header line in data file", __FILE__, __LINE__, "", errno);
}

string birdseed::readWord(FILE *fp)
{
    int ret = fscanf(fp, "%s", buf);
    if (ret == EOF) {
        if (ferror(fp)) {
            throw BroadException("Error scanning data file", __FILE__, __LINE__, "", errno);
        } else {
            return "";
        }
    }
    return buf;
}

bool birdseed::readDouble(FILE *fp, double *val)
{
    int ret = fscanf(fp, "%lf", val);
    if (ret == 1) {
        return true;
    } else if (ret == 0) {
        return false;
    } else if (ret == EOF) {
        if (ferror(fp)) {
            throw BroadException("Error scanning data file", __FILE__, __LINE__, "", errno);
        } else {
            return false;
        }
    } else {
            throw BroadException("Should never happen", __FILE__, __LINE__);
    }
}

bool birdseed::readUnsignedInt(FILE *fp, unsigned int *val)
{
    int ret = fscanf(fp, "%u", val);
    if (ret == 1) {
        return true;
    } else if (ret == 0) {
        return false;
    } else if (ret == EOF) {
        if (ferror(fp)) {
            throw BroadException("Error scanning data file", __FILE__, __LINE__, "", errno);
        } else {
            return false;
        }
    } else {
            throw BroadException("Should never happen", __FILE__, __LINE__);
    }
}



/******************************************************************/
/**************************[END OF FileUtil.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
