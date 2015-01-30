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

/*
 * FILE BroadUtil.h
 */

#ifndef _BROADUTIL_H_
#define _BROADUTIL_H_

//
#include "broadutil/BroadException.h"
//
#include <cerrno>
#include <cstdio>
#include <cstring>
//

#define ARRAY_LEN(ary) (sizeof(ary)/sizeof(ary[0]))

void safe_strncpy(char *to, const char *from, size_t len);

bool startswith(const char *target, const char *test);

bool endswith(const char *str, const char * suffix);

FILE *fopen_check(const char *fname, const char *mode);

void fclose_check(FILE *fp);

void fwrite_check(const void *data, size_t size, size_t nitems, FILE *fp);

void fread_check(void *data, size_t size, size_t nitems, FILE *fp);

void fseek_check(FILE *stream, long offset, int whence);

long ftell_check(FILE *stream);

#endif /* _BROADUTIL_H_ */

/******************************************************************/
/**************************[END OF BroadUtil.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
