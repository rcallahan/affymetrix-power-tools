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

//
#include "broadutil/BroadUtil.h"
//
#include <cerrno>
#include <cstring>
//

void safe_strncpy(char *to, const char *from, size_t len)
{
  strncpy(to, from, len);
  to[len] = '\0';
}

bool startswith(const char *target, const char *test)
{
  return (strncmp(target, test, strlen(test)) == 0);
}

bool endswith(const char *str, const char * suffix)
{
  size_t len = strlen(str);
  size_t suffixLen = strlen(suffix);
  if (len < suffixLen) {
    return false;
  }
  return strcmp(str + len - suffixLen, suffix) == 0;
}

FILE *fopen_check(const char *fname, const char *mode)
{
  FILE *strm = fopen(fname, mode);
  if (strm == NULL) {
    throw BroadException("Could not open file", __FILE__, __LINE__, fname, errno);
  }
  return strm;
}

void fclose_check(FILE *fp)
{
  ;
  if (fclose(fp) != 0) {
    throw BroadException("Could not close file", __FILE__, __LINE__, "", errno);
  }
}

void fwrite_check(const void *data, size_t size, size_t nitems, FILE *fp)
{
  size_t ret = fwrite(data, size, nitems, fp);
  if (ret != nitems) {
    throw BroadException("Problem writing file", __FILE__, __LINE__, "", errno);
  }
}

void fread_check(void *data, size_t size, size_t nitems, FILE *fp)
{
  size_t ret = fread(data, size, nitems, fp);
  if (ret != nitems) {
    throw BroadException("Problem writing file", __FILE__, __LINE__, "", errno);
  }
}

void fseek_check(FILE *stream, long offset, int whence)
{
  int ret = fseek(stream, offset, whence);
  if (ret != 0) {
    throw BroadException("Error fseeking", __FILE__, __LINE__, "", errno);
  }
}

long ftell_check(FILE *stream)
{
  long ret = ftell(stream);
  if (ret == -1) {
    throw BroadException("Error ftelling", __FILE__, __LINE__, "", errno);
  }
  return ret;
}

/******************************************************************/
/**************************[END OF BroadUtil.cpp]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
