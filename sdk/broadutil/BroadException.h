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
 * FILE BroadException.h
 */

#ifndef _BROADEXCEPTION_H
#define _BROADEXCEPTION_H

#include <climits>
#include <cstdio>
#include <exception>
//

#ifdef _MSC_VER
#include <windows.h>
#define PATH_MAX MAX_PATH
#endif

/*
 * general-purpose exception class.  Can be subclassed, but need not be.
 */
class BroadException: public std::exception 
{
    enum {kMessageMax = 1023,
          kFormattedMessageMax = PATH_MAX * 3};
    
    
  public:
    char m_msg[kMessageMax + 1];
    char m_sourcefile[PATH_MAX + 1];
    const size_t m_sourceline;
    char m_filename[PATH_MAX + 1];
    const int m_errno;
    char m_formattedMsg[kFormattedMessageMax + 1];

    BroadException(const char *msg, const char *sourcefile, size_t sourceline);
    BroadException(const char *msg, const char *sourcefile, size_t sourceline, const char *filename);
    BroadException(const char *msg, const char *sourcefile, size_t sourceline, const char *filename, int theErrno);
    void report(FILE *fp);

  private:
    void formatException();
};

#endif /* _BROADEXCEPTION_H */

/******************************************************************/
/**************************[END OF BroadException.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
