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

#include "broadutil/BroadException.h"
//
#include "broadutil/BroadUtil.h"
//
#include <sstream>
//


BroadException::BroadException(const char *msg, const char *sourcefile, size_t sourceline):
    m_msg(),
    m_sourcefile(),
    m_sourceline(sourceline),
    m_filename(),
    m_errno(0)
{
    m_filename[0] = '\0';
    safe_strncpy(m_msg, msg, sizeof(m_msg));
    safe_strncpy(m_sourcefile, sourcefile, sizeof(m_sourcefile));
    formatException();
}

BroadException::BroadException(const char *msg, const char *sourcefile, size_t sourceline, const char *filename):
    m_msg(),
    m_sourcefile(),
    m_sourceline(sourceline),
    m_filename(),
    m_errno(0)
{
    safe_strncpy(m_msg, msg, sizeof(m_msg));
    safe_strncpy(m_sourcefile, sourcefile, sizeof(m_sourcefile));
    safe_strncpy(m_filename, filename, sizeof(m_filename));
    formatException();
}

BroadException::BroadException(const char *msg, const char *sourcefile, size_t sourceline, const char *filename, int theErrno):
    m_msg(),
    m_sourcefile(),
    m_sourceline(sourceline),
    m_filename(),
    m_errno(theErrno)
{
    safe_strncpy(m_msg, msg, sizeof(m_msg));
    safe_strncpy(m_sourcefile, sourcefile, sizeof(m_sourcefile));
    safe_strncpy(m_filename, filename, sizeof(m_filename));
    formatException();
}

void BroadException::report(FILE *fp) 
{
    fprintf(fp, "%s\n", m_formattedMsg);
}

void BroadException::formatException() 
{
    const char *errnoStr = "";
    if (m_errno != 0) {
        errnoStr = strerror(m_errno);
    }
    std::ostringstream strm;
    strm << "ERROR: " << m_sourcefile << ":" << m_sourceline << " " << m_msg << ", file " << m_filename << "; " << errnoStr;
    safe_strncpy(m_formattedMsg,  strm.str().c_str(), sizeof(m_formattedMsg));
}
    

/******************************************************************/
/**************************[END OF BroadException.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
