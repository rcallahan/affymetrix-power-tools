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
 * FILE CelFileStreamReader.h
 */

#ifndef _CELFILESTREAMREADER_H
#define _CELFILESTREAMREADER_H

#include "broadutil/APTUtil.h"
//
#include "file/CELFileData.h"
//
#include <cstring>
#include <fstream>
#include <string>
//

// Class for reading a cel file without loading the whole thing into RAM at once.
class CelFileStreamReader
{
  public:
    CelFileStreamReader(const std::string &filename);

    // The copied or assigned-to object starts its iteration at the beginning.
    CelFileStreamReader(const CelFileStreamReader& other);
    CelFileStreamReader& operator=(const CelFileStreamReader &other);
    
    ~CelFileStreamReader();

    
  CONSTHACK affxcel::CCELFileHeaderData &getHeader() CONSTHACK;
    

    // Each of the methods below advances to the next element in the cel file
    // and returns the requested value.
    // If you want to get more than one attribute from a single element, you
    // must use getNextEntry, and then pick the attributes you want from the
    // struct that gets filled out.
    float getNextIntensity();
    float getNextStdv();
    short getNextPixels();
    void getNextEntry(affxcel::CELFileEntryType *val);

    // Returns the index of the next element to be read from the cel file.
  size_t getNextIndex();

    static size_t determineHeaderSize(const std::string &filename);

  private:
    void startIterating();
  affxcel::CELFileEntryType getNextEntryUnswapped();
    
  private:
    std::string m_filename;
    affxcel::CCELFileHeaderData m_header;
    size_t m_nextIndex;
    std::ifstream m_instream;
};



#endif /* _CELFILESTREAMREADER_H */

/******************************************************************/
/**************************[END OF CelFileStreamReader.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
