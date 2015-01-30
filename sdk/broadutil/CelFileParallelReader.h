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
 * FILE CelFileParallelReader.h
 */

#ifndef _CELFILEPARALLELREADER_H
#define _CELFILEPARALLELREADER_H

#include "broadutil/BroadException.h"
#include "broadutil/CelFileStreamReader.h"
//
#include <cassert>
#include <cstring>
#include <string>
#include <vector>
//

// Used only to select what an instantiation of CelFileParallelReader will read.
enum ParallelReaderDesiredAttribute {
    ParallelReader_Intensity,
    ParallelReader_Stdv,
    ParallelReader_Pixels,
    ParallelReader_Entry
};
    

// Low-memory use tool for iterating through a set of cel files in parallel.
// Instead of reading all files into memory, it uses CelFileStreamReader to
// read through the files as needed.
// This class currently has two limitations:
// * Only supports reading binary cel files;
// * Can process only about 1000 files at a time (OS file handle limitation).
// Both of these limitations could be removed if necessary.
// This class is a template (possibly foolishly) for two reasons:
// 1) in the hope that the switch statement will be optimized out by the compiler
//    because the value being switched on is constant;
// 2) because eventually this class will be changed so that it reads ahead in each
//    file (if there are more input files than available file handles), and it will
//    read ahead only the desired attribute, instead of the entire entry, in order
//    to maximize memory use.
template<class C, ParallelReaderDesiredAttribute ATTRIBUTE>
class CelFileParallelReader
{
  public:
    CelFileParallelReader(const std::vector<std::string> &filenames)
    {
        if (filenames.empty()) {
            throw BroadException("It doesn't make sense to create a CelFileParallelReader with an empty list of files",
                                 __FILE__, __LINE__);
        }
        m_readers.reserve(filenames.size());
        for (std::vector<std::string>::const_iterator it = filenames.begin();
             it != filenames.end(); ++it) {
            m_readers.push_back(CelFileStreamReader(*it));
        }
    }

    // Get the next value (which value to get is determined by the
    // way the template is instantiated), from each file, and return
    // the values in the arg values, in the same order as the filenames
    // passed to the ctor.
    // Any previous contents of values vector are discarded.
    // Advances each file so that a subsequence call will return the next
    // value.
    void getNext(std::vector<C> *values)
    {
        values->clear();
        values->resize(m_readers.size());
        for (size_t i = 0; i < m_readers.size(); ++i) {
            switch (ATTRIBUTE) {
            case ParallelReader_Intensity:
                (*values)[i] = m_readers[i].getNextIntensity();
                break;
            case ParallelReader_Stdv:
                (*values)[i] = m_readers[i].getNextStdv();
                break;
            case ParallelReader_Pixels:
                (*values)[i] = m_readers[i].getNextPixels();
                break;
            case ParallelReader_Entry:
                // The reason the cast below is necessary is for the situation
                // in which ATTRIBUTE != ParallelReader_Entry, in which this
                // line will never be executed, but it must be compilable.
                m_readers[i].getNextEntry((affxcel::CELFileEntryType*)&((*values)[i]));
                break;
            default:
                assert(false);
                break;
            }
        }
    }

    // Returns the index of the next element to be read from all the cel files.
    size_t getNextIndex()
    {
        return m_readers[0].getNextIndex();
    }

    // Returns the file header for the file at the given index
    CONSTHACK affxcel::CCELFileHeaderData &getHeader(size_t whichFile) CONSTHACK
    {
        return m_readers[0].getHeader();
    }
    
  private:
    // Disallow copy ctor and assignment of this class
    CelFileParallelReader<C, ATTRIBUTE>(const CelFileParallelReader<C, ATTRIBUTE> &);
    CelFileParallelReader<C, ATTRIBUTE> &operator=(const CelFileParallelReader<C, ATTRIBUTE> &);

  private:
    std::vector<CelFileStreamReader> m_readers;
};

#endif /* _CELFILEPARALLELREADER_H */

/******************************************************************/
/**************************[END OF CelFileParallelReader.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
