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

#include "broadutil/CelFileStreamReader.h"
//
#include "broadutil/BroadException.h"
//
#include "file/FileIO.h"
//
#include <cerrno>
#include <cstdio>
#include <iostream>
#include <sstream>
//

using namespace std;
using namespace affxcel;

CONSTHACK affxcel::CCELFileHeaderData &CelFileStreamReader::getHeader() CONSTHACK  {
     return m_header;
}

// Returns the index of the next element to be read from the cel file.
size_t CelFileStreamReader::getNextIndex()  {
  return m_nextIndex;
}

CelFileStreamReader::CelFileStreamReader(const string &filename):
    m_filename(filename),
    m_header(),
    m_nextIndex(0),
    m_instream()
{
    CCELFileData cel;
    cel.SetFileName(m_filename.c_str());
    if (!cel.ReadHeader()) {
        throw BroadException("Error reading cel header", __FILE__, __LINE__, m_filename.c_str());
    }
    if (cel.GetFileFormat() != CCELFileData::XDA_BCEL) {
        ostringstream errmsg("Wrong file format: ");
        errmsg <<  cel.GetFileFormat();
        throw BroadException(errmsg.str().c_str(), __FILE__, __LINE__, m_filename.c_str());
    }
    // Fortunately, default copy ctor should work since no ptrs.
    m_header = cel.GetHeader();

}

// A copy starts a new iteration, because we can't copy the file pointer.
CelFileStreamReader::CelFileStreamReader(const CelFileStreamReader& other):
    m_filename(other.m_filename),
    m_header(other.m_header),
    m_nextIndex(0),
    m_instream()
{
}

// A copy starts a new iteration, because we can't copy the file pointer.
CelFileStreamReader& CelFileStreamReader::CelFileStreamReader::operator=(const CelFileStreamReader &other)
{
    m_filename = other.m_filename;
    m_header = other.m_header;
    m_nextIndex = 0;
    if (m_instream.is_open()) {
        m_instream.close();
    }
    return *this;
}



CelFileStreamReader::~CelFileStreamReader()
{
    if (m_instream.is_open()) {
        m_instream.close();
    }
}

// Open the file for iterating
void CelFileStreamReader::startIterating()
{
    size_t headerSize = determineHeaderSize(m_filename);

    m_instream.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit);
	m_instream.open(m_filename.c_str(), std::ios::in | std::ios::binary);

    if (!m_instream) {
        throw BroadException("Error opening file", __FILE__, __LINE__, m_filename.c_str(), errno);
    }
    m_instream.seekg(headerSize);
}

static void seek_throw(FILE *stream, long offset, int whence, const char *filename)
{
    if (fseek(stream, offset, whence) != 0) {
        throw BroadException("Seek error", __FILE__, __LINE__, filename, errno);
    }
}

static size_t tell_throw(FILE *stream, const char *filename)
{
    int ret = ftell(stream);
    if (ret == -1) {
        throw BroadException("Tell error", __FILE__, __LINE__, filename, errno);
    }
    return size_t(ret);
}


static void read_throw(void *ptr, size_t size, FILE *stream, const char *filename)
{
    if (fread(ptr, size, 1, stream) != 1) {
        throw BroadException("Read error", __FILE__, __LINE__, filename, errno);
    }
}

static void skip_string(FILE *stream, const char*filename)
{
    uint32_t slen;
    read_throw(&slen, sizeof(slen), stream, filename);
    seek_throw(stream, slen, SEEK_CUR, filename);
}


size_t CelFileStreamReader::determineHeaderSize(const string &filename)
{
    // It would be nice if this were incorporated into CCelFileData

    // I tried using ifstream, but got hard-to-debug problem after seeking.
    // It turns out the problem was bad data, not bad code, but I'm leaving this using
    // stdio since it should work fine either way.
	FILE *headerstrm = fopen(filename.c_str(), "rb");
	if (headerstrm == NULL)	{
        throw BroadException("Error opening file", __FILE__, __LINE__, filename.c_str(), errno);
    }
    seek_throw(headerstrm, INT_SIZE + // magic number
               INT_SIZE + // version
               INT_SIZE + // num rows
               INT_SIZE + // num cols
               INT_SIZE, // num cells
               SEEK_SET,
               filename.c_str());

    // Header string
    skip_string(headerstrm, filename.c_str());
    // Algorithm name
    skip_string(headerstrm, filename.c_str());
    // Algorithm params
    skip_string(headerstrm, filename.c_str());

    // Skip over cell margin, # of outlier cells, # of masked cells, # of sub-grids
    size_t ret = tell_throw(headerstrm, filename.c_str()) + 2 * INT_SIZE + 2 * UINT32_SIZE;
    fclose(headerstrm);
    return ret;
}



float CelFileStreamReader::getNextIntensity()
{
    return MmGetFloat_I(&(getNextEntryUnswapped().Intensity));
}

float CelFileStreamReader::getNextStdv()
{
    return MmGetFloat_I(&(getNextEntryUnswapped().Stdv));
}

short CelFileStreamReader::getNextPixels()
{
    return MmGetInt16_I(&(getNextEntryUnswapped().Pixels));
}

void CelFileStreamReader::getNextEntry(CELFileEntryType *val)
{
    CELFileEntryType unswapped = getNextEntryUnswapped();
    val->Intensity = MmGetFloat_I(&unswapped.Intensity);
    val->Stdv = MmGetFloat_I(&unswapped.Stdv);
    val->Pixels = MmGetInt16_I(&unswapped.Pixels);
}

affxcel::CELFileEntryType CelFileStreamReader::getNextEntryUnswapped()
{
  if (!m_instream.is_open()) {
    startIterating();
  }
  ++m_nextIndex;
  affxcel::CELFileEntryType entry;
  m_instream.read((char *)&entry, sizeof(affxcel::CELFileEntryType));
  return entry;
}

/******************************************************************/
/**************************[END OF CelFileStreamReader.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
