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
 * FILE CelUtil.h
 */

#ifndef _CELUTIL_H_
#define _CELUTIL_H_

//
#include "broadutil/BroadException.h"
//
#include "file/CELFileWriter.h"
#include "file/TsvFile/TsvFile.h"
//
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>
//

std::vector<size_t> getDimensionsFromCelFile(const std::string &fileName);

bool writeCelFile(affxcel::CCELFileWriter &cel, const std::string &format);

// Read 1-based probe #s from file into probeSet, converting to 0-based
// JHG: NOTE: size_t changes size from host to host. fix its size at uint32_t so
// it doesnt break TsvFile with the template.
// void readProbeSetFromFile(std::vector<size_t> *probeSet, const char *probesetFilename);
void readProbeSetFromFile(std::vector<uint32_t> *probeSet, const char *probesetFilename);

void readQuantilesFromFile(std::vector<double> *quantiles, const char *quantilesFilename);

void readCelFilesFromFile(std::vector<std::string> *celfiles, const char *celfilesFilename);

template<class T>
void readSingleColumnFile(std::vector<T> *outputVec, const std::string& filename, const std::string& columnName)
{
    outputVec->clear();
    affx::TsvFile tsv;
    if (columnName == "") {
        tsv.m_optHasColumnHeader = false;
        tsv.m_optAutoColumns = true;
    }
    if(tsv.open(filename) != affx::TSV_OK) {
      throw BroadException("Error opening file", __FILE__, __LINE__, filename.c_str(), errno);
    }
    T val;
    if (columnName != "") {
        tsv.bind(0, columnName, &val, affx::TSV_BIND_REQUIRED);
    }
    while(tsv.nextLevel(0) == affx::TSV_OK) {
        if (columnName == "") {
            tsv.get(0, 0, val);
        }
        outputVec->push_back(val);
    }
    tsv.close();
}
    

#endif /* _CELUTIL_H_ */

/******************************************************************/
/**************************[END OF CelUtil.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
