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

#include "broadutil/CelUtil.h"
//
#include "broadutil/BroadException.h"
//
#include "file/CELFileData.h"
//

using namespace affx;
using namespace affxcel;
using namespace std;

vector<size_t> getDimensionsFromCelFile(const string &fileName)
{
    CCELFileData cel;
    cel.SetFileName(fileName.c_str());
    if (!cel.ReadHeader()) {
        throw BroadException("Error reading header", __FILE__, __LINE__,  fileName.c_str());
    }
    size_t numCols = cel.GetCols();
    size_t numRows = cel.GetRows();
    cel.Close();
    vector<size_t> dimensions(2);
    dimensions[0] = numRows;
    dimensions[1] = numCols;
    return dimensions;
}

// Copied from apt-cel-transformer.cpp
bool writeCelFile(affxcel::CCELFileWriter &cel, const string &format) {
  bool success = false;
  if(format == "text") 
    success = cel.WriteTextCel();
  else if(format == "xda")
    success = cel.WriteXDABCel();
  else if(format == "bcel")
    success = cel.WriteTranscriptomeBCel();
  else if(format == "ccel") 
    success = cel.WriteCompactBCel();
  else 
    throw BroadException("Don't recognize output format type", __FILE__, __LINE__, format.c_str());
  return success;
}

// Read 1-based probe #s from file into probeSet, converting to 0-based
void readProbeSetFromFile(std::vector<uint32_t> *probeSet, const char *probesetFilename)
{
    readSingleColumnFile(probeSet, probesetFilename, "probe_id");
    // Convert to 0-based
    for (size_t i = 0; i < probeSet->size(); ++i) {
        --((*probeSet)[i]);
    }
}

void readQuantilesFromFile(std::vector<double> *quantiles, const char *quantilesFilename)
{
    readSingleColumnFile(quantiles, quantilesFilename, "intensities");
}

void readCelFilesFromFile(std::vector<string> *celfiles, const char *celfilesFilename)
{
    readSingleColumnFile(celfiles, celfilesFilename, "cel_files");
}

/******************************************************************/
/**************************[END OF CelUtil.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
