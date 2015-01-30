////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#ifndef SPFREADER_H
#define SPFREADER_H

//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "chipstream/CumulativeStats.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/ProbeSet.h"
#include "file/TsvFile/ClfFile.h"
#include "file/TsvFile/PgfFile.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <map>
#include <set>
#include <string>
#include <vector>

/**
 * Class for parsing the Simple Probe Format (SPF) text files. Reads
 * one record at a time to ease memory requirements
 */
class SpfReader {

public:

/** Open a text file for reading. */
void openSpf(const std::string &fileName);

/** Destructor to cleanup open files. */
~SpfReader() {
m_IncrementalIn.close();
}

/** Read the next record. Returned probeset must be deleted by
    caller */
ProbeSet *readNextProbeSet();

bool isFileOpen() {
return m_IncrementalIn.is_open();
}

int getNumProbeSets() {
  std::string numPSetsString;
  m_IncrementalIn.getHeader("num-probesets", numPSetsString);
  return Convert::toInt(numPSetsString.c_str());
}


private:

/* Read probeset from our buffers */
ProbeSet *readProbeSetRecord(std::string &name, int type, int numBlocks, std::string &blockSizesStr,
                             std::string &blockAnnsStr, int numMatch, int numProbes, std::string &probesStr,
                             int lineNumber);

/** Utility function to convert string to vector. */
void chopToIntVector(const std::string &s, char delim, std::vector<int> &vec);
/// Datafile to read in
affx::TsvFile m_IncrementalIn;
/// Buffers to read into
std::string m_IncName, m_IncBlockSizesStr, m_IncBlockAnnsStr, m_IncProbesStr;
int m_IncType, m_IncNumBlocks, m_IncNumMatch, m_IncNumProbes;

};

#endif
