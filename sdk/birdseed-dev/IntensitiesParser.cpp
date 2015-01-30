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

//
#include "birdseed-dev/IntensitiesParser.h"
//
#include "broadutil/BroadException.h"
#include "broadutil/BroadUtil.h"
#include "broadutil/FileUtil.h"
//
#include <cerrno>
#include <cstdio>

using namespace std;
using namespace birdseed::dev;

IntensitiesParser::~IntensitiesParser()
{
    if (fp != NULL) {
        fclose(fp);
    }
}

IntensitiesParser::IntensitiesParser(const std::string& path, bool calculateCorrectionFactor):
    dataPath(path),
    fp(fopen(path.c_str(), "r")),
    correctionFactor(0.0),
    header(),
    snpName(),
    intensities()
{
    if (fp == NULL) {
      throw BroadException("Could not open file", __FILE__, __LINE__, path.c_str(), errno);
    }
    header = readHeader(fp, "probeset_id\t");
    if (calculateCorrectionFactor) {
        long dataStartOffset = ftell(fp);
        if (dataStartOffset == -1) {
          throw BroadException("ftell error", __FILE__, __LINE__, path.c_str(), errno);
        }
        double total = 0.0;
        size_t numIntensities = 0;
        while (true) {
            // Skip over SNP name
            snpName = readWord(fp);
            if (snpName.empty()) {
                break;
            }
            double val;
            while (readDouble(fp, &val)) {
                total += val;
                ++numIntensities;
            }
        }
        correctionFactor = total/numIntensities;

        if (fseek(fp, dataStartOffset, SEEK_SET) == -1) {
          throw BroadException("fseek error", __FILE__, __LINE__, path.c_str(), errno);
        }
    }
}

bool IntensitiesParser::advanceSNP()
{
    // Skip over first SNP name
    string snpNameA = readWord(fp);
    if (snpNameA.empty()) {
        snpName = "";
        intensities.reset(NULL);
        return false;
    }
    if (!endswith(snpNameA.c_str(), "-A")) {
        std::stringstream strm;
        strm << "Did not find -A suffix for SNP name: " << snpNameA << ".";
        throw BroadException(strm.str().c_str(), __FILE__, __LINE__, dataPath.c_str());
    }
    snpName = snpNameA.substr(0, snpNameA.size() - 2);
    
    vector<double> aValues;
    double val;
    while (readDouble(fp, &val)) {
        aValues.push_back(val);
    }

    // Skip over next SNP name
    string snpNameB = readWord(fp);
    if (!endswith(snpNameB.c_str(), "-B")) {
        std::stringstream strm;
        strm << "Did not find -B suffix for SNP name: " << snpNameB << ".";
        throw BroadException(strm.str().c_str(), __FILE__, __LINE__, dataPath.c_str());
    }
    if (snpName != snpNameB.substr(0, snpNameB.size() - 2)) {
        std::stringstream strm;
        strm << "B intensities do not have same SNP name as A.  SNP-A: " << snpNameA << "; SNP-B: " << snpNameB << " .";
        throw BroadException(strm.str().c_str(), __FILE__, __LINE__, dataPath.c_str());
    }
        
    intensities.reset(new IntensityMatrix(aValues.size()));

    // Syntactic sugar
    IntensityMatrix &intensitiesRef = *(intensities.get());
    
    for (size_t i = 0; i < aValues.size(); ++i) {
        intensitiesRef[i][0] = aValues[i];
        if (!readDouble(fp, &val)) {
            throw BroadException("Number of B values is less than number of A values", __FILE__, __LINE__, dataPath.c_str());
        }
        intensitiesRef[i][1] = val;
    }
    return true;
}

/******************************************************************/
/**************************[END OF IntensitiesParser.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
