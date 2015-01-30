////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

#include "exp_report/src/ExpressionControlsParameterExtraction.h"
//
#include "calvin_files/parameter/src/ParameterFileData.h"
#include "calvin_files/parsers/src/ParameterFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
//

using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_io;
using namespace ExpressionReport;
using namespace std;

/*! A map between a parameter name and control pointer. */
typedef map<wstring, ExpressionControl *> ExpressionControlMap;

/*
 * Extract the parameters from the parameter file.
 */
bool ExpressionControlsParameterExtraction::ExtractParameters
(
    const char *fileName, 
    int &probePairThreshold,
    ExpressionControls &controls
)
{
    // Clear the controls
    controls.Clear();

    // Read the file
    ParameterFileData paramData;
    ParameterFileReader reader;
    if (reader.Read(fileName, paramData) == false)
        return false;

    // Map the parameters from the file to the controls object.
    ExpressionControlMap spike;
    ExpressionControlMap housekeeping;
    for (ParameterTypeList::iterator it=paramData.Parameters().begin(); it!=paramData.Parameters().end(); it++)
    {
        bool isSpike = false;
        bool isHouse = false;

        // Get the probe array type
        if (it->name == L"ArrayType")
        {
            controls.ProbeArrayType() = StringUtils::ConvertWCSToMBS(it->currentValue);
        }

        else if (it->name == L"ProbePairThreshold")
        {
            string mbs = StringUtils::ConvertWCSToMBS(it->currentValue);
            probePairThreshold = strtol(mbs.c_str(),NULL,10);
        }

        // Check if spike control
        else if (it->name.find(L"spike") != -1)
        {
            isSpike = true;
        }

        // Check if housekeeping control
        else if (it->name.find(L"housekeeping") != -1)
        {
            isHouse = true;
        }

        // Extract the probe set index and type from the parameter.
        if (isSpike == true || isHouse == true)
        {
            ExpressionControl *control;
            ExpressionControlMap *controlMap;
            if (isSpike == true)
                controlMap = &spike;
            else
                controlMap = &housekeeping;

            ExpressionControlMap::iterator mapIt = controlMap->find(it->displayName);
            if (mapIt != controlMap->end())
            {
                control = mapIt->second;
            }
            else
            {
                control = new ExpressionControl();
                control->Name() = StringUtils::ConvertWCSToMBS(it->displayName);
                (*controlMap)[it->displayName] = control;
            }
            string mbs = StringUtils::ConvertWCSToMBS(it->index);
            int index = strtol(mbs.c_str(),NULL,10);
            if (it->category == L"3")
                control->SetProbeSetIndex(ExpressionControl::THREE_PRIME_PROBE_SET, index);
            else if (it->category == L"M")
                control->SetProbeSetIndex(ExpressionControl::MIDDLE_PROBE_SET, index);
            else if (it->category == L"5")
                control->SetProbeSetIndex(ExpressionControl::FIVE_PRIME_PROBE_SET, index);
        }
    }

    // Add the controls to the object
    for(ExpressionControlMap::iterator it=spike.begin(); it!=spike.end(); ++it)
    {
        controls.SpikeControls().push_back(*it->second);
    }
    for(ExpressionControlMap::iterator it=housekeeping.begin(); it!=housekeeping.end(); ++it)
    {
        controls.HousekeepingControls().push_back(*it->second);
    }

    return true;
}

