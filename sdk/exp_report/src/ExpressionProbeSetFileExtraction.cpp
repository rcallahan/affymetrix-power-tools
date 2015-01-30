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

#include "exp_report/src/ExpressionProbeSetFileExtraction.h"
//
#include "calvin_files/parameter/src/ParameterFileData.h"
#include "calvin_files/parsers/src/ParameterFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "file/TsvFile/TsvFile.h"
//

using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_io;
using namespace affx;
using namespace std;

/*
 * Extract the parameters from the parameter file.
 */
bool ExpressionProbeSetFileExtraction::ExtractParameters
(
    const char *fileName, 
    ProbeSetFileEntryMap &probeSets
)
{
    // Clear the map
    probeSets.clear();

    // Open the file.
    TsvFile tsv;
    try
    {
        if (tsv.open(fileName) == TSV_OK)
        {
            // Bind the columns and read the data.
            string id;
            string name;
            string category;
            string signal_in_header;
            tsv.bind(0,"probeset_id",&id);
            tsv.bind(0,"probeset_name",&name);
            tsv.bind(0,"group_name",&category);
            tsv.bind(0,"quantification_in_header", &signal_in_header);
            while(tsv.nextLevel(0) == TSV_OK)
            {
                if (signal_in_header == "0")
                    continue;

                if (name.length() == 0)
                    continue;

                if (category.length() == 0)
                    probeSets[id] = name;
                else
                    probeSets[id] = category + string("-") + name;
            }
            tsv.close();
            return true;
        }
    }
    catch(...) {}
    return false;
}

