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

#include "mas5-stat/src/IntensityFileType.h"
//
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;

/*
 * Check the DAT header for the HP style.
 */
static bool IsDatHeaderStringFromHP(const string &datHeader)
{
    int index = (int) datHeader.find(":CLS=");
    string strTypeAndID = datHeader.substr(index+1);
    strTypeAndID = strTypeAndID.substr(67);

    // There are several sub-fields in this field. The 
    // first sub field is the scanner ID, followed by the scanner type, 
    // followed by three spaces. If the scanner ID is absent, 
    // the field consists of four spaces
    // Find the 3 spaces that are after the scanner type
    int iStartScanType = (int) strTypeAndID.find(0x14);

    // If the start of the 3 spaces is at the begining, there is 
    // no scannertype and scanner ID
    string strScannerType;
    if (iStartScanType > 0)
    {
        strScannerType = strTypeAndID.substr(0, iStartScanType - 3);
        iStartScanType = (int) strScannerType.rfind(" ");
        if(iStartScanType <= 0)
            return true;
        
        strScannerType = strScannerType.substr(iStartScanType + 1, strScannerType.length() - iStartScanType - 1);
        return (strScannerType == "HP");
    }
    return true;
}

/*
 * Return a flag indicating if the CEL file is from an HP scanner.
 */
bool FromHP(FusionCELData &cel)
{ 
    // Check if from the HT scanner. This had a upper left grid of 1,1
    affymetrix_fusion_io::FGridCoords grid = cel.GetGridCorners();
    if (grid.upperleft.x < 1.001f && grid.upperleft.y < 1.001f)
        return false;

    GenericData *gdata = cel.GetGenericData();
    if (gdata != NULL)
    {
        // Check for the scanne type parameter.
        int nparents = gdata->Header().GetGenericDataHdr()->GetParentCnt();
        for (int iparent=0; iparent<nparents; ++iparent)
        {
            GenericDataHeader pheader = gdata->Header().GetGenericDataHdr()->GetParent(iparent);
            ParameterNameValueType nvt;
            if (pheader.FindNameValParam(SCANNER_TYPE_PARAM_NAME, nvt) == true)
            {
                wstring val = nvt.ToString();
                if (val.size() == 0 || val == L"HP")
                    return true;
            }
            else if (pheader.FindNameValParam(DAT_HEADER_PARAM_NAME, nvt) == true)
            {
                wstring val = nvt.ToString();
                return IsDatHeaderStringFromHP(StringUtils::ConvertWCSToMBS(val));
            }
        }
        return false;
    }
    else
    {
        // Check the DAT header. It it contains HP or blank for the scanner type
        // then it came from the HP scanner.
        string datHeader = StringUtils::ConvertWCSToMBS(cel.GetDatHeader());
        return IsDatHeaderStringFromHP(datHeader);
    }
}

