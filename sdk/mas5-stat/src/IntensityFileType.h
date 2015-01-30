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

/*! \file IntensityFileType.h This file contains functions to check if a CEL file is from an HP scanner. */

#ifndef _IntensityFileType_HEADER_
#define _IntensityFileType_HEADER_

#include "calvin_files/fusion/src/FusionCELData.h"
//

/*! Checks if the CEL file is from an HP scanner.
 * @param cel The CEL file object.
 * @return True if from a HP scanner.
 */
bool FromHP(affymetrix_fusion_io::FusionCELData &cel);

#endif

