////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////


#ifndef _GCOSParameterNames_HEADER_
#define _GCOSParameterNames_HEADER_

namespace affymetrix_fusion_io
{

// Should we move these to their own file?
/*! The parameter name for DATIO7 ImageData::fArcRadius (float) */
#define ARCRADIUS_PARAM_NAME L"affymetrix-fusion-arcradius"

/*! The parameter name for DATIO7 ImageData::fLaserSpotSize (float) */
#define LASER_SPOTSIZE_PARAM_NAME L"affymetrix-fusion-laser-spotsize"

/*! The parameter name for DATIO7 ImageData::SetExperimentName() (text) */
#define EXPERIMENT_NAME_PARAM_NAME L"affymetrix-fusion-experiment-name"

/*! The parameter name format for DATIO7 ImageData::GetChipInfo() (text)
 *	This is a simulated parameter.  It's value is derived from many other 
 *	parameter values.  This parameter will not appear in the iteration
 *	returned by GetParameterIts.
 */
#define GET_CHIP_INFO_PARAM_NAME L"affymetrix-fusion-get-chipinfo"
/*! The parameter name for DATIO7 ImageData::GetScanInfo() (text)
 *	This is a simulated parameter.  It's value is derived from many other 
 *	parameter values.  This parameter will not appear in the iteration
 *	returned by GetParameterIts.
 */
#define GET_SCAN_INFO_PARAM_NAME L"affymetrix-fusion-get-scaninfo"

/*! The parameter name format for DATIO7 ImageData::SetChipInfo(int) and ImageData::GetChipInfo(int) (text)
 *	The value for the variable part of the name may be between 0-9 inclusive
 *	Not sure that this is beig used.
 */
#define CHIPINFO_NAME_FORMAT L"affymetrix-fusion-set-chipinfo-%d"

}

#endif	//_GCOSParameterNames_HEADER_
