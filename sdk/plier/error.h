////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The PLIER (Probe Logarithmic Error Intensity Estimate) method produces
an improved signal by accounting for experimentally observed patterns 
in probe behavior and handling error at the appropriately at low and 
high signal values.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/

 /*
 * error.h: error code for PLIER
 *
 */

#if !defined(PLIER_ERROR_HEADER)
#define PLIER_ERROR_HEADER

#define NO_PLIER_ERROR          0
#define PLIER_ERROR             8000
#define NO_DATAMEM              PLIER_ERROR + 1
#define INVALID_NUM_EXP         PLIER_ERROR + 2
#define INVALID_NUM_FEATURE       PLIER_ERROR + 3
#define INVALID_PM              PLIER_ERROR + 4
#define INVALID_MM              PLIER_ERROR + 5
#define INVALID_TARGETRESPONSE   PLIER_ERROR + 6
#define INVALID_FEATURERESPONSE        PLIER_ERROR + 7
#define INVALID_AUGMENTATION    PLIER_ERROR + 8
#define INVALID_GMCUTOFF        PLIER_ERROR + 9
#define INVALID_DROPMAX         PLIER_ERROR + 10
#define INVALID_DIFFERENTIALTARGETPENALTY     PLIER_ERROR + 11
#define INVALID_DIFFERENTIALFEATUREPENALTY    PLIER_ERROR + 12
#define INVALID_OPTIMIZATION    PLIER_ERROR + 13
#define INVALID_SEAITERATION    PLIER_ERROR + 14
#define INVALID_PLIERITERATION  PLIER_ERROR + 15
#define MAXIT_SEA_REACHED       PLIER_ERROR + 16
#define MAXIT_PLIER_REACHED     PLIER_ERROR + 17
#define INVALID_SAFETYZERO	PLIER_ERROR + 18
#define INVALID_NUMERICALCONVERGENCE PLIER_ERROR + 19

#define PLIER_ERROR_START       NO_DATAMEM
#define PLIER_ERROR_END         MAXIT_PLIER_REACHED

#endif /* !defined(PLIER_ERROR_HEADER) */
