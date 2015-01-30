////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
/**
 * @file   ATDRegression.h
 * @author Mybrid Spalding
 * @date   Wed Jun  4 10:20:47 PDT 2008
 * @brief  DMET3 Allele translation regression header.
 */

#ifndef TRANSLATION_REGRESSION_H
#define TRANSLATION_REGRESSION_H

#include <translation/RunTimeEnvironment.h>
//
#include <cstring>
#include <string>
//

// CONSTANTS

const std::string INPUT_DIR           = TEST_DATA_REGRESSION_DIR;
const std::string TRANSLATION_PROGRAM = "apt-dmet-translation";
  
const std::string DMET2_REGRESSION_MARKER_FILE_EXT = "dmet2_marker.reg";
const std::string DMET2_REGRESSION_HAPLOTYPE_FILE_EXT = "dmet2_haplotype.reg";


#endif /* TRANSLATION_REGRESSION__H */
