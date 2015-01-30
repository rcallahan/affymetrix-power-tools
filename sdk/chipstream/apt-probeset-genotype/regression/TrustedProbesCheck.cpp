////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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
 * @file   TrustedProbesCheck.cpp
 * @brief  Extend MatrixCheck to allow nested columns
 * where the outer columns are tab delimited. 
 */

#include "chipstream/apt-probeset-genotype/regression/TrustedProbesCheck.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/AffxByteArray.h"

#include "util/Fs.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <string>
//

bool PosteriorMatrixCheck::check(std::string &msg) {

  Verbose::out(1,m_Generated + " post processing file by changing commas to tabs.");

  affx::TsvFile::replaceCharInFile(m_Generated, ',' , '\t');
  
  return MatrixCheck::check(msg);
}



