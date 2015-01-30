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

#ifndef CELLISTENER_H
#define CELLISTENER_H

#include "chipstream/ChipSummary.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
//

/** 
 * Class for objects that need to actually access real cel file data.
 */
class CelListener {

public:

  /** 
   * Cel file of data.
   * 
   * @param cel - Handle to an open cel file.
   */
  virtual void newChip(affymetrix_fusion_io::FusionCELData *cel) = 0;

  /**
   * Virtual destructor
   */
  virtual ~CelListener(){}

};

#endif /* CELLISTENER_H */
