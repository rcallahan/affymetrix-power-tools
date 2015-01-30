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

#ifndef _CNAnalysisMethodSegment_H_
#define _CNAnalysisMethodSegment_H_
/**
 * @file CNAnalysisMethodSegment.h
 *
 * @brief This header contains the CNAnalysisMethodSegment class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

/**
 * @brief  The Segment analysis method.
 *
 */
class CNAnalysisMethodSegment : public CNAnalysisMethod
{
public:
    CNAnalysisMethodSegment();
    virtual ~CNAnalysisMethodSegment() {}

    virtual bool isSegmentTypeAnalysis() {return true;}

    virtual AffxString getName() = 0;

    virtual void run() = 0;

    virtual int newSegments(	int iSegmentType, 
                                CNProbeSetArray*,
                                int iMinSegSeparation=1000000000 );
};

#endif


