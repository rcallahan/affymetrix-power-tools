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

#ifndef _CNReferenceMethodAdditionalWaves_h_
#define _CNReferenceMethodAdditionalWaves_h_

/**
 * @file CNWaveCorrection.h
 *
 * @brief This header contains the CNWaveCorrection class definition.
 */

#define NEWMAT_USING_FLOAT
#include "copynumber/CNReferenceMethodWaveCorrection.h"
//
#include "util/AffxArray.h"
//
#include "../external/newmat/newmatap.h"
//
#include <vector>
//

/**
 * @brief  A class for calulating WaveCorrection.
 *
 */
class CNReferenceMethodAdditionalWaves : public CNReferenceMethodWaveCorrection
{
public:
    static std::string getType() { return "additional-waves-reference-method"; }
    static std::string getDescription() { return "CopyNumber AdditionalWaves"; }
    static std::string getVersion() { return "1.0"; }

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate * newObject(std::map<std::string, std::string> & param);
};

#endif // _CNReferenceModuleAdditionalWaves_h_


