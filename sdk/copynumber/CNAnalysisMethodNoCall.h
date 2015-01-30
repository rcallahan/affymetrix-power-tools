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

#ifndef _CNAnalysisMethodNoCall_H_
#define _CNAnalysisMethodNoCall_H_
/**
 * @file CNAnalysisMethodNoCall.h
 *
 * @brief This header contains the CNAnalysisMethodNoCall class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

/**
 * @brief  The NoCall analysis method.
 *
 */
class CNAnalysisMethodNoCall : public CNAnalysisMethod
{
public:
    static std::string getType() {return "no-call";}
    static std::string getDescription() {return "Copynumber NoCall";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodNoCall();
    virtual ~CNAnalysisMethodNoCall() {}

    virtual bool isSegmentTypeAnalysis() {return true;}

    virtual AffxString getName() {return getType();}

    virtual void run();
};

#endif


