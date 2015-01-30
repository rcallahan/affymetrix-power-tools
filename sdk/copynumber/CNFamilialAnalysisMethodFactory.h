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

#ifndef _CNFamilialAnalysisMethodFactory_H_
#define _CNFamilialAnalysisMethodFactory_H_
/**
 * @file CNFamilialAnalysisMethodFactory.h
 *
 * @brief This header contains the CNFamilialAnalysisMethodFactory class definition.
 */

#include "copynumber/CNFamilialAnalysisMethod.h"
#include "copynumber/CNFamilialReporter.h"
//
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//

/**
 * @brief  A factory class for making copy number analysis methods based on a string representation.
 *
 */
class CNFamilialAnalysisMethodFactory
{
protected:
    std::vector<SelfDoc> m_Docs;
    std::vector<SelfCreate::selfCreator> m_Creators;

public:
    CNFamilialAnalysisMethodFactory();
    virtual ~CNFamilialAnalysisMethodFactory();

    /**
     * @brief  Register the method
     * @param const SelfDoc& - The SelfDoc object
     * @param SelfCreate::selfCreator - The newObject function
     */
    void registerMethod(const SelfDoc& doc, SelfCreate::selfCreator maker)
    {
        m_Docs.push_back(doc);
        m_Creators.push_back(maker);
    }

    CNFamilialAnalysisMethod* CNFamilialAnalysisMethodForString(const std::string &spec);
    CNFamilialReporter* CNFamilialReporterForString(const std::string &spec);

    std::vector<SelfDoc> getDocs() {return m_Docs;}
};

#endif


