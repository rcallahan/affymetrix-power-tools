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

#ifndef _CNReporterCnchp_H_
#define _CNReporterCnchp_H_
/**
 * @file CNReporterCnchp.h
 *
 * @brief This header contains the CNReporterCnchp class definition.
 */

#include "copynumber/CNReporter.h"

using namespace std;

class BaseEngine;
/**
 * @brief  The SNP data CN Reporter.
 *
 */
class CNReporterCnchp : public CNReporter
{
public:

    CNReporterCnchp();
    virtual ~CNReporterCnchp() {}

    virtual void defineOptions(BaseEngine& e);
    virtual void checkOptions(BaseEngine& e);
    virtual void run();

private:
    set<string> excludedOptionNames;
    bool isExcluded(const string& name);
};

#endif


