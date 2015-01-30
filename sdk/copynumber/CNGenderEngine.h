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

#ifndef _CNGenderEngine_H_
#define _CNGenderEngine_H_
/**
 * @file CNGenderEngine.h
 *
 * @brief This header contains the CNGenderEngine class definition.
 */

//#include "copynumber/CNLog2RatioData.h"
//
//#include "chipstream/apt-geno-qc/GenoQC.h"
#include "util/BaseEngine.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include<string>



using namespace std;

/**
 * @brief  The engine for controlling the copynumber Cytos.
 */
class CNGenderEngine : public BaseEngine
{
private:

    /// List of aliases for standardized methods
    void checkOptionsImp();
    void runImp();
    void defineOptions();
    void extraHelp(){};
    void explain();

    void error(const std::string& strMessage) {
        Err::errAbort(strMessage);
    }
public:
        virtual std::string getEngineName() { return CNGenderEngine::EngineName(); }
        static const std::string EngineName() { return "CNGenderEngine"; }

        /**
         * Constructor
         */
        CNGenderEngine();

        /**
         * Destructor
         */
        ~CNGenderEngine();


protected:

private:
    bool checkForInvalidCelFiles(const std::vector<std::string>& celFiles);
};


#endif

