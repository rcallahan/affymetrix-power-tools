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

#ifndef _CNCorrelationEngine_H_
#define _CNCorrelationEngine_H_
/**
 * @file CNCorrelationEngine.h
 *
 * @brief This header contains the CNCorrelationEngine class definition.
 */

#include "util/AffxString.h"
#include "util/BaseEngine.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//

using namespace std;

/**
 * @brief  The engine for controlling the copynumber Cytos.
 */
class CNCorrelationEngine : public BaseEngine
{
private:
    /// List of aliases for standardized methods
    std::map<std::string,std::string>* m_pstdMethods;

    int* m_pGenderCalls;

    void checkOptionsImp();
    void runImp();
    void defineOptions();
    void defineStates();
    void printStandardMethods(std::ostream &out);
    void extraHelp();
    void explain();

    /**
     * Print out a warning message
     */
    void warning(const std::string& strMessage) {
        Verbose::out(1, "WARNING: " + strMessage);
    }

    /**
     * Error out by calling clear() and Err::errAbort()
     */
    void error(const std::string& strMessage) {
        clear();
        Err::errAbort(strMessage);
    }

public:
        void clear();
        void checkParameters();
        virtual std::string getEngineName() { return CNCorrelationEngine::EngineName(); }
        static const std::string EngineName() { return "CNCorrelationEngine"; }
        /**
         * Constructor
         */
        CNCorrelationEngine();

        /**
         * Destructor
         */
        ~CNCorrelationEngine();

        /*! A class to register the engine. */
        class Reg : public EngineReg
        {
        public:
            /*! Constructor - register the engine. */
            Reg() : EngineReg(CNCorrelationEngine::EngineName())
            {
            }

            /*! Creates an object.
             * @return The object.
             */
            BaseEngine *MakeObject() { return new CNCorrelationEngine; }
        };

        /*! The one and only registration object. */
        static Reg reg;

        /*! Converts the type to the engine type.
         * @param chip The pointer to the base engine object.
         * @return The engine type or NULL if not compatible.
         */
        static CNCorrelationEngine * FromBase(BaseEngine *engine);

    void addCychp(const std::string& strCelFileName);

protected:
    void calculateCorrelation();
    void defineStdMethods();
};


#endif


