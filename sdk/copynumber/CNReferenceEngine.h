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

#ifndef _CNReferenceEngine_H_
#define _CNReferenceEngine_H_
/**
 * @file CNReferenceEngine.h
 *
 * @brief This header contains the Log2RatioEngine class definition.
 */

#include "copynumber/CNLog2RatioData.h"
//
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//

using namespace std;

/**
 * @brief  The engine for calculating the log2 ratios.
 */
class CNReferenceEngine : public BaseEngine
{
private:
    BaseEngine* m_pEngine;
    CNLog2RatioData m_data;
    int m_iProbeSetCount;
    int m_iExperimentCount;
    int m_iProbeSetsToProcess;
    int m_iExperimentsToProcess;
    void defineOptions();
    void defineStates();
    void checkOptionsImp();
    void runImp();

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
        virtual std::string getEngineName() { return CNReferenceEngine::EngineName(); }
        static const std::string EngineName() { return "CNReferenceEngine"; }
        /**
         * Constructor
         */
        CNReferenceEngine();

        /**
         * Destructor
         */
        ~CNReferenceEngine();

        /*! A class to register the engine. */
        class Reg : public EngineReg
        {
        public:
            /*! Constructor - register the engine. */
            Reg() : EngineReg(CNReferenceEngine::EngineName())
            {
            }

            /*! Creates an object.
             * @return The object.
             */
            BaseEngine *MakeObject() { return new CNReferenceEngine; }
        };

        /*! The one and only registration object. */
        static Reg reg;

        /*! Converts the type to the engine type.
         * @param chip The pointer to the base engine object.
         * @return The engine type or NULL if not compatible.
         */
        static CNReferenceEngine * FromBase(BaseEngine *engine);


    void defineSharedOptions();

    void defineSharedStates();

    /**
     * @brief Set the engine pointer. If NULL then set to this.
     * @param BaseEngine* - The engine pointer to set from
     */
    void setEngine(BaseEngine* p)
    {
        if (p == NULL) {p = this;}
        m_pEngine = p;
    }

protected:
    bool determineMemoryUsage();
    bool processData();
    void transferProbeSetDataToSnps();
    bool combineOutput();
};


#endif


