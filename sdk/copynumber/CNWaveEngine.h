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

#ifndef _CNWaveEngine_H_
#define _CNWaveEngine_H_
/**
 * @file CNWaveEngine.h
 *
 * @brief This header contains the CNWaveEngine class definition.
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
class CNWaveEngine : public BaseEngine
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
        virtual std::string getEngineName() { return CNWaveEngine::EngineName(); }
        static const std::string EngineName() { return "CNWaveEngine"; }
        /**
         * Constructor
         */
        CNWaveEngine();

        /**
         * Destructor
         */
        ~CNWaveEngine();

        /*! A class to register the engine. */
        class Reg : public EngineReg
        {
        public:
            /*! Constructor - register the engine. */
            Reg() : EngineReg(CNWaveEngine::EngineName())
            {
            }

            /*! Creates an object.
             * @return The object.
             */
            BaseEngine *MakeObject() { return new CNWaveEngine; }
        };

        /*! The one and only registration object. */
        static Reg reg;

        /*! Converts the type to the engine type.
         * @param chip The pointer to the base engine object.
         * @return The engine type or NULL if not compatible.
         */
        static CNWaveEngine * FromBase(BaseEngine *engine);

    void addCychp(const std::string& strCelFileName);

    static bool isCN5Reference(const AffxString& strReferenceFileName);
    static bool isCopyNumberReference(const AffxString& strReferenceFileName);
    // This method is called by ChAS without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static bool isCyto2Reference(const AffxString& strReferenceFileName);
    // This method is called by ChAS without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationParameter(const AffxString& strFileName, const AffxString& strParamterName);
    // This method is called by ChAS without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationArrayType(const AffxString& strFileName) {return getAnnotationParameter(strFileName, "affymetrix-array-type");}
    // This method is called by ChAS without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationVersion(const AffxString& strFileName) {return getAnnotationParameter(strFileName, "affymetrix-algorithm-param-netaffx-build");}

protected:
    void modifyReference();
    void defineStdMethods();
};


#endif


