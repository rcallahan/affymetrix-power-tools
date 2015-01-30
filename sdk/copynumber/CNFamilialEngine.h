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

/**
 * @file   CNFamilialEngine.h
 *
 * @brief Analyzes trios and duos or ccychp files.
 */
#ifndef _CNFamilialEngine_H_
#define _CNFamilialEngine_H_

#include "copynumber/CNCychp.h"
#include "copynumber/CNFamilialAnalysisMethod.h"
#include "copynumber/CNFamilialReporter.h"
//
#include "util/BaseEngine.h"
//
#include <cstring>
#include <ctime>
#include <ostream>
#include <string>
#include <set>
#include <map>
//

typedef std::map<std::string, std::set<std::string> > groupDatasets_t;

class SegmentOverlap;
class CNFamilialEngine : public BaseEngine {

    public:

        virtual std::string getEngineName() { return CNFamilialEngine::EngineName(); }
        static const std::string EngineName() { return "CNFamilialEngine"; }
        /**
         * Constructor
         */
        CNFamilialEngine();

        /**
         * Destructor
         */
        ~CNFamilialEngine();

        /*! A class to register the engine. */
        class Reg : public EngineReg
        {
        public:
            /*! Constructor - register the engine. */
            Reg() : EngineReg(CNFamilialEngine::EngineName())
            {
            }

            /*! Creates an object.
             * @return The object.
             */
            BaseEngine *MakeObject() { return new CNFamilialEngine; }
        };

        /*! The one and only registration object. */
        static Reg reg;

        /*! Converts the type to the engine type.
         * @param chip The pointer to the base engine object.
         * @return The engine type or NULL if not compatible.
         */
        static CNFamilialEngine * FromBase(BaseEngine *engine);

        CNCychp& getCychpIndex() {return m_cychpIndex;}
        CNCychp& getCychpMother() {return m_cychpMother;}
        CNCychp& getCychpFather() {return m_cychpFather;}

        void clear();

    private:
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

        void checkOptionsImp();
        void runImp();

        void defineOptions();
        void defineStates();
        void defineStdMethods();
        void reportBasics();

        void createAnalysis();
        void prepareGroupsDatasets();

        void loadCychpFiles();
        void familialAnalysis();
        void loadSegmentOverlaps(AffxArray<SegmentOverlap>& vSegmentOverlaps);
        void loadParam(const AffxString& strName, PgOpt::PgOptType type, const AffxString& strValue, affymetrix_calvin_parameter::ParameterNameValueType& param);
        void checkFamilialConsistency();
        void checkFamilialConsistencyPair(CNCychp& cychp1, CNCychp& cychp2);

    private:
        AffxArray<CNFamilialAnalysisMethod> m_vCNFamilialAnalysisMethods;
        AffxArray<CNFamilialReporter> m_vCNFamilialReporters;

        std::map<std::string,std::string>* m_pstdMethods;

        // Calvin cychp groups/datasets required by analyses stored in m_vCNFamilialAnalysisMethods
        groupDatasets_t m_groupDatasetsAllFamilial;

        CNCychp m_cychpIndex;
        CNCychp m_cychpMother;
        CNCychp m_cychpFather;
};

#endif /* _CNFamilialEngine_H_ */

