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

#ifndef _CNAnalysisEngine_H_
#define _CNAnalysisEngine_H_
/**
 * @file CNAnalysisEngine.h
 *
 * @brief This header contains the CNAnalysisEngine class definition.
 */

#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNExperiment.h"
#include "copynumber/CNReporter.h"
#include "copynumber/CNReporterCnchp.h"
//
#include "util/AffxMultiDimensionalArray.h"
#include "util/BaseEngine.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//

using namespace std;

struct CalibrationParams
{
    double alphaCNCalibration;
    double alphaCNCalibration_X;
    double alphaCNCalibration_Y;
    double betaCNCalibration;
    double betaCNCalibration_X;
    double betaCNCalibration_Y;
};

/**
 * @brief  The engine for copy number calculations post log2 ratio calculation.
 * This engine works on one experiment at a time.
 *
 */
class CNAnalysisEngine : public BaseEngine
{
private:
    BaseEngine* m_pEngine;

    AffxArray<CNAnalysisMethod> m_vCNAnalysisMethods;
    AffxArray<CNReporter> m_vCNReporters;
    bool m_bComponentRun;
    void defineOptions();
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
        virtual std::string getEngineName() { return CNAnalysisEngine::EngineName(); }
        static const std::string EngineName() { return "CNAnalysisEngine"; }
        /**
         * Constructor
         */
        CNAnalysisEngine();

        /**
         * Destructor
         */
        ~CNAnalysisEngine();

        /*! A class to register the engine. */
        class Reg : public EngineReg
        {
        public:
            /*! Constructor - register the engine. */
            Reg() : EngineReg(CNAnalysisEngine::EngineName())
            {
            }

            /*! Creates an object.
             * @return The object.
             */
            BaseEngine *MakeObject() { return new CNAnalysisEngine; }
        };

        /*! The one and only registration object. */
        static Reg reg;

        /*! Converts the type to the engine type.
         * @param chip The pointer to the base engine object.
         * @return The engine type or NULL if not compatible.
         */
        static CNAnalysisEngine * FromBase(BaseEngine *engine);

    void defineSharedOptions();
    void defineSharedStates();

    void createAnalysis();
    void setEngine(BaseEngine* p)
    {
        if (p == NULL) {p = this;}
        m_pEngine = p;
    }
    void process(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets, CNProbeArray* vProbes=NULL);
    void postProcessing(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets);
    void calculateCalibrationParameters();
    bool needToCalculateCalibrationParameters();
    CalibrationParams retrieveCalibrationParams();
};


/*
 * Auxiliary predicate for searches using equal_range on probe set names
 */
struct ProbeSetNameCompare : std::binary_function<const CNProbeSet*, const CNProbeSet*, bool>
{
    bool operator()(const CNProbeSet* lhs, const CNProbeSet* rhs) const {
        return lhs->getProbeSetName() < rhs->getProbeSetName();
    }
};



#endif


