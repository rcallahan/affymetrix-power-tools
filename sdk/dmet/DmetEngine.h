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
 * @file   DmetEngine.h
 * @author Alan Williams
 * @date   Mon Jun 23 14:57:34 PDT 2008
 * 
 * @brief Analysis engine for DMET 3.0
 */
#ifndef _DMETENGINE_H_
#define _DMETENGINE_H_

#include "chipstream/AnalysisInfo.h"
#include "util/BaseEngine.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @brief Class to run the DMET analysis
 */
class DmetEngine : public BaseEngine {

    public:


		virtual std::string getEngineName() { return DmetEngine::EngineName(); }
		static const std::string EngineName() { return "DmetEngine"; }

		/**
		* Constructor
		*/
		DmetEngine();

		/**
		* Destructor
		*/
		~DmetEngine();

		/*! A class to register the summary engine. */
		class Reg : public EngineReg
		{
		public:
			/*! Constructor - register the summary engine. */
			Reg() : EngineReg(DmetEngine::EngineName())
			{
			}

			/*! Creates an object.
			 * @return The object.
			 */
			BaseEngine *MakeObject() { return new DmetEngine; }
		};

		/*! The one and only registration object. */
		static Reg reg;

		/*! Converts the type to the summary engine type.
		 * @param chip The pointer to the base engine object.
		 * @return The summary engine type or NULL if not compatible.
		 */
		static DmetEngine * FromBase(BaseEngine *engine);


        virtual int parseArgv( const char * const * const argv, int start = 1 );

    private:
        /** 
         * Do the analysis
         * Will call Verbose::out() and Err::errAbort().
         */
        void runImp();

        std::string makeVal(std::string name, int optionIndex = 0);
        void defineOptions();
        void defineStates();
        void checkOptionsImp();
        bool computeCnSummaries();
        void computeCnState();
        void computeGenotypes();
        void generateChpFiles();
        void fillInAnalysisInfo(AnalysisInfo &info);
        void writeProbesetList(std::string filename, std::vector<std::string> &probesets);
        bool generateApsProbesetList(std::string requestedProbesets, std::string probesetsInRegion, std::string output, bool allowProbesetConsent = false);
        int m_ArgvPosAPS;
        int m_ArgvPosCN;
        int m_ArgvPosAPG;
        const char * const * m_argv;
};

#endif /* _DMETENGINE_H_ */
