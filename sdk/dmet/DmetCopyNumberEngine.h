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
 * @file   DmetCopyNumberEngine.h
 * @author Chuck Sugnet
 * @date   Wed Mar 22 16:01:09 2006
 * 
 * @brief Copy Number State Calling for DMET3
 */
#ifndef _COPYNUMBERENGINEDMET_H_
#define _COPYNUMBERENGINEDMET_H_

#include "dmet/DmetCopyNumberData.h"
//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/TsvReport.h"
#include "util/BaseEngine.h"
//
#include <cstring>
#include <set>
#include <string>
//

class DmetCopyNumberEngine : public BaseEngine {

    public:

		virtual std::string getEngineName() { return DmetCopyNumberEngine::EngineName(); }
		static const std::string EngineName() { return "DmetCopyNumberEngine"; }

		/**
		* Constructor
		*/
		DmetCopyNumberEngine();

		/**
		* Destructor
		*/
		~DmetCopyNumberEngine();

		/*! A class to register the summary engine. */
		class Reg : public EngineReg
		{
		public:
			/*! Constructor - register the summary engine. */
			Reg() : EngineReg(DmetCopyNumberEngine::EngineName())
			{
			}

			/*! Creates an object.
			 * @return The object.
			 */
			BaseEngine *MakeObject() { return new DmetCopyNumberEngine; }
		};

		/*! The one and only registration object. */
		static Reg reg;

		/*! Converts the type to the summary engine type.
		 * @param chip The pointer to the base engine object.
		 * @return The summary engine type or NULL if not compatible.
		 */
		static DmetCopyNumberEngine * FromBase(BaseEngine *engine);


    private:

        void defineOptions();
        void defineStates();
        void defineStdMethods();
        void checkOptionsImp();

        /** 
         * Do the analysis
         * Will call Verbose::out() and Err::errAbort().
         */
        void runImp();


        double minDouble(double a, double b);
        double maxDouble(double a, double b);
        double Log2(double input); 
        double tDensity(double ww, double degFreedom);
        void parseSummaryDataFile(std::string summaryDataFileName, DmetCopyNumberData &copyNumberData);
        void parseSummaryDataFile5(std::string summaryDataFileName, DmetCopyNumberData &copyNumberData);
        void parseRegionDataFile(std::string regionDataFileName, DmetCopyNumberData &copyNumberData);
        void parseProbeSetDataFile(std::string probeSetDataFileName, DmetCopyNumberData &copyNumberData);
	void parseCNRegionGenotypeProbesetFile(std::string cn_region_gt_marker_filename, DmetCopyNumberData &copyNumberData);
        double calculate01CopyNumberEstimate(	std::string region, 
					std::set<std::string> probeSetNames,
					std::string sampleName,
					const DmetCopyNumberData &copyNumberData); 
        double calculate12CopyNumberEstimate(	std::string region, 
					std::set<std::string> probeSetNames,
					std::string sampleName,
					const DmetCopyNumberData &copyNumberData); 
        double calculateCNDensity(
                double inputValue,
				double mean,
				double stDeviation,
				double degFreedom);
        int findMaxIndex(const std::vector<double> inputVector);
        int tailQuantileLower(const std::vector<double> inputVector, double pValue);
        int tailQuantileUpper(const std::vector<double> inputVector, double pValue);
        void configureReport(affx::TsvReport *tsv, affx::TsvReport::TsvReportFmt_t format, AnalysisInfo &info, const std::string &guid);
};

#endif /* _COPYNUMBERENGINEDMET_H_ */

