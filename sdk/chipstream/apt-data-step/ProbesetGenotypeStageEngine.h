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
 * @file   ProbesetGenotypeStageEngine.h
 * @author Chuck Sugnet
 * @date   Wed Mar 22 16:01:09 2006
 * 
 * @brief Core routines for probeset-genotype binaries. By separating the
 * command line parsing form the computation we allow a GUI application to share
 * the core computation once setting up the option class. To run the algorithm
 * setup the Options class appropriately for a particular run (i.e. specify
 * files and parameters) then call probesetGenotype() to perform analysis.
 */

#ifndef _PROBESETGENOTYPEENGINE_H_
#define _PROBESETGENOTYPEENGINE_H_

//
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/AptTypes.h"
#include "chipstream/CelStatListener.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/GenderCalls.h"
#include "chipstream/GenoSeed.h"
#include "chipstream/InbredStatus.h"
#include "chipstream/MetaProbeset.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/QuantBRLMM.h"
#include "chipstream/QuantLabelZ.h"
#include "chipstream/SpecialSnps.h"
#include "util/BaseEngine.h"
//
#include <cstring>
#include <ctime>
#include <ostream>
#include <string>
#include <vector>

class ProbesetGenotypeStageEngine : public BaseEngine {

    public:

        virtual std::string getEngineName() { return ProbesetGenotypeStageEngine::EngineName(); }
        static const std::string EngineName() { return "ProbesetGenotypeStageEngine"; }
        /**
         * Constructor
         */
        ProbesetGenotypeStageEngine(
                ChipLayout *chipLayout = NULL, 
                std::vector<const char *> *probesetNames = NULL);

        /**
         * Destructor
         */
        ~ProbesetGenotypeStageEngine();

        /*! A class to register the engine. */
        class Reg : public EngineReg {
            public:
                /*! Constructor - register the engine. */
                Reg() : EngineReg(ProbesetGenotypeStageEngine::EngineName()) { }
        
                /*! Creates an object.
                 * @return The object.
                 */
                BaseEngine *MakeObject() { return new ProbesetGenotypeStageEngine; }
        };
    
        /*! The one and only registration object. */
        static Reg reg;
        
        /*! Converts the type to the engine type.
         * @param chip The pointer to the base engine object.
         * @return The engine type or NULL if not compatible.
         */
        static ProbesetGenotypeStageEngine * FromBase(BaseEngine *engine);
        
        /**
         * Compare needed disk space to available disk space
         */


    private:

        void defineOptions();
        void defineStates();
        void defineStdMethods();
        void checkOptionsImp();
        void checkDiskSpaceImp();
        void checkBirdseedTsvPriorFile(const std::string &path, const std::string &specialSnpType);
        void checkSnpLists(const std::map<std::string,bool> &haploidSnps, ChipLayout *layout);
        int DetectNeedForSeeds(std::vector<AnalysisStream *> analysis);
        void runImp();
        void reportBasics();
        void extraHelp();
        void explain();
        void printStandardMethods(std::ostream &out);

        void loadPSFile(const std::string &fileName,
                   std::vector<const char *> *vNames,
                   std::set<const char *, Util::ltstr> *sNames);
        void loadLayout(ChipLayout &layout,
                   std::vector<const char *> &probesetNames,
                   probeidmap_t &killList,
                   bool justStats,
                   std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad);
        void fillInSpecialSnps(SpecialSnpMap &SpecialSnps,
                   ChipLayout &layout, bool permissive);
        void setPrior(std::string &outFile, std::vector<std::string> &sample,
                   ChipLayout &layout, IntensityMart &iMart,
                   std::vector<ChipStream *> &chipStream, PmAdjuster &pmAdjuster,
                   QuantBRLMM *qBRLMM);
        void setNormLabelZ(std::string &outFile, std::vector<std::string> sample,
                   ChipLayout &layout, IntensityMart &iMart,
                   std::vector<ChipStream *> &chipStream, PmAdjuster &pmAdjuster,
                           QuantLabelZ *qLabelZ);
        AnalysisInfo makeAnalysisInfo(ChipLayout &layout,
                                      const std::vector<const char*> &probesetNames,
                                      GenderCalls *genders,
                                      const std::string &analysisName);

        AnalysisInfo makeAnalysisInfo(AnalysisStream *as,
                   ChipLayout &layout,
                   const std::vector<const char*> &probesetNames,
                   GenderCalls *genders);
        void createAnalysisFromSpec(
                   std::vector<AnalysisStream *> &analysis,
                   AnalysisStreamFactory &asFactory,
                   ChipLayout &layout,
                   std::vector<const char*> &probesetNames,
                   unsigned int numProbeSets);
        void setInitialCalls(ChipLayout &layout, GenoSeed &seeds,
                   std::map<const char *,std::vector<affx::GType>, Util::ltstr > *knownGenotypes,
                   std::map<std::string,bool> &haploidSnps,
		   std::vector<AnalysisStream *> &analysis,
                   std::vector<const char *> toRunProbesets);
        void makeToSampleVector(std::vector<const char *> &toSample, 
                    map<std::string, int> &sampleMap,
                   ChipLayout *layout, 
                   const std::vector<const char *> &probesetNames,
                   const std::vector<const char *> &toRunProbesetNames,
                   std::set<affxcdf::GeneChipProbeSetType> &psTypesToSample);
        void ObtainSnpSample(
                   ChipLayout *layout,
                   std::vector<const char *> &selectFrom,
                   std::vector<std::string> &SnpSampleNames);
        void sortProbesetsByMapOrder(
                   vector<std::string>::iterator begin, 
                   vector<std::string>::iterator end, 
                   map<std::string, int> &mapOrder);
        void addReporters_Expression(affx::TsvReport::TsvReportFmt_t report_format,
                   AnalysisInfo& info,
                   QuantGTypeMethod *qMethod,
                   AnalysisStream *as,
                   std::vector<ChipSummary *> &chipSummaries,
                   GenderCalls *gender,
                   std::set<const char *, Util::ltstr> *probeSetsToReport,
                                     bool doCompact=false);
        void addReporters_Expression(affx::TsvReport::TsvReportFmt_t report_format,
                                     AnalysisInfo& info,
                                     std::vector<QuantMethodReport *> &exprReporters,
                                     vector<ChipSummary *> &chipSummaries,
                                     GenderCalls *gender,
                                     std::set<const char *, Util::ltstr> *probeSetsToReport,
                                     bool doCompact);

        void addReporters(
                   QuantGTypeMethod *qMethod,
                   AnalysisStream *as,
                   std::vector<ChipSummary *> &chipSummaries,
                   GenderCalls *gender,
                   std::set<const char *, Util::ltstr> *probeSetsToReport,
                   SpecialSnpMap &SpecialSnps);
        
        void makeReporters(
                           const std::string &qMethodName,
                           std::vector<QuantMethodReport *> &reporters,
                           std::vector<QuantMethodReport *> &exprReporters,
                           AnalysisInfo &info,
                           vector<ChipSummary *> &chipSummaries,
                           GenderCalls *gender,
                           std::set<const char *, Util::ltstr> *probeSetsToReport,
                           SpecialSnpMap &SpecialSnps);
        void DoSetup(std::vector<AnalysisStream *> &analysis,
                   ChipLayout &localLayout,
                   IntensityMart &iMart,
                   std::vector<std::string> &priorSnpNames,
                   std::vector<std::string> &normSnpNames,
                   std::vector<std::string> &fileNames,
                   GenderCalls *genders,
        InbredStatus *inbred,
                   std::map<std::string,bool> &haploidSnps,
                   SpecialSnpMap &SpecialSnps,
                   CopyNumberMap &copyNumberMap,
                   std::vector<ChipSummary *> &chipSummaries,
                   const std::vector<const char *> &probesetNames,
                   std::set<const char *, Util::ltstr> *probeSetsToReport);
        void setIterationData(ChipLayout &localLayout, 
                   std::vector<AnalysisStream *> &analysis,
                   std::vector<std::string> &priorSnpNames, std::vector<std::string> &normSnpNames,
                   std::vector<const char *> &toRunProbesets,
                   SpecialSnpMap &specialSnps);
        void fillInAnalysisInfo(AnalysisInfo &info, QuantMethod *qMethod);
        void fillInAnalysisInfo(AnalysisInfo &info, std::string prefix);

        void fillInAnalysisInfo(
                   AnalysisInfo &info, 
                   AnalysisStream *as, 
                   std::string prefix = "apt-");
        void closeGlobalA5();

    private:
        void fillInDesiredOrder(std::vector<int> &desiredOrder);

        void fillInSelectedGenders(GenderCalls **selectedGenders, 
                                   std::string &analysisName, 
                                   std::map<std::string, GenderCalls *> &genderCallsMap,
                                   NoGenderCalls *noGenderCalls);

        void fillInChipLayout();

        /// List of aliases for standardized methods
        map<string,string> m_stdMethods; 
        /// Global A5 Input File
        affx::File5_File* m_a5_global_input_file;
        affx::File5_Group* m_a5_global_input_group;
        /// Global A5 Output File
        affx::File5_File* m_a5_global_output_file;
        affx::File5_Group* m_a5_global_output_group;
        /// counts
        int m_celfile_count;
        int m_probeset_count;
        int m_channel_count;
        /// Chip Layout Object
        ChipLayout *m_ChipLayout;
        bool m_FreeChipLayout;
        std::vector<const char *> m_ProbesetNames;  /// names of all probesets in cdf order.
};

#endif /* _PROBESETGENOTYPEENGINE_H_ */

