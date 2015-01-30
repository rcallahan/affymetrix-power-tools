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
 * @file   CanaryEngine.h
 * @author HZ
 * @date   Wed Mar 22 16:01:09 2006
 * 
 * @brief Give a brief description of Canary
 */
#ifndef _CANARYENGINE_H_
#define _CANARYENGINE_H_

//
#include "canary/CanaryOptions.h"
//
#include "chipstream/AnalysisStream.h"
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/CelReader.h"
#include "file/TsvFile/TsvFile.h"
#include "portability/apt-win-dll.h"
#include "util/BaseEngine.h"
//
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//

class APTLIB_API CanaryEngine : public BaseEngine {

    public:

        virtual std::string getEngineName() { return CanaryEngine::EngineName(); }
        static const std::string EngineName() { return "CanaryEngine"; }

        /**
         * Constructor
         */
        CanaryEngine();
        
        /**
         * Destructor
         */
        ~CanaryEngine();
        
        /*! A class to register the summary engine. */
        class Reg : public EngineReg
        {
            public:
            /*! Constructor - register the summary engine. */
            Reg() : EngineReg(CanaryEngine::EngineName())
            {
            }
            
            /*! Creates an object.
             * @return The object.
             */
            BaseEngine *MakeObject() { return new CanaryEngine; }
        };
        
        /*! The one and only registration object. */
        static Reg reg;
        
        /*! Converts the type to the summary engine type.
         * @param chip The pointer to the base engine object.
         * @return The summary engine type or NULL if not compatible.
         */
        static CanaryEngine * FromBase(BaseEngine *engine);


    private:

        void extraHelp();
        void defineOptions();
        void defineStates();
        /** 
         * Make sure that our options are sane.
         * 
         * @param opts - Options to be checked.
         */
        void checkOptionsImp();

        void checkDiskSpaceImp();

        void runImp();

        void fillInOptions(CanaryOptions &o);

        /** 
         * This is where all of the computation begins. To kick off a run set up the
         * options desired and pass the options object to this function.
         * 
         * @param opts - Options specifying files, parameters, etc.
         */
        void CanaryCompute(CanaryOptions &opts);

        /*
         *  Take the output from pre-canary steps and average the A and B alleles.
         */
        void averageAlleles(CanaryOptions &opts);

        /*
         *  Compute region summaries
         */
        void writeCnvSummaries(CanaryOptions &opts, map<string,vector<string> > &cnv_region_map);

        /*
         *  The actual Canary part
         */
        void mainCanary(CanaryOptions &opts, map<string,vector<string> > &cnv_region_map, 
                vector<string> &cnv_region_names);

        set<string> get_cnv_probes(string ifname);
        set<string> get_cnv_normalization_probes(string ifname);
        string get_genome_version(string cnvRegionFile);
        void read_region_file(string cnvRegionFile, map<string,vector<string> > &cnv_region_map, vector<string> &cnv_region_names) ;
        void check_input_file_headers(vector<string> chipTypes, 
                vector<string> files, string &mapName, string &mapVersion, bool check = true);

        vector<string> get_available_psets(set<string> & cnv_probes,
                set<string> & normalization_probes, set<string> & all_probes,
                CanaryOptions & opts);

        void loadCdfLayout_canary(ChipLayout & layout, CanaryOptions & opts, set<string> & libpsets);
        set<string> get_cnv_psets(map<string,vector<string> > &cnv_region_map);
        set<string> get_normalization_psets(string ifname);
        void get_meta_tags(string file, vector<string> &keys, vector<string> &vals);
};

#endif
