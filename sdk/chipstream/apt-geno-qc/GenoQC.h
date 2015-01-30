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

#ifndef _GENOQC_H_
#define _GENOQC_H_

//
#include "chipstream/CelListener.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/QCAnalysisOptions.h"
#include "chipstream/QuantMethodRunReport.h"
#include "chipstream/KitAODb.h"
//
#include "util/BaseEngine.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

/**
 * @brief Object for generating chip qc statistics from a 
 * Encapsulates all the work of generating chip qc stats 
 * from multiple cel files with different statistical algorithms
 */
class GenoQC : public BaseEngine {
    public:
			
		virtual std::string getEngineName() { return GenoQC::EngineName(); }
		static const std::string EngineName() { return "GenoQC"; }

		/**
		* Constructor
		*/
		GenoQC(ChipLayout *chipLayout = NULL);

		/**
		* Destructor
		*/
		~GenoQC();

		/*! A class to register the summary engine. */
		class Reg : public EngineReg
		{
		public:
			/*! Constructor - register the summary engine. */
			Reg() : EngineReg(GenoQC::EngineName())
			{
			}

			/*! Creates an object.
			 * @return The object.
			 */
			BaseEngine *MakeObject() { return new GenoQC; }
		};

		/*! The one and only registration object. */
		static Reg reg;

		/*! Converts the type to the summary engine type.
		 * @param chip The pointer to the base engine object.
		 * @return The summary engine type or NULL if not compatible.
		 */
		static GenoQC * FromBase(BaseEngine *engine);


        /** 
         * @brief Compute the stats given the input options
         */
        void ComputeStats();

        /**
         * @brief Get CEL Listeners which have stat stuff
         */
        std::vector<CelListener *> getCelListeners() {
            return m_CelListeners;
        }

    private:
        void checkOptionsImp();
        void runImp();
        void defineOptions();
        void defineStates();
        void readQCAFile(std::string qcaFile);
        void readQCCFile(std::string qccFile);
  
        void clear();

        /// QCA File
        QCAnalysisOptions qcAnalysisOpts;
        /// QCC File
        QCProbesetOptions qcProbesetOpts;
        /**
         * @brief Check chip type on all cel files, all cel files 
         *        are assumed to have the same chip type
         */
        void CheckChipType();

        /// Vector of CelListeners which compute/hold the various stats
        std::vector<CelListener *> m_CelListeners;
        /// Reporter
        QuantMethodRunReport *m_Report;
        /// Chip Layout Object
        ChipLayout *m_ChipLayout;
        bool m_FreeChipLayout;

  //
  std::string m_kitao_filename;
  KitAODb* m_kitao_db;
};

#endif // _GENOQC_H_
