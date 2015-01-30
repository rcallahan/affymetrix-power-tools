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

#ifndef _SNPSUMMARYENGINE_H_
#define _SNPSUMMARYENGINE_H_

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
class SnpSummaryEngine : public BaseEngine {
    public:
			
		virtual std::string getEngineName() { return SnpSummaryEngine::EngineName(); }
		static const std::string EngineName() { return "SnpSummaryEngine"; }

		/**
		* Constructor
		*/
		SnpSummaryEngine();
		

		/**
		* Destructor
		*/
		~SnpSummaryEngine();
		

		/*! A class to register the summary engine. */
		class Reg : public EngineReg
		{
		public:
			/*! Constructor - register the summary engine. */
			Reg() : EngineReg(SnpSummaryEngine::EngineName())
			{
			}

			/*! Creates an object.
			 * @return The object.
			 */
			BaseEngine *MakeObject() { return new SnpSummaryEngine; }
		};

		/*! The one and only registration object. */
		static Reg reg;

		/*! Converts the type to the summary engine type.
		 * @param chip The pointer to the base engine object.
		 * @return The summary engine type or NULL if not compatible.
		 */
		static SnpSummaryEngine * FromBase(BaseEngine *engine);
		std::vector<std::string> chpFiles;


        /** 
         * @brief Compute the stats given the input options
         */
       //$$ void ComputeSnpSummary();

       

    private:
        void checkOptionsImp();
        void runImp();
        void defineOptions();
        void defineStates();
		void CheckChipType();
        void clear(); 
    
};

#endif // _SNPSUMMARY_H_
