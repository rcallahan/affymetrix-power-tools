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
 * @file   MAS5Engine.h
 * @author Alan Williams
 * @date   Wed Apr 15 15:29:23 PDT 2009
 * 
 * @brief Core routines for running MAS5
 */
#ifndef _MAS5ENGINE_H_
#define _MAS5ENGINE_H_

#include "mas5-stat/workflow/MAS5Workflow.h"
//
#include "util/BaseEngine.h"
//
#include <cstring>
#include <string>
#include <vector>
//

class MAS5Engine : public BaseEngine {

    public:

		virtual std::string getEngineName() { return MAS5Engine::EngineName(); }
		static const std::string EngineName() { return "MAS5Engine"; }

        /**
         * Constructor
         */
        MAS5Engine();

        /**
         * Destructor
         */
        ~MAS5Engine();

		/*! A class to register the summary engine. */
		class Reg : public EngineReg
		{
		public:
			/*! Constructor - register the summary engine. */
			Reg() : EngineReg(MAS5Engine::EngineName())
			{
			}

			/*! Creates an object.
			 * @return The object.
			 */
			BaseEngine *MakeObject() { return new MAS5Engine; }
		};

		/*! The one and only registration object. */
		static Reg reg;

		/*! Converts the type to the summary engine type.
		 * @param chip The pointer to the base engine object.
		 * @return The summary engine type or NULL if not compatible.
		 */
		static MAS5Engine * FromBase(BaseEngine *engine);

    private:

        /** 
         * Do the analysis
         * Will call Verbose::out() and Err::errAbort().
         */
        void runImp();

        void defineOptions();
        void defineStates();
        void checkOptionsImp();

    private:
	    MAS5Workflow m_Mas5;
};

#endif /* _MAS5ENGINE_H_ */

