////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
* @file   OnePointNormalizationEngine.h
* @author Luis Jevons
*/
#ifndef _OnePointNormalizationEngine_H_
#define _OnePointNormalizationEngine_H_

#include "util/BaseEngine.h"
//
#include <cstring>
#include <fstream>
#include <list>
#include <map>
#include <string>
//

class OnePointNormalizationEngine : public BaseEngine {

public:

	virtual std::string getEngineName() { return OnePointNormalizationEngine::EngineName(); }
	static const std::string EngineName() { return "OnePointNormalizationEngine"; }

	/**
	* Constructor
	*/
	OnePointNormalizationEngine();

	/**
	* Destructor
	*/
	~OnePointNormalizationEngine();

	/*! A class to register the summary engine. */
	class Reg : public EngineReg
	{
	public:
		/*! Constructor - register the summary engine. */
		Reg() : EngineReg(OnePointNormalizationEngine::EngineName())
		{
		}

		/*! Creates an object.
		 * @return The object.
		 */
		BaseEngine *MakeObject() { return new OnePointNormalizationEngine; }
	};

	/*! The one and only registration object. */
	static Reg reg;

	/*! Converts the type to the summary engine type.
	 * @param chip The pointer to the base engine object.
	 * @return The summary engine type or NULL if not compatible.
	 */
	static OnePointNormalizationEngine * FromBase(BaseEngine *engine);

private:
	void defineOptions();
	void checkOptionsImp();
	void defineStates();
	void runImp();
	void writeHeader(std::ofstream &out);
};

#endif

