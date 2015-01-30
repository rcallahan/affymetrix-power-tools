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

#ifndef _ExportDataFilesEngine_HEADER_
#define _ExportDataFilesEngine_HEADER_

#include <string>
#include <vector>
#include "util/BaseEngine.h"

class ExportDataFilesEngine : public BaseEngine
{
public:
	
	/** The engine name. */
	virtual std::string getEngineName() { return ExportDataFilesEngine::EngineName(); }

	/** The engine name. */
	static const std::string EngineName() { return "ExportDataFilesEngine"; }

	/**
	* Constructor
	*/
	ExportDataFilesEngine();

	/**
	* Destructor
	*/
	~ExportDataFilesEngine();

	/*! A class to register the merge engine. */
	class Reg : public EngineReg
	{
	public:
		/*! Constructor - register the merge engine. */
		Reg() : EngineReg(ExportDataFilesEngine::EngineName())
		{
		}

		/*! Creates an object.
		 * @return The object.
		 */
		BaseEngine *MakeObject() { return new ExportDataFilesEngine; }
	};

	/*! The one and only registration object. */
	static Reg reg;

	/*! Converts the type to the summary engine type.
	 * @param chip The pointer to the base engine object.
	 * @return The summary engine type or NULL if not compatible.
	 */
	static ExportDataFilesEngine * FromBase(BaseEngine *engine);

public:
	void defineOptions();
	void checkOptionsImp();
	void defineStates();
	void runImp();
};

#endif
