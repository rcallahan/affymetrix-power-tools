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
* @file   GenotypeExtendedReportEngine.h
* @author Luis Jevons
* @date   Wed Map 13 16:01:09 2009
*/
#ifndef _GenotypeExtendedReportEngine_H_
#define _GenotypeExtendedReportEngine_H_

#include "util/BaseEngine.h"
//
#include <cstring>
#include <fstream>
#include <string>
//

class GenotypeExtendedReportEngine : public BaseEngine {

public:

	virtual std::string getEngineName() { return GenotypeExtendedReportEngine::EngineName(); }
	static const std::string EngineName() { return "GenotypeExtendedReportEngine"; }

	/**
	* Constructor
	*/
	GenotypeExtendedReportEngine();

	/**
	* Destructor
	*/
	~GenotypeExtendedReportEngine();

	/*! A class to register the summary engine. */
	class Reg : public EngineReg
	{
	public:
		/*! Constructor - register the summary engine. */
		Reg() : EngineReg(GenotypeExtendedReportEngine::EngineName())
		{
		}

		/*! Creates an object.
		 * @return The object.
		 */
		BaseEngine *MakeObject() { return new GenotypeExtendedReportEngine; }
	};

	/*! The one and only registration object. */
	static Reg reg;

	/*! Converts the type to the summary engine type.
	 * @param chip The pointer to the base engine object.
	 * @return The summary engine type or NULL if not compatible.
	 */
	static GenotypeExtendedReportEngine * FromBase(BaseEngine *engine);

private:

	virtual void defineOptions();
	virtual void checkOptionsImp();
	virtual void defineStates();
	virtual void runImp();
	void writeHeader(std::ofstream &out);
	void processSingleChannelFiles();
	void processMultiChannelFiles();
	void readMixFile();
	void readAnnotationFile();
private:
	std::map<std::string, bool> controlMap;
	bool annotationFileHasOneCode;
};

#endif

