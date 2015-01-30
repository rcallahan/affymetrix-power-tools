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
* @file   ListSummaryEngine.h
* @author Luis Jevons
* @date   Wed April 17 16:01:09 2009
* 
* @brief Core routines for list-summary binaries. By separating the
* command line parsing form the computation we allow a GUI application to share
* the core computation once setting up the option class. 
*/
#ifndef _LISTSUMMARYENGINE_H_
#define _LISTSUMMARYENGINE_H_

#include "util/BaseEngine.h"
//
#include <cstring>
#include <fstream>
#include <list>
#include <map>
#include <string>
//

class ListSummaryEngine : public BaseEngine {

public:

	virtual std::string getEngineName() { return ListSummaryEngine::EngineName(); }
	static const std::string EngineName() { return "ListSummaryEngine"; }

	/**
	* Constructor
	*/
	ListSummaryEngine();

	/**
	* Destructor
	*/
	~ListSummaryEngine();

	/*! A class to register the summary engine. */
	class Reg : public EngineReg
	{
	public:
		/*! Constructor - register the summary engine. */
		Reg() : EngineReg(ListSummaryEngine::EngineName())
		{
		}

		/*! Creates an object.
		 * @return The object.
		 */
		BaseEngine *MakeObject() { return new ListSummaryEngine; }
	};

	/*! The one and only registration object. */
	static Reg reg;

	/*! Converts the type to the summary engine type.
	 * @param chip The pointer to the base engine object.
	 * @return The summary engine type or NULL if not compatible.
	 */
	static ListSummaryEngine * FromBase(BaseEngine *engine);

private:
	std::map<std::string, std::pair<std::string, int> > mixMap;
	std::map<std::string, std::pair<std::string, std::string> > annotMap;
	std::list<std::string> markerOrder;
	bool annotationFileHasOneCode;
	//std::map<std::string, std::pair<uint64_t,int> > computeMedianData(const std::string& listFile, float &unknownCodeRate);
	std::map<std::string, std::pair<std::vector<uint64_t>, int> > computeMedianData(const std::string& listFile, float &unknownCodeRate);
	void defineOptions();
	void checkOptionsImp();
	void defineStates();
	void runImp();
	void readMixFile();
	void readAnnotationFile();
	void writeHeader(std::ofstream &out, const std::string &signalColumn, std::list<std::pair<std::string, float> > *unknownCodeRates = NULL);
};

#endif

