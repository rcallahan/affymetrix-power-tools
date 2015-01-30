////////////////////////////////////////////////////////////////
//
// Copyright (C) 2013 Affymetrix, Inc.
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

#ifndef _CELTransformer_HEADER_
#define _CELTransformer_HEADER_

#include <string>
#include <map>
#include <vector>
#include <list>
#include "chipstream/ChipStream.h"
#include "util/BaseEngine.h"

namespace affx
{

/*! Program to generically transform cel file data using
 * chipstreams and write the data to a new cel file. Also serves as relatively simple
 * code sample to illustrate usage of chipstream code.
 */
class CELTransformerEngine : public BaseEngine
{
private:
	
	std::vector<std::string> celFiles;

	/*! Define all options. */
	void defineOptions();

	/*! Run the workflow to compute the phenotypes */
	void runImp();

	/*! Check the options. */
	void checkOptionsImp();

	void transformOld(int celIx, std::vector<ChipStream *> &cStreamVec);
	void transformNew(int celIx, std::vector<ChipStream *> &cStreamVec);

public:
	
	/*! Get the engine name */
	virtual std::string getEngineName() { return CELTransformerEngine::EngineName(); }

	/*! Get the engine name */
	static const std::string EngineName() { return "CELTransformerEngine"; }

	/*! A class to register the summary engine. */
	class Reg : public EngineReg
	{
	public:
		/*! Constructor - register the summary engine. */
		Reg() : EngineReg(CELTransformerEngine::EngineName())
		{
		}

		/*! Creates an object.
		 * @return The object.
		 */
		BaseEngine *MakeObject() { return new CELTransformerEngine; }
	};

	/*! The one and only registration object. */
	static Reg reg;

	/*! Converts the type to the summary engine type.
	 * @param chip The pointer to the base engine object.
	 * @return The summary engine type or NULL if not compatible.
	 */
	static CELTransformerEngine * FromBase(BaseEngine *engine);

	/*! Constructor */
	CELTransformerEngine();

	/*! Destructor */
	~CELTransformerEngine();
};

}

#endif
