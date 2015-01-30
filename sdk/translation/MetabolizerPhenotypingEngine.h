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

#ifndef _MetabolizerPhenotyping_HEADER_
#define _MetabolizerPhenotyping_HEADER_

#include <string>
#include <map>
#include <vector>
#include <list>
#include "util/BaseEngine.h"

namespace affymetrix
{

/*! The phenotype and associated activities. */
typedef struct _PhenotypeCallType
{
public:
	std::string phenotype;
	std::string activity1;
	std::string activity2;

	_PhenotypeCallType() { };

	_PhenotypeCallType(const std::string &p, const std::string &a1, const std::string &a2)
	{
		phenotype = p;
		activity1 = a1;
		activity2 = a2;
	};
} PhenotypeCallType;

/*! A map of allele pairs to phenotype */
typedef std::map<std::string /*allele pair*/, PhenotypeCallType /*phenotype*/> GeneMetabolizerPhenotypeMap;

/*! A map of phenotypes for all genes */
typedef std::map<std::string /*gene*/, GeneMetabolizerPhenotypeMap> MetabolizerPhenotypeMap;

/*! Holds the results of the metabolizer analysis for a single gene. */
typedef struct _MetabolizerResultType
{	
	std::string index;
	std::string chpFile;
	std::string gene;
	std::string knownCall;
	std::string unknownCall;
	std::string code;
	std::vector<std::string> phenotypes;
	std::vector<std::string> activities;
	std::string finalPhenotype;
	std::string finalActivity;
} MetabolizerResultType;

/*! A list of results */
typedef std::list<MetabolizerResultType> MetabolizerResultsType;

/*! Implements a method to calculate metabolizer phenotypes */
class MetabolizerPhenotypingEngine : public BaseEngine
{
private:

	/*! The map of alleles to phenotype */
	MetabolizerPhenotypeMap metabolizerPhenotypeMap;

	/*! The results of the metabolizer analysis */
	MetabolizerResultsType metabolizerResults;

	/*! The description of the phenotype calls. */
	std::vector<std::string> phenotypeCallDescs;

	/*! Information in the metabolizer file to pass along to the user in the output translation file. */
	std::vector<std::string> userInformation;

	/*! Read the data from the metabolizer bin file.
	 * @param fileName The name of the bin file.
	 */
	void ReadMetabolizerBinFile();

	/*! Extract the contents of the allele translation table file and compute
	 * the metabolizer phenotypes for each gene of each CHP file.
	 */
	void ProcessAlleleTranslationFile();

	/*! Create MD5's of the input and output files. */
	void CreateMd5File();

	/*! Initialize the members. */
	void Init();

	/*! Define all options. */
	void defineOptions();

	/*! Run the workflow to compute the phenotypes */
	void runImp();

	/*! Check the options. */
	void checkOptionsImp();

public:
	
	/*! Get the engine name */
	virtual std::string getEngineName() { return MetabolizerPhenotypingEngine::EngineName(); }

	/*! Get the engine name */
	static const std::string EngineName() { return "MetabolizerPhenotypingEngine"; }

	/*! A class to register the summary engine. */
	class Reg : public EngineReg
	{
	public:
		/*! Constructor - register the summary engine. */
		Reg() : EngineReg(MetabolizerPhenotypingEngine::EngineName())
		{
		}

		/*! Creates an object.
		 * @return The object.
		 */
		BaseEngine *MakeObject() { return new MetabolizerPhenotypingEngine; }
	};

	/*! The one and only registration object. */
	static Reg reg;

	/*! Converts the type to the summary engine type.
	 * @param chip The pointer to the base engine object.
	 * @return The summary engine type or NULL if not compatible.
	 */
	static MetabolizerPhenotypingEngine * FromBase(BaseEngine *engine);

	/*! Constructor */
	MetabolizerPhenotypingEngine();

	/*! Destructor */
	~MetabolizerPhenotypingEngine();
	
	/*! The results of the metabolizer analysis */
	const MetabolizerResultsType& MetabolizerResults() const { return metabolizerResults; }
	
	/*! The description of the phenotype calls. */
	const std::vector<std::string>& PhenotypeCallDescs() const { return phenotypeCallDescs; }
	
	/*! Information in the metabolizer file to pass along to the user in the output translation file. */
	const std::vector<std::string>& UserInformation() const { return userInformation; }
};

}

#endif
