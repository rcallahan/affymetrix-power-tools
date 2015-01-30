////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

/*
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The PLIER (Probe Logarithmic Error Intensity Estimate) method produces
an improved signal by accounting for experimentally observed patterns 
in probe behavior and handling error at the appropriately at low and 
high signal values.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/

/*
 * test.cpp : defines the entry point for the plier interfaces test program
 */

//
#include "plier/error.h"
#include "plier/iaffyplier.h"
//
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string.h>
#include <string>
#include <vector>
//

#ifdef _MSC_VER
#pragma warning(disable: 4996) // don't show deprecated warnings.
#endif

using namespace std;

#define BUF_LEN	1024

/*
 * a plier wrapper function to use the iaffyplier interfaces
 */

void runPlierWrapper(
	iaffyplier* plier,
	long numExp,
	long numProbe,
	long* replicate,
	double** pm,
	double** mm,
	double* targetResponse,
	double* featureResponse,
	long* errorCode)
{
	*errorCode = 0;
	if (plier)
	{
		plier->setNumExp(numExp);
		plier->setNumFeature(numProbe);
		plier->setReplicate(replicate);
		plier->setPM(pm);
		plier->setMM(mm);
		plier->setTargetResponse(targetResponse);
		plier->setFeatureResponse(featureResponse);
		plier->run(errorCode);
	}
}

/*
 * a plier algorithm testing class
 */

class plierAlg
{
public:
	plierAlg();
	~plierAlg();

	void readFile(const char* file, bool& status);
	void run();

protected:

	long numExp;
	long numProbes;
	long* replicate;
	double** pm;
	double** mm;

	bool useHgu133aDefault;
	double augmentation;
	double gmCutoff;
	double differentialFeaturePenalty;
	double differentialTargetPenalty;
	double defaultFeatureResponse;
	double defaultTargetResponse;
	double attenuation;
	double seaConvergence;
	long seaIteration;
	double plierConvergence;
	long plierIteration;
	bool useMM;
	bool useModel;
	bool fitFeatureResponse;
	double dropMax;
	double lambdaLimit;
	long optimization;
	double safetyZero;
	double numericalTolerance;
	bool fixPrecomputed;
	bool fixFeatureEffect;

	vector<vector<double> > qPM;
	vector<vector<double> > qMM;

	affy_ptr<iaffyplier> sp;

private:

	void interpretLine(const char* line);
	void readData(const char* line, const long dataLineCount);
	bool checkData();
	void prepareData();
	void clearData();
	void setParams();

};

/*
 * constructor
 */

plierAlg::plierAlg()
{
	numExp = 0;
	numProbes = 0;
	replicate = 0;
	pm = 0;
	mm = 0;
	useHgu133aDefault =  false;
}

/*
 * destructor
 */

plierAlg::~plierAlg()
{
	if (replicate)
		delete[] replicate;
	clearData();
}

/*
 * set the plier parameters
 * if use_hgu133a_default is true, call the set_to_hgu133a_default method
 * and then set the parameters one by one.
 */

void plierAlg::setParams()
{
	if (sp)
	{
		if (useHgu133aDefault)
			sp->setDefault();

		// one should really check that these values were
		// provided in the input file
		sp->setInitAugmentation(augmentation);
		sp->setInitDefaultFeatureResponse(defaultFeatureResponse);
		sp->setInitDefaultTargetResponse(defaultTargetResponse);

		sp->setSeaAttenuation((float) attenuation);
		sp->setSeaOptConvergence(seaConvergence);
		sp->setSeaOptIteration(seaIteration);

		sp->setPlierGmCutoff((float) gmCutoff);
		sp->setPlierDifferentialFeaturePenalty((float) differentialFeaturePenalty);
		sp->setPlierDifferentialTargetPenalty((float) differentialTargetPenalty);
		sp->setPlierUseMMLikelihood(useMM);
		sp->setPlierUseInputModel(useModel);
		sp->setPlierFitFeatureResponse(fitFeatureResponse);
		
		sp->setPlierOptConvergence(plierConvergence);
		sp->setPlierOptIteration(plierIteration);
		sp->setPlierOptDropMax(dropMax);
		sp->setPlierOptLambdaLimit(lambdaLimit);
		sp->setPlierOptOptimizationMethod(optimization);

		sp->setFixPrecomputed(fixPrecomputed);
		sp->setNumericalTolerance(numericalTolerance);
		sp->setSafetyZero(safetyZero);
		sp->setFixFeatureEffect(false);

	}
}

/*
 * read the data inputs
 * each line contains either pm/mm probe intensities (separated by space or tab)
 * number of columns = number of probes
 * expected a line of MMs to go with the PMs
 */

void plierAlg::readData(const char* line, const long dataLineCount)
{
	char *token;
	char seps[] = " \t\n";
	char src[BUF_LEN];
	strcpy(src, line);
	token = strtok(src, seps);

	long curSize;
	if (dataLineCount % 2)
	{
		curSize = qMM.size();
		qMM.resize(curSize+1);
	}
	else
	{
		curSize = qPM.size();
		qPM.resize(curSize+1);
	}
	while (token != 0)
	{
		if (dataLineCount % 2)
			qMM[curSize].push_back(atof(token));
		else
			qPM[curSize].push_back(atof(token));
		token = strtok(0, seps);
	}
}

/*
 * check the input data is valid in their sizes
 */

bool plierAlg::checkData()
{
	if (qPM.size() != qMM.size())
		return false;

	if (qPM.size() == 0)
		return false;

	int i;
	int curSize = qPM[0].size();
	for (i=0; i<(int)qPM.size(); i++)
		if (qPM[i].size() != curSize || qMM[i].size() != curSize)
			return false;
	return true;
}

/*
 * prepare data for plier interfaces
 * convert from the vector to the pointer structures
 */

void plierAlg::prepareData()
{
	numExp = qPM.size();
	numProbes = qPM[0].size();
	pm = new double*[numExp];
	mm = new double*[numExp];
	long i,j;
	for (i=0; i<numExp; i++)
	{
		pm[i] = new double[numProbes];
		mm[i] = new double[numProbes];
		for (j=0; j<numProbes; j++)
		{
			pm[i][j] = qPM[i][j];
			mm[i][j] = qMM[i][j];
		}
	}
}

/*
 * clear data structures
 */

void plierAlg::clearData()
{
	long i;
	if (pm)
	{
		for (i=0; i<numExp; i++)
			delete[] pm[i];
		delete[] pm;
		pm = 0;
	}
	if (mm)
	{
		for (i=0; i<numExp; i++)
			delete[] mm[i];
		delete[] mm;
		mm = 0;
	}
	qPM.clear();
	qMM.clear();
}

/*
 * interpret the input line
 */

void plierAlg::interpretLine(const char* line)
{
	static bool startReadData = false;
	static long dataLineCount = 0;

	if (line)
	{
		if (line[0] == 0 || line[0] == ';' || line[0] == ' ' || line[0] == '\r' || line[0] == '\n') return;

		char command[BUF_LEN];
		char value[BUF_LEN];

		printf("%s", line);

		if (startReadData == true)
		{
			readData(line, dataLineCount);
			dataLineCount++;
			return;
		}

		const char* ptr = strchr(line, ':');
		if (ptr) 
		{
			int pos = ptr-line;
			strncpy(command, line, pos);
			command[pos] = 0;
			strcpy(value, ptr+2);
		}
		else
			return;

		if (strcmp(command, "exp") == 0)
			numExp = atoi(value);
		else if (strcmp(command, "probes") == 0)
			numProbes = atoi(value);
		else if (strcmp(command, "use_hgu133a_default") == 0)
			useHgu133aDefault = atoi(value) == 1? true: false;
		else if (strcmp(command, "InitAugmentation") == 0)
			augmentation = atof(value);
		else if (strcmp(command, "InitDefaultFeatureResponse") == 0)
			defaultFeatureResponse = atof(value);
		else if (strcmp(command, "InitDefaultTargetResponse") == 0)
			defaultTargetResponse = atof(value);
		else if (strcmp(command, "SeaAttenuation") == 0)
			attenuation = atof(value);
		else if (strcmp(command, "SeaOptConvergence") == 0)
			seaConvergence = atof(value);
		else if (strcmp(command, "SeaOptIteration") == 0)
			seaIteration = atoi(value);
		else if (strcmp(command, "PlierGmCutoff") == 0)
			gmCutoff = atof(value);
		else if (strcmp(command, "PlierDifferentialFeaturePenalty") == 0)
			differentialFeaturePenalty = atof(value);
		else if (strcmp(command, "PlierDifferentialTargetPenalty") == 0)
			differentialTargetPenalty = atof(value);
		else if (strcmp(command, "PlierUseMMLikelihood") == 0)
			useMM = atoi(value) == 1? true: false;
		else if (strcmp(command, "PlierUseInputModel") == 0)
			useModel = atoi(value) == 1? true: false;
		else if (strcmp(command, "PlierFitFeatureResponse") == 0)
			fitFeatureResponse = atoi(value) == 1? true: false;
		else if (strcmp(command, "PlierOptConvergence") == 0)
			plierConvergence = atof(value);
		else if (strcmp(command, "PlierOptIteration") == 0)
			plierIteration = atoi(value);
		else if (strcmp(command, "PlierOptDropMax") == 0)
			dropMax = atof(value);
		else if (strcmp(command, "PlierOptLambdaLimit") == 0)
			lambdaLimit = atof(value);
		else if (strcmp(command, "PlierOptOptimizationMethod") == 0)
			optimization = atoi(value);
		else if (strcmp(command, "data") == 0)
			startReadData = true;
		else if (strcmp(command, "SafetyZero") == 0)
			safetyZero = atof(value);
		else if (strcmp(command, "FixPrecomputed") == 0)
			fixPrecomputed = atof(value) == 1 ? true : false;
		else if (strcmp(command, "FixFeatureEffect") == 0)
			fixFeatureEffect = atof(value) == 1 ? true : false;
		else if (strcmp(command, "NumericalTolerance") == 0)
			numericalTolerance = atof(value);
		else
			printf("invalid input line %s\n", command);
	}
}

/*
 * read in the input file
 */

void plierAlg::readFile(const char* file, bool& status)
{
	status = true;
	char line[BUF_LEN];
	FILE* fd = fopen(file, "rt");
	if (fd)
	{
		while (fgets(line, BUF_LEN, fd))
			interpretLine(line);
		fclose(fd);
	}
	else
		status = false;
}

/*
 * run the plier algorithm
 */

void plierAlg::run()
{
	try
	{
		if (checkData())
		{
			prepareData();

			double* targetResponse=0;
			double* featureResponse=0;
			long errorCode=0;
			targetResponse = new double[numExp];
			featureResponse = new double[numProbes];

			createPlierObject(0, &sp);
			if (sp)
			{
				setParams();
				runPlierWrapper(sp, numExp, numProbes, replicate, pm, mm, targetResponse, featureResponse, &errorCode);
				if (errorCode != 0)
				{
					char errorStr[MAX_ERROR_LENGTH];
					getPlierError(errorCode, errorStr);
					printf("Error in running plier: %s\n", errorStr);
				}
				else
				{
					printf("\n============= PLIER results ===============\n");
					printf("TargetResponse:\n");
					long i;
					for (i=0; i<numExp; i++)
						printf("%f\t", targetResponse[i]);
					printf("\n\n");
					printf("FeatureResponse:\n");
					for (i=0; i<numProbes; i++)
						printf("%f\t", featureResponse[i]);
					printf("\n");
				}
			}
			delete[] targetResponse;
			delete[] featureResponse;
		}
		else
			printf("invalid input data\n");
	}
	catch(...)
	{
	}
}

/*
 * main program to get input file and run plier algorithm 
 */

int main(int argc, char* argv[])
{
	char filename[BUF_LEN];
	filename[0] = 0;
	if (argc >= 2)
		strcpy(filename, argv[1]);
	else
	{
		while (strcmp(filename, "") == 0)
		{
			printf("Please enter input file name: ");
			fgets(filename,BUF_LEN,stdin);
      // chop string at first '\n'
      char *lf=strchr(filename,'\n');
      if (lf!=NULL) {
        *lf=0;
      }
		}
	}

	plierAlg myplier;
	bool status;
	myplier.readFile(filename, status);
	if (status)
		myplier.run();
	else
		printf("error in reading file %s\n", filename);

	return 0;
}

