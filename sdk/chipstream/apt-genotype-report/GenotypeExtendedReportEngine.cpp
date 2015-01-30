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

//
#include "chipstream/apt-genotype-report/GenotypeExtendedReportEngine.h"
//
#include "chipstream/EngineUtil.h"
#include "chipstream/GenoUtility.h"
//
#include "file/TsvFile/TsvFile.h"
#include "stats/statfun.h"
#include "util/PgOptions.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//

//
using namespace std;
using namespace affx;

GenotypeExtendedReportEngine::Reg GenotypeExtendedReportEngine::reg;

GenotypeExtendedReportEngine * GenotypeExtendedReportEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == GenotypeExtendedReportEngine::EngineName())
		return (GenotypeExtendedReportEngine *)engine;
	return NULL;
}

GenotypeExtendedReportEngine::GenotypeExtendedReportEngine() : annotationFileHasOneCode(true) {
	defineOptions();
}

GenotypeExtendedReportEngine::~GenotypeExtendedReportEngine() {
}

void GenotypeExtendedReportEngine::defineOptions() {
	defineOption("", "calls-file", PgOpt::STRING_OPT,
		"Text file specifying genotype calls.",
		"");
	defineOption("", "forced-calls-file", PgOpt::STRING_OPT,
		"Text file specifying forced genotype calls.",
		"");
	defineOption("", "count-file", PgOpt::STRING_OPT,
		"Text file specifying counts.",
		"");
	defineOption("", "cv-file", PgOpt::STRING_OPT,
		"Text file specifying cvs.",
		"");
	defineOption("", "confidence-file", PgOpt::STRING_OPT,
		"Text file specifying the genotype confidence values.",
		"");
	defineOption("m", "summaries", PgOpt::STRING_OPT,
		"Text file specifying the signal summaries.",
		"");
	defineOption("", "out-file", PgOpt::STRING_OPT,
		"The output file.",
		"");
	defineOption("", "output-sample-col", PgOpt::STRING_OPT,
		"The column name to use for the sample column.",
		"Sample");
	defineOption("", "output-marker-col", PgOpt::STRING_OPT,
		"The column name to use for the marker column.",
		"Marker");
	defineOption("", "marker-content-file", PgOpt::STRING_OPT,
		"Text file specifying the mix file.",
		"");
	defineOption("", "annotation-file", PgOpt::STRING_OPT,
		"The annotation file.",
		"");
	defineOption("", "note-col", PgOpt::STRING_OPT,
		"The column name for the note.",
		"Note");
	defineOption("", "probe-id-col", PgOpt::STRING_OPT,
		"The column name for the probe id.",
		"Probe Id");
	defineOption("", "type-col", PgOpt::STRING_OPT,
		"The column name for the type.",
		"Type");
	defineOption("", "summaries-norm", PgOpt::STRING_OPT,
		"The normalized signal summary file.",
		"");
	defineOption("", "stretch-constant", PgOpt::DOUBLE_OPT,
		"The stretch constant for the cluster transformation.",
		"1.5");
	defineOption("", "cluster-transformation", PgOpt::STRING_OPT,
		"The type of cluster transformation. Can be one of CES, CCS, MvA, RvT.",
		"CES");
}

void GenotypeExtendedReportEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void GenotypeExtendedReportEngine::checkOptionsImp()
{
    defineStates();

	string str = getOpt("out-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an output file name.");

	str = getOpt("confidence-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an input confidence file.");

	str = getOpt("cv-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an input cv file.");

	str = getOpt("count-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an input count file.");

	str = getOpt("summaries");
	if (str.empty() == true)
		Err::errAbort("Must specify an input summaries file.");

	str = getOpt("calls-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an input call file.");

	str = getOpt("marker-content-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an input marker file.");

	str = getOpt("annotation-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an annotation file.");

	str = getOpt("summaries-norm");
	if (str.empty() == true)
		Err::errAbort("Must specify an input normalized signal summary file.");
}

void GenotypeExtendedReportEngine::writeHeader(std::ofstream &out)
{
	std::vector<std::string> option_names;
	getOptionNames(option_names);
	std::vector<PgOpt::PgOptType> option_types;
	getOptionTypes(option_types);
	for (int i = 0; i < option_names.size(); i++) 
	{
		switch (option_types[i]) 
		{
		case PgOpt::BOOL_OPT:
			out << "#%" << option_names[i] << "=" << getOptBool(option_names[i]) << endl;
			break;
		case PgOpt::DOUBLE_OPT:
			out << "#%" << option_names[i] << "=" << getOptDouble(option_names[i]) << endl;
			break;
		case PgOpt::INT_OPT:
			out << "#%" << option_names[i] << "=" << getOptInt(option_names[i]) << endl;
			break;
		case PgOpt::STRING_OPT:
			if (getOpt(option_names[i]).empty() == false)
				out << "#%" << option_names[i] << "=" << getOpt(option_names[i]) << endl;
		}
	}
}

void GenotypeExtendedReportEngine::processSingleChannelFiles()
{
	string callFile = getOpt("calls-file");
	string countFile = getOpt("count-file");
	string cvFile = getOpt("cv-file");
	string forcedFile = getOpt("forced-calls-file");
	string confFile = getOpt("confidence-file");
	string sumFile = getOpt("summaries");
	string outFile = getOpt("out-file");
	string sampleCol = getOpt("output-sample-col");
	string markerCol = getOpt("output-marker-col");
	string normFile = getOpt("summaries-norm");
	double k = getOptDouble("stretch-constant");
	enum Transformation transformation = transformationForString(getOpt("cluster-transformation"));

	// Open the input files.
	affx::TsvFile callTsv;
	if(callTsv.open(callFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + callFile);
	callTsv.rewind();
	affx::TsvFile forcedTsv;
	if (forcedFile.empty() == false)
	{
		if(forcedTsv.open(forcedFile) != TSV_OK)
			Err::errAbort("Couldn't open the file: " + forcedFile);
		forcedTsv.rewind();
	}
	affx::TsvFile confTsv;
	if(confTsv.open(confFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + confFile);
	confTsv.rewind();
	affx::TsvFile countTsv;
	if(countTsv.open(countFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + countFile);
	countTsv.rewind();
	affx::TsvFile cvTsv;
	if(cvTsv.open(cvFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + cvFile);
	cvTsv.rewind();
	affx::TsvFile sumTsv;
	if(sumTsv.open(sumFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + sumFile);
	sumTsv.rewind();
	affx::TsvFile normTsv;
	if(normTsv.open(normFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + normFile);
	normTsv.rewind();

	// Check the number of columns.
	if ((callTsv.getColumnCount(0) == confTsv.getColumnCount(0) &&
		 confTsv.getColumnCount(0) == sumTsv.getColumnCount(0) &&
		 sumTsv.getColumnCount(0) == countTsv.getColumnCount(0) &&
		 countTsv.getColumnCount(0) == cvTsv.getColumnCount(0) &&
		 cvTsv.getColumnCount(0) == normTsv.getColumnCount(0)) == false)
		 Err::errAbort("The column counts in the input files do not match");

	// Create an output file
	ofstream out(outFile.c_str(), ios::out);
	writeHeader(out);
	out << getOpt("output-sample-col") << "\t" << getOpt("output-marker-col") << "\t" << "Call\t";
	if (forcedFile.empty() == false)
		out << "Forced Call\t";
	out << "Confidence\tRaw_Signal_A\tRaw_Signal_B\tCount_A\tCount_B\tCV_A\tCV_B\tContrast\tStrength" << endl;

	// The first column are the marker names. The remaining columns are the data.
	// The call codes are: -1=NN, 0=AA, 1=AB, 2=BB
	// Output the data by sample then by marker. One row per marker per sample.
	int ncols = callTsv.getColumnCount(0);
	string sample;
	string marker;
	int callData;
	int forcedData;
	string name;
	string confData;
	string signalData;
	string aSignal;
	string bSignal;
	int countData;
	int aCount;
	int bCount;
	string cvData;
	string aCV;
	string bCV;
	double aNormSignal;
	double bNormSignal;
	double normData;
	double strength;
	double contrast;
	bool dataAvailable;
	string callCodes[] = {"NoCall", "AA", "AB", "BB"};

	// Work on one sample at a time.
	for (int i=1; i<ncols; i++)
	{
		sample = callTsv.getColumnName(0, i);

		callTsv.unbindAll();
		callTsv.rewind();
		callTsv.bind(0, 0, &marker, TSV_BIND_REQUIRED);
		callTsv.bind(0, i, &callData, TSV_BIND_REQUIRED);

		if (forcedFile.empty() == false)
		{
			forcedTsv.unbindAll();
			forcedTsv.rewind();
			forcedTsv.bind(0, i, &forcedData, TSV_BIND_REQUIRED);
		}

		confTsv.unbindAll();
		confTsv.rewind();
		confTsv.bind(0, i, &confData, TSV_BIND_REQUIRED);

		countTsv.unbindAll();
		countTsv.rewind();
		countTsv.bind(0, i, &countData, TSV_BIND_REQUIRED);

		cvTsv.unbindAll();
		cvTsv.rewind();
		cvTsv.bind(0, i, &cvData, TSV_BIND_REQUIRED);

		sumTsv.unbindAll();
		sumTsv.rewind();
		sumTsv.bind(0, 0, &name, TSV_BIND_REQUIRED);
		sumTsv.bind(0, i, &signalData, TSV_BIND_REQUIRED);

		normTsv.unbindAll();
		normTsv.rewind();
		normTsv.bind(0, i, &normData, TSV_BIND_REQUIRED);

		// Read all of the output files for just the one sample's data.
		while (true)
		{
			// If we are at the end of the file then break out of the while loop.
			if (sumTsv.nextLevel(0) != TSV_OK)
				break;

			// Check if the probe set is a control and skip it if its a control.
			if (controlMap.find(name) != controlMap.end())
			{
				countTsv.nextLevel(0);
				cvTsv.nextLevel(0);
				continue;
			}

			countTsv.nextLevel(0);
			aCount = countData;
			countTsv.nextLevel(0);
			bCount = countData;
			dataAvailable = (aCount > 0 && bCount > 0);

			normTsv.nextLevel(0);
			aNormSignal = normData;
			normTsv.nextLevel(0);
			bNormSignal = normData;
			if (dataAvailable == true)
				valueForTransformation(aNormSignal, bNormSignal, transformation, k, contrast, strength);

			aSignal = signalData;
			sumTsv.nextLevel(0);
			bSignal = signalData;

			cvTsv.nextLevel(0);
			aCV = cvData;
			cvTsv.nextLevel(0);
			bCV = cvData;

			callTsv.nextLevel(0);
			confTsv.nextLevel(0);

			out << sample << "\t" << marker << "\t" << callCodes[callData+1] << "\t";
			if (forcedFile.empty() == false)
			{
				forcedTsv.nextLevel(0);
				out << callCodes[forcedData+1] << "\t";
			}

			out << (dataAvailable == true ? confData : "NaN") << "\t";
			out << (aCount > 0 ? aSignal : "NaN") << "\t";
			out << (bCount > 0 ? bSignal : "NaN") << "\t";
			out << aCount << "\t";
			out << bCount << "\t";
			out << (aCount > 0 ? aCV : "NaN") << "\t";
			out << (bCount > 0 ? bCV : "NaN") << "\t";
			if (dataAvailable == true)
				out << contrast << "\t" << strength << endl;
			else
				out << "NaN\tNaN" << endl;
 		}
	}

	// Close the files.
	out.close();
	callTsv.close();
	confTsv.close();
	sumTsv.close();
	if (forcedFile.empty() == false)
		forcedTsv.close();
}

static vector<string> foo(string data[2], string names[2], bool oneFile)
{
	vector<string> result(2);

	result[0] = "NaN";
	result[1] = "NaN";

	int count = oneFile == true ? 1 : 2;
	for (int i = 0; i < count; ++i)
	{
		string::size_type n = names[i].find_last_of("-");
		if (n != string::npos)
		{
			string x = names[i].substr(n + 1, 1);

			if (x == "A")
				result[0] = data[i];
			else
				result[1] = data[i];
		}
	}

	return result;
}

void GenotypeExtendedReportEngine::processMultiChannelFiles()
{
	string callFile = getOpt("calls-file");
	string countFile = getOpt("count-file");
	string cvFile = getOpt("cv-file");
	string forcedFile = getOpt("forced-calls-file");
	string confFile = getOpt("confidence-file");
	string sumFile = getOpt("summaries");
	string outFile = getOpt("out-file");
	string sampleCol = getOpt("output-sample-col");
	string markerCol = getOpt("output-marker-col");
	string normFile = getOpt("summaries-norm");
	double k = getOptDouble("stretch-constant");
	enum Transformation transformation = transformationForString(getOpt("cluster-transformation"));

	// Check that there is two summary and cv files.
	vector<string> sumFiles = StringUtils::Split(sumFile, "?");
	bool oneSumFile = sumFiles.size() == 1;
	if (sumFiles.size() > 2 )
	{
		Err::errAbort("Too many input summaries files.");
	}

	vector<string> cvFiles = StringUtils::Split(cvFile, "?");
	bool oneCVFile = cvFiles.size() == 1;
	if (cvFiles.size() > 2)
	{
		Err::errAbort("Two many input cv files.");
	}

	// Open the input files.
	affx::TsvFile callTsv;
	if(callTsv.open(callFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + callFile);
	callTsv.rewind();
	affx::TsvFile forcedTsv;
	if (forcedFile.empty() == false)
	{
		if(forcedTsv.open(forcedFile) != TSV_OK)
			Err::errAbort("Couldn't open the file: " + forcedFile);
		forcedTsv.rewind();
	}
	affx::TsvFile confTsv;
	if(confTsv.open(confFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + confFile);
	confTsv.rewind();
	affx::TsvFile countTsv;
	if(countTsv.open(countFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + countFile);
	countTsv.rewind();
	affx::TsvFile cvTsv[2];
	if(cvTsv[0].open(cvFiles[0]) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + cvFiles[0]);
	cvTsv[0].rewind();
	if(oneCVFile == false)
	{
		if (cvTsv[1].open(cvFiles[1]) != TSV_OK)
		{
			Err::errAbort("Couldn't open the file: " + cvFiles[1]);
		}
		cvTsv[1].rewind();
	}
	affx::TsvFile sumTsv[2];
	if(sumTsv[0].open(sumFiles[0]) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + sumFiles[0]);
	sumTsv[0].rewind();
	if (oneSumFile == false)
	{
		if (sumTsv[1].open(sumFiles[1]) != TSV_OK)
		{
			Err::errAbort("Couldn't open the file: " + sumFiles[1]);
		}
		sumTsv[1].rewind();
	}
	affx::TsvFile normTsv;
	if(normTsv.open(normFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + normFile);
	normTsv.rewind();

	// Check the number of columns.
	if ((callTsv.getColumnCount(0) == confTsv.getColumnCount(0) &&
		 confTsv.getColumnCount(0) == sumTsv[0].getColumnCount(0) &&
		 sumTsv[0].getColumnCount(0) == countTsv.getColumnCount(0) &&
		 countTsv.getColumnCount(0) == cvTsv[0].getColumnCount(0) &&
		 cvTsv[0].getColumnCount(0) == normTsv.getColumnCount(0)) == false &&
		 (oneCVFile == true ||  cvTsv[1].getColumnCount(0) == countTsv.getColumnCount(0)) &&
		 (oneSumFile == true ||  sumTsv[1].getColumnCount(0) == countTsv.getColumnCount(0)))
		 Err::errAbort("The column counts in the input files do not match");

	// Create an output file
	ofstream out(outFile.c_str(), ios::out);
	writeHeader(out);
	out << getOpt("output-sample-col") << "\t" << getOpt("output-marker-col") << "\t" << "Call\t";
	if (forcedFile.empty() == false)
		out << "Forced Call\t";
	out << "Confidence\tRaw_Signal_A\tRaw_Signal_B\tCount\tCV_A\tCV_B\tContrast\tStrength" << endl;

	// The first column are the marker names. The remaining columns are the data.
	// The call codes are: -1=NN, 0=AA, 1=AB, 2=BB
	// Output the data by sample then by marker. One row per marker per sample.
	int ncols = callTsv.getColumnCount(0);
	string sample;
	string marker;
	int callData;
	int forcedData;
	string sumNames[2];
	string cvNames[2];
	string confData;
	string signalData[2];
	vector<string> signalInAlleleOrder;
	//string aSignal;
	//string bSignal;
	//int countData;
	int count;
	string cvData[2];
	vector<string> cvInAlleleOrder;
	//string aCV;
	//string bCV;
	double aNormSignal;
	double bNormSignal;
	double normData;
	double strength;
	double contrast;
	bool dataAvailable;
	string callCodes[] = {"NoCall", "AA", "AB", "BB"};

	// Work on one sample at a time.
	for (int i=1; i<ncols; i++)
	{
		sample = callTsv.getColumnName(0, i);

		callTsv.unbindAll();
		callTsv.rewind();
		callTsv.bind(0, 0, &marker, TSV_BIND_REQUIRED);
		callTsv.bind(0, i, &callData, TSV_BIND_REQUIRED);

		if (forcedFile.empty() == false)
		{
			forcedTsv.unbindAll();
			forcedTsv.rewind();
			forcedTsv.bind(0, i, &forcedData, TSV_BIND_REQUIRED);
		}

		confTsv.unbindAll();
		confTsv.rewind();
		confTsv.bind(0, i, &confData, TSV_BIND_REQUIRED);

		countTsv.unbindAll();
		countTsv.rewind();
		countTsv.bind(0, i, &count, TSV_BIND_REQUIRED);

		cvTsv[0].unbindAll();
		cvTsv[0].rewind();
		cvTsv[0].bind(0, 0, &cvNames[0], TSV_BIND_REQUIRED);
		cvTsv[0].bind(0, i, &cvData[0], TSV_BIND_REQUIRED);

		if (oneCVFile == false)
		{
			cvTsv[1].unbindAll();
			cvTsv[1].rewind();
			cvTsv[1].bind(0, 0, &cvNames[1], TSV_BIND_REQUIRED);
			cvTsv[1].bind(0, i, &cvData[1], TSV_BIND_REQUIRED);
		}

		sumTsv[0].unbindAll();
		sumTsv[0].rewind();
		sumTsv[0].bind(0, 0, &sumNames[0], TSV_BIND_REQUIRED);
		sumTsv[0].bind(0, i, &signalData[0], TSV_BIND_REQUIRED);

		if (oneSumFile == false)
		{
			sumTsv[1].unbindAll();
			sumTsv[1].rewind();
			sumTsv[1].bind(0, 0, &sumNames[1], TSV_BIND_REQUIRED);
			sumTsv[1].bind(0, i, &signalData[1], TSV_BIND_REQUIRED);
		}

		normTsv.unbindAll();
		normTsv.rewind();
		normTsv.bind(0, i, &normData, TSV_BIND_REQUIRED);

		// Read all of the output files for just the one sample's data.
		while (true)
		{
			// If we are at the end of the file then break out of the while loop.
			if (sumTsv[0].nextLevel(0) != TSV_OK)
				break;

			if (oneSumFile == false && sumTsv[1].nextLevel(0) != TSV_OK)
				break;

			// Check if the probe set is a control and skip it if its a control.
			if (controlMap.find(sumNames[0]) != controlMap.end())
			{
				countTsv.nextLevel(0);
				cvTsv[0].nextLevel(0);
				if (oneCVFile == false)
					cvTsv[1].nextLevel(0);
				continue;
			}

			countTsv.nextLevel(0);
			dataAvailable = (count > 0);

			normTsv.nextLevel(0);
			aNormSignal = normData;
			normTsv.nextLevel(0);
			bNormSignal = normData;
			if (dataAvailable == true)
				valueForTransformation(aNormSignal, bNormSignal, transformation, k, contrast, strength);

			signalInAlleleOrder = foo(signalData, sumNames, oneSumFile);

			cvTsv[0].nextLevel(0);
			if (oneCVFile == false)
				cvTsv[1].nextLevel(0);

			cvInAlleleOrder = foo(cvData, cvNames, oneCVFile);

			callTsv.nextLevel(0);
			confTsv.nextLevel(0);

			out << sample << "\t" << marker << "\t" << callCodes[callData+1] << "\t";
			if (forcedFile.empty() == false)
			{
				forcedTsv.nextLevel(0);
				out << callCodes[forcedData+1] << "\t";
			}

			out << (dataAvailable == true ? confData : "NaN") << "\t";
			out << (count > 0 ? signalInAlleleOrder[0] : "NaN") << "\t";
			out << (count > 0 ? signalInAlleleOrder[1] : "NaN") << "\t";
			out << count << "\t";
			out << (count > 0 ? cvInAlleleOrder[0] : "NaN") << "\t";
			out << (count > 0 ? cvInAlleleOrder[1] : "NaN") << "\t";
			if (dataAvailable == true)
				out << contrast << "\t" << strength << endl;
			else
				out << "NaN\tNaN" << endl;
 		}
	}

	// Close the files.
	out.close();
	callTsv.close();
	confTsv.close();
	sumTsv[0].close();
	if (oneSumFile == false)
		sumTsv[1].close();
	cvTsv[0].close();
	if (oneCVFile == false)
		cvTsv[1].close();
	if (forcedFile.empty() == false)
		forcedTsv.close();
}

void GenotypeExtendedReportEngine::readMixFile()
{
	string mixFile = getOpt("marker-content-file");

	// Open the mix file and store those negative type entries.
	affx::TsvFile tsv;
	if(tsv.open(mixFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + mixFile);
	tsv.rewind();
	int typeValue;
	string note;
	string id;
	tsv.bind(0, getOpt("type-col"), &typeValue, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("note-col"), &note, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("probe-id-col"), &id, TSV_BIND_REQUIRED);
	while(tsv.nextLevel(0) == TSV_OK)
	{
		if (typeValue < 0)
		{
			if (id.empty() == true)
				id = note;
			controlMap[id] = true;
		}
	}
	tsv.clear();
}

void GenotypeExtendedReportEngine::readAnnotationFile()
{
	string annotFile = getOpt("annotation-file");

	// Open the annotation file to determine how to read the count file.
	affx::TsvFile tsv;
	if(tsv.open(annotFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + annotFile);
	int numChannels = 1;
	if (tsv.getHeader("num-channels", numChannels) == TSV_OK)
	{
		annotationFileHasOneCode = (numChannels == 1);
	}
	tsv.clear();
}

/**
   This is the "main()" equivalent of the engine.
*/
void GenotypeExtendedReportEngine::runImp()
{
	readMixFile();
	readAnnotationFile();

	if (annotationFileHasOneCode == true)
	{
		processSingleChannelFiles();
	}
	else
	{
		processMultiChannelFiles();
	}
}


