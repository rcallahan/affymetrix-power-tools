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

#include "chipstream/apt-genotype-report/GenotypeReportEngine.h"
//
#include "chipstream/EngineUtil.h"
//
#include "file/TsvFile/TsvFile.h"
#include "stats/statfun.h"
#include "util/PgOptions.h"
#include "calvin_files/utils/src/StringUtils.h"

//
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//
using namespace std;
using namespace affx;

// Too bad we don't have the boost::tuple class.
template <class T1, class T2, class T3> struct triplet
{
	typedef T1 first_type;
	typedef T2 second_type;
	typedef T3 third_type;

	T1 first;
	T2 second;
	T3 third;
	triplet() : first(T1()), second(T2()), third(T3()) {}
	triplet(const T1& x, const T2& y, const T3& z) : first(x), second(y), third(z) {}
	template <class U, class V, class W>
		triplet (const triplet<U,V,W> &p) : first(p.first), second(p.second), third(p.third) { }
};

static float percentile(const vector<float> &data, int percentile)
{
	if (data.size() == 0)
		return numeric_limits<float>::quiet_NaN();
	double dindex = percentile * ((data.size()-1) / 100.0f);
	int index = (int) (dindex + 0.5f);
	if (abs(dindex - index) <= numeric_limits<float>::epsilon())
		return data[index];
	else
		return (data[index-1] + data[index])/2;
}

GenotypeReportEngine::Reg GenotypeReportEngine::reg;

GenotypeReportEngine * GenotypeReportEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == GenotypeReportEngine::EngineName())
		return (GenotypeReportEngine *)engine;
	return NULL;
}

GenotypeReportEngine::GenotypeReportEngine() {
	defineOptions();
}

GenotypeReportEngine::~GenotypeReportEngine() {
}

void GenotypeReportEngine::defineOptions() {
	defineOption("i", "calls-file", PgOpt::STRING_OPT,
		"Text file specifying genotype calls process.",
		"");
	defineOption("s", "sample-out-file", PgOpt::STRING_OPT,
		"File to write sample results into.",
		"");
	defineOption("m", "marker-out-file", PgOpt::STRING_OPT,
		"File to write marker results into.",
		"");
	defineOption("", "call-rate-col", PgOpt::STRING_OPT,
		"The call rate column header.",
		"Call Rate");
	defineOption("", "sample-col", PgOpt::STRING_OPT,
		"The sample column header.",
		"Sample");
	defineOption("", "sample-status-col", PgOpt::STRING_OPT,
		"The sample status column header.",
		"Sample Status");
	defineOption("", "hetfreq-col", PgOpt::STRING_OPT,
		"The het frequency column header.",
		"Het Freq");
	defineOption("", "count-file", PgOpt::STRING_OPT,
		"Text file specifying counts.",
		"");
	defineOption("", "cv-file", PgOpt::STRING_OPT,
		"Text files specifying cvs. Multiple files should be delimited by ? character.",
		"");
	defineOption("", "summaries", PgOpt::STRING_OPT,
		"Text files specifying the signal summaries. Multiple files should be delimited by ? character.",
		"");
	defineOption("", "cv-percentile", PgOpt::INT_OPT,
		"The percentile to report for the CV data.",
		"75");
	defineOption("", "signal-percentile", PgOpt::INT_OPT,
		"The percentile to report for the signal data.",
		"50");
	defineOption("", "max-marker-het-freq", PgOpt::DOUBLE_OPT,
		"The maximum marker het frequency.",
		"0");
	defineOption("", "min-marker-call-rate", PgOpt::DOUBLE_OPT,
		"The minimum marker call rate.",
		"0");
	defineOption("", "min-marker-hwp", PgOpt::DOUBLE_OPT,
		"The minimum marker hw p-value.",
		"0");
	defineOption("", "min-marker-maf", PgOpt::DOUBLE_OPT,
		"The minimum marker minor allele frequency.",
		"0");
	defineOption("", "marker-content-file", PgOpt::STRING_OPT,
		"Text file specifying the mix file.",
		"");
	defineOption("", "type-col", PgOpt::STRING_OPT,
		"The column name for the type.",
		"Type");
	defineOption("", "note-col", PgOpt::STRING_OPT,
		"The column name for the note.",
		"Note");
	defineOption("", "probe-id-col", PgOpt::STRING_OPT,
		"The column name for the probe id.",
		"Probe Id");
	defineOption("", "pass-text", PgOpt::STRING_OPT,
		"The text to show for markers that pass the threshold test.",
		"Pass");
	defineOption("", "fail-text", PgOpt::STRING_OPT,
		"The text to show for markers that fail the threshold test.",
		"Fail");
	defineOption("", "cr-thr-value", PgOpt::DOUBLE_OPT,
		"The call rate threshold.",
		"95.0");
	defineOption("", "cr-thr-method", PgOpt::STRING_OPT,
		"The method to compare the call rate. Can be one of \"lt\", \"gt\", \"le\", \"ge\", \"eq\", \"ne\".",
		"gt");
	defineOption("", "hetfreq-thr-value", PgOpt::DOUBLE_OPT,
		"The call rate threshold.",
		"70.0");
	defineOption("", "hetfreq-thr-method", PgOpt::STRING_OPT,
		"The method to compare the call rate. Can be one of \"lt\", \"gt\", \"le\", \"ge\", \"eq\", \"ne\".",
		"lt");
	defineOption("", "min-particle-count", PgOpt::INT_OPT,
		"The minimum particle count.",
		"10");
}

void GenotypeReportEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void GenotypeReportEngine::checkOptionsImp()
{
	defineStates();
	string str = getOpt("sample-out-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an output sample summary file.");

	str = getOpt("marker-out-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an output marker summary file.");

	str = getOpt("calls-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an input calls file.");

	str = getOpt("count-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an input count file.");

	str = getOpt("cv-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an input cv file.");

	str = getOpt("summaries");
	if (str.empty() == true)
		Err::errAbort("Must specify an input summary file.");

	str = getOpt("marker-content-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an input marker file.");

	str = getOpt("cr-thr-method");
	if (str != "lt" && str != "gt" && str != "le" && str != "ge" && str != "eq" && str != "ne")
		Err::errAbort("The method must be one of {\"lt\", \"gt\", \"le\", \"ge\", \"eq\", \"ne\"}.");

	str = getOpt("hetfreq-thr-method");
	if (str != "lt" && str != "gt" && str != "le" && str != "ge" && str != "eq" && str != "ne")
		Err::errAbort("The method must be one of {\"lt\", \"gt\", \"le\", \"ge\", \"eq\", \"ne\"}.");
}

void GenotypeReportEngine::writeHeader(std::ofstream &out)
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

static bool PassThreshold(double callRate, double thr, const string &thrMethod)
{
	if (thrMethod == "lt" && callRate < thr)
		return true;
	
	else if (thrMethod == "gt" && callRate > thr)
		return true;

	else if (thrMethod == "le" && callRate <= thr)
		return true;
	
	else if (thrMethod == "ge" && callRate >= thr)
		return true;
	
	else if (thrMethod == "eq" && abs(callRate - thr) <= numeric_limits<double>::epsilon())
		return true;
	
	else if (thrMethod == "ne" && abs(callRate - thr) > numeric_limits<double>::epsilon())
		return true;

	return false;
}

map<string, bool> GenotypeReportEngine::CalculateMarkerReport()
{
	map<string, bool> passedMarkers;
	double callRate;
	double hetFreq;
	double hwp;
	double maf;
	double maxMarkerHetFreq = getOptDouble("max-marker-het-freq");
	double minMarkerCallRate = getOptDouble("min-marker-call-rate");
	double minMarkerHWP = getOptDouble("min-marker-hwp");
	double minMarkerMAF = getOptDouble("min-marker-maf");
	string pass = getOpt("pass-text");
	string fail = getOpt("fail-text");

	// Open the input call file.
	affx::TsvFile tsv;
	string in_file = getOpt("calls-file");
	if(tsv.open(in_file) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + in_file);
	tsv.rewind();

	// The first column are the marker names. The
	// remaining columns are the genotype calls.
	int nsamples = tsv.getColumnCount(0) - 1;
	vector<int> calls(nsamples);
	string marker;
	tsv.bind(0, 0, &marker, TSV_BIND_REQUIRED);
	for (int i=0; i<nsamples; i++)
		tsv.bind(0, i+1, &calls[i], TSV_BIND_REQUIRED);

	// Create the output file for the per marker results
	string out_file = getOpt("marker-out-file");
	ofstream out(out_file.c_str(), ios::out);
	writeHeader(out);
	out.setf(ios::fixed, ios::floatfield);
	out <<
		"Marker" << "\t" <<
		"Marker Status" << "\t" <<
		"Call Rate" << "\t" <<
		"Het Freq" << "\t" <<
		"MAF" << "\t" <<
		"HW p-value" <<
		endl;

	// Read the data and calculate the per snp metrics.
	// The call codes are: -1=NN, 0=AA, 1=AB, 2=BB
	while(tsv.nextLevel(0) == TSV_OK)
	{
		bool passMaxMarkerHetFreq = true;
		bool passMinMarkerCallRate = true;
		bool passMinMarkerHWP = true;
		bool passMinMarkerMAF = true;
		int aa=0;
		int ab=0;
		int bb=0;
		int nc=0;
		for (int i=0; i<nsamples; i++)
		{
			if (calls[i] == 0)
				++aa;
			else if (calls[i] == 1)
				++ab;
			else if (calls[i] == 2)
				++bb;
			else
				++nc;
		}

		callRate = 100.0f * (aa+ab+bb) / nsamples;
		passMinMarkerCallRate = (callRate >= minMarkerCallRate);

		if (aa + ab + bb > 0)
		{
			hetFreq = (float) ab / (aa+ab+bb);
			passMaxMarkerHetFreq = (hetFreq <= maxMarkerHetFreq);

			float freqA = (float)(2*aa + ab) / (2*aa + 2*ab + 2*bb);
			float freqB = (float)(2*bb + ab) / (2*aa + 2*ab + 2*bb);
			maf = min(freqA, freqB);
			passMinMarkerMAF = (maf >= minMarkerMAF);
		}

		hwp = affxstat::CalcHWEqPValue(aa, ab, bb);
		if (hwp != std::numeric_limits<double>::signaling_NaN() && hwp != std::numeric_limits<double>::quiet_NaN() && hwp != std::numeric_limits<double>::infinity())
			passMinMarkerHWP = (hwp >= minMarkerHWP);



		// Output the marker name.
		out << marker << "\t";

		// Output the pass/fail status
		if (passMaxMarkerHetFreq && passMinMarkerCallRate && passMinMarkerHWP && passMinMarkerMAF)
		{
			passedMarkers[marker] = true;
			out << pass;
		}
		else
		{
			out << fail + ":";
			if (passMaxMarkerHetFreq == false)
				out << " HF";
			if (passMinMarkerCallRate == false)
				out << " CR";
			if (passMinMarkerHWP == false)
				out << " HWP";
			if (passMinMarkerMAF == false)
				out << " MAF";
		}
		out << "\t";

		// Output the call rate and determine if it passes the threshold test
		out.precision(1);
		out << callRate << "\t";

		// Output the het freq and MAF and determine if it passes the threshold test
		if (aa + ab + bb == 0)
		{
			out << "NaN" << "\t";
			out << "NaN" << "\t";
		}
		else
		{
			out.precision(1);
			out << hetFreq * 100. << "\t";
			out << maf * 100. << "\t";
		}
		
		// Output the HW p-value and determine if it passes the threshold test
		out.precision(6);
		if (hwp == std::numeric_limits<double>::signaling_NaN() || hwp == std::numeric_limits<double>::quiet_NaN() || hwp == std::numeric_limits<double>::infinity())
			out << "NaN";
		else
			out << hwp;

		out << endl;
	}
	out.close();
	tsv.close();
	return passedMarkers;
}

static vector<string> GetSampleNames(const string& fileName)
{
	// Get the sample names from the input file. These are the column names (except the first column).
	affx::TsvFile tsv;
	if(tsv.open(fileName) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + fileName);
	tsv.rewind();
	int nsamples = tsv.getColumnCount(0) - 1;
	vector<string> samples(nsamples);
	for (int i=0; i<nsamples; i++)
		samples[i] = tsv.getColumnName(0, i+1);
	tsv.close();
	return samples;
}

static vector<triplet<float, float, float> > CalculateCallRates(const string &fileName, int nsamples, const map<string, bool> &passedMarkers)
{
	// Bind to the call columns.
	affx::TsvFile tsv;
	if(tsv.open(fileName) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + fileName);
	tsv.rewind();
	string marker;
	vector<int> calls(nsamples);
	for (int i=0; i<nsamples; i++)
		tsv.bind(0, i+1, &calls[i], TSV_BIND_REQUIRED);
	tsv.bind(0, 0, &marker, TSV_BIND_REQUIRED);

	// Read the data and calculate call counts for each sample.
	int rows = 0;
	int passedRows = 0;
	vector<triplet<float,float,float> > callRates(nsamples);
	callRates.assign(nsamples, triplet<float,float,float>(0.0f, 0.0f, 0.0f));
	while(tsv.nextLevel(0) == TSV_OK)
	{
		++rows;
		for (int i=0; i<nsamples; i++)
		{
			if (calls[i] != -1)
				++callRates[i].first;
			if (calls[i] == 1)
				++callRates[i].third;
		}

		if (passedMarkers.find(marker) != passedMarkers.end())
		{
			++passedRows;
			for (int i=0; i<nsamples; i++)
			{
				if (calls[i] != -1)
					++callRates[i].second;
			}
		}
	}
	tsv.close();

	// Calculate the call rates
	for (int i=0; i<nsamples; i++)
	{
		callRates[i].first = (callRates[i].first * 100.0f) / rows;
		callRates[i].second = (passedRows > 0 ? (callRates[i].second * 100.0f) / passedRows : 0.0f);
		callRates[i].third = (callRates[i].third * 100.0f) / rows;
	}

	return callRates;
}

static vector<float> CalculatePercentileReportedCV(const string &fileName, const string &countFile, int nsamples, int cvPercentile, vector<string>& signalColumns)
{
	// Open the count file
	affx::TsvFile countTsv;
	if(countTsv.open(countFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + countFile);
	countTsv.rewind();
	vector<int> count(nsamples);
	for (int i=0; i<nsamples; i++)
		countTsv.bind(0, i+1, &count[i], TSV_BIND_REQUIRED);

	vector<float> sampleCV(nsamples);
	vector<vector<float> > sampleCVs(nsamples);

	vector<string> cvfiles = StringUtils::Split(fileName, "?");

	for (vector<string>::iterator ii = cvfiles.begin(); ii != cvfiles.end(); ++ii)
	{
		// Read all of the CV's into a vector, one per sample
		affx::TsvFile tsv;
		if(tsv.open(*ii) != TSV_OK)
			Err::errAbort("Couldn't open the file: " + *ii);

		string val;
		tsv.headersBegin();
		if (tsv.headersFindNext("signal-column", val) == TSV_OK)
		{
			signalColumns.push_back(val);
		}
		else
		{
			signalColumns.push_back("");
		}
		tsv.rewind();
		for (int i=0; i<nsamples; i++)
			tsv.bind(0, i+1, &sampleCV[i], TSV_BIND_REQUIRED);
		int rows = 0;

		countTsv.rewind();

		// Collect the CV's into a vector, one per sample.
		// This assumes that the marker order in the cv file is the same as the counts file.
		while(tsv.nextLevel(0) == TSV_OK)
		{
			countTsv.nextLevel(0);
			++rows;
			for (int i=0; i<nsamples; i++)
			{
				if (count[i] > 0)
					sampleCVs[i].push_back(sampleCV[i]);
			}
 		}
		tsv.close();
	}

	countTsv.close();

	// Calculate the CV's for each sample.
	sampleCV.assign(nsamples, 0.0f);
	for (int i=0; i<nsamples; i++)
	{
		if (sampleCVs[i].empty() == false)
		{
			std::sort(sampleCVs[i].begin(), sampleCVs[i].end());
			sampleCV[i] = percentile(sampleCVs[i], cvPercentile);
		}
	}

	return sampleCV;
}

static void CalculatePercentileSignal(const string &signalFileName, const string &countFileName, map<string, bool> controls, int nsamples, int signalPercentile, int minCount, vector<vector<float> > &signalData, vector<vector<float> > &wavg)
{
	// Bind the count data
	affx::TsvFile countTsv;
	if(countTsv.open(countFileName) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + countFileName);
	countTsv.rewind();
	vector<int> sampleCount(nsamples);
	for (int i=0; i<nsamples; i++)
		countTsv.bind(0, i+1, &sampleCount[i], TSV_BIND_REQUIRED);

	bool firstSignalFile = true;
	vector<string> signalfiles = StringUtils::Split(signalFileName, "?");

	vector<float> boundSignalData(nsamples);
	signalData.resize(signalfiles.size());
	wavg.resize(signalfiles.size());
	vector<vector<vector<float> > > sampleSignals(signalfiles.size());
	vector<int> wcount(nsamples, 0);

	for (int ichannel = 0; ichannel < signalfiles.size(); ++ichannel)
	{
		// Bind the signal data
		affx::TsvFile signalTsv;
		if(signalTsv.open(signalfiles[ichannel]) != TSV_OK)
			Err::errAbort("Couldn't open the file: " + signalfiles[ichannel]);
		signalTsv.rewind();
		string name;
		signalTsv.bind(0, 0, &name, TSV_BIND_REQUIRED);
		for (int i=0; i<nsamples; i++)
			signalTsv.bind(0, i+1, &boundSignalData[i], TSV_BIND_REQUIRED);

		countTsv.rewind();

		sampleSignals[ichannel].resize(nsamples);
		wavg[ichannel].resize(nsamples, 0);

		// Read the signal and count files and store the signal values and
		// compute the weighted control signals.
		while(signalTsv.nextLevel(0) == TSV_OK)
		{
			countTsv.nextLevel(0);	// moved per Luis 1-22-2010
			for (int i=0; i<nsamples; i++)
			{
				if (sampleCount[i] >= minCount)
					sampleSignals[ichannel][i].push_back(boundSignalData[i]);
			}
			if (controls.find(name) != controls.end())
			{
				for (int i=0; i<nsamples; i++)
				{
					if (firstSignalFile == true)
						wcount[i] += sampleCount[i];
					wavg[ichannel][i] += (boundSignalData[i] * sampleCount[i]);
				}
			}
 		}
		firstSignalFile = false;
		signalTsv.close();
	}

	countTsv.close();

	// Compute the controls weighted average
	for (int ichannel = 0; ichannel < signalfiles.size(); ++ichannel)
	{
		for (int i=0; i<nsamples; i++)
			wavg[ichannel][i] = (wcount[i] > 0 ? wavg[ichannel][i] / wcount[i] : std::numeric_limits<float>::quiet_NaN());
	}

	// Compute the percentile value for each sample
	for (int ichannel = 0; ichannel < signalfiles.size(); ++ichannel)
	{
		signalData[ichannel].assign(nsamples, 0.0f);

		for (int i=0; i<nsamples; i++)
		{
			if (sampleSignals[ichannel][i].empty() == false)
			{
				std::sort(sampleSignals[ichannel][i].begin(), sampleSignals[ichannel][i].end());
				signalData[ichannel][i] = percentile(sampleSignals[ichannel][i], signalPercentile);
			}
		}
	}
}

static void CalculateSampleCounts(const string &fileName, int nsamples, vector<float> &medianCodeCount, vector<float> &knownCodeCount, vector<float> &minCodeCount, vector<float> &unknownCodeRate)
{
	// Read all of the count's into a vector, one per sample
	affx::TsvFile tsv;
	if(tsv.open(fileName) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + fileName);
	tsv.rewind();

	// Allocate memory for the results
	medianCodeCount.resize(nsamples);
	medianCodeCount.assign(nsamples, 0);
	knownCodeCount.resize(nsamples);
	knownCodeCount.assign(nsamples, 0);
	minCodeCount.resize(nsamples);
	minCodeCount.assign(nsamples, 0);
	unknownCodeRate.resize(nsamples);
	unknownCodeRate.assign(nsamples, 0);

	// Read the header line for the unknown rates
	TsvFileHeaderLine* hdr;
	int sampleIdx = 0;
	while ((hdr = tsv.nextHeaderPtr()) != NULL)
	{
		if (hdr->m_key.find("Unknown Code Rate.") != string::npos)
		{
			unknownCodeRate[sampleIdx++] = atof(hdr->m_value.c_str());
		}
	}

	// Get the counts
	vector<float> sampleCount(nsamples);
	vector<vector<float> > sampleCounts(nsamples);
	for (int i=0; i<nsamples; i++)
		tsv.bind(0, i+1, &sampleCount[i], TSV_BIND_REQUIRED);
	int rows = 0;
	while(tsv.nextLevel(0) == TSV_OK)
	{
		++rows;
		for (int i=0; i<nsamples; i++)
		{
			sampleCounts[i].push_back(sampleCount[i]);
		}
	}
	tsv.close();

	// Compute the median, min, and known code counts
	for (int i=0; i<nsamples; i++)
	{
		std::sort(sampleCounts[i].begin(), sampleCounts[i].end());
		medianCodeCount[i] = percentile(sampleCounts[i], 50);
		minCodeCount[i] = sampleCounts[i][0];
		for (int j=0; j<(int)sampleCounts[i].size(); j++)
			knownCodeCount[i] += sampleCounts[i][j];
	}
}

void GenotypeReportEngine::CalculateSampleReport(const map<string, bool> &passedMarkers)
{
	// Get the sample names and metrics
	double crThr = getOptDouble("cr-thr-value");
	string crThrMethod = getOpt("cr-thr-method");
	double hetfreqThr = getOptDouble("hetfreq-thr-value");
	string hetfreqThrMethod = getOpt("hetfreq-thr-method");
	vector<string> samples = GetSampleNames(getOpt("calls-file"));
	int nsamples = (int) samples.size();
	vector<triplet<float,float,float> > callRates = CalculateCallRates(getOpt("calls-file"), nsamples, passedMarkers);
	vector<string> signalColumns;
	map<string, bool> controls = GetControlsFromMixFile();
	vector<float> cvData = CalculatePercentileReportedCV(getOpt("cv-file"), getOpt("count-file"), nsamples, getOptInt("cv-percentile"), signalColumns);
	vector<vector<float> > signalData;
	vector<vector<float> > wavg;
	CalculatePercentileSignal(getOpt("summaries"), getOpt("count-file"), controls, nsamples, getOptInt("signal-percentile"), getOptInt("min-particle-count"), signalData, wavg);
	vector<float> medianCodeCount;
	vector<float> knownCodeCount;
	vector<float> minCodeCount;
	vector<float> unknownCodeRate;
	CalculateSampleCounts(getOpt("count-file"), nsamples, medianCodeCount, knownCodeCount, minCodeCount, unknownCodeRate);

	// Output the metrics to the report file.
	// Med Code Signal1: for each sample, median signal value in signals_report (raw, not normalized)
	// Known Code Count: for each sample, sum of counts that made it into the count report
	// Med Code Count: for each sample, what is the median count across all codes that made it into the count report?
	// Min Code Count: for each sample, what is the minimum count across all codes that made it into the count report?
	// UQuart Code CV1: for each sample, what is the 75%-ile rank of reported CV values of the first fluorescence channel
	string out_file = getOpt("sample-out-file");
	ofstream out(out_file.c_str(), ios::out);
	writeHeader(out);
	out.setf(ios::fixed, ios::floatfield);
	out << getOpt("sample-col") <<
		"\t" << "Sample Status" <<
		"\t" << getOpt("call-rate-col") <<
		"\t" << getOpt("call-rate-col") << " (passing markers)" <<
		"\t" << getOpt("hetfreq-col");
	
	int nchannel = wavg.size();
	for (int ichannel = 0; ichannel < nchannel; ++ichannel)
	{
		out << "\t" << "Median Sig - Bkgd";
		if (nchannel > 1 && signalColumns[ichannel].empty() == false)
			out << " (" << signalColumns[ichannel] << ")";
	}
	for (int ichannel = 0; ichannel < nchannel; ++ichannel)
	{
		out << "\t" << "Median Code Signal";
		if (nchannel > 1 && signalColumns[ichannel].empty() == false)
			out << " (" << signalColumns[ichannel] << ")";
	}
	for (int ichannel = 0; ichannel < nchannel; ++ichannel)
	{
		out << "\t" << "NegCtrl Signal";
		if (nchannel > 1 && signalColumns[ichannel].empty() == false)
			out << " (" << signalColumns[ichannel] << ")";
	}

	out << "\t" << "Known Code Count" <<
		"\t" << "Median Code Count" <<
		"\t" << "Minimum Code Count" <<
		//"\t" << "Upper Quartile Code CV" <<
		"\t" << "Unknown Code Rate" <<
		endl;

	for (int i=0; i<nsamples; i++)
	{
		out.precision(1);
		out <<
			samples[i] <<
			"\t" << (PassThreshold(callRates[i].first, crThr, crThrMethod) == true && 
								PassThreshold(callRates[i].third, hetfreqThr, hetfreqThrMethod) == true ? "Pass" : "Fail") <<
			"\t" << callRates[i].first <<
			"\t" << callRates[i].second <<
			"\t" << callRates[i].third;

		for (int ichannel = 0; ichannel < wavg.size(); ++ichannel)
		{
			if (wavg[ichannel][i] != wavg[ichannel][i])
				out << "\t" << "NaN";
			else
				out << "\t" << signalData[ichannel][i] - wavg[ichannel][i];
		}
		for (int ichannel = 0; ichannel < nchannel; ++ichannel)
		{
			out << "\t" << signalData[ichannel][i];
		}
		for (int ichannel = 0; ichannel < nchannel; ++ichannel)
		{
			if (wavg[ichannel][i] != wavg[ichannel][i])
				out << "\t" << "NaN";
			else
				out << "\t" << wavg[ichannel][i];
		}
		out.precision(0);
		out <<
			"\t" << knownCodeCount[i] <<
			"\t" << medianCodeCount[i] <<
			"\t" << minCodeCount[i];
		//out.precision(2);
		//out << "\t" << cvData[i];
		out.precision(1);
		out << 
			"\t" << unknownCodeRate[i] <<
			endl;
	}
	out.close();
}

/**
* This is the "main()" equivalent of the engine.
*/
void GenotypeReportEngine::runImp()
{
	map<string, bool> passedMarkers = CalculateMarkerReport();
	CalculateSampleReport(passedMarkers);
}

map<string, bool> GenotypeReportEngine::GetControlsFromMixFile()
{
	// Open the mix file and store those negative type entries.
	string mixFile = getOpt("marker-content-file");
	affx::TsvFile mixTsv;
	if(mixTsv.open(mixFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + mixFile);
	mixTsv.rewind();
	map<string, bool> controlMap;
	int typeValue;
	string note;
	string id;
	mixTsv.bind(0, getOpt("type-col"), &typeValue, TSV_BIND_REQUIRED);
	mixTsv.bind(0, getOpt("note-col"), &note, TSV_BIND_REQUIRED);
	mixTsv.bind(0, getOpt("probe-id-col"), &id, TSV_BIND_REQUIRED);
	while(mixTsv.nextLevel(0) == TSV_OK)
	{
		if (typeValue < 0)
		{
			if (id.empty() == true)
				id = note;
			controlMap[id] = true;
		}
	}
	mixTsv.clear();
	return controlMap;
}
