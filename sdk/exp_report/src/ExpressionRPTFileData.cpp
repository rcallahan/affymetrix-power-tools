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

#include "exp_report/src/ExpressionRPTFileData.h"
//
#include "file/FileIO.h"
//
#include <cstring>
#include <fstream>
#include <istream>
#include <string.h>
#include <string>
#include <sys/stat.h>
//

using namespace ExpressionReport;
using namespace std;

#ifdef WIN32
#pragma warning(disable: 4996)
#endif

/*
 * Check for the file existance
 */
bool ExpressionRPTFileData::Exists()
{
	struct stat st;
	return (stat(fileName.c_str(), &st) == 0);
}

/*
 * Get the ith token from the line
 */
static string GetStringToken(char *line, const char *seps, int index)
{
	char *t = strtok(line, seps);
	for (int i=0; i<index; i++)
	{
		t = strtok(NULL, seps);
	}
	return string(t);
}

/*
 * Get the ith token from the line
 */
static int GetIntToken(char *line, const char *seps, int index)
{
	string t = GetStringToken(line, seps, index);
	if (t.length() > 0)
		return atoi(t.c_str());
	return 0;
}

/*
 * Get the ith token from the line
 */
static float GetFloatToken(char *line, const char *seps, int index)
{
	string t = GetStringToken(line, seps, index);
	if (t.length() > 0)
		return (float)atof(t.c_str());
	return 0.0f;
}

/*
 * Get all of the tokens from the line
 */
static list<string> GetTokens(char *line, const char *seps)
{
	list<string> tokens;
	string token;

	char *t = strtok(line, seps);
	while (t)
	{
		token = t;
		tokens.push_back(token);
		t = strtok(NULL, seps);
	}
	return tokens;
}

/*
 * Remove the colon from the end of the string.
 */
string RemoveColon(string s)
{
	if (s[s.length()-1] == ':')
		s = s.substr(0,s.length()-1);
	return s;
}

/*
 * Read the header of the file.
 */
static void ReadHeader(ifstream &instr, ExpressionReportData &data)
{
	const int max_line_len = 1024;
	char line[max_line_len];

	// Get the date
	ReadNextLine(instr, line, max_line_len);
	data.Date() = GetStringToken(line, "\t", 1);

	// Get the next line (spacer)
	ReadNextLine(instr, line, max_line_len);

	// Get the file name, array type, algorithm, threshold and control flag.
	ReadNextLine(instr, line, max_line_len);
	data.CHPFileName() = GetStringToken(line, "\t", 1);
	ReadNextLine(instr, line, max_line_len);
	data.ArrayType() = GetStringToken(line, "\t", 1);
	ReadNextLine(instr, line, max_line_len);
	data.AlgName() = GetStringToken(line, "\t", 1);
	ReadNextLine(instr, line, max_line_len);
	data.ProbePairThreshold() = GetIntToken(line, "\t", 1);
	ReadNextLine(instr, line, max_line_len);
	data.AntiSenseControls() = (GetStringToken(line, "\t", 1) == "Antisense" ? true : false);

	// Get the next line (spacer)
	ReadNextLine(instr, line, max_line_len);
}

/*
 * Read the parameters.
 */
static void ReadParameters(ifstream &instr, ExpressionReportData &data)
{
	const int max_line_len = 1024;
	char line[max_line_len];

	// Get the parameters
	const char *sep = "__________";
	NameValuePairList params;
	NameValuePair param;
	while (1)
	{
		ReadNextLine(instr, line, max_line_len);
		if (strncmp(line, sep, strlen(sep)) == 0)
			break;

		list<string> nameval = GetTokens(line, "\t");
		if (nameval.size() == 2)
		{
			list<string>::iterator it = nameval.begin();
			param.name = RemoveColon(*it);
			++it;
			param.value = *it;
			params.push_back(param);
		}
	}
	data.AlgParams() = params;
}

/*
 * Read the background, noise and controls information.
 */
static void ReadBackgroundNoiseControls(ifstream &instr, ExpressionReportData &data)
{
	const char *sep = "__________";
	const int max_line_len = 1024;
	char line[max_line_len];
	while (1)
	{
		ReadNextLine(instr, line, max_line_len);
		if (strncmp(line, sep, strlen(sep)) == 0)
			break;

		// Background.
		if (strcmp(line, "Background:") == 0)
		{
			ReadNextLine(instr, line, max_line_len);
			list<string> nameval = GetTokens(line, "\t ");
			list<string>::iterator it = nameval.begin();
			++it;
			data.BackgroundStats().avg = (float) atof(it->c_str());
			++it;
			++it;
			data.BackgroundStats().std = (float) atof(it->c_str());
			++it;
			++it;
			data.BackgroundStats().min = (float) atof(it->c_str());
			++it;
			++it;
			data.BackgroundStats().max = (float) atof(it->c_str());
		}

		// Noise
		else if (strcmp(line, "Noise:") == 0)
		{
			ReadNextLine(instr, line, max_line_len);
			list<string> nameval = GetTokens(line, "\t ");
			list<string>::iterator it = nameval.begin();
			++it;
			data.NoiseStats().avg = (float) atof(it->c_str());
			++it;
			++it;
			data.NoiseStats().std = (float) atof(it->c_str());
			++it;
			++it;
			data.NoiseStats().min = (float) atof(it->c_str());
			++it;
			++it;
			data.NoiseStats().max = (float) atof(it->c_str());
		}
		
		// Corner+
		else if (strcmp(line, "Corner+") == 0 || strcmp(line, "Corner-") == 0 ||
			strcmp(line, "Central+") == 0 || strcmp(line, "Central-") == 0)
		{
			NameAvgCount value;
			value.name = line;
			ReadNextLine(instr, line, max_line_len);
			list<string> nameval = GetTokens(line, "\t ");
			list<string>::iterator it = nameval.begin();
			++it;
			value.avg = (float) atof(it->c_str());
			++it;
			++it;
			value.count = atoi(it->c_str());
			data.ControlStats().push_back(value);
		}
	}
}

/*
 * Read the probe set stats.
 */
static void ReadProbeSetStats(ifstream &instr, ExpressionReportData &data)
{
	const int max_line_len = 1024;
	char line[max_line_len];

	// Skip until "Total Probe Sets" is found.
	const char *tps = "Total Probe Sets";
	while (1)
	{
		ReadNextLine(instr, line, max_line_len);
		if (strncmp(line, tps, strlen(tps)) == 0)
			break;
	}

	data.ProbeSetResults().NumSets() = GetIntToken(line, "\t", 1);
	ReadNextLine(instr, line, max_line_len);
	data.ProbeSetResults().PresentCalls().Count() = GetIntToken(line, "\t", 1);
	ReadNextLine(instr, line, max_line_len);
	data.ProbeSetResults().AbsentCalls().Count() = GetIntToken(line, "\t", 1);
	ReadNextLine(instr, line, max_line_len);
	data.ProbeSetResults().MarginalCalls().Count() = GetIntToken(line, "\t", 1);

	ReadNextLine(instr, line, max_line_len);
	data.ProbeSetResults().PresentCalls().Signal() =
		data.ProbeSetResults().PresentCalls().Count() * GetFloatToken(line, "\t", 1);

	ReadNextLine(instr, line, max_line_len);
	data.ProbeSetResults().AbsentCalls().Signal() =
		data.ProbeSetResults().AbsentCalls().Count() * GetFloatToken(line, "\t", 1);

	ReadNextLine(instr, line, max_line_len);
	data.ProbeSetResults().MarginalCalls().Signal() =
		data.ProbeSetResults().MarginalCalls().Count() * GetFloatToken(line, "\t", 1);

	// Skip the next separator.
	ReadNextLine(instr, line, max_line_len);
}

/*
 * Convert the string to a detection value.
 */
ReportDataAccessor::DetectionCall DetectionCallFromString(const string &det)
{
	if (det == "p" || det == "P")
		return ReportDataAccessor::DetectionPresent;
	else if (det == "a" || det == "A")
		return ReportDataAccessor::DetectionAbsent;
	else if (det == "m" || det == "M")
		return ReportDataAccessor::DetectionMarginal;
	else
		return ReportDataAccessor::DetectionNoCall;
}

/*
 * Read the controls
 */
static void ReadControls(ifstream &instr, ExpressionControlResultList &controls)
{
	const int max_line_len = 1024;
	char line[max_line_len];

	// Skip the header
	ReadNextLine(instr, line, max_line_len);

	// Read until a separator is found.
	const char *sep = "__________";
	const int COLS_PER_CONTROL_LINE = 9;
	while (1)
	{
		ReadNextLine(instr, line, max_line_len);
		if (strncmp(line, sep, strlen(sep)) == 0)
			break;

		list<string> tokens = GetTokens(line, "\t");
		if (tokens.size() == COLS_PER_CONTROL_LINE)
		{
			list<string>::iterator it = tokens.begin();
			
			ExpressionControlResult control;
			control.SetName(*it);
			++it;
			control.SetControlSignalResult(ExpressionControl::FIVE_PRIME_PROBE_SET, (float)atof(it->c_str()));
			++it;
			control.SetControlDetectionResult(ExpressionControl::FIVE_PRIME_PROBE_SET, DetectionCallFromString(*it));
			++it;
			control.SetControlSignalResult(ExpressionControl::MIDDLE_PROBE_SET, (float)atof(it->c_str()));
			++it;
			control.SetControlDetectionResult(ExpressionControl::MIDDLE_PROBE_SET, DetectionCallFromString(*it));
			++it;
			control.SetControlSignalResult(ExpressionControl::THREE_PRIME_PROBE_SET, (float)atof(it->c_str()));
			++it;
			control.SetControlDetectionResult(ExpressionControl::THREE_PRIME_PROBE_SET, DetectionCallFromString(*it));
			++it; // ignore the sig(all) - this is just the average of the other signal values.
			++it;
			control.SetThreeFiveRatio((float)atof(it->c_str()));
			controls.push_back(control);
		}
	}
}

/*
 * Read the controls
 */
static void ReadControls(ifstream &instr, ExpressionReportData &data)
{
	const int max_line_len = 1024;
	char line[max_line_len];
	while (instr.eof() == false)
	{
		ReadNextLine(instr, line, max_line_len);
		if (strcmp(line, "Housekeeping Controls:") == 0)
			ReadControls(instr, data.HousekeepingStats());

		else if (strcmp(line, "Spike Controls:") == 0)
			ReadControls(instr, data.SpikeStats());
	}
}

/*
 * Read the file contents.
 */
bool ExpressionRPTFileData::Read(ExpressionReportData &data)
{
	data.Clear();

	// Open the file.
	ifstream instr(fileName.c_str(), ios::in);
	if (!instr)
		return false;

	const int max_line_len = 1024;
	char line[max_line_len];

	// Get the first line and check if the file is an expression report.
	ReadNextLine(instr, line, max_line_len);
	if (strcmp(line, "Report Type:	Expression Report") != 0)
		return false;

	// Read the data.
	ReadHeader(instr, data);
	ReadParameters(instr, data);
	ReadBackgroundNoiseControls(instr, data);
	ReadProbeSetStats(instr, data);
	ReadControls(instr, data);

	return true;
}
