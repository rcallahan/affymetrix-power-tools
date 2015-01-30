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
* @file   SnpSummaryReporter.h
* @author Luis Jevons
* @date   Wed April 17 16:01:09 2009
* 
* @brief Core routines for snp-summary. By separating the
* command line parsing form the computation we allow a GUI application to share
* the core computation once setting up the option class. 
*/
#ifndef _SNPSUMMARYREPORTER_H_
#define _SNPSUMMARYREPORTER_H_

//
#include "file/SimpleBinaryFile.h"
#include "calvin_files/portability/src/AffymetrixBaseTypes.h"
#include "util/SQLite.h"
//
#include <string>
#include <vector>

using namespace std;

class SnpSummaryReporter
{
public:
	SnpSummaryReporter(void);
	~SnpSummaryReporter(void);

bool CreateSnpSummaryReport(const string &outputFile,  const vector<string> &chpFiles, 
		const std::string &outputFormat, int &bufferSize, const std::string &callFile = "");

private:
	 SimpleBinaryFile binFile;
	 std::ofstream textFile;
	 SQLiteDatabase db; 

	 int ntotalsets;
	 std::wstring arrayType;	

	 std::vector<std::string> GetMetricNames();
	 void InitializeOutputFile(std::vector<std::string> &colNames, const std::string &outputFormat, const std::string &outputFile);
	 bool GetArrayInformation(const std::vector<std::string> &chps);
	 void CollectData(const std::vector<std::string> &chps, int startIndex, int numSets, std::vector<std::vector<u_int8_t> > &calls, std::vector<std::string> &snpids);
	 void AllocateDataArrays(int nFiles, int nSets, vector<std::vector<u_int8_t> > &calls, std::vector<std::vector<float> > &metrics, int nMetrics, std::vector<string> &snpids);
	 void InitializeCallsAndMetrics(vector<vector<u_int8_t> > &calls, vector<vector<float> > &metrics);
	 void CloseReportFile(const std::string& outputFormat);
	 char* ToStringSQLFormat(char* szBuffer, const float val, bool bCommas);
	 char* ToString(char* szBuffer, const float val);


	 bool isNaN(const float& f) {return f != f;}
	
	


};

#endif

