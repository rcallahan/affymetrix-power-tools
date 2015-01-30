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

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//
#include "chipstream/apt-snp-summary/SnpSummaryReporter.h"
//
#include "SnpSummaryReporter.h"
//
#include "SnpSummaryStats.h"
#include "chipstream/BioTypes.h"
#include "calvin_files/fusion/src/FusionCHPData.h"
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/BioTypes.h"
#include "chipstream/apt-snp-summary/SnpSummaryStats.h"
#include "file/CHPFileData.h"
#include "file/SimpleBinaryFile.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/SQLite.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <iomanip>
#include <limits>

using namespace std;
using namespace affxchp;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;


SnpSummaryReporter::SnpSummaryReporter(void)
{
  ntotalsets = 0;
}

SnpSummaryReporter::~SnpSummaryReporter(void)
{
}

/*
 * Convert a float to a string for text output.
 */
char* SnpSummaryReporter::ToString(char* szBuffer, const float val)
{
	if(isNaN(val))
        strcpy(szBuffer, "");
    else
        sprintf(szBuffer, "%.7g", val);
	return szBuffer;
}

/*
 * Convert a float to a string with SQL formating.
 */
char* SnpSummaryReporter::ToStringSQLFormat(char* szBuffer, const float val, bool bCommas)
{
	if (bCommas) {
	    if(isNaN(val))
            strcpy(szBuffer, "null,");
        else
            sprintf(szBuffer, "%.7g,", val);
	} else {
	    if(isNaN(val))
            strcpy(szBuffer, "null)");
        else
            sprintf(szBuffer, "%.7g)", val);
	}
	return szBuffer;
}


/*
 * Allocate memory for the call buffer
 */
void SnpSummaryReporter::AllocateDataArrays(int nFiles, int nSets, vector<vector<u_int8_t> > &calls, vector<vector<float> > &metrics, int nMetrics, vector<string> &snpids)
{
	try
	{
		snpids.resize(nSets);
		calls.resize(nSets);
		metrics.resize(nSets);
		for (int iset=0; iset<nSets; iset++)
			calls[iset].resize(nFiles);
		for (vector<vector<float> >::iterator metricIt=metrics.begin(); metricIt!=metrics.end(); metricIt++)
			metricIt->resize(nMetrics);

	}
	catch(...)
	{
		throw (string) "FATAL ERROR: Insufficient memory.";
	}
}

/*
 * Set all of the calls and metrics to 0
 */
void SnpSummaryReporter::InitializeCallsAndMetrics(vector<vector<u_int8_t> > &calls, vector<vector<float> > &metrics)
{
    for (vector<vector<u_int8_t> >::iterator setIt=calls.begin(); setIt!=calls.end(); ++setIt)
    {
        for (vector<u_int8_t>::iterator fileIt=setIt->begin();  fileIt!=setIt->end(); ++fileIt)
            *fileIt=0;
    }
    for (vector<vector<float> >::iterator metricIt=metrics.begin(); metricIt!=metrics.end(); metricIt++)
    {
        for (vector<float>::iterator setIt=metricIt->begin(); setIt!=metricIt->end(); setIt++)
            *setIt = 0.0f;
    }
}

/*
 * Extract the data from the CHP file.
 */
void SnpSummaryReporter::CollectData(const vector<string> &chps, int startIndex, int numSets, vector<vector<u_int8_t> > &calls, vector<string> &snpids)
{
	FusionCHPData *chp = NULL;
	FusionCHPMultiDataData *dchp = NULL;

	for (int ichp=0; ichp<(int)chps.size(); ichp++)
    {
      std::string chp_filename=Fs::convertToUncPath(chps[ichp]);
		try
		{
			chp = FusionCHPDataReg::ReadHeader(chp_filename);
			dchp = FusionCHPMultiDataData::FromBase(chp);
		}
		catch(...)
		{
			if(dchp != NULL) {
			    Freez(dchp);
                chp = NULL;
            }
            else
                Freez(chp);
      std::string message = "FATAL ERROR: Unable to read the file: "+FS_QUOTE_PATH(chp_filename);
			throw (string) message;
		}

        if (dchp == NULL)
		{
			Freez(chp);
      std::string message = "Unable to read the file: "+FS_QUOTE_PATH(chp_filename);
			throw (string) message;
		}
		if (startIndex == 0)
		{
			if (arrayType != dchp->GetArrayType())
				throw (string) "FATAL ERROR: There was more than one array type included in the list of chp files.";
		}

        for (int iset=0; iset<numSets; iset++)
        {
            calls[iset][ichp] = dchp->GetGenoCall(GenotypeMultiDataType, iset+startIndex);
            if (ichp == 0)
                snpids[iset] = dchp->GetProbeSetName(GenotypeMultiDataType, iset+startIndex);
        }

        Freez(dchp); // this will free chp
        chp = NULL;
    }
}

/*
 * Get the metric names
 */
vector<string> SnpSummaryReporter::GetMetricNames()
{
    vector<string> names;
    names.push_back("SNP Call Rate");
    names.push_back("SNP %AA");
    names.push_back("SNP %AB");
    names.push_back("SNP %BB");
    names.push_back("Minor Allele Frequency");
    names.push_back("H.W. p-Value");
    return names;
}

/*
 * Get the number of probe sets in the files. Only from the first
 * file of each array.
 */
bool SnpSummaryReporter::GetArrayInformation(const std::vector<std::string> &chps)
{
	FusionCHPData *chp = NULL;
	FusionCHPMultiDataData *dchp = NULL;
	bool status = true;
  std::string chp_filename=Fs::convertToUncPath(chps[0]);
	try
	{
		chp = FusionCHPDataReg::ReadHeader(chp_filename);
		dchp = FusionCHPMultiDataData::FromBase(chp);
		if (dchp != NULL)
		{
			ntotalsets = dchp->GetEntryCount(GenotypeMultiDataType);
			arrayType = dchp->GetArrayType();
            Freez(dchp); // this will free chp
            chp = NULL;
		}
		else
		{
			status = false;
            Freez(chp); // this will free chp
		}
	}
	catch(...)
	{
		if(dchp != NULL) {
		    Freez(dchp);
            chp = NULL;
        }
        else
            Freez(chp);
		status = false;
	}
	return status;
}

/*
 * Open the output file and write the header
 */
void SnpSummaryReporter::InitializeOutputFile(std::vector<std::string> &colNames, 
                                              const std::string& outputFormat,
                                              const std::string &outputFile)
{
  std::string tmp_unc_path;
  tmp_unc_path=Fs::convertToUncPath(outputFile);

	try
	{
		if (outputFormat == "db")
		{
			db.open(tmp_unc_path);
			db.beginTransaction();
			db.execute("delete from SnpSummary", true);
			db.commitTransaction();
			db.beginTransaction();
		}
		else if (outputFormat == "txt")
		{
      // @todo use TsvFile.
			Fs::mustOpenToWrite(textFile, tmp_unc_path);
			textFile << "# Array Type" << "\t" << StringUtils::ConvertWCSToMBS(arrayType) << endl;;
			textFile <<  "SNPID" << "\t";
			for (size_t i = 0; i<colNames.size(); i++) {
				textFile << colNames[i] << "\t";
      }
			textFile << endl;
			textFile.setf(ios::fixed, ios::floatfield);
		}
		else
		{
 			if (binFile.CreateNewFile(tmp_unc_path) == false)
				Err::errAbort("Couldn't open file: " + tmp_unc_path + " to write.");
			int maxSnpIdLength = 24;
			binFile.WriteHeader(ntotalsets, (int) colNames.size(), "SNPID", maxSnpIdLength, colNames, colNames);
		}
	}
	catch(...)
	{
		Err::errAbort("Couldn't open file: "+FS_QUOTE_PATH(tmp_unc_path)+" to write.");
	}
}

/*
 * Create the summary report
 */
bool SnpSummaryReporter::CreateSnpSummaryReport(const string &outputFile,
                                                const vector<string> &chpFiles,
                                                const std::string &outputFormat, 
                                                int &bufferSize,
                                                const string &callFile)
{

	bool status = false;
	vector<string> colNames = GetMetricNames();
	const vector<string> &chps = chpFiles;
	ntotalsets = 0;

	// get number of probesets and array type
	if (callFile != "") {
		if (outputFormat != "db" && outputFormat != "txt") {
            Err::errAbort("Cannot use bin format will table file input.");
        }
    } else {
        if((GetArrayInformation(chps) == false) || (ntotalsets == 0))
	    {
        std::string msg = "Fatal Error:  Unable to read " + Fs::convertToUncPath(chpFiles[0]);
		    Verbose::out(1, msg);
        Verbose::out(1,"Note: SNP Summary reports can only be created for CHP files in AGCC format.");
		    return false;
	    }
    }

	 // Open and Write the header of the file
	InitializeOutputFile(colNames, outputFormat, outputFile);
	char szBuffer[64];
	if(callFile!="") {
        affx::TsvFile tsv;
#ifdef WIN32
        tsv.m_optEscapeOk = false;
#endif
        std::string snpId;
        tsv.bind(0, "probeset_id", &snpId, affx::TSV_BIND_REQUIRED);
        if (tsv.open(callFile) != affx::TSV_OK) {
            Err::errAbort("Couldn't open call-file: " + callFile);
        }
        tsv.rewind();
        int numCol = tsv.getColumnCount(0);
        while (tsv.nextLevel(0) == affx::TSV_OK) {
	        SnpSummaryStats stats;
            vector<vector<u_int8_t> > calls(1);
            vector<vector<float> > metrics(1);
            metrics[0].resize(colNames.size());

            for(int i=1; i<numCol; i++) {
                int val;
                if(tsv.get(0,i,val)==affx::TSV_OK)
                    calls[0].push_back((u_int8_t)affx::GTypeCallToChpValue(affx::GType_from_int(val)));
            }

            stats.CalculateMetrics(calls, 1, metrics);

		    if (outputFormat == "txt")
		    {
			    textFile  << snpId << "\t";
			    for (size_t i = 0; i< metrics[0].size(); i++)
				    textFile << ToString(szBuffer, metrics[0][i]) << "\t";
			    textFile << endl;
		    }
		    else
			    binFile.WriteRow(snpId, metrics[0]);
        }
        tsv.close();
        status = true;
    }
	 else
	 {
		// Start a progress meter
		int nCachedSets = bufferSize;
		int progressLimit = (ntotalsets / nCachedSets) + 1;
		Verbose::progressBegin(1, "Creating SNP summary file", progressLimit, 0, progressLimit);

		try
		{
			// Allocate memory for the results.
			vector<string> snpids;
			vector<vector<u_int8_t> > calls;
			vector<vector<float> > metrics;

			// compute the stats.
			int nfiles = (int) chpFiles.size();
			AllocateDataArrays(nfiles, nCachedSets, calls, metrics, (int) colNames.size(), snpids);

			// Calculate the metrics and write them to file
			int iset = 0;

			SnpSummaryStats stats;
			int size = metrics[0].size() -1;
			int stringSize = (metrics[0].size() * 8) + 80;


			while (iset < ntotalsets)
			{
				Verbose::progressStep(1);
				int nProcessSets = min(nCachedSets, ntotalsets -iset);
				InitializeCallsAndMetrics(calls, metrics);
				CollectData(chps, iset, nProcessSets, calls, snpids);
				stats.CalculateMetrics(calls, nProcessSets, metrics);
				std::string strSQL;

				if (outputFormat == "db")
				{
					for (int imetric=0; imetric<nProcessSets; imetric++)
					{
						strSQL.resize(stringSize);
						strSQL.clear();
						strSQL = "insert into SNPSummary values ('" + snpids[imetric] + "',";
						for (int i = 0; i<size; i++)
						{
							strSQL.append(ToStringSQLFormat(szBuffer, metrics[imetric][i], true));
						}
						strSQL.append(ToStringSQLFormat(szBuffer, metrics[imetric][size], false));
						db.execute(strSQL);
					}
				}
				else if (outputFormat == "txt")
				{
					for (int imetric=0; imetric<nProcessSets; imetric++)
					{
						textFile  << snpids[imetric] << "\t";
						for (size_t i = 0; i< metrics[imetric].size(); i++)
						{
							textFile << ToString(szBuffer, metrics[imetric][i]) << "\t";
						}
						textFile << endl;
					}
				}
				else
				{
					for (int imetric=0; imetric<nProcessSets; imetric++)
					{
						binFile.WriteRow(snpids[imetric], metrics[imetric]);
					}
				}
			  iset += nProcessSets;
			}

		status = true;
		}
		catch(string err)
		{
			Verbose::out(1, "\n" + err);
			status = false;
		}
		catch(...)
		{
			Verbose::out(1,"FATAL ERROR: An error was detected with writing the file.");
			status = false;
		}

	CloseReportFile(outputFormat);

    if (status == false)
	{
        FileUtils::RemoveFile(outputFile.c_str());
		Verbose::progressEnd(1, "Summary file was not created successfully");
	}
	else
		Verbose::progressEnd(1, "Summary file created");
	}

    return status;
}

//close the report file
void SnpSummaryReporter::CloseReportFile(const std::string& outputFormat)
{
	try
	{
		if (outputFormat == "db")
		{
			db.commitTransaction();
			db.close();
		}
		else if (outputFormat == "txt")
			textFile.close();
		else
			binFile.Close();
	}
	catch(...)
	{
	}
}


