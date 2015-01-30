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
#include "chipstream/apt-export-impute/ExportDataFilesEngine.h"
//
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "file/CHPFileData.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <fstream>
#include <list>
#include <memory>
#include <string>

using namespace std;
using namespace affx;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

#define PROBE_SET_ID_COL_NAME "probeset_id"

ExportDataFilesEngine::Reg ExportDataFilesEngine::reg;

/*
 * Convert pointer
 */
ExportDataFilesEngine * ExportDataFilesEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == ExportDataFilesEngine::EngineName())
		return (ExportDataFilesEngine *)engine;
	return NULL;
}

/*
 * Constructor
 */
ExportDataFilesEngine::ExportDataFilesEngine()
{
    defineOptions();
}

/*
 * Destructor
 */
ExportDataFilesEngine::~ExportDataFilesEngine()
{
}

/*
 * Define the program options.
 */
void ExportDataFilesEngine::defineOptions() {
	defineOption("", "chp-files", PgOpt::STRING_OPT,
		"Text file specifying input CHP files to process.",
		"");
	defOptMult("", "chps", PgOpt::STRING_OPT,
		"chp files to process.",
        "");
	defOptMult("", "cols", PgOpt::STRING_OPT,
		"Columns to process.",
        "");
	defineOption("", "output-file", PgOpt::STRING_OPT,
		"The output file name.",
		"");
	defineOption("", "output-codes", PgOpt::BOOL_OPT,
		"Output the calls as codes (-1=NN, 0=AA, 1=AB, 2=BB).",
		"true");
}

/*
 * Define the states
 */
void ExportDataFilesEngine::defineStates() { }


/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void ExportDataFilesEngine::checkOptionsImp()
{
    defineStates();
	string file = getOpt("output-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an output file.");

	if (getOptVector("cols").size() == 0)
		Err::errAbort("Must specify columns to output.");

	file = getOpt("chp-files");
	vector<string> files;
	if (file.empty() == false)
	{
		affx::TsvFile tsv;
#ifdef WIN32
		tsv.m_optEscapeOk = false;
#endif
		string chpFile;
		tsv.bind(0, "file", &chpFile, TSV_BIND_REQUIRED);
		try 
		{
			if(tsv.open(file) != TSV_OK)
			{
				Err::errAbort("Couldn't open file: " + file);
			}
		}
		catch (Except&)
		{
			tsv.close();	// need to close the file
			throw;
		}
		tsv.rewind();
		while(tsv.nextLevel(0) == TSV_OK) {
			files.push_back(chpFile);
		}
		tsv.close();
		Verbose::out(1, "Read " + ToStr(files.size()) + " files from: " + Fs::basename(file));
	}
	else
	{
		files = getOptVector("chps");
	}
	setOpt("chps",files);
}

/*
 * Write the header of the output file.
 */
static void OutputHeader(ofstream &out, const vector<string> &fnames, const vector<string> &cnames)
{
	out << PROBE_SET_ID_COL_NAME;
	for (int ifile=0; ifile<(int)fnames.size(); ifile++)
	{
		for (int icol=1; icol<(int)cnames.size(); icol++)
		{
			out << "\t" << Fs::basename(fnames[ifile]) << "-" << cnames[icol];
		}	
	}
	out << endl;
}

/*
 * Create a map of column name to index value.
 */
static void CreateColumnIndexMap(DataSet *set, const vector<string> &cnames, map<int, int> &colIdx, map<int, DataSetColumnTypes> &colTypes)
{
	// Create a map of requested columns to actual columns.
	int nc = set->Header().GetColumnCnt();
	for (int ic=0; ic<nc; ic++)
	{
		ColumnInfo ci = set->Header().GetColumnInfo(ic);
		string cname = StringUtils::ConvertWCSToMBS(ci.GetName());
		for (int icn=0; icn<(int)cnames.size(); icn++)
		{
			if (cnames[icn] == cname)
			{
				colIdx[icn] = ic;
				colTypes[icn] = ci.GetColumnType();
				break;
			}
		}
	}
}

/*
 * Get the data from the current row to a buffer.
 */
static vector<string> GetRowData(int nr, int currentRow, int ic, map<int, int> &colIdx, map<int, DataSetColumnTypes> &colTypes, DataSet *set)
{
	vector<string> rowBuffer(nr);
	for (int ir=0; ir<nr; ir++)
	{
		switch (colTypes[ic])
		{
		case ByteColType:
		{
			int8_t bval;
			set->GetData(ir+currentRow, colIdx[ic], bval);
			rowBuffer[ir] = ToStr(bval);
			break;
		}

		case UByteColType:
		{
			u_int8_t ubval;
			set->GetData(ir+currentRow, colIdx[ic], ubval);
			rowBuffer[ir] = ToStr(ubval);
			break;
		}

		case ShortColType:
		{
			int16_t sval;
			set->GetData(ir+currentRow, colIdx[ic], sval);
			rowBuffer[ir] = ToStr(sval);
			break;
		}

		case UShortColType:
		{
			u_int16_t usval;
			set->GetData(ir+currentRow, colIdx[ic], usval);
			rowBuffer[ir] = ToStr(usval);
			break;
		}

		case IntColType:
		{
			int32_t ival;
			set->GetData(ir+currentRow, colIdx[ic], ival);
			rowBuffer[ir] = ToStr(ival);
			break;
		}

		case UIntColType:
		{
			u_int32_t uival;
			set->GetData(ir+currentRow, colIdx[ic], uival);
			rowBuffer[ir] = ToStr(uival);
			break;
		}

		case FloatColType:
		{
			float fval;
			set->GetData(ir+currentRow, colIdx[ic], fval);
			rowBuffer[ir] = ToStr(fval);
			break;
		}

		case ASCIICharColType:
		{
			string sval;
			set->GetData(ir+currentRow, colIdx[ic], sval);
			rowBuffer[ir] = sval;
			break;
		}

		default:
			rowBuffer[ir] = "";
			break;
		}
	}
	return rowBuffer;
}

/*
 * Convert detection codes to string values.
 */
static string ExpressionDetection(char c)
{
	switch (c)
	{
	case ABS_PRESENT_CALL:
		return "P";
	case ABS_MARGINAL_CALL:
		return "M";
	case ABS_ABSENT_CALL:
		return "A";
	default:
		return "No Call";
	}
}

/*
 * Convert genotype call codes to string values.
 */
static string GenotypeCall(char c)
{
	switch (c)
	{
	case ALLELE_A_CALL:
		return "AA";
	case ALLELE_B_CALL:
		return "BB";
	case ALLELE_AB_CALL:
		return "AB";
	default:
		return "No Call";
	}
}

/*
 * Convert genotype call codes (in CHP file) to APT codes.
 */
static string GenotypeCode(char c)
{
	switch (c)
	{
	case ALLELE_A_CALL:
		return "0";
	case ALLELE_AB_CALL:
		return "1";
	case ALLELE_B_CALL:
		return "2";
	default:
		return "-1";
	}
}

/*
 * Write the buffer to the output file.
 */
static void WriteBuffer(ofstream &out, const vector<vector<string> > &buffer, const vector<string> &dataColumnNames, bool outputAsCodes)
{
	// Output the buffer
	if (buffer.size() == 0)
		return;
	int nc = buffer.size();
	int nr = buffer[0].size();
	int np = (int)((nr / 10.f) + 0.5f);
	for (int ir=0; ir<nr; ir++)
	{
		if (ir % np == 0)
			Verbose::progressStep(1);

		// Write the data
		for (int ic=0; ic<nc; ic++)
		{
			const string &colName = (ic == 0 ? "" : dataColumnNames[((ic-1) % (dataColumnNames.size()-1)) + 1]);
			if (colName == "Detection")
			{
				out << ExpressionDetection(buffer[ic][ir][0]);
			}
			else if (colName == "Call" || colName == "Forced Call")
			{
				if (outputAsCodes == false)
					out << GenotypeCall(buffer[ic][ir][0]);
				else
					out << GenotypeCode(buffer[ic][ir][0]);
			}
			else
				out << buffer[ic][ir];
			if (ic < nc-1)
				out << "\t";
		}
		out << endl;
	}
}

/*
 * Determine the number of lines to read in the cache (20MB)
 */
static int DetermineCacheSize(int nfiles, const map<int, DataSetColumnTypes> &colTypes)
{
	uint64_t bytesPerRow = 0;
	uint64_t maxBufferSize = 20 * 1024 * 1024; // 20MB
	uint64_t freeRam = 0, totalRam = 0, swapAvail = 0, memAvail = 0;
	Util::memInfo(freeRam, totalRam, swapAvail, memAvail, false);
	for (map<int, DataSetColumnTypes>::const_iterator it=colTypes.begin(); it!=colTypes.end(); it++)
	{
		switch (it->second)
		{
		case ByteColType:
		case UByteColType:
			bytesPerRow += 1;
			break;

		case ShortColType:
		case UShortColType:
			bytesPerRow += 2;
			break;

		case IntColType:
		case UIntColType:
		case FloatColType:
			bytesPerRow += 4;
			break;

		case ASCIICharColType:
			bytesPerRow += (4 + 32); // assuming a 32 character limit.
			break;

		default:
			break;
		}
	}
	return (int) (min(memAvail, maxBufferSize) / (nfiles * bytesPerRow));
}

/*
 * This is the "main" equivalent
 */
void ExportDataFilesEngine::runImp()
{
	// Create the output file.
	bool outputAsCodes = getOptBool("output-codes");
	string outputFile = getOpt("output-file");
	ofstream out(outputFile.c_str(), ios::out);
	if (!out)
		Err::errAbort("Unable to open the output file");

	// Output the header
	vector<string> dataFileNames = getOptVector("chps");
	vector<string> dataColumnNames = getOptVector("cols");
	OutputHeader(out, dataFileNames, dataColumnNames);

	try
	{
		// Output the data
		vector<vector<string> > buffer;
		map<int, int> colIdx;
		map<int, DataSetColumnTypes> colTypes;
		int currentRow = 0;
		int numRows = 0;
		int rowCacheSize = -1;
		int nr = -1;
		bool continueParsing=true;
		while (continueParsing)
		{
			buffer.clear();
			buffer.resize(1 + ((dataColumnNames.size()-1) * dataFileNames.size()));
			for (int ifile=0; ifile<(int)dataFileNames.size(); ifile++)
			{
				// Open the file.
                std::auto_ptr<FusionCHPData> chpData;
				chpData.reset(FusionCHPDataReg::Read(dataFileNames[ifile]));
				if (chpData.get()  == NULL)
					Err::errAbort("Failed to open the data file: " + dataFileNames[ifile]);
				DataSet *set = chpData->GetGenericData()->DataSet(0, 0);
				set->Open();

				// Create a map of requested columns to actual columns.
				if (colIdx.empty() == true)
					CreateColumnIndexMap(set, dataColumnNames, colIdx, colTypes);

				// Determine the buffer row count
				numRows = set->Header().GetRowCnt();
				if (rowCacheSize == -1)
				{
					rowCacheSize = DetermineCacheSize(dataFileNames.size(), colTypes);
					nr = min(rowCacheSize, numRows - currentRow);
					int nloops = (numRows / nr) + 1;
					int steps = dataFileNames.size()*nloops + 10*nloops;
					Verbose::progressBegin(1, "Exporting data and annotations", steps, 0, steps);
				}
				nr = min(rowCacheSize, numRows - currentRow);
				Verbose::progressStep(1);

				// Read the columns.
				if (ifile == 0)
				{
					set->GetData(colIdx[0], currentRow, nr, buffer[0]);
				}
				for (int ic=1; ic<(int)dataColumnNames.size(); ic++)
				{
					int bidx = (dataColumnNames.size()-1)*ifile + ic;
					buffer[bidx] = GetRowData(nr, currentRow, ic, colIdx, colTypes, set);
				}

				// Destry the data set (which will also close it).
				set->Delete();

				// Increment the row index when on the last file. Also check if
				// all rows have been processed.
				if (ifile == (int)dataFileNames.size() - 1)
				{
					currentRow += rowCacheSize;
					if (currentRow >= numRows)
					{
						continueParsing = false;
						break;
					}
				}
			}

			// Output the buffer
			WriteBuffer(out, buffer, dataColumnNames, outputAsCodes);
		}
		Verbose::progressEnd(1, "The data and annotations have been exported to: " + outputFile);
		out.close();
	}
    catch(const Except &e)
	{
		out.close();
		FileUtils::RemoveFile(outputFile.c_str());
		Err::errAbort(e.what());
    }
    catch(const std::bad_alloc &e)
	{
		out.close();
		FileUtils::RemoveFile(outputFile.c_str());
		Err::errAbort(e.what());
    }
    catch(const std::exception &e)
	{
		out.close();
		FileUtils::RemoveFile(outputFile.c_str());
		Err::errAbort(e.what());
    }
	catch (...)
	{
		out.close();
		FileUtils::RemoveFile(outputFile.c_str());
		Err::errAbort("Unknown error.");
	}
}
