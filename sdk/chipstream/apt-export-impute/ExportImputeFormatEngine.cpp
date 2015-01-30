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
#include "chipstream/apt-export-impute/ExportImputeFormatEngine.h"
//
#include "calvin_files/portability/src/AffymetrixBaseTypes.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "file/TsvFile/TsvFile.h"
#include "util/AffxConv.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <fstream>
#include <list>
#include <stdio.h>
#include <string>

using namespace std;
using namespace affx;
using namespace affymetrix_calvin_utilities;

#define PROBE_SET_ID_COL_NAME "probeset_id"

ExportImputeFormatEngine::Reg ExportImputeFormatEngine::reg;

/*
 * The information for each marker in the annotation file.
 */
typedef struct _AnnotationEntryType
{
	string name;
	string rsId;
	int pos;
	char chr;
	char alleleA;
	char alleleB;
	int operator<(const _AnnotationEntryType &rhs) const
	{
		if (chr == rhs.chr)
			return (pos < rhs.pos);
		else
			return (chr < rhs.chr);
	}
} AnnotationEntryType;

/*
 * Compliments the base call
 */
static char Compliment(char call)
{
	if (call == 'a')
		return 't';
	else if (call == 'A')
		return 'T';
	else if (call == 'c')
		return 'g';
	else if (call == 'C')
		return 'G';
	else if (call == 'g')
		return 'c';
	else if (call == 'G')
		return 'C';
	else if (call == 't')
		return 'a';
	else if (call == 'T')
		return 'A';
	return 'N';
}

/*
 * Convert pointer
 */
ExportImputeFormatEngine * ExportImputeFormatEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == ExportImputeFormatEngine::EngineName())
		return (ExportImputeFormatEngine *)engine;
	return NULL;
}

/*
 * Constructor
 */
ExportImputeFormatEngine::ExportImputeFormatEngine()
{
    defineOptions();
}

/*
 * Destructor
 */
ExportImputeFormatEngine::~ExportImputeFormatEngine()
{
}

/*
 * Defines the options for the engine.
 */
void ExportImputeFormatEngine::defineOptions()
{
	defineOption("t", "txt-files", PgOpt::STRING_OPT,
		"Text file specifying input TXT files to process.",
		"");
	defOptMult("", "files", PgOpt::STRING_OPT,
		"The text files to process.",
        "");
	defineOption("a", "annotation-file", PgOpt::STRING_OPT,
		"The annotation file name.",
		"");
}

/*
 * Defines the states for the engine.
 */
void ExportImputeFormatEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void ExportImputeFormatEngine::checkOptionsImp()
{
    defineStates();

	string file = getOpt("annotation-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an annotation file.");

	file = getOpt("txt-files");
	vector<string> files;
	if (file.empty() == false)
	{
		affx::TsvFile tsv;
#ifdef WIN32
		tsv.m_optEscapeOk = false;
#endif
		string txtFile;
		tsv.bind(0, "file", &txtFile, TSV_BIND_REQUIRED);
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
			files.push_back(txtFile);
		}
		tsv.close();
		Verbose::out(1, "Read " + ToStr(files.size()) + " files from: " + Fs::basename(file));
	}
	else
	{
		files = getOptVector("files");
	}
	setOpt("files",files);
}

/*
 * Convert the chromosome string to an integer.
 */
static char ChrToInt(const string &chr)
{
	if (chr == "X")
		return 125;
	else if (chr == "Y")
		return 126;
	else if (chr == "MT")
		return 127;
	else
		return atoi(chr.c_str());
}

static string IntToChr(char c)
{
	if (c == 125)
		return "X";
	else if (c == 126)
		return "Y";
	else if (c == 127)
		return "MT";
	else
	{
		char buf[10];
		sprintf(buf, "%d", (int) c);
		return buf;
	}
}

/*
 * Cache the annotations in a sorted list.
 */
static list<AnnotationEntryType> CacheAnnotations(const string &annotationFileName)
{
	TsvFile tsv;
	if (tsv.open(annotationFileName) != TSV_OK)
		Err::errAbort("Unable to open the TSV file");

	list<AnnotationEntryType> annotations;
	string chr;
	string strand;
	string alleleA;
	string alleleB;
	string position;
	AnnotationEntryType entry;
	tsv.bind(0, "Probe Set ID", &entry.name, TSV_BIND_REQUIRED);
	tsv.bind(0, "dbSNP RS ID", &entry.rsId, TSV_BIND_REQUIRED);
	tsv.bind(0, "Chromosome", &chr, TSV_BIND_REQUIRED);
	tsv.bind(0, "Physical Position", &position, TSV_BIND_REQUIRED);
	tsv.bind(0, "Strand", &strand, TSV_BIND_REQUIRED);
	tsv.bind(0, "Allele A", &alleleA, TSV_BIND_REQUIRED);
	tsv.bind(0, "Allele B", &alleleB, TSV_BIND_REQUIRED);
	while (tsv.nextLine() == TSV_OK)
	{
		if (position.empty() == true)
			continue;
		if (position[0] == '-')
			continue;

		entry.pos = getInt(position);
		entry.chr = ChrToInt(chr);
		entry.alleleA = alleleA[0];
		entry.alleleB = alleleB[0];
		if (strand == "-")
		{
			entry.alleleA = Compliment(entry.alleleA);
			entry.alleleB = Compliment(entry.alleleB);
		}
		annotations.push_back(entry);
	}
	tsv.close();
	annotations.sort();
	return annotations;
}

/*
 * Create a mapping of probe set id to line number.
 */
static map<string, linenum_t> CacheMarkerPositions(TsvFile &tsv)
{
	map<string, linenum_t> positions;
	string name;
	tsv.bind(0, PROBE_SET_ID_COL_NAME, &name, TSV_BIND_REQUIRED);
	while (tsv.nextLine() == TSV_OK)
	{
		if (positions.find(name) == positions.end())
			positions[name] = tsv.lineNum();
	}
	tsv.unbindAll();
	tsv.rewind();
	return positions;
}

/*
 * Merge the annotation and data in the IMPUTE format.
 */
static void MergeAndExportData(const string &outputFileBaseName, TsvFile &callTsv, TsvFile &confTsv, map<string, linenum_t> &namePositions, const list<AnnotationEntryType> &annotations)
{
	string name;
	u_int8_t lastChr = -1;
	ofstream out;
	int nDataCols;
	vector<string> calls;
	vector<float> confidences;
	bool includeConfScores = (confTsv.is_open() == 1);

	// Bind the columns in the data file
	nDataCols = callTsv.getColumnCount(0) - 1;
	calls.resize(nDataCols);
	callTsv.bind(0, 0, &name, TSV_BIND_REQUIRED);
	for (int i=0; i<nDataCols; i++)
		callTsv.bind(0, i+1, &calls[i],  TSV_BIND_REQUIRED);
	if (includeConfScores == true)
	{
		if (nDataCols != confTsv.getColumnCount(0) - 1)
			Err::errAbort("There is a mismatch in the number of columns in the confidence and call files.");
		confidences.resize(nDataCols);
		for (int i=0; i<nDataCols; i++)
			confTsv.bind(0, i+1, &confidences[i],  TSV_BIND_REQUIRED);
	}

	// Go through the annotations, find the data and export the results
	for (list<AnnotationEntryType>::const_iterator annotIt=annotations.begin(); annotIt!=annotations.end(); annotIt++)
	{
		// Get the next annotation (this is in sorted order)
		const AnnotationEntryType &annot = *annotIt;
		if (namePositions.find(annot.name) == namePositions.end())
			continue;

		// Go to the associated data line
		linenum_t line = namePositions[annot.name];
		callTsv.seekLine(line - 1);
		callTsv.nextLine();
		assert(name == annot.name);
		if (includeConfScores == true)
		{
			confTsv.seekLine(line - 1);
			confTsv.nextLine();
		}

		// Create a new file for each chromosome
		if (out.is_open() == false || lastChr != annot.chr)
		{
			if (out.is_open() == true)
				out.close();
			string fileName = outputFileBaseName + string("chr") + IntToChr(annot.chr) + string(".txt");
			if (FileUtils::Exists(fileName.c_str()) == true)
				out.open(fileName.c_str(), ios::out | ios::app);
			else
				out.open(fileName.c_str(), ios::out);
			if (!out)
				Err::errAbort("Unable to create the output file: " + fileName);
			lastChr = annot.chr;
		}

		// Output the annotation information.
		out << name << " " // probe set name
			<< annot.rsId << " " // rs ID
			<< annot.pos << " " // position
			<< annot.alleleA << " " // Allele A
			<< annot.alleleB << " "; // Allele B

		// Output the genotype calls
		for (int i=0; i<nDataCols; i++)
		{
			if (calls[i] == "AA" || calls[i] == "A" || calls[i] == "0")
			{
				if (includeConfScores == false)
					out << "1 0 0";
				else
				{
					float conf = confidences[i];
					out << 1 - conf << " " << conf << " " << "0";
				}
			}
			else if (calls[i] == "AB" || calls[i] == "1")
			{
				if (includeConfScores == false)
					out << "0 1 0";
				else
				{
					float conf = confidences[i];
					out << 0.5f*conf << " " << 1-conf << " " << 0.5f*conf;
				}
			}
			else if (calls[i] == "BB" || calls[i] == "B" || calls[i] == "2")
			{
				if (includeConfScores == false)
					out << "0 0 1";
				else
				{
					float conf = confidences[i];
					out << "0 " << conf << " " << 1-conf;
				}
			}
			else
			{
				out << "0 0 0";
			}
			if (i < nDataCols - 1)
				out << " ";
		}
		out << endl;
	}
	out.close();
}

/*
 * Determine if confidence files are part of the input files.
 */
static bool IncludeConfidences(const vector<string> &dataFileNames)
{
	bool found = false;
	for (vector<string>::const_iterator fileIt=dataFileNames.begin(); fileIt!=dataFileNames.end() && found == false; fileIt++)
	{
		string baseName = Fs::basename(*fileIt);
		found = (baseName.find(".confidences.") != -1);
	}
	return found;
}

/*
 * Check the order of the input files to make sure calls and confidences are together.
 */
static void CheckOrder(const vector<string> &dataFileNames)
{
	if (dataFileNames.size() % 2 != 0)
		Err::errAbort("There must be the same number of call and confidence files.");

	string callBase;
	string confBase;
	for (int i=0; i<(int)dataFileNames.size(); i++)
	{
		string baseName = Fs::basename(dataFileNames[i]);
		int idx = baseName.find(".confidences.");
		if (idx != -1)
			confBase = baseName.substr(0, idx);
		else
		{
			idx = baseName.find(".calls.");
			if (idx != -1)
				callBase = baseName.substr(0, idx);
			else
				Err::errAbort("The input files must have either 'calls' or 'confidences' as part of the name.");
		}
		if (i%2 == 1)
		{
			if (callBase != confBase)
				Err::errAbort("The call and confidence input files must be paired together.");
			callBase = "";
			confBase = "";
		}
	}
}

/*
 * This is the "main" equivalent
 */
void ExportImputeFormatEngine::runImp()
{
	try
	{
		// Determine if any of the files are confidence files.
		vector<string> dataFileNames = getOptVector("files");
		int nloops = (int) dataFileNames.size();
		bool includeConfidenceValues = IncludeConfidences(dataFileNames);
		if (includeConfidenceValues == true)
		{
			CheckOrder(dataFileNames);
			nloops /= 2;
		}
		Verbose::progressBegin(1, "Exporting data to the IMPUTE format", 2*nloops + 1, 0, 2*nloops + 1);

		// Create a mapping of the annotations.
		string annotationFileName = getOpt("annotation-file");
		list<AnnotationEntryType> annotations = CacheAnnotations(annotationFileName);

		// Update the progress meter
		Verbose::progressStep(1);

		// Export the genotypes to the IMPUTE format.
		for (vector<string>::iterator fileIt=dataFileNames.begin(); fileIt!=dataFileNames.end(); fileIt++)
		{
			// Determine the output file base name
			string callFile;
			string confFile;
			string outputFileBaseName;
			if (includeConfidenceValues == false)
			{
				callFile = *fileIt;
				outputFileBaseName = (*fileIt).substr(0, fileIt->length() - strlen(".calls.txt")) + ".impute.";
			}
			else
			{
				if (fileIt->find(".calls.") != -1)
				{
					outputFileBaseName = (*fileIt).substr(0, fileIt->length() - strlen(".calls.txt")) + ".impute.";
					callFile = *fileIt;
					confFile = *(++fileIt);
				}
				else
				{
					outputFileBaseName = (*fileIt).substr(0, fileIt->length() - strlen(".confidences.txt")) + ".impute.";
					confFile = *fileIt;
					callFile = *(++fileIt);
				}
			}
			

			// Open the call file and create a map of name to line number
			TsvFile callTsv;
			if (callTsv.open(callFile) != TSV_OK)
				Err::errAbort("Unable to open the file: " + callFile);
			map<string, linenum_t> namePositions = CacheMarkerPositions(callTsv);

			// Update the progress meter
			Verbose::progressStep(1);

			// Open the confidence file
			TsvFile confTsv;
			if (includeConfidenceValues == true)
			{
				if (confTsv.open(confFile) != TSV_OK)
					Err::errAbort("Unable to open the file: " + confFile);
			}

			// Output the data
			MergeAndExportData(outputFileBaseName, callTsv, confTsv, namePositions, annotations);

			// Update the progress meter
			Verbose::progressStep(1);
		}
		Verbose::progressEnd(1, "The data have been exported.");
	}
    catch(const Except &e)
	{
		Err::errAbort(e.what());
    }
    catch(const std::bad_alloc &e)
	{
		Err::errAbort(e.what());
    }
    catch(const std::exception &e)
	{
		Err::errAbort(e.what());
    }
	catch (...)
	{
		Err::errAbort("Unknown error.");
	}
}
