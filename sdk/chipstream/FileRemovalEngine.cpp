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
#include "chipstream/FileRemovalEngine.h"
//
#include "chipstream/EngineUtil.h"
#include "util/Fs.h"
#include "util/PgOptions.h"
//
#include <cstring>
#include <string>
#include <vector>

using namespace std;
using namespace affx;

FileRemovalEngine::Reg FileRemovalEngine::reg;

FileRemovalEngine * FileRemovalEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == FileRemovalEngine::EngineName())
		return (FileRemovalEngine *)engine;
	return NULL;
}

FileRemovalEngine::FileRemovalEngine() {
    defineOptions();
}

FileRemovalEngine::~FileRemovalEngine() {
}

void FileRemovalEngine::defineOptions() {
	defineOption("", "files-list", PgOpt::STRING_OPT,
		"Text file specifying files to process, one per line with the first line being 'files'.",
		"");
	defOptMult("", "files", PgOpt::STRING_OPT,
		"The files to process.",
		"");
	defineOption("", "delete-listed-files", PgOpt::STRING_OPT,
		"Should the engine delete the listed files.  This option allows the user to choose.",
		"yes");
}

void FileRemovalEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void FileRemovalEngine::checkOptionsImp()
{
    defineStates();

	// Read in list files list from other file if specified.
	vector<string> files;
	string fileName = getOpt("files-list");
	if (fileName.empty() == false)
	{
		affx::TsvFile tsv;
#ifdef WIN32
		tsv.m_optEscapeOk = false;
#endif
		string file;
		tsv.bind(0, "file", &file, TSV_BIND_REQUIRED);
		if(tsv.open(fileName) != TSV_OK) {
			Err::errAbort("Couldn't open files file: " + fileName);
		}
		tsv.rewind();
		while(tsv.nextLevel(0) == TSV_OK) {
			files.push_back(file);
		}
		tsv.close();
		Verbose::out(1, "Read " + ToStr(files.size()) + " files from: " + Fs::basename(fileName));
	}
	else
	{
		files = getOptVector("files");
	}
	setOpt("files", files);
}

static vector<string> split(const string &inputString)
{
	vector<string> tokens;
	size_t substrBegin = 0;
	for (;;)
	{
		size_t substrEnd = inputString.find ("?", substrBegin);
		if (substrEnd == string::npos)
		{
			// No more '?' - save what's left, quit.
			string subString = inputString.substr (substrBegin);
			// Avoid returning a null string from a terminating '?' or an empty inputString.
			if (! subString.empty())
				tokens.push_back (subString);
			break;
		}
		// Avoid null strings from an initial '?' or '??'.
		if (substrEnd != substrBegin)
			tokens.push_back (inputString.substr (substrBegin, substrEnd - substrBegin) );
		// Continue following the '?'
		substrBegin = substrEnd + 1;
	}
	return tokens;
}

/**
   This is the "main()" equivalent of the engine.
*/
void FileRemovalEngine::runImp()
{
	string delete_listed_files = getOpt("delete-listed-files");

	if (delete_listed_files.empty() == false && (delete_listed_files[0] == 'y' || delete_listed_files[0] == 'Y'))
	{
		// Read the input files and delete them.
		vector<string> files = getOptVector("files");
		for (vector<string>::iterator it=files.begin(); it!=files.end(); it++)
		{
			vector<string> tokens = split(*it);
			for (vector<string>::iterator tokenIt=tokens.begin(); tokenIt!=tokens.end(); tokenIt++)
				Fs::rm(tokenIt->c_str(), false);
		}
	}
}
