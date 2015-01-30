////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

#include "util/AptVersionInfo.h"
#include "util/md5sum.h"
#include "util/PgOptions.h"
#include "util/Verbose.h"
#include "file/TsvFile/TsvFile.h"
#include "util/FsPath.h"
#include "util/Fs.h"

using namespace std;
using namespace affx;

/**
 * Define the options for the program.
 * @param opts The command line options.
 */
static void define_options(PgOptions* opts)
{
	opts->setUsage("apt-md5-check - check the integrity of files given an MD5 value file.\n"
		"The MD5 value file is a tab delimited text file with two columns: \"File\" and \"MD5\"\n\n"
		"usage:\n"
		"   apt-md5-check -a <application> -i <md5-file.rpt>");
	opts->defineOption("h", "help", PgOpt::BOOL_OPT,
		"This message.",
		"false");
	opts->defineOption("i", "", PgOpt::STRING_OPT,
		"The MD5 value file.",
		"");
	opts->defineOption("a", "", PgOpt::STRING_OPT,
		"The name of the application that generated the MD5 value file.",
		"");
	opts->defineOption("s", "", PgOpt::STRING_OPT,
		"Text to add to the file contents when calculating MD5 values.",
		"");
	opts->defineOption("v", "", PgOpt::INT_OPT,
		"The verbosity: 0 = normal, 2 = detailed messages",
		"0");
	opts->defineOption("", "version", PgOpt::BOOL_OPT,
		"Display version information.",
		"false");
}

/** 
 * Determine the salt value based on the command line inputs and known application salts
 * @param opts The command line options
 * @return The salt value to use when computing the MD5.
 */
static string determineSalt(PgOptions *opts)
{
	string salt = opts->get("s");
	if (salt.empty() == true)
	{
		string app = opts->get("a");
		if (app == "DMET Console")
			salt = "salt";
	}
	return salt;
}

/**
 * Check the MD5 values stored in the input file versus a newly computed value.
 * @param opts The command line options
 * @return The number of files that do not match the MD5.
 */
static int checkFiles(PgOptions* opts)
{
	int nerrors = 0;
	string salt=determineSalt(opts);
	string file;
	string md5value;
	md5sum md5;
	string fileMd5;

	string md5File = opts->get("i");
	if (Fs::exists(md5File) == false)
		throw string("The input file does not exist");
	FsPath base(md5File);
	string basePath = base.getDirName();

	Verbose::setLevel(opts->getInt("v"));

	TsvFile tsv;
	if (tsv.open(md5File) != affx::TSV_OK) {
		Err::errAbort("Unable to open the MD5 file.");
	}
	tsv.bind(0, "File", &file, affx::TSV_BIND_REQUIRED);
	tsv.bind(0, "MD5", &md5value, affx::TSV_BIND_REQUIRED);

	while (tsv.nextLine() == TSV_OK)
	{
		Verbose::out(2, "Checking: " + file);
		Verbose::out(2, "MD5: " + md5value);

		md5.ofFile(Fs::join(basePath, file), fileMd5, salt);
		Verbose::out(2, "Calculated MD5: " + fileMd5);
		if (fileMd5 != md5value)
		{
			++nerrors;
			Verbose::out(0, "The MD5 value does not match for: " + file);
		}
		else
		{
			Verbose::out(0, "The MD5 value matches for: " + file);
		}

		Verbose::out(2, "");
	}

	tsv.close();

	return nerrors;
}

/** 
 * The main entry point.
 * @return The number of files that do not match the input MD5 values.
 */
int main(int argc, char* argv[])
{
	int nerrors = 0;
	try
	{
		const string version = AptVersionInfo::versionToReport();

		PgOptions *opts = new PgOptions();
		define_options(opts);
		opts->parseArgv(argv);

		// Do we need help?
		if(opts->getBool("help") || argc == 1)
		{
			opts->usage();
			cout << "version: " << version << endl;
		}
		// Show version number
		else if(opts->getBool("version"))
		{
			cout << "version: " << version << endl;
		}
		// Check the MD5 values.
		else
		{
			nerrors = checkFiles(opts);
		}
		delete opts;
	}
	catch(string err)
	{
		Verbose::out(1, "ERROR: " + err);
		nerrors = -1;
	}
	catch(...)
	{
		Verbose::out(1,"Unexpected Error: uncaught exception.");
		nerrors = -1;
	}
	return nerrors;
}
