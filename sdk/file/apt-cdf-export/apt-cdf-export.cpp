////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

//
#include "file/CDFFileData.h"
//
#include "util/AptVersionInfo.h"
#include "util/Err.h"
#include "util/PgOptions.h"
#include "util/md5sum.h"
//
#include <iostream>
#include <map>

// @todo: APT-586: rename this to "apt-cdf-util".

using namespace std;
using namespace affxcdf;


/*! File not specified exception */
#define FILE_NOT_SPECIFIED 0

/*! File not found exception */
#define FILE_NOT_FOUND 1

/*! Unable to read the file. */
#define READ_ERROR 2

/*! A string message corresponding to the error codes. */
static const char *error_string[] =
{
	"The CDF file was not specified",
	"The CDF file does not exist",
	"Unable to read the CDF file"
};

void define_aptcdfexport_options(PgOptions* opts)
{
  opts->setUsage("This program will output the probe set information stored in a CDF file.\n"
                 "usage:\n"
                 "   apt-cdf-export -c myCdfFile.cdf -ps probeset1 -ps probeset2 -ps probeset3 ...\n"
                 "     or\n"
                 "   apt-cdf-export -c myCdfFile.cdf\n"
                 "     for all probesets.");
  // generic options
  opts->defOpt("h", "help", PgOpt::BOOL_OPT,
               "This help message.",
               "false");
  opts->defOpt("", "version", PgOpt::BOOL_OPT,
               "Display version information.",
               "false");

  // Action options
  opts->defOpt("","export",PgOpt::BOOL_OPT,
               "export a CDF file.",
               "false");

  opts->defOpt("","file-md5",PgOpt::BOOL_OPT,
               "compute the MD5 of the entire file.",
               "false");

  // 
  opts->defOpt("c", "cdf-file", PgOpt::STRING_OPT,
               "The full path to the CDF file.",
               "");
  opts->defOptMult("", "ps", PgOpt::STRING_OPT,
                   "The name of a probe set. "
                   "Repeat for each probe set to output, or omit for all probe sets.",
                   "");
  
}

/*! The options for the program. */
typedef struct _Options
{
	string file;			/*! The name of the CDF file. */
	map<string, bool> sets;	/*! The name of probe sets to extract. */
} Options; /*! The options for the program. */


/*! Gets the name of the input CDF file from the command line.
 * @param opts PgOptions object describing the command line arguments.
 * @param options The program options.
 */
static void get_options(PgOptions& opts, Options &options)
{
	// get the required arguments
	options.file = opts.get("cdf-file");
  //
  PgOpt* psOpt = opts.mustFindOpt("ps");
  for (int i=0;i<psOpt->getValueCount();i++) {
    options.sets[psOpt->getValue(i)] = true;
  }
}

/*! Reads the CDF file.
 * @param options The program options.
 * @param cdf The CDF file object to fill.
 * @exception FILE_NOT_FOUND, READ_ERROR
 */
static void read_file(const Options &options, CCDFFileData &cdf)
{
	cdf.SetFileName(options.file.c_str());
	if (cdf.Exists() == false)
		throw FILE_NOT_FOUND;

	if (cdf.Read() == false)
		throw READ_ERROR;
}

/*! Output the probe set information.
 * @param name The probe set name.
 * @param set The probe set information.
 */
static void output_set(const string &name, CCDFProbeSetInformation &set)
{
	CCDFProbeGroupInformation group;
	CCDFProbeInformation cel;

	cout << "#probe_set_name\t" << name << endl;
	cout << "#probe_set_type\t";
	switch (set.GetProbeSetType())
	{
	case affxcdf::ExpressionProbeSetType:
		cout << "expression" << endl;
		break;

	case affxcdf::GenotypingProbeSetType:
		cout << "genotyping" << endl;
		break;

	case affxcdf::ResequencingProbeSetType:
		cout << "resequencing" << endl;
		break;

	case affxcdf::TagProbeSetType:
		cout << "tag" << endl;
		break;

	default:
		cout << "unknown" << endl;
		break;
	}
	cout << "#num_groups\t" << set.GetNumGroups() << endl;
	cout << endl;

	int ng = set.GetNumGroups();
	for (int ig=0; ig<ng; ig++)
	{
		set.GetGroupInformation(ig, group);
		cout << "#group_index\t" << ig << endl;
		cout << "#group_name\t" << group.GetName() << endl;
		cout << "#direction\t" << (group.GetDirection() == affxcdf::AntiSenseDirection ? "Antisense" : "Sense") << endl;
		cout << "#num_cells\t" << group.GetNumCells() << endl;
		
		int nc = group.GetNumCells();
		cout << "#cell_header\t" << "x\t" << "y\t" << "expos\t" << "pbase\t" << "tbase\t" << "atom" << endl;
		for (int ic=0; ic<nc; ic++)
		{
			group.GetCell(ic, cel);
			cout <<
				cel.GetX() << "\t" <<
				cel.GetY() << "\t" <<
				cel.GetExpos() << "\t" <<
				cel.GetPBase() << "\t" <<
				cel.GetTBase() << "\t" <<
				cel.GetListIndex() << endl;
		}
		cout << endl;
	}
}

/*! Determines if the probe set should be output.
 * @param name The probe set name.
 * @param options The program options.
 * @return true if the probe set should be displayed.
 */
static bool output_set(const string &name, const Options &options)
{
	if (options.sets.empty())
		return true;

	map<string, bool>::const_iterator pos = options.sets.find(name);
	return (pos != options.sets.end());
}

/*! Outputs the data to the command line.
 * @param options The program options.
 * @param cdf The CDF file object.
 */
static void output_data(const Options &options, CCDFFileData &cdf)
{
	CCDFProbeSetInformation set;

	// Loop over all probe sets and output the data.
	int ns = cdf.GetHeader().GetNumProbeSets();
	cout << "#num_probe_sets\t" << ns << endl << endl;
	for (int is=0; is<ns; is++)
	{
		string name = cdf.GetProbeSetName(is);

		// Output the set.
		if (output_set(name, options) == true)
		{
			cdf.GetProbeSetInformation(is, set);
			output_set(name, set);
		}
	}
}

/*! The main entry point. This application will extract information from a CDF file
 *  and print the contents to the command line.
 * @param argc The number of command line arguments
 * @param argv The command line arguments
 * @return 0 if successful.
 */
int main(int argc, char ** argv)
{
  try {
    const string version = AptVersionInfo::versionToReport();

	try
	{
		// Get the program options.
		PgOptions opts;
    define_aptcdfexport_options(&opts);
		opts.parseArgv(argv);

		if (opts.getBool("help")||(argc == 1)) {
			opts.usage();
			cout << "version: " << version << endl;
			return 0;
		}

    if (opts.getBool("version")) {
			cout << "version: " << version << endl;
			return 0;
		}

    if (opts.getBool("file-md5")) {
      affx::md5sum md5sum;
      std::string file_path;
      std::string file_md5;
      int err_cnt=0;

      for (int i=0;i<opts.getArgCount();i++) {
        file_path=opts.getArg(i);
        int rv=md5sum.ofFile(file_path,file_md5);
        {
          cout << file_path << ": ";
          if (rv==0) {
            cout << file_md5;
          }
          else {
            cout << " error=" << rv;
            err_cnt++;
          }
          cout << endl;
        }
      }
      return ((err_cnt==0)?0:-1);
    }


    // "--export" option.
		Options options;
		get_options(opts, options);

		// Read the file
		CCDFFileData cdf;
		read_file(options, cdf);

		// Output the data to the command line.
		output_data(options, cdf);
	}
	catch (int err)
	{
        Err::errAbort(error_string[err]);
	}

	return 0;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
