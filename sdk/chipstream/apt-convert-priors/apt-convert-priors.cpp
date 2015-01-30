////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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

// @file   apt-convert-priors
// @author Harley Gorrell
// @brief  Convert priors from v1 to v2 (and back)

#include "chipstream/QuantLabelZIO.h"
//
#include "util/PgOptions.h"
//

//////////

void define_options(PgOptions* opts)
{
  opts->setUsage("apt-convert-priors --\n"
                 "\n"
                 );
  //
  opts->defOpt("h","help",PgOpt::BOOL_OPT,
               "Print this message.",
               "false");
  opts->defOpt("v","verbose",PgOpt::INT_OPT,
               "Print this message.",
               "1");
  //
  opts->defOpt("f","format-ver",PgOpt::INT_OPT,
               "output format version (1 or 2).",
               "2");
  opts->defOpt("o","output",PgOpt::STRING_OPT,
               "Output file.",
               "");
  opts->defOpt("","tsv-name",PgOpt::STRING_OPT,
               "Name of tsv5 object.",
               "");

}

int main(int argc,const char* argv[]) {
  try {
    PgOptions* opts;

    // Create our option parser, define the options and parse.
    opts = new PgOptions;
    define_options(opts);
    opts->parseArgv(argv);
    
    // Print our help message if necessary.
    if (opts->getBool("help")) {
        opts->usage();
    }
    
    //
    int verbose=opts->getInt("verbose");
    int format_ver=opts->getInt("format-ver");
    std::string file_in;
    std::string tsv_name=opts->get("tsv-name");
    std::string file_out=opts->get("output");
    
    //
    if (!((format_ver==1)||(format_ver==2))) {
        Err::errAbort("format-ver: format must be 1 or 2.");
    }
    
    //
    if (tsv_name!="") {
        file_in=opts->getArg(0);
        if (file_out=="") {
            file_out=file_in+".v"+opts->get("format-ver");
        }
        QuantLabelZ__convert_priors_tsv5(file_in,tsv_name,file_out,format_ver,verbose);
    }
    else {
        if ((file_out!="") && (opts->getArgCount()>1)) {
            Err::errAbort("cant use output with more than one input.");
        }
        //
        for (int i=0;i<opts->getArgCount();i++) {
            std::string file_in=opts->getArg(i);
            if (file_out=="") {
                file_out=file_in+".v"+opts->get("format-ver");
            }
            QuantLabelZ__convert_priors_tsv(file_in,file_out,format_ver,verbose);
            file_out="";
        }
    }
    
    //
    delete opts;
    return 0;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
