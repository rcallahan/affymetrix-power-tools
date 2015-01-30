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

//
#include "chipstream/QCProbesetOptions.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Convert.h"
#include "util/Util.h"
//
#include <utility>

QCProbesetOptions::QCProbesetOptions() {
  m_exprAnalysisString="quant-norm.target=1000.sketch=50000,pm-only,median";
  m_hasLigationBase=false;
  m_tempdir=".";
  useNames = false;
}

// Private utility function to parse groups
bool QCProbesetOptions::SplitGroups(std::string& str,std::string& val,const char sep)
{
  int index = str.find_first_of(sep);
  if(index == -1) {
    val = str;
    return false;
  }
  val = str.substr(0,index);  
  str = str.substr(index+1,str.size()-index);
  return true;
}

// Static method to read in a QCCFile and populate QCProbesetOptions
void QCProbesetOptions::readQCCFile(const std::string &qccFileName, QCProbesetOptions &opts) {
  affx::TsvFile tsv;
  //int lines_out=0;

  if ((tsv.open(qccFileName)) != affx::TSV_OK) {
    Err::errAbort("Failed to open '" + qccFileName);
  }

  // vars to bind
  int probeset_id = -1;
  int quantification_in_header = -1;
  std::string group_names="";
  std::string probeset_name="";
  std::string ligation_base="?";

  //
  tsv.bind(0, "probeset_id",              &probeset_id,              affx::TSV_BIND_OPTIONAL);
  tsv.bind(0, "group_name",               &group_names,              affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "probeset_name",            &probeset_name,            affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "quantification_in_header", &quantification_in_header, affx::TSV_BIND_OPTIONAL);
  tsv.bind(0, "ligation_base",            &ligation_base,            affx::TSV_BIND_OPTIONAL);

  // Is the primary key probeset_id or probeset_name?
  std::string primaryKey;
  if(tsv.getHeader("primary_key",primaryKey) != affx::TSV_OK) {
    opts.useNames = false;
    ///@todo check that probeset_id was bound
  }
  else if (primaryKey == "probeset_name") {
    opts.useNames = true;
  } 
  else {
    Err::errAbort("Invalid primary key '" + primaryKey + "' found. Expecting 'probeset_id' or 'probeset_name'.");
  }

  if (tsv.cname2cidx(0,"ligation_base")>=0) {
    opts.m_hasLigationBase=true;
  }
  else {
    opts.m_hasLigationBase=false;
  }
    
  // get the column values for each method
  int rv;
  while ((rv = tsv.nextLevel(0))==affx::TSV_OK) {

    if (opts.useNames && probeset_name == "") {
      Err::errAbort("Invalid probeset_name (empty)");
    }
    if (!opts.useNames && probeset_id < 1) {
      Err::errAbort("Invalid probeset_id '" + ToStr(probeset_id) + "'. Expecting integer greater than 0.");
    }
    if (group_names == "") {
      Err::errAbort("Invalid group_name (empty)");
    }

    // figure out what groups this probeset belongs to
    std::vector<std::string> group_name_vec;

    ///@todo AW: This block should go, or otherwise only be executed if the parsing of ligation bases from probeset names was requested.
    if (!opts.m_hasLigationBase) {
      APT_ERR_ASSERT(probeset_name!="","name required.");
      // If it has any of these, then it isnt the correct type.
      if (Util::endsWithStr(probeset_name,"AG_",1) ||
          Util::endsWithStr(probeset_name,"AC_",1) ||
          Util::endsWithStr(probeset_name,"GT_",1) ||
          Util::endsWithStr(probeset_name,"CT_",1)) {
        APT_ERR_ABORT("QCAnalysisOptions::CreateMultiChannelNonOverlapCelListener: "
                      "Found SNPs when expecting non-polymorphic probesets."
                      "'"+probeset_name+"'");
      }
      //
      if (Util::endsWithStr(probeset_name,"A_",1) || Util::endsWithStr(probeset_name,"A")) {
        ligation_base="A";
      }
      else if (Util::endsWithStr(probeset_name,"C_",1) || Util::endsWithStr(probeset_name,"C")) {
        ligation_base="C";
      }
      else if (Util::endsWithStr(probeset_name,"G_",1) || Util::endsWithStr(probeset_name,"G")) {
        ligation_base="G";
      }
      else if (Util::endsWithStr(probeset_name,"T_",1) || Util::endsWithStr(probeset_name,"T")) {
        ligation_base="T";
      }
      else {
        // Non falcon, eg SNP6 chips hit this point...
        //APT_ERR_ABORT("Ending not recognized.");
      }
    }

    // put 
    Util::chopString(group_names,' ',group_name_vec);
    for(int i=0; i<group_name_vec.size(); i++) {
      if (!opts.useNames) {
        probeset_name=ToStr(probeset_id);
      }
      //
      opts.probesetGroups[group_name_vec[i]].push_back(std::pair<std::string,std::string>(probeset_name,ligation_base));
    }

    //
    //if (lines_out++<100) {
    //Verbose::out(1,"QCC: "+ToStr(lines_out)+":"+probeset_name+":"+ligation_base);
    //}

    // Quantification in CHP header
    if(quantification_in_header > 0) {
      if(opts.useNames) 
        opts.headerProbesets.push_back(probeset_name);
      else
        opts.headerProbesets.push_back(ToStr(probeset_id));
    }
    //
    ligation_base="?";
  } // while

  if (rv != affx::TSV_OK && rv != affx::TSV_ERR_EOF && rv != affx::TSV_ERR_FILEIO) {
    Err::errAbort("Failed to read line: " + tsv.getError());
  }
  tsv.close();
}
