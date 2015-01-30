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

/// @file   apt-convert-spf.cpp
/// @brief  Converts an spf file from one SPF format to another.  Reformats DMET SPF files.

//
#include "chipstream/ProbeListFactory.h"
#include "chipstream/ProbeListStl.h"
//
#include "file/CDFFileData.h"
#include "file/TsvFile/SpfFile.h"
#include "util/Err.h"
#include "util/PgOptions.h"
#include "util/Verbose.h"
//
#include "affy-pcre.h"
#include "pcrecpp.h"
//
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
//



/// @brief     Convert an spf from one format to another.
/// @param     spf_i_name
/// @param     spf_o_name
/// @param     spf_o_format
/// @param     extra_columns
void apt_convert_spf(const std::string& spf_i_name,
                     const std::string& spf_o_name,
                     int spf_o_format,
                     int extra_columns
                     )
{
  printf("apt-convert-spf: '%s' => '%s' (fmt=%d)\n",
         spf_i_name.c_str(),
         spf_o_name.c_str(),
         spf_o_format);
  //
  ProbeListFactory* plf=new ProbeListFactory();
  affx::SpfFile spf_i;
  affx::SpfFile spf_o;

  spf_i.openSpf(spf_i_name);
  plf->readSpfFile(spf_i);
  spf_i.close();

  //
  spf_o.addHeadersFrom(spf_i,affx::TSV_ADD_ALL);
  //
  spf_o.m_has_allele_info=extra_columns;
  spf_o.m_has_context_info=extra_columns;
  spf_o.m_has_channel_info=extra_columns;
  //
  spf_o.writeSpf(spf_o_name,spf_o_format);
  plf->writeToSpfFile(spf_o);
  //
  delete plf;
}

//////////

/// A handy macro.
#define PRINT_PS_CNT(_lbl,_plf) {                                             \
    std::string tmp_str=_lbl;                                           \
    printf("%s psCount=%d\n",tmp_str.c_str(),_plf->getProbeSetCount()); \
  }

ProbeListFactory*
apt_falcon_merge_spf2plf(const std::string& spf_i_name)
{
  ProbeListFactory* plf=new ProbeListFactory();
  plf->readSpfFile(spf_i_name);
  PRINT_PS_CNT("read probes from '"+spf_i_name+"'",plf);
  //
  // printf("dump:start: '%s'\n",spf_i_name.c_str());
  // plf->dump();
  // printf("dump:end: '%s'\n",spf_i_name.c_str());
  //
  return plf;
}

void
apt_falcon_merge_chiptype(ProbeListFactory* plf_out,
                          ProbeListFactory* plf_from)
{
  plf_out->m_chipTypes   = plf_from->m_chipTypes;
  plf_out->m_numCols     = plf_from->m_numCols;
  plf_out->m_numRows     = plf_from->m_numRows;
  plf_out->m_numChannels = plf_from->m_numChannels;
}

// A simple concat of the input Factory to the output Factory.
void
apt_falcon_merge_plf2plf(ProbeListFactory* plf_to,
                         ProbeListFactory* plf_from,
                         int chan)
{
  ProbeListStl pl_tmp;
  //
  for (int i=0;i<plf_from->getProbeSetCount();i++) {
    pl_tmp=plf_from->getProbeListAtIndex(i);
    // set the channel
    for (int b=0;b<pl_tmp.block_cnt();b++) {
      pl_tmp.set_blockChannel(b,chan);
    }
    //
    plf_to->add_ProbeList(pl_tmp);
  }
}

// How to do the mapping.
//
// * NonAS
//
//   5298423AC_r 1 1 5 0 1 5 1502317,1502569,1502066,1503094,1502830
//     This is mono-morphic. (we think)
//     G will ask Therea
//     
//   snp_rs2820041_0w_r 1 1 6 0 1 6 689782,688627,688019,678473,313521,687601
//     Goes on both channels. (0 and 1)
//
//   snp_rs27908_0w_f 1 1 6 0 1 6 1255986,1255761,882153,1256244,1255478,1255186
//   snp_rs27908_1wa0_r 1 1 6 0 1 6 492244,867320,867782,867007,867514,866742
//   snp_rs27908_1wa1_r 1 1 6 0 1 6 492243,866741,867781,867319,867513,867008
//     Do not merge (a0 and a1) into a single probeset.
//     Instead treat as above.  Put on both channels.
//     (This makes for 3 probesets.)
//   
//   indel_rs1182182_0w_f
//     These may only be in the NonAS file
//
// * AS (AT/GC)
// 
//   snp_rs10000235_0w_r-A	1	1	3	0	1	3	84962,97117,66491
//   snp_rs10000235_0w_r-B	1	1	3	0	1	3	83474,95629,67979
//     merge into one probeset "snp_rs10000235_0w_r" 
//     A and B go to same channel based on filename
//   
//   The "-A" or "-B" is required.
//

// Constuct the RE only once.
pcrecpp::RE getTypeFromName_Marker_re("^(snp|indel|snpind)_");

int getTypeFromName(const std::string& ps_name)
{
  // Does this name look like a Marker name?
  if (getTypeFromName_Marker_re.PartialMatch(ps_name)) {
    // yep.
    return affxcdf::MarkerProbeSetType;
  }
  // nope.
  return affxcdf::UnknownProbeSetType;
}

/// put the nonas probes on both channels.
void
apt_falcon_merge_plf2plf_nonas(ProbeListFactory* plf_out,
                               ProbeListFactory* plf_from)
{
  ProbeListStl pls_orig;
  ProbeListStl pls_out;
  std::string match_name;

  for (int idx=0;idx<plf_from->getProbeSetCount();idx++) {
    pls_orig.clear();
    pls_out.clear();

    //
    pls_orig=plf_from->getProbeListAtIndex(idx);
    APT_ERR_ASSERT(pls_orig.block_cnt()==1,"Only one block allowed per probeset in input. "+
                   ToStr(pls_orig.get_name())+" has "+ToStr(pls_orig.block_cnt())+".");

    // init output ProbeList based on input.
    pls_out.set_name(pls_orig.get_name());
    // type
    pls_out.set_type(getTypeFromName(pls_out.get_name()));
    // force to "1"
    pls_out.set_numMatch(1);

    // make a copy for each channel.
    // block 0 always A=0,C=1
    pls_orig.set_blockAllele(0,0);
    pls_orig.set_blockContext(0,0); // always 0
    pls_orig.set_blockChannel(0,1);
    pls_out.push_probelist(pls_orig); // push modified block as first block on output.
    // block 1 always A=1,C=0
    pls_orig.set_blockAllele(0,1);
    pls_orig.set_blockContext(0,0);
    pls_orig.set_blockChannel(0,0);
    pls_out.push_probelist(pls_orig); // push modified block as second block on output
    //
    plf_out->add_ProbeList(pls_out);
  }
}

void
apt_falcon_merge_plf2plf_AB(ProbeListFactory* plf_to,
                            ProbeListFactory* plf_from,
                            const int channel_to)
{
  //const int warn_max=50;
  //int warn_count=0;
  //
  std::string match_name;
  std::string match_allele_ab;
  int match_allele;
  //int match_context;
  std::string match_dir;
  //
  std::string out_name;
  //
  ProbeListStl pls;
  //
  std::map<std::string,ProbeListStl> name2pls_map;
  std::map<std::string,ProbeListStl>::iterator name2pl_map_iter;
  std::map<std::string,ProbeListStl>::iterator iter;

  // A  "^...$" is implied by FullMatch()
  pcrecpp::RE allele_AB_re("(.*)-([AB])");

  //
  for (int idx=0;idx<plf_from->getProbeSetCount();idx++) {
    //
    match_name="";
    match_allele=-1;
    match_dir="";
    //
    pls=plf_from->getProbeListAtIndex(idx);
    //
    if (allele_AB_re.FullMatch(pls.get_name_string(),
                               &match_name,
                               &match_allele_ab)) {
      out_name=match_name;
      if (match_allele_ab=="A") {
        match_allele=0;
      }
      else if (match_allele_ab=="B") {
        match_allele=1;
      }
      else {
        APT_ERR_ABORT("Probeset name does not end with 'A' or 'B'.");
      }
    }
    else {
      APT_ERR_ABORT("Probeset name does match RE.");
    }

    // look it up.
    iter=name2pls_map.find(out_name);
    if (iter==name2pls_map.end()) {
      std::pair<std::map<std::string,ProbeListStl>::iterator,bool> rv;
      rv=name2pls_map.insert(std::pair<std::string,ProbeListStl>(out_name,ProbeListStl()));
      iter=rv.first;
      iter->second.set_name(out_name);
      iter->second.set_type(getTypeFromName(out_name));
      // force to "1"
      iter->second.set_numMatch(1);
    }
    // the probeset we are about to merge may only have one block.
    // (We could do more, but that is all we want to for now.)
    APT_ERR_ASSERT(pls.block_cnt()==1,"May only have one block. Index="+ToStr(idx));
    // it should also have the same numMatch as what it is being merged with.
    // (No merging a pm/mm with a pm-only probeset!)
    APT_ERR_ASSERT(pls.get_numMatch()==iter->second.get_numMatch(),"mismatch in numMatch!");
    //
    // pls is what is about to be pushed.
    pls.set_blockAllele(0,match_allele);
    pls.set_blockContext(0,0);
    pls.set_blockChannel(0,channel_to);
    //
    iter->second.push_probelist(pls);
  }

 // put the merged probesets into the output ProbeListFactory
  for (std::map<std::string,ProbeListStl>::iterator iter=name2pls_map.begin();iter!=name2pls_map.end();iter++) {
    ProbeListStl* pls=&iter->second;
    if (pls->block_cnt()!=2) {
      printf("WARNING: probeset '%s' has %d blocks.\n",pls->get_name_cstr(),pls->block_cnt());
      pls->dump();
    }
    // make sure the blocks are in the correct order.  (A then B)
    APT_ERR_ASSERT(pls->get_blockAllele(0)==0,"Block 0 is not allele 0.");
    APT_ERR_ASSERT(pls->get_blockAllele(1)==1,"Block 1 is not allele 1.");
    //
    plf_to->add_ProbeList(*pls);
  }
}

void
apt_falcon_merge(const std::vector<std::string>& spf_i_nonas_name_vec,
                 const std::vector<std::string>& spf_i_at_name_vec,
                 const std::vector<std::string>& spf_i_gc_name_vec,
                 int change_probe_type,
                 bool do_sort,
                 //
                 const std::string& spf_o_name)
{
  //
  //
  ProbeListFactory* plf_out=new ProbeListFactory();
  //
  std::vector<int> channel_vec;

  // nonas is on both channels.
  for (int i=0;i<spf_i_nonas_name_vec.size();i++) {
    ProbeListFactory* plf_nonas=apt_falcon_merge_spf2plf(spf_i_nonas_name_vec[i]);
    apt_falcon_merge_plf2plf_nonas(plf_out,plf_nonas);
    PRINT_PS_CNT("After adding nonas:",plf_out);
    //
    apt_falcon_merge_chiptype(plf_out,plf_nonas);
    delete plf_nonas;
  }

  // Map AT to channel 1
  for (int i=0;i<spf_i_at_name_vec.size();i++) {
    ProbeListFactory* plf_at=apt_falcon_merge_spf2plf(spf_i_at_name_vec[i]);
    apt_falcon_merge_plf2plf_AB(plf_out,plf_at,1);
    PRINT_PS_CNT("After adding AT:   ",plf_out);
    //
    apt_falcon_merge_chiptype(plf_out,plf_at);
    delete plf_at;
  }

  // Map GC to channel 0
  for (int i=0;i<spf_i_gc_name_vec.size();i++) {
    ProbeListFactory* plf_gc=apt_falcon_merge_spf2plf(spf_i_gc_name_vec[i]);
    apt_falcon_merge_plf2plf_AB(plf_out,plf_gc,0);
    PRINT_PS_CNT("After adding GT:   ",plf_out);
    //
    apt_falcon_merge_chiptype(plf_out,plf_gc);
    delete plf_gc;
  }

  // make it ordered for searching...
  if (do_sort) {
    printf("Sorting probesets by name...\n");
    plf_out->sortFactoryByName();
  }

  // do we want to change the probe type?
  if (change_probe_type!=-1) {
    for (int psi=0;psi<plf_out->getProbeSetCount();psi++) {
      plf_out->getProbeListAtIndex(psi).set_type(change_probe_type);
    }
  }

  // force the channel count to 2 as this is the whole point.
  plf_out->m_numChannels=2;
  // ...and write it to our new SpfFile
  printf("Writing to '%s'...\n",spf_o_name.c_str());
  plf_out->writeSpfFile(spf_o_name,4,1,1);

  //
  delete plf_out;
}

//////////

void define_apt_convert_spf_options(PgOptions* opts)
{
  opts->setUsage("apt-convert-spf - Program for converting spf files.\n"
"\n"
"EXAMPLES: \n"
"\n"
"* To merge falcon spf files and assign the probes to channels:\n"
"    apt-convert-spf \\\n"
"	     --output falcon-merged.spf \\\n"
"	     --merge-nonas HTFalconScreen01-nonAS_geno.spf \\\n"
"	     --merge-at HTFalconScreen01-norm_AS_HM_AT.spf \\\n"
"	     --merge-gc HTFalconScreen01-norm_AS_HM_GC.spf \\\n"
"\n"
"* To convert the format of an spf file:\n"
"    apt-convert-spf -f v3 --output foo-v3.spf GenomeWideSNP_6.spf.small\n"

"The Row x Col size is taken from the input Spfs.\n"
"When --merge is used, num-channels is set to '2'.\n"
"\n"
);
  
  opts->defOpt("h", "help", PgOpt::BOOL_OPT,
                     "Display program options short blurb on usage.",
                     "false");
  //
  opts->defOpt("f", "output-format", PgOpt::INT_OPT,
                     "Output format: [2|3|4] default= 4",
                     "4");
  //
  opts->defOpt("o","output",PgOpt::STRING_OPT,
                     "The output file to generate.",
                     "");
  opts->defOpt("","extra-columns",PgOpt::BOOL_OPT,
                     "Add allele,context and channel columns.",
                     "true");
  //
  opts->defOpt("d","dump",PgOpt::BOOL_OPT,
               "Dump the probes to stdout.",
               "false");
  opts->defOpt("s","sort",PgOpt::BOOL_OPT,
               "Sort the probes by name before dumping.",
               "true");
  //
  opts->defOptMult("","merge-nonas",PgOpt::STRING_OPT,
                   "The SPF files to merge onto both channels. (May be repeated.)",
                   "");
  opts->defOptMult("","merge-at",PgOpt::STRING_OPT,
                   "The files to merge into the AT channel. (May be repeated.)",
                   "");
  opts->defOptMult("","merge-gc",PgOpt::STRING_OPT,
                   "The files to merge into the GC channel.  (May be repeated.)",
                   "");
  //
  opts->defOpt("cpt","change-probe-type",PgOpt::INT_OPT,
               "Change the probe type.  (-1 means no change)",
               "-1");
  
}

int main(int argc, char *argv[]) {
  try {
    PgOptions *opts = NULL;
    opts = new PgOptions();
    define_apt_convert_spf_options(opts);
    opts->parseArgv(argv);
    
    int extra_columns;
    int output_format;
    std::string output_file;
    std::string arg;
    
    if(opts->getBool("help") || argc == 1) {
        opts->usage(true);
        return 0;
    }
    
    if (opts->getBool("dump")) {
        for (int i=0;i<opts->getArgCount();i++) {
        arg=opts->getArg(i);
        //
        ProbeListFactory* plf=new ProbeListFactory();
        plf->readSpfFile(arg);
        if (opts->getBool("sort-dump")) {
            plf->sortFactoryByName();
        }
        plf->dump();
        delete plf;
        }
        return 0;
    }
    
    extra_columns=opts->getBool("extra-columns")?1:0;
    output_format=opts->getInt("output-format");
    output_file=opts->get("output");
    
    // if any are set we are being asked to merge.
    if ((opts->get("merge-nonas")!="")||
        (opts->get("merge-at")!="") ||
        (opts->get("merge-gc")!="")) {
        //
        apt_falcon_merge(opts->mustFindOpt("merge-nonas")->getValueVector(),
                        opts->mustFindOpt("merge-at")->getValueVector(),
                        opts->mustFindOpt("merge-gc")->getValueVector(),
                        opts->getInt("change-probe-type"),
                        opts->getBool("sort"),
                        //
                        opts->get("output"));
        //
        return 0;
    }
    
    // a single file
    if (output_file!="") {
        if (opts->getArgCount()!=1) {
        APT_ERR_ABORT("Can only give one input arg when using --output.");
        }
        apt_convert_spf(opts->getArg(0),output_file,output_format,extra_columns);
        return 0;
    }
    
    if (opts->getArgCount()>0) {
        for (int i=0;i<opts->getArgCount();i++) {
        arg=opts->getArg(i);
        output_file=arg+"-v"+ToStr(output_format);
        //
        apt_convert_spf(arg,output_file,output_format,extra_columns);
        }
    }
    
    //
    printf("Need to give an input file.\n");
    opts->usage(true);
    return 0;
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
