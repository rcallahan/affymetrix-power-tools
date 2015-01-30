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
#include "chipstream/QCAnalysisOptions.h"
//
#include "chipstream/QCProbesetOptions.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Convert.h"
#include "util/Util.h"
//
#include <cctype>
//

using namespace std;

CelListener*
QCMethod::createCelListener(ChipLayout *layout,
                            std::vector<string>& celFiles,
                            QCProbesetOptions &psOpts)
{
  /// @todo err checking that probeset list for groupName exists and is not empty

  // single channel methods
  if (analysis == QCMethod::DM) {
    return createDmListener(layout, psOpts.probesetGroups[groupName], celFiles);
  }
  if (analysis == QCMethod::GENDER) {
    return createEmGenderCelListener(layout, psOpts.probesetGroups[groupName], celFiles);
  }
  if (analysis == QCMethod::HOMHILO) {
    return createHomHiLoCelListener(layout, psOpts.probesetGroups[groupName], celFiles);
  }
  if (analysis == QCMethod::ADJHOMHILO) {
    return createAdjHomHiLoCelListener(layout, psOpts, celFiles);
  }

  // @todo Add the rest of the methods here -- jhg
  // the multichannel methods.
  if (analysis == QCMethod::MULTICHANNEL_DISHQC) {
    return createMultiChannelNonOverlapCelListener(layout, psOpts.probesetGroups[groupName], psOpts);
    return NULL;
  }
  if (analysis == QCMethod::MULTICHANNEL_HOMHILO) {
    return createMultiChannelHomHiLoCelListener(layout, psOpts.probesetGroups[groupName], psOpts);
    return NULL;
  }
  if (analysis == QCMethod::MULTICHANNEL_SIGNALBACKGROUND) {
    return createSignalBackgroundCelListener(layout, psOpts.probesetGroups[groupName], psOpts);
  }
  if (analysis == QCMethod::MULTICHANNEL_VARIABILITY) {
    return createVariabilityScoreCelListener(layout, psOpts.probesetGroups[groupName], psOpts);
  }

  // opps!
  APT_ERR_ABORT("Unable to create CelListener for " + ToStr(analysis));

  // keep g++ quiet
  return NULL;
}

/////

SignalBackgroundCelListener*
QCMethod::createSignalBackgroundCelListener(ChipLayout *layout,
                                            const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                            QCProbesetOptions &psOpts)
{
  SignalBackgroundCelListener *sbCelListener = NULL;

  // SignalBackground takes a vector of ProbeLists
  vector<ProbeListPacked> probeSets;
  fillProbeListVector(probeSets, probesetNames, layout);
  sbCelListener = new SignalBackgroundCelListener(probeSets,psOpts,layout,analysisName);

  map<string, vector<bool> > masks;
  masks["AT"].resize(layout->getProbeCount());
  masks["GC"].resize(layout->getProbeCount());
  masks["A"].resize(layout->getProbeCount());
  masks["C"].resize(layout->getProbeCount());
  masks["G"].resize(layout->getProbeCount());
  masks["T"].resize(layout->getProbeCount());

  // Verbose::out(4, "processing probesetname: " + ToStr(probesetname));
  std::string probeclass;
  std::string ATGCclass; //A probe, G pobe, T probe or C probe?
  map<string, vector<bool> >::iterator maskIter;
  std::string ps_name;

  int ps_idx=0;
  std::string psn_name;
  std::string psn_ligation_base;
  //int line_cnt=0;

  // layout->m_PlFactory.dump();

  // @todo This is so disgusting... - jhg
  for(int psn_i = 0; psn_i < probesetNames.size(); psn_i++) {
    psn_name=probesetNames[psn_i].first;
    psn_ligation_base=probesetNames[psn_i].second;
    //
    //if (line_cnt++<20) {
    //  printf("1:createSignalBackgroundCelListener: '%s':'%s'\n",psn_name.c_str(),psn_ligation_base.c_str());
    //}
    ProbeListPacked ps = layout->getProbeListByName(psn_name);
    if (ps.isNull()) {
      continue;
    }
    //if (line_cnt++<20) {
    //  printf("2:createSignalBackgroundCelListener: '%s':'%s'\n",psn_name.c_str(),psn_ligation_base.c_str());
    ///

    //
    if ((psn_ligation_base=="A")||(psn_ligation_base=="T")) {
      probeclass="AT";
    }
    else {
      probeclass="GC";
    }

    // turn on all the masks.
    for (int pIx = 0; pIx<probeSets[ps_idx].probe_cnt(); pIx++) {
      int probeid=probeSets[ps_idx].get_probeId(pIx);
      masks[probeclass].at(probeid) = true;
      masks[psn_ligation_base].at(probeid) = true;
    }
    //
    ps_idx++;
  } // for each probelist in each probeset
  
  for(maskIter = masks.begin(); maskIter != masks.end(); maskIter++) {
    // count the number of probes set.
    int count = 0;
    for(vector<bool>::iterator mIx = maskIter->second.begin(); mIx != maskIter->second.end(); mIx++) {
      if(*mIx) {
        count++;
      }
    }
    //
    Verbose::out(2,"Number of " + maskIter->first + " probes is: " + ToStr(count));
    sbCelListener->addProbeMask(maskIter->first, maskIter->second);
  }

  Verbose::out(1, "ready to go to SignalBackgroundCelListener...");

  return sbCelListener;
}

/*
  map<string, vector<bool> > masks;
  masks["AT"].resize(layout->getProbeCount());
  masks["GC"].resize(layout->getProbeCount());
  masks["G"].resize(layout->getProbeCount());
  masks["C"].resize(layout->getProbeCount());
  masks["A"].resize(layout->getProbeCount());
  masks["T"].resize(layout->getProbeCount());

  std::vector<std::string> ATlist; //list of AT probeset endings
  ATlist.push_back("A");
  ATlist.push_back("A_f");
  ATlist.push_back("A_r");
  ATlist.push_back("T");
  ATlist.push_back("T_f");
  ATlist.push_back("T_r");

  std::vector<std::string> GClist; //list of GC probeset endings
  GClist.push_back("C");
  GClist.push_back("C_f");
  GClist.push_back("C_r");
  GClist.push_back("G");
  GClist.push_back("G_f");
  GClist.push_back("G_r");

  // Verbose::out(4, "processing probesetname: " + ToStr(probesetname));
  std::string probeclass;
  std::string ATGCclass; //A probe, G pobe, T probe or C probe?
  map<string, vector<bool> >::iterator maskIter;
  std::string ps_name;
  
  for (int i = 0; i< probeSets.size(); i++){
    ps_name=probeSets[i].get_name_string();
    probeclass="";
    ATGCclass="";
    
    // If it has any of these, then it isnt the correct type.
    if (Util::endsWithStr(ps_name,"AG_",1) ||
        Util::endsWithStr(ps_name,"AC_",1) ||
        Util::endsWithStr(ps_name,"GT_",1) ||
        Util::endsWithStr(ps_name,"CT_",1)) {
      APT_ERR_ABORT("QCAnalysisOptions::CreateMultiChannelNonOverlapCelListener: "
                    "Found SNPs when expecting non-polymorphic probesets."
                    "'"+ps_name+"'");
    }
    
    // Verbose::out(4, "Found NonPolymorphic Probesetname: '" + ps_name+"'");
    for (int sIx = 0; sIx < ATlist.size(); sIx++) {
      // @todo This should only occur at the end.
      if (Util::endsWithStr(ps_name,ATlist[sIx])) {
        probeclass="AT";
        if (Util::endsWithStr(ps_name,"A_f") ||
            Util::endsWithStr(ps_name,"A_r") ||
            Util::endsWithStr(ps_name,"A")) {
          ATGCclass = "A";
        }
        else {
          ATGCclass = "T";
        }
        break;
      }
    }
    
    // @todo should skip if above hit.
    for (int sIx = 0; sIx < GClist.size(); sIx++) {
      if (Util::endsWithStr(ps_name,GClist[sIx])) {
        probeclass="GC";
        if (Util::endsWithStr(ps_name,"C_f") ||
            Util::endsWithStr(ps_name,"C_r") ||
            Util::endsWithStr(ps_name,"C")) {
          ATGCclass = "C";
        }
        else {
          ATGCclass = "G";
        }
        break;
      }
    }
    
    APT_ERR_ASSERT(probeclass!="","internal error");
    APT_ERR_ASSERT(ATGCclass!="","internal error");
    
    // turn on all the masks.
    for (int pIx = 0; pIx<probeSets[i].probe_cnt(); pIx++) {
      masks[probeclass].at(probeSets[i].get_probeId(pIx)) = true;
      masks[ATGCclass].at(probeSets[i].get_probeId(pIx)) = true;
    }
  } // for each probelist in each probeset

  for(maskIter = masks.begin(); maskIter != masks.end(); maskIter++) {
    int count = 0;
    for(vector<bool>::iterator mIx = maskIter->second.begin(); mIx != maskIter->second.end(); mIx++) {
      if(*mIx) {
        count++;
      }
    }

    Verbose::out(2,"Number of " + maskIter->first + " probes is: " + ToStr(count));
    sbCelListener->addProbeMask(maskIter->first, maskIter->second);
  }

  Verbose::out(1, "ready to go to SignalBackgroundCelListener...");

  return sbCelListener;
}
*/

/////

VariabilityScoreCelListener*
QCMethod::createVariabilityScoreCelListener(ChipLayout *layout,
                                            const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                            QCProbesetOptions &psOpts)
{
  VariabilityScoreCelListener *VSCelListener = NULL;
  // SignalBackground takes a vector of Probe Lists
  vector<ProbeListPacked> probeSets;

  fillProbeListVector(probeSets, probesetNames, layout);
  VSCelListener = new VariabilityScoreCelListener(probeSets,psOpts,layout,analysisName);
  Verbose::out(1, "ready to go to VariabilityScoreCelListener...");

  return VSCelListener;
}

/////
/*
void
generateMasksFromNames(MultiChannelNonOverlapCelListener* nonoverlapCelListener,
                       vector<ProbeListPacked>& probeSets)
{
  map<string, vector<bool> > masks;
  masks["AT"].resize(probeSets.size());
  masks["GC"].resize(probeSets.size());
  
  map<string, vector<bool> >::iterator maskIter;
  
  std::string probeclass;
  std::string ps_name;

  for (int i = 0; i< probeSets.size(); i++) {
    ps_name=probeSets[i].get_name_string();
  
    //Verbose::out(4, "processing probesetname: " + ToStr(probesetname));
    if (Util::endsWithStr(ps_name,"AG_",1) ||
        Util::endsWithStr(ps_name,"AC_",1) ||
        Util::endsWithStr(ps_name,"GT_",1) ||
        Util::endsWithStr(ps_name,"CT_",1)) {
      APT_ERR_ABORT("QCAnalysisOptions::CreateMultiChannelNonOverlapCelListener: "
                    "Found SNPs when expecting non-polymorphic probesets." 
                    "'"+ps_name+"'");
    }
    
    for (int sIx = 0; sIx < ATlist.size(); sIx++) {
      if (Util::endsWithStr(ps_name,ATlist[sIx])) {
        probeclass="AT";
        break;
      }
    }
    for (int sIx = 0; sIx < GClist.size(); sIx++) {
      if (Util::endsWithStr(ps_name,GClist[sIx])) { 
        probeclass="GC";
        break;
      }
    }

    masks[probeclass].at(i) = true; 
  } // for each probelist in each probeset
  
  for (maskIter = masks.begin(); maskIter != masks.end(); maskIter++) {
    int count = 0;
    for(vector<bool>::iterator mIx = maskIter->second.begin(); mIx != maskIter->second.end(); mIx++) {
      if(*mIx) {
        count++; 
      }
    }
    Verbose::out(2,"Number of " + maskIter->first + " probesets is: " + ToStr(count));
    nonoverlapCelListener->addProbeMask(maskIter->first, maskIter->second); //mask name, and vector of bools
  }
}
*/
  /*
  std::vector<std::string> ATlist; //list of AT probeset endings
  ATlist.push_back("A");
  ATlist.push_back("A_f");
  ATlist.push_back("A_r");
  ATlist.push_back("T");
  ATlist.push_back("T_f");
  ATlist.push_back("T_r");

  std::vector<std::string> GClist; //list of GC probeset endings
  GClist.push_back("C");
  GClist.push_back("C_f");
  GClist.push_back("C_r");
  GClist.push_back("G");
  GClist.push_back("G_f");
  GClist.push_back("G_r");
  
  map<string, vector<bool> >::iterator maskIter;
  
  std::string probeclass;
  std::string ps_name;

  for (int i = 0; i< probeSets.size(); i++) {
    ps_name=probeSets[i].get_name_string();
  
    //Verbose::out(4, "processing probesetname: " + ToStr(probesetname));
    if (Util::endsWithStr(ps_name,"AG_",1) ||
        Util::endsWithStr(ps_name,"AC_",1) ||
        Util::endsWithStr(ps_name,"GT_",1) ||
        Util::endsWithStr(ps_name,"CT_",1)) {
      APT_ERR_ABORT("QCAnalysisOptions::CreateMultiChannelNonOverlapCelListener: "
                    "Found SNPs when expecting non-polymorphic probesets." 
                    "'"+ps_name+"'");
    }
    
    for (int sIx = 0; sIx < ATlist.size(); sIx++) {
      if (Util::endsWithStr(ps_name,ATlist[sIx])) {
        probeclass="AT";
        break;
      }
    }
    for (int sIx = 0; sIx < GClist.size(); sIx++) {
      if (Util::endsWithStr(ps_name,GClist[sIx])) { 
        probeclass="GC";
        break;
      }
    }

    masks[probeclass].at(i) = true; 
  } // for each probelist in each probeset
  
  for (maskIter = masks.begin(); maskIter != masks.end(); maskIter++) {
    int count = 0;
    for(vector<bool>::iterator mIx = maskIter->second.begin(); mIx != maskIter->second.end(); mIx++) {
      if(*mIx) {
        count++; 
      }
    }
    Verbose::out(2,"Number of " + maskIter->first + " probesets is: " + ToStr(count));
    nonoverlapCelListener->addProbeMask(maskIter->first, maskIter->second); //mask name, and vector of bools
  }
}
  */
/*
void
generateMasksFromChannels(MultiChannelNonOverlapCelListener* nonoverlapCelListener,
                          vector<ProbeListPacked>& probeSets)
{
  map<string, vector<bool> > masks;
  masks["AT"].resize(probeSets.size());
  masks["GC"].resize(probeSets.size());
  
  map<string, vector<bool> >::iterator maskIter;
  
  std::string probeclass;
  std::string ps_name;

  ProbeListPacked pl;
  for (int i = 0; i< probeSets.size(); i++) {
    pl=probeSets[i];
    
    if (pl.block_cnt()!=1) {
      APT_ERR_ABORT("internal error.");
    }

    if (pl.get_blockChannel(0)==0) {
      probeclass="GC";
    }
    else if (pl.get_blockChannel(0)==1) {
      probeclass="AT";
    }
    else {
      APT_ERR_ABORT("internal error.");
    }

    masks[probeclass].at(i) = true; 
  } // for each probelist in each probeset
  
  for (maskIter = masks.begin(); maskIter != masks.end(); maskIter++) {
    int count = 0;
    for(vector<bool>::iterator mIx = maskIter->second.begin(); mIx != maskIter->second.end(); mIx++) {
      if(*mIx) {
        count++; 
      }
    }
    Verbose::out(2,"Number of " + maskIter->first + " probesets is: " + ToStr(count));
    nonoverlapCelListener->addProbeMask(maskIter->first, maskIter->second); //mask name, and vector of bools
  }
}
*/
MultiChannelNonOverlapCelListener*
QCMethod::createMultiChannelNonOverlapCelListener(ChipLayout *layout,
                                                  const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                                  QCProbesetOptions &psOpts)
{
  MultiChannelNonOverlapCelListener *nonoverlapCelListener = NULL;
  // SignalBackground takes a vector of Probe Lists
  vector<ProbeListPacked> probeSets;
  fillProbeListVector(probeSets, probesetNames, layout);
  nonoverlapCelListener = new MultiChannelNonOverlapCelListener(probeSets,psOpts,layout,analysisName);

  // generateMasksFromNames(nonoverlapCelListener,probeSets);
  // generateMasksFromChannels(nonoverlapCelListener,probeSets);

  map<string, vector<bool> > masks;
  masks["AT"].resize(probeSets.size());
  masks["GC"].resize(probeSets.size());
  
  map<string, vector<bool> >::iterator maskIter;
  
  std::string probeclass;
  std::string psn_name;
  std::string psn_ligation_base;
  
  int ps_i=0;
  for(int psn_i = 0; psn_i < probesetNames.size(); psn_i++) {
    psn_name=probesetNames[psn_i].first;
    psn_ligation_base=probesetNames[psn_i].second;

    ProbeListPacked ps = layout->getProbeListByName(psn_name);
    if (ps.isNull()) {
      continue;
    }

    // this has already been done.
    //Verbose::out(4, "processing probesetname: " + ToStr(probesetname));
    //if (Util::endsWithStr(ps_name,"AG_",1) ||
    //    Util::endsWithStr(ps_name,"AC_",1) ||
    //    Util::endsWithStr(ps_name,"GT_",1) ||
    //    Util::endsWithStr(ps_name,"CT_",1)) {
    //  APT_ERR_ABORT("QCAnalysisOptions::CreateMultiChannelNonOverlapCelListener: "
    //                "Found SNPs when expecting non-polymorphic probesets." 
    //                "'"+ps_name+"'");
    //}
    
    if ((psn_ligation_base=="A") || (psn_ligation_base=="T")) {
      probeclass="AT";
    }
    else {
      probeclass="GC";
    }
    
    //
    masks[probeclass].at(ps_i) = true; 
    //
    ps_i++;
  }
  
  for (maskIter = masks.begin(); maskIter != masks.end(); maskIter++) {
    int count = 0;
    for(vector<bool>::iterator mIx = maskIter->second.begin(); mIx != maskIter->second.end(); mIx++) {
      if(*mIx) {
        count++; 
      }
    }
    Verbose::out(2,"Number of " + maskIter->first + " probesets is: " + ToStr(count));
    nonoverlapCelListener->addProbeMask(maskIter->first, maskIter->second); //mask name, and vector of bools
  }

  //
  Verbose::out(1, "ready to go to SignalBackgroundCelListener...");
  return nonoverlapCelListener;
}

/////

MultiChannelHomHiLoCelListener*
QCMethod::createMultiChannelHomHiLoCelListener(ChipLayout* layout,
                                               const std::vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                               QCProbesetOptions& psOpts)
{
  MultiChannelHomHiLoCelListener *hiLoCelListener = NULL;
  // Figure out options
  double k = 2.0f;
  double emThresh = 0.05f;
  double binSize = 0.02f;

  if (options.find("em.k") != options.end()) {
    k = QCAnalysisOptions::ToFloat(options["em.k"]);
  }
  if (options.find("em.thresh") != options.end()) {
    emThresh = QCAnalysisOptions::ToFloat(options["em.thresh"]);
  }
  if (options.find("binsize") != options.end()) {
    binSize = QCAnalysisOptions::ToFloat(options["binsize"]);
  }

  // HomHiLoCelListener takes a vector of Probe Lists
  vector<ProbeListPacked> probeSets;
  fillProbeListVector(probeSets, probesetNames, layout);
  Verbose::out(1, "ready to go to MultiChannelHomHiLoCelListener...");
  hiLoCelListener = new MultiChannelHomHiLoCelListener(probeSets,psOpts,layout,k,emThresh,binSize,analysisName);
  return hiLoCelListener;
}

/////

DmListener *
QCMethod::createDmListener(ChipLayout *layout,
                           const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                           const vector<string>& celFiles)
{
  DmListener *dmListener = NULL;

  // Figure out options
  double dmCutOff = 0.33f;
  double dmHetMult = 1.25f;
  if(options.find("dm.cutoff") != options.end())
    dmCutOff = QCAnalysisOptions::ToFloat(options["dm.cutoff"]);
  if(options.find("dm.hetMult") != options.end())
    dmHetMult = QCAnalysisOptions::ToFloat(options["dm.hetMult"]);

  // DmListener needs a map, not a vector
  vector<ProbeListPacked> probesets;
  for(int i = 0; i < probesetNames.size(); i++) {
    ProbeListPacked ps = layout->getProbeListByName(probesetNames[i].first);
    if (!ps.isNull()) {
      probesets.push_back(ps);
    }
  }
  map<string, bool> emptyChrXSnps;

  // Setup thresholds to report
  vector<float> thresholds(2);
  thresholds[0] = dmCutOff;
  thresholds[1] = 0.33f; // for gender, we are hard coded to 0.33

  // Set the column label for call rate
  // We use a "hack" here to make the old labels in QCC file
  // look like the default types of labels
  string label;
  for(int i=0;i<analysisName.size(); i++) {
    // make lower case
    char c = std::tolower(analysisName[i]);
    // replace space and "/" with "-"
    if(c=='/' || c==' ')
      c = '-';
    // remove "(",")"
    if(c!='(' && c!=')')
      label.push_back(c);
  }

  // Create the listener
  dmListener = new DmListener(probesets, dmHetMult, celFiles.size(),
                              thresholds, emptyChrXSnps, label);

  // Set the default index for genotype calling
  dmListener->setDefaultCallIndex(0);
  // Set the default index for gender calling. Currently
  dmListener->setDefaultChrXCallIndex(1);

  return dmListener;
}

EmGenderCelListener*
QCMethod::createEmGenderCelListener(ChipLayout* layout,
                                    const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                    const vector<string>& celFiles)
{
  EmGenderCelListener *emGenderCelListener = NULL;

  // Figure out options
  double k = 4.0f;
  double emThresh = 0.05f;
  double emCutOff = 0.5f;
  double genderCutOff = 0.1f;

  if(options.find("em.k") != options.end())
    k = QCAnalysisOptions::ToFloat(options["em.k"]);
  if(options.find("em.thresh") != options.end())
    emThresh = QCAnalysisOptions::ToFloat(options["em.thresh"]);
  if(options.find("em.cutoff") != options.end())
    emCutOff = QCAnalysisOptions::ToFloat(options["em.cutoff"]);
  if(options.find("gender.cutoff") != options.end())
    genderCutOff = QCAnalysisOptions::ToFloat(options["gender.cutoff"]);

  // EmGenderCelListenertakes a vector of Probe Lists
  vector<ProbeListPacked> probeSets;
  fillProbeListVector(probeSets, probesetNames, layout);

  emGenderCelListener = new EmGenderCelListener(probeSets,k,emThresh,emCutOff,genderCutOff);

  return emGenderCelListener;
}

HomHiLoCelListener*
QCMethod::createHomHiLoCelListener(ChipLayout *layout,
                                   const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                   const vector<string> &celFiles)
{
  HomHiLoCelListener *hiLoCelListener = NULL;

  // Figure out options
  double k = 2.0f;
  double emThresh = 0.05f;
  double binSize = 0.02f;
  if(options.find("em.k") != options.end())
    k = QCAnalysisOptions::ToFloat(options["em.k"]);
  if(options.find("em.thresh") != options.end())
    emThresh = QCAnalysisOptions::ToFloat(options["em.thresh"]);
  if(options.find("binsize") != options.end())
    binSize = QCAnalysisOptions::ToFloat(options["binsize"]);

  // HomHiLoCelListener takes a vector of Probe Lists
  vector<ProbeListPacked> probeSets;
  fillProbeListVector(probeSets, probesetNames, layout);

  hiLoCelListener = new HomHiLoCelListener(probeSets,analysisName,k,emThresh,binSize);
  return hiLoCelListener;
}

AdjHomHiLoCelListener*
QCMethod::createAdjHomHiLoCelListener(ChipLayout *layout,
                                      QCProbesetOptions &psOpts,
                                      vector<string> &celFiles)
{
  AdjHomHiLoCelListener *adjHiLoCelListener = NULL;

  // Figure out options
  double k = 2.0f;
  double emThresh = 0.05f;
  double binSize = 0.02f;
  double hiLoDiffCutOff = 2.0f;
  double hiLoFailedValue = 0.0f;
  if(options.find("em.k") != options.end())
    k = QCAnalysisOptions::ToFloat(options["em.k"]);
  if(options.find("em.thresh") != options.end())
    emThresh = QCAnalysisOptions::ToFloat(options["em.thresh"]);
  if(options.find("binsize") != options.end())
    binSize = QCAnalysisOptions::ToFloat(options["binsize"]);
  if(options.find("adj.cutoff") != options.end())
    hiLoDiffCutOff = QCAnalysisOptions::ToFloat(options["adj.cutoff"]);
  if(options.find("adj.failed") != options.end())
    hiLoFailedValue = QCAnalysisOptions::ToFloat(options["adj.failed"]);

  // AdjHomHiLoCelListener takes a vector of Probe Lists
  vector<ProbeListPacked> randProbesets, set1Probesets, set2Probesets;
  fillProbeListVector(randProbesets, psOpts.probesetGroups[options["ps.rand"]], layout);
  fillProbeListVector(set1Probesets, psOpts.probesetGroups[options["ps.set1"]], layout);
  fillProbeListVector(set2Probesets, psOpts.probesetGroups[options["ps.set2"]], layout);

  adjHiLoCelListener = new AdjHomHiLoCelListener(randProbesets, set1Probesets, set2Probesets,
                                                 options["label.rand"], options["label.set1"], options["label.set2"],
                                                 k, emThresh, binSize, analysisName, hiLoDiffCutOff, hiLoFailedValue);

  return adjHiLoCelListener;
}

//////////

// private method to fill in a probelist vector
void QCMethod::fillProbeListVector(vector<ProbeListPacked>& probeSets,
                                   const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                   ChipLayout *layout) {
  for(int i = 0; i < probesetNames.size(); i++) {
    ProbeListPacked ps = layout->getProbeListByName(probesetNames[i].first);
    if (!ps.isNull()) {
      probeSets.push_back(ps);
    }
  }
}

// private method to fill in a probelist vector
void QCMethod::fillProbeListVector(vector<ProbeListPacked>& probeSets,
                                   const vector<string>& probesetNames,
                                   ChipLayout *layout) {
  for(int i = 0; i < probesetNames.size(); i++) {
    ProbeListPacked ps = layout->getProbeListByName(probesetNames[i]);
    if (!ps.isNull()) {
      probeSets.push_back(ps);
    }
  }
}

// Private utility function to convert string to float
float QCAnalysisOptions::ToFloat(string valString) {
  string format = "%f";
  float val = 0.0;
  if(sscanf(valString.c_str(),format.c_str(),&val) == 0)
    Err::errAbort("Unable to convert analysis option to a number.");
  return val;
}

// Private utility function to parse options
bool QCAnalysisOptions::SplitOptions(string& str,pair<string,string>& val,const char sep1,const char sep2)
{
  int index = str.find_first_of(sep1);
  val.first = str.substr(0,index);
  int stop = str.find_first_of(sep2);
  if(stop == -1) {
    val.second = str.substr(index+1,str.length()-index);

    // done parsing
    return false;
  } else {
    val.second = str.substr(index+1,stop-index-1);

    // return unparsed part of string
    str = str.substr(stop+1,str.size()-index);
    return true;
  }

}

// Static method to read in a QCAFile and populate QCAnalysisOptions
void QCAnalysisOptions::readQCAFile(const string &qcaFileName, QCAnalysisOptions &opts) {

  APT_ERR_ASSERT((qcaFileName!=""),"qcaFilename is blank.");

  affx::TsvFile tsv;

  if ((tsv.open(qcaFileName)) != affx::TSV_OK)
    Err::errAbort("Failed to open '" + qcaFileName);

  // vars to bind
  string analysis, analysis_name, group_name, options;
  tsv.bind(0, "analysis_name", &analysis_name, affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "group_name",    &group_name,    affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "analysis",      &analysis,      affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "options",       &options,       affx::TSV_BIND_OPTIONAL);

  // get the default analysis name from the file header
  if(tsv.getHeader("default_analysis_name",opts.defaultAnalysisName) != affx::TSV_OK)
    Err::errAbort("Failed to find 'default_analysis_name' in the file header.");

  // get the column values for each method
  int rv;
  std::string analysis_lc;
  std::string analysis_name_lc;

  while ((rv = tsv.nextLevel(0))==affx::TSV_OK) {
    QCMethod method;

    //
    analysis_lc=Util::downcaseString(analysis);
    analysis_name_lc=Util::downcaseString(analysis_name);

    APT_ERR_ASSERT((analysis_name_lc!=""),"Invalid analysis name (empty)");
    APT_ERR_ASSERT((analysis!=""),"Invalid analysis (empty)");

    method.analysisName = analysis_name;
    method.groupName = group_name;

    if (analysis_lc=="dm") {
      method.analysis = QCMethod::DM;
    }
    else if (analysis_lc=="gender") {
      method.analysis = QCMethod::GENDER;
    }
    else if (analysis_lc=="homhilo") {
      method.analysis = QCMethod::HOMHILO;
    }
    else if (analysis_lc=="adj-homhilo") {
      method.analysis = QCMethod::ADJHOMHILO;
    }
    //
    else if ((analysis_lc=="mc-dishqc")) {
      method.analysis = QCMethod::MULTICHANNEL_DISHQC;
    }
    else if ((analysis_lc=="mc-homhilo")) {
      method.analysis = QCMethod::MULTICHANNEL_HOMHILO;
    }
    else if ((analysis_lc=="mc-signalbackground")) {
      method.analysis = QCMethod::MULTICHANNEL_SIGNALBACKGROUND;
    }
    else if ((analysis_lc=="mc-variability")) {
      method.analysis = QCMethod::MULTICHANNEL_VARIABILITY;
    }
    else {
      APT_ERR_ABORT("Invalid analysis requested '" + analysis_lc + "'");
    }

    if (options != "" && options != "NA") {
      bool stop = false;;
      do {
        pair<string,string> val;
        stop = SplitOptions(options,val);
        method.options.insert(val);
      } while (stop && options != "");
    }
    opts.methods.push_back(method);
  }

  if (rv != affx::TSV_OK && rv != affx::TSV_ERR_EOF && rv != affx::TSV_ERR_FILEIO) {
    /// @todo add line number
    /// @todo have tsv.getError add the file and linenumber
    APT_ERR_ABORT("Failed to read line: " + tsv.getError());
  }
  tsv.close();
}
