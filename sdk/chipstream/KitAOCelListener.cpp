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
#include "chipstream/KitAOCelListener.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "chipstream/CelReader.h"
#include "chipstream/ChipStream.h"
#include "chipstream/ChipStreamFactory.h"
#include "chipstream/SparseMart.h"
#include "util/Convert.h"
#include "util/Err.h"
//
#include <algorithm>

KitAOCelListener::KitAOCelListener() {
  init();
  // explain the output
  addHeaderComment("reagent_version legend:");
  addHeaderComment("  Unknown==0");
  addHeaderComment("  Kit_A==1");
  addHeaderComment("  Kit_O==2");
  addHeaderComment("");
  //
  declareMetric("reagent_version",ChipSummary::Metric::Integer);
  declareMetric("reagent_discrimination_value",ChipSummary::Metric::Double,6);
  m_layout = NULL;
  m_return_kit_id = 0;
  m_return_kit_value = 0;
}

KitAOCelListener::~KitAOCelListener() {
  clear();
}

//////////

void KitAOCelListener::init() {
  m_db=NULL;
  //
  m_exprAnalysisString="quant-norm.sketch=50000";
}

void KitAOCelListener::clear() {
  //m_db.clear();
}

void KitAOCelListener::setDb(KitAODb* db) {
  m_db=db;
}
void KitAOCelListener::setLayout(ChipLayout* layout) {
  m_layout=layout;
}

void KitAOCelListener::newChip(affymetrix_fusion_io::FusionCELData* cel) {
  int chipIx=getNextChipIdx();

  //
  m_celFileName=cel->GetFileName();
  computeKitAO();

  //
  setMetric(chipIx,"reagent_version",m_return_kit_id);
  setMetric(chipIx,"reagent_discrimination_value",m_return_kit_value);

  //
  setValid(true);
}

double KitAOCelListener::getClassiferAValue() const {
  APT_ERR_ASSERT(m_db!=NULL,"m_db is NULL");
  return m_db->getClassiferAValue();
}
double KitAOCelListener::getClassiferOValue() const {
  APT_ERR_ASSERT(m_db!=NULL,"m_db is NULL");
  return m_db->getClassiferOValue();
}

// Again, like the  we only have the celfile and no
void KitAOCelListener::computeKitAO() {
  //
  //printf("KitAOCelListener::computeKitAO()...\n");

  // check our pointers are valid.
  APT_ERR_ASSERT(m_db!=NULL,"internal error");
  APT_ERR_ASSERT(m_layout!=NULL,"internal error");

  APT_ERR_ASSERT(m_sketchFileName!="","a sketch is required.");
  //////////

  ChipStreamFactory chipStreamFactory;
  if (m_sketchFileName!="") {
    chipStreamFactory.readTargetSketchFromFile(m_sketchFileName);
  }
  ChipStream* chipStream=NULL;
  std::string junk;
  chipStream=chipStreamFactory.chipStreamForString(m_exprAnalysisString,*m_layout,junk);

  //
  int probeCount=m_layout->getProbeCount();
  std::vector<int> dummy_order(probeCount);
  for (int i = 0; i < probeCount; i++) {
    dummy_order[i] = i;
  }

  std::vector<std::string> celFileNames;
  celFileNames.push_back(m_celFileName);
  //
  SparseMart sparseMart = SparseMart(dummy_order,celFileNames,true);
  //
  CelReader celReader;
  celReader.setFiles(celFileNames);
  celReader.registerStream(chipStream);
  celReader.registerIntensityMart(&sparseMart);

  celReader.readFiles();

  //////////

  // reset the return info to initial values.
  m_return_kit_id=-1;
  m_return_kit_value=0.0;

  KitAODbEntry* kitao_entry;

  //
  double total_signal_size=0.0;
  ProbeListPacked probeset;
  std::string probeset_name;
  double probeset_signal_size;
  std::vector<int> probeset_probeids;
  //
  double probe_intensity_ch0;
  double probe_intensity_ch1;
  int probe_id;

  //
  total_signal_size=0.0;
  std::vector<double> log2vals;

  double log2=log(2.0);
  int missing_probe_count=0;

  // loop over the probesets on this chip.
  for (KitAODb::name2entry_iter_t entry_i=m_db->m_name2entry.begin();
       entry_i!=m_db->m_name2entry.end();
       entry_i++) {
    //
    probeset_name=entry_i->first;
    kitao_entry=entry_i->second;

    //printf("probeset: '%s': \n",probeset_name.c_str());

    //
    ProbeListPacked probeset = m_layout->getProbeListByName(probeset_name);
    if (probeset.isNull()) {
      missing_probe_count++;
      Verbose::out(1,"KitAOCelListener::computeKitAO(): didnt find probeset '"
                   +probeset_name+"' in layout.");
      continue;
    }

    // We want just the first block of probesets.
    // the second channel should be on the second block.
    probeset_probeids.clear();
    probeset.get_probeIdsForBlock(0,probeset_probeids);

    probeset_signal_size=0.0;
    log2vals.clear();
    for (int p_idx=0;p_idx<probeset_probeids.size();p_idx++) {
      probe_id=probeset_probeids[p_idx];
      probe_intensity_ch0=chipStream->getTransformedIntensity(probe_id,0,0);
      probe_intensity_ch1=chipStream->getTransformedIntensity(probe_id,0,1);
      //
      //printf("%d=%.10f/%.10f \n",probe_id,probe_intensity_ch0,probe_intensity_ch1);
      double log2val=(log(probe_intensity_ch0+probe_intensity_ch1)/log2/2.0);
      log2vals.push_back(log2val);

      // probeset_id probe_id chan0 chan1 log2.A.B
//      printf("%s\t%d\t%f\t%f\t%f\n",
//             probeset_name.c_str(),probe_id,
//             probe_intensity_ch0,probe_intensity_ch1,
//             log2val);
    }

    // pick out the median
    if (log2vals.size()==0) {
      APT_ERR_ABORT("Probeset '"+probeset_name+"' has 0 probeids.");
    }
    std::sort(log2vals.begin(),log2vals.end());
    double median_log2val;
    int median_idx=(log2vals.size()-1)/2;
    if (log2vals.size()==1) {
      median_log2val=log2vals[0];
    }
    else if ((log2vals.size()%2)==1) {
      median_log2val=log2vals[median_idx];
    }
    else {
      median_log2val=(log2vals[median_idx]+log2vals[median_idx+1])/2.0;
    }

    // a Format to compare with Johnathans R values.
    //
    // printf(" raw med_log=%.10f",median_log2val);
    // printf(" means=%.10f ",kitao_entry->m_means);
    // printf(" PC1=%.10f ",kitao_entry->m_pc1);
    // probeset_signal_size has the sum of probes log2(A+B)/2.0,
    // Now subtract the mean and multiply by PC1.
    probeset_signal_size=(median_log2val-kitao_entry->m_means)*kitao_entry->m_pc1;
    // printf(" adj ss=%f",probeset_signal_size);
    // Add this signal_size to the total for this chip.
    total_signal_size+=probeset_signal_size;
    // printf(" total=%.10f",total_signal_size);
    // printf("\n");
  }

  // now do the abort.
  if (missing_probe_count!=0) {
    APT_ERR_ABORT("KitAOCelListener.cpp: "+ToStr(missing_probe_count)+" probes missing.");
  }
  
  // now classify this cel file.
  m_return_kit_value=total_signal_size;
  //
  double classifierAValue=getClassiferAValue();
  Verbose::out(1,"KitAOCelListener: getClassiferCutoffAValue()=="+ToStr(classifierAValue));
  double classifierOValue=getClassiferOValue();
  Verbose::out(1,"KitAOCelListener: getClassiferCutoffOValue()=="+ToStr(classifierOValue));

  APT_ERR_ASSERT(classifierAValue<classifierOValue,"A-Value must be less than O-Value.");
  
  if (total_signal_size<=classifierAValue) {
    m_return_kit_id=KIT_A;
  }
  else if (total_signal_size>=classifierOValue) {
    m_return_kit_id=KIT_O;
  }
  else {
    m_return_kit_id=KIT_UNKNOWN;
  }

  // clean up resources.
  delete chipStream;

  //printf("KitAOCelListener: id=%d val=%0.6f\n",m_return_kit_id,m_return_kit_value);
}
