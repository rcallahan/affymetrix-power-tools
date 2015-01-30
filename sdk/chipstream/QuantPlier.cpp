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

/**
 * @file   QuantPlier.cpp
 * @author Chuck Sugnet
 * @date   Mon Oct 24 10:56:04 2005
 * 
 * @brief  Class for doing probe set quantification using the PLIER (Probe
 * Logarithmic Error Intensity Estimate) method.
 */

//
#include "chipstream/QuantPlier.h"

using namespace std;

/** Constructor. */
QuantPlier::QuantPlier(QuantPlierParams &plierParams)
  : QuantPlierBase(plierParams) {
  setupSelfDoc(*this);
  m_Type = getDocName();
  // set the pointer to NULL
  zeroMemoryVars();
  //
  setUsePrecompFeatureEffects(false);
}

/** 
 * Destructor
 */
QuantPlier::~QuantPlier() {
  freeMemory();
}

void QuantPlier::zeroMemoryVars() {
  m_MaxChips=0;
  m_MaxProbes=0;
  m_ProbeEffects=NULL;
  m_ChipEffects=NULL;
  m_Residuals=NULL;
  m_PM=NULL;
  m_MM=NULL;
}

/** 
 * @brief Free up all the memory that has been allocated.
 */
void QuantPlier::freeMemory() {
  unsigned int i = 0;
  delete [] m_ProbeEffects;
  delete [] m_ChipEffects;
  for(i = 0; i < m_MaxChips; i++) {
    delete [] m_PM[i];
    delete [] m_MM[i];
    delete [] m_Residuals[i];
  }
  delete [] m_Residuals;
  delete [] m_PM;
  delete [] m_MM;
  // Now zero the vars in case freeMemory is called again.
  zeroMemoryVars();
}

void QuantPlier::allocMemory(unsigned int maxChips, unsigned int maxProbes) {
  freeMemory();

  unsigned int i = 0;
  m_MaxChips = maxChips;
  m_MaxProbes = maxProbes;

  m_ProbeEffects = new double[m_MaxProbes];
  m_ChipEffects = new double[m_MaxChips];

  m_PM = new double *[m_MaxChips];
  m_MM = new double *[m_MaxChips];
  m_Residuals = new double *[m_MaxChips];
  for(i = 0; i < m_MaxChips; i++) {
    m_PM[i] = new double[m_MaxProbes];
    m_MM[i] = new double[m_MaxProbes];
    m_Residuals[i] = new double[m_MaxProbes];
  }
}

/** 
 * @brief Do the heavy lifting of estimation.
 */
void QuantPlier::computeEstimate() {
  long errorCode = 0;
  m_Plier.setNumExp(m_ChipCount);
  m_Plier.setNumFeature(m_ProbeCount);
  m_Plier.setPM(m_PM);
  m_Plier.setMM(m_MM);
  m_Plier.setResiduals(m_Residuals);
  m_Plier.setTargetResponse(m_ChipEffects);
  m_Plier.setFeatureResponse(m_ProbeEffects);
  m_Plier.run(&errorCode);
  if(errorCode != 0) 
    Err::errAbort("Problem running plier. Error code: " + 
                    ToStr(errorCode));
}

/** 
 * @brief Return number of features (probes).
 * @return feature count.
 */
inline unsigned int QuantPlier::getNumFeatures()  { return m_Probes.size(); }

/** 
 * Fill in the information for Self documentation.
 * @param doc - Self documenter to be filled in.
 */
void QuantPlier::setupSelfDoc(SelfDoc &doc) {
  doc.setDocName(QUANTPLIERSTR);
  doc.setDocDescription("The PLIER (Probe Logarithmic Error Intensity Estimate) method produces an improved signal by accounting for experimentally observed patterns in feature behavior and handling error at the appropriately at low and high signal values. This version of PLIER differs from the previous version by the addition of a SafteyZero, NumericalTolerance, and FixPrecomputed. These options are intended to improve the stability of PLIER results when using precomputed feature reponse values. To get the older PLIER behavior set SafetyZero to 0.0, NumericalTolerance to 0.0, and FixPrecomputed to false.");
  /* Set common parameters and documentation used by both QuantPlier and QuantIterPlier. */
  doc.setDocOptions(getDefaultDocOptions());
}

/** 
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.  The default
 * parameters are from the inputs.txt file in the sdk example.
 * @return SelfDoc
 */
SelfDoc QuantPlier::explainSelf() { 
  SelfDoc doc;
  setupSelfDoc(doc);
  return doc;
}

/** 
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 * The default values used for the plier parameters are from the
 * inputs.txt file in the sdk example.
 * 
 * @param param - Map of key/value pairs to initialize the object.
 * 
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate *QuantPlier::newObject(std::map<std::string,std::string> &param) {
  QuantPlierParams plierParams;
  SelfDoc doc = explainSelf();
  fillPlierParams(param, plierParams, doc);
  QuantPlier *plier = new QuantPlier(plierParams);
  plier->setDocOptions(doc.getDocOptions());
  return plier;
}

/** 
 * @brief Specify precomputed feature effects to be used as map from
 * probe id index to the feature effect.  This make plier much
 * faster to do as doesn't have to iterate through anything, but
 * have to make sure that the feature effects being used match the
 * analysis that is being done. For example the feature effects
 * might be different using pm-only, pm-gcbg and pm-mm adjustments.
 * 
 * @param featureEffects - mapping from probe idex to feature effects.
 */
void QuantPlier::setFeaturePriorEffects(double* featureEffects, int iSize) {
  Verbose::out(3, "QuantPlier::setFeaturePriorEffects()");
  if (m_FeaturePriorEffect != NULL) {delete[] m_FeaturePriorEffect;}
  m_FeaturePriorEffect = new double[iSize];
  for (int iIndex = 0; (iIndex < iSize); iIndex++)
    {
      m_FeaturePriorEffect[iIndex] = featureEffects[iIndex];
    }
  setUsePrecompFeatureEffects(true);
}
