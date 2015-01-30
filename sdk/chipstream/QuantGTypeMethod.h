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
 * @file   QuantGTypeMethod.h
 * @author Chuck Sugnet
 * @date   Fri Feb 24 13:27:20 2006
 * 
 * @brief  Interface for methods computing SNP genotyping calls.
 */
#ifndef QUANTGTYPEMETHOD_H
#define QUANTGTYPEMETHOD_H

//
#include "chipstream/BioTypes.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodReport.h"
//
#include <cstring>
#include <string>
//

/**
 * Quantification methods used for making genotyping calls implement (currently
 * just brlmm) this interface.
 */
class QuantGTypeMethod : public QuantMethod {

public:
  virtual ~QuantGTypeMethod();

  /** 
   * How many genotyping calls do we have? Also indicates how many
   * samples there are.
   *
   * @return Number of genotyping calls made.
   */
  virtual size_t getNumCalls() = 0;
  
  /** 
   * Get the genotype call at specified index (sample).
   * @param index - sample of interest.
   * @return - Genotyping call made.
   */

  virtual affx::GType getGTypeCall(unsigned int index);

  virtual affx::GType getGTypeForcedCall(unsigned int index);

  /** 
   * Get the genotype call at specified index (sample).
   * @param index - sample of interest.
   * @return - Genotyping call made.
   */

  virtual int getForcedCall(unsigned int index);

  /** 
   * Get the genotype call at specified index (sample).
   * @param index - sample of interest.
   * @return - Genotyping call made.
   */

  virtual int getCall(unsigned int index);

  /** 
   * Get our confidence value for a particular call in a particular sample.
   * @param index - sample of interest.
   * @return - Confidence in the predicted call.
   */
  virtual double getConfidence(unsigned int index) = 0;

  /** 
   * Get our confidence value for a particular call in a particular sample.
   * @param genoType -  Genotype of interest.
   * @param index - sample of interest.
   * @return - Confidence in the predicted call.
   */
  virtual double getConfGenotype(affx::GType  genoType, unsigned int index) = 0;

  /** 
   * Get the summary value after transformation for the A and B alleles respectively.
   * 
   * @param index - Which chip to get values for.
   * @param aValue - Filled in with transformed a value.
   * @param bValue - Filled in with transformed b value.
   */
  virtual void getAlleleValues(unsigned int index, double &aValue, double &bValue) = 0;

  /** 
   * Get the summary value names for the A and B alleles respectively.
   * 
   * @param aName - Filled in with a value name.
   * @param bName - Filled in with b value name.
   */
  virtual void getAlleleValueNames(std::string &aName, std::string &bName) = 0;

  /** 
   * Get the name of the probeset that these calls are being made for. 
   * @return name of probeset.
   */
  virtual const std::string &getProbeSetName() = 0;

  /** 
   * Get a string representation for a particular SNP.
   * @return - a string version of a SNP
   */
  virtual std::string getModelString();

  /**
   * Get the max value for a call
   */
  virtual double getMaxThresh() = 0;

  /**
   * Get the min value for a call
   */
  virtual double getMinThresh() = 0;

  /** 
   * Add an expression reporter to output the summary values for each
   * allele, residuals, etc.
   * 
   * @param reporter - the object that will output all summary values.
   */
  virtual void addExprReporter(QuantMethodReport *reporter) = 0;

  virtual std::string getModelFile() {
      APT_ERR_ABORT("Model file not suppored.");
      return "";
  }

protected:
  /** clear a probe set*/
  void clearProbeSet(ProbeSet &ps); 

  /** fill in the probe set with the allele values */
  bool fillInAlleleProbeSets(const ProbeSet &gtPs, ProbeSet &aAllele, ProbeSet &bAllele); 

  /** summarize this allele */
  bool summarizeAllele(ProbeSet *pSet,
                       std::vector<double> &summaryValues,
                       const IntensityMart &iMart,
                       std::vector<ChipStream *> &iTrans, 
                       PmAdjuster &pmAdjust, 
                       QuantExprMethod *quantMethod,
                       bool lowPrecision,
                       bool doReport, 
                       std::vector<QuantMethodReport *> reporters);

  /** summarize this allele - no reports */
  bool summarizeAllele(ProbeSet *pSet,
                       std::vector<double> &summaryValues,
                       const IntensityMart &iMart,
                       std::vector<ChipStream *> &iTrans, 
                       PmAdjuster &pmAdjust, 
                       QuantExprMethod *quantMethod,
                       bool lowPrecision);

  bool extractProbes(std::vector <std::vector <double> > &probeIntensity, 
                     ProbeSet *pSet,
                     std::vector< unsigned int > &probeIds,
                     const IntensityMart &iMart, 
                     std::vector<ChipStream *> &iTrans, 
                     PmAdjuster &pmAdjust, 
                     bool lowPrecision, 
                     QuantExprMethod *quantMethod); 
  
  bool canSetUpProbeSet(const ProbeSet* gtPs);

protected:
  std::vector<affx::GType> m_Calls;
  const ProbeSet *m_GtProbeSet; ///< Genotyping probeset currently using.
  std::string m_ProbesetName; ///< Name of current probeset.
  ProbeSet m_Aallele; ///< Probeset to detect the A allele
  ProbeSet m_Ballele; ///< Probeset to detect the B allele
  std::vector<double> m_AValues;  ///< Our summarized values for A allele intensity
  std::vector<double> m_BValues;  ///< Our summarized values for B allele intensity
};

#endif /* QUANTGTYPEMETHOD_H */
