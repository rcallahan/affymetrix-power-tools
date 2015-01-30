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
 * @file   AnalysisStreamExpression.h
 * @author Chuck Sugnet
 * @date   Mon Feb 13 10:18:56 2006
 * 
 * @brief AnalysisStream class specialized for quantification methods that are
 * used for mRNA expression genechips like the ubiquitous U133 human chip. 
 */

#ifndef _ANALYSISSTREAMEXPRESSION_H_
#define _ANALYSISSTREAMEXPRESSION_H_

#include "chipstream/AnalysisStream.h"
#include "chipstream/QuantExprMethod.h"
//
#include "util/Util.h"
//
#include <set>
//


class AnalysisStreamExpression : public AnalysisStream {

public:
  /** 
   * Constructor.
   * @param doGenoTypes - Do we want to process genotyping probesets.
   * @param doStrand - Do we want to break out genotyping probesets by
   * strand in addition to allele?
   */
  AnalysisStreamExpression(bool doGenoTypes=false, bool doStrand=false, bool aOnly=false);

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions();

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc);

  QuantExprMethod *getQuantExprMethod();

  /** 
   * Set the object for summarizing an individual probe set.
   * @param qMethod - Quantification object.
   */
  void setQuantMethod(QuantMethod *qMethod);

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf();
  
  /** 
   * @brief This static function should be overridden by child classes
   * to return an object of the correct type initialized correctly
   * with the parameters in the string, string map. All objects
   * created this way should be deleted when finished using.
   * 
   * @param param - Map of key/value pairs to initialize the object.
   * 
   * @return Pointer toCreate object, this should be sub casted as necessary.
   */
  static SelfCreate *newObject(std::map<std::string,std::string> &param);

  /** 
   * Fill in psNewGroup with probes from psGroup based on the block
   * indexes in the blocks vector. To divide a genotyping probeset
   * based on Allele (A or B) you would use the 1,3 blocks for A and
   * 2,4 blocks for B. If you just want one strand of a particular
   * allele then it would be 1 for A+, 2 for B+, 3 for A- and 4 for
   * B-. This is sensitive to the vagaries of the genotyping probeset
   * structure which is poorly documented and ever
   * "innovating"... Eventually this should be moved to the ChipLayout
   * class or ProbeSet class.
   * 
   * @param psNewGroup 
   * @param psGroup 
   * @param blocks 
   */
  void fillInProbeSetFromBlocks(ProbeSetGroup &psNewGroup, ProbeSetGroup &psGroup,
                                const std::vector<short> &blocks);

  void fillInGroup(ProbeSetGroup &psAllele, ProbeSetGroup &psGroup, int blockIx, const char *basename, const char *suffix);

  /** 
   * Do the analysis for a particular group of probe sets. If we see a
   * genotype probeset and user wishes it try to break it up into
   * expression probesets. The code for dividing up a genotyping
   * probeset should probably live in ChipLayout eventually, but for
   * now the initial implementation is here.
   * 
   * @param psGroup - Collection of probe sets to get probes from.
   * @param iMart - Object containing raw data values for all chips.
   * @param doReport - Should the quantification report object be called?
   * 
   * @return true if success, false otherwise.
   */
  virtual bool doAnalysis(ProbeSetGroup &psGroup, IntensityMart &iMart, bool doReport, bool alleleSummaryOnly = false);

  /** 
   * Make a new probeset group which is a subset of the original based on the
   * probeset ids contained in the goodIds set.
   * 
   * @param selectGroup - New probeset group to be filled in.
   * @param goodIds - Set with ids of probes to keep.
   * @param psGroup - Original probeset to get probes and structure from.
   */
  static void makeSelectProbesetGroup(ProbeSetGroup &selectGroup, 
                                      std::set<probeid_t> &goodIds,
                                      ProbeSetGroup &psGroup);
private:
  QuantExprMethod *m_QExprMethod;
  bool m_DoGenotypes; ///< Should we break genotyping probesets into different alleles and summarize them?
  bool m_DoStrand; ///< Should we separate genotyping probesets by strand in addition to block
  bool m_Aonly; ///< Should only the A allele be processed for
                ///genotyping probesets? Useful for pm adjuster pm-sum
                ///where A and B allele give exact same results.
};

#endif /* _ANALYSISSTREAMEXPRESSION_H_ */
