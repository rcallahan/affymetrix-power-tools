////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
 * @file   QuantLabelZMulti.h
 * @author Martin Gilchrist 
 * @date   Mon Jun 23 09:17:00 2008
 * 
 * @brief  Interface for methods computing SNP genotyping calls.
 */


#ifndef QUANTLABELZMULTI_H
#define QUANTLABELZMULTI_H

//#include "chipstream/SelfCreate.h" 
/*
#include "algorithm/covarnorm/covarnorm.h"
#include "algorithm/em/PrimeEM.h"
#include "algorithm/selector/ProbeSelector.h"
#include "chipstream/BioTypes.h"
#include "chipstream/ClusterZ.h"
#include "chipstream/GenoUtility.h"

#include "chipstream/ProbeSet.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantGTypeMethod.h"

#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodExprReport.h"
#include "file5/File5.h"
#include "label/snp.label.h"
#include "newmat.h"
#include "stats/stats.h"
#include "util/AffxArray.h"
*/
#include "chipstream/QuantLabelZ.h"
//
#include "calvin_files/utils/src/GenoCallCoder.h"
//
#include <cfloat>
#include <cstring>
#include <map>
#include <string>
#include <vector>
//
#define QUANTBRLMMPMULTI "brlmm-p-multi"

typedef vector<double> SampleSummaries;
typedef map<int, SampleSummaries> ContextMap;
typedef map<int, ContextMap> AlleleMap;

typedef pair<int, int> AlleleContext;
typedef pair<AlleleContext, double> ASummary;
typedef vector<ASummary> AlleleVector;
typedef vector<ASummary>  MaxContextVector;
typedef vector<MaxContextVector> SampleMaxs;

class LessThan : public binary_function<ASummary, ASummary, bool> {
  public:
    bool operator()( const ASummary &summary1, const ASummary &summary2)
      {return summary1.second > summary2.second;}
};


class AlleleLessThan : public binary_function<ASummary, ASummary, bool> {
  public:
    bool operator()( const ASummary &summary1, const ASummary &summary2)
      {return summary1.first.first < summary2.first.first;}
};

class QuantLabelZMulti : public QuantLabelZ {

  public:

    // Constructor:  We pass the information into the base class QuantLabelZ via the initialization list.
    QuantLabelZMulti(enum Transformation transform, double K, bool lowPrecision, int praThresh, int ccMaxAlleles, string ccType, string ccVersion) 
        : QuantLabelZ (transform, K, lowPrecision) {

        m_coder = new GenoCallCoder(ccMaxAlleles,ccType,ccVersion, '\0');  

        // Setup self doc
        setupSelfDoc(*this);
        setOptValue("transform", transform);
        setOptValue("K",ToStr(K));
        setOptValue("lowprecision", lowPrecision);
        setOptValue("pra-thresh", ToStr(praThresh));
        setOptValue("cc-alleles", ToStr(ccMaxAlleles));
        setOptValue("cc-type", ccType);
        setOptValue("cc-version", ccVersion);

        m_PraThresh = praThresh;
    }  

    virtual ~QuantLabelZMulti() {
      delete m_coder;
    } 

    static SelfCreate *newObject(std::map<std::string,std::string> &param); 

    /** 
    * @brief Supply a little how/what/why about the algorithms this
     * class performs and what parameters it takes.
     * @return SelfDoc
     */
    static SelfDoc explainSelf();

    void getSequentialPrior(std::string &name, snp_param &tsp);
   
    /** Fill in our self documentation from current state. */
    static void setupSelfDoc(SelfDoc &doc);


    /*  The functions setUp and computeEstimate are two functions that are required of a QuantMethod.

        In QuantLabelZ setUp performs 3 things: 
            1) creates individual probeGroups for the alleles (only 2 alleles) 
            2) summarizes the probes  for each allele.
            3) creates the contrast/strength values.

        These constrast/strength values are passed to the computeEstimate function which does the genotyping.

        In QuantGTypeMultiMethod 1 and 2 will be performed in setUp and the 3rd will be done in computeEstimate
        along with the genotyping.  The reason for this is that we only need to perform contrast/strength computations
        on the summary values which correspond to maximum contexts, as these are the only alleles genotyped. These
        max contexts cannot be determined until all allele/contexts have been summarized and thus to avoid storing
        all contrast/strength values we calculate only those needed after all summarization has been done.

        As well, we will not create member data structures to hold the probeSets for each allele/context since
        the actual probeSets are only needed to compute Summary values, and thus they can be created in memory
        and destroyed after summarization.  (The actual probeSets in QuantLabelZ ARE stored but I don't think this
        is necessary.) 
    */ 


    bool setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
               std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust);

    void getSequentialPrior(std::string &name, snp_param &tsp, snp_labeled_distribution &sDist);

    void getPrior(std::string &TmpName, int presentCopyNumber, snp_param &tsp, snp_labeled_distribution &sDist);

    void computeEstimate();

    bool SetUpProbeSet(const ProbeSet *gtPs,
                        const IntensityMart &iMart, std::vector<ChipStream *> &iTrans,
                        PmAdjuster &pmAdjust, bool doReport);
  

  void setCopyNumberMap(std::map<std::string, std::vector<int>  > &copyNumberMap){
      ///@todo should we just store a pointer to keep memory use down?
      m_CopyNumberMap = copyNumberMap;
  }

    /** 
     * @brief What is the name of the quantification method?
     * @return name of adjuster.
     */  
    std::string getType() {
      return std::string(QUANTBRLMMPMULTI);
    }


  /**
    @brief  Gets the integer valued (not affx::GType) forced call for the probeset. 
   */
  virtual int getForcedCall(unsigned int index);

  /**
    @brief  Gets the integer valued (not affx::GType) call for the probeset. 
   */
  virtual int getCall(unsigned int index);

  /** 
   * Get the genotype forced call at specified index (sample).
   * @param index - sample of interest.
   * @return - Genotyping call made.
   */
  virtual affx::GType getGTypeForcedCall(unsigned int index);

  /** 
   * Get the genotype call at specified index (sample).
   * @param index - sample of interest.
   * @return - Genotyping call made.
  */
  virtual affx::GType getGTypeCall(unsigned int index);

   /** 
    * How many genotyping calls do we have? Also indicates how many
    * samples there are.
    *
    * @return Number of genotyping calls made.
    */
    virtual inline size_t getNumCalls() { return m_CallsMulti.size(); }

    virtual string getContext(unsigned int index) {
      if(index >= m_Context.size()) {
        Err::errAbort("Asking for context at sample index " + ToStr(index) + " when Probeset " + 
                      m_ProbesetName + " has only " + ToStr(m_Context.size()) + " samples.");
      }
      return m_Context[index];
    }


  public:
    GenoCallCoder*                              m_coder;  

  protected:



  private:

    int translateBrlmmpOutputToCall(int copyNumber, int brlmmpCall,
            int firstAllele, int secondAllele );
 
    void determineCopyNumber();
 
    SampleMaxs findMaxContextForEachSample(int numberOfSamples); 
  
    void amalgamateAtoms(const ProbeSet* gtPs, ProbeSet &amalgamatedProbeset); 

  private:
    AlleleMap                                   m_Summaries; 
    std::map<std::string, std::vector<int> >    m_CopyNumberMap; 
    int                                         m_numberOfSamples;    
    vector<int>                                 m_CallsMulti;    
    vector<int>                                 m_CopyNumbers;    
    std::vector<std::string>                    m_Context;
    vector<double>                              m_PriorObs;
    int                                         m_PraThresh;

}; 


#endif /* QUANTLABELZMULTI_H */

