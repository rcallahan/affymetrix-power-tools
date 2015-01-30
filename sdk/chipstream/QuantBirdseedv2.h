////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

/*  
 * FILE QuantBirdseedv2.h
 */

#ifndef _QUANTBIRDSEEDV2_H_
#define _QUANTBIRDSEEDV2_H_


//
#include "birdseed-dev/PriorsReader.h"
#include "chipstream/BioTypes.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantBirdseed.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodExprReport.h"
#include "stats/stats.h"
#include "util/Fs.h"
//
#include "newmat.h"
//
#include <cfloat>
#include <cstring>
#include <string>
#include <vector>
//

/// String describing quant method
#define QUANTBIRDSEEDV2 "birdseed-v2"

class QuantBirdseedv2 : public QuantBirdseed {

  public:

    /** 
     * What version of the BRLMM algorithm are we implementing?
     * @return version string.
     */
    std::string getVersion() {
        return "2.4";
    }

    /** Constructor, currently creates prior estimates from pre-calculated data from R. */
    QuantBirdseedv2(double confidenceThreshold, double correctionFactor): QuantBirdseed(confidenceThreshold, correctionFactor)
    {
        setupSelfDoc(*this);

        setOptValue("correction-factor", ToStr(m_CorrectionFactor));
        setOptValue("conf-threshold", ToStr(m_ConfidenceThreshold));

    }

    ~QuantBirdseedv2()
    {
    }
        
    /** Fill in our self documentation from current state. */
    static void setupSelfDoc(SelfDoc &doc) {
        doc.setDocName(QUANTBIRDSEEDV2);
        doc.setDocDescription("Do genotyping calls using the Birdseed v2 algorithm.");
        doc.setDocOptions(getDefaultDocOptions());
    }
  
    /** 
     * @brief Do the heavy lifting of estimation.
     */
    void computeEstimate();


    /** 
     * @brief What is the name of the quantification method?
     * @return name of adjuster.
     */  
    std::string getType() {
        return std::string(QUANTBIRDSEEDV2);
    }


    /** 
     * @brief Default Getter method for parameters and their documentation.
     * @return map of parameters and their descriptions.
     */
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();

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
     * How many genotyping calls do we have? Also indicates how many
     * samples there are.
     *
     * @return Number of genotyping calls made.
     */
    inline size_t getNumCalls() { return m_Calls.size(); }

    /** 
     * Get our confidence value for a particular call in a particular sample.
     * @param index - sample of interest.
     * @return - Confidence in the predicted call.
     */
    double getConfidence(unsigned int index) { 
        if(index >= m_Confidences.size()) {
            Err::errAbort("Asking for call at index " + ToStr(index) + " when Probeset " + 
                          m_ProbesetName + " has only " + ToStr(m_Confidences.size()) + " calls.");
        }
        return m_Confidences[index];
    }

    /** 
     * Get our confidence value for a particular call in a particular sample.
     * @param genoType -  Genotype of interest.
     * @param index - sample of interest.
     * @return - Confidence in the predicted call.
     */
    double getConfGenotype(affx::GType genoType, unsigned int index) {
        assert(index < m_Distances.size());
        if(genoType == affx::NN) 
            return FLT_MAX;
        return m_Distances[index][genoType];
    }

    /** 
     * Get the summary value after transformation for the A and B alleles respectively.
     * 
     * @param index - Which chip to get values for.
     * @param aValue - Filled in with transformed a value.
     * @param bValue - Filled in with transformed b value.
     */
    void getAlleleValues(unsigned int index, double &aValue, double &bValue) {
        assert(index < m_AValues.size());
        assert(m_AValues.size() == m_BValues.size());
        aValue = m_AValues[index];
        bValue = m_BValues[index];
    }

    /** 
    * Get the summary value names for the A and B alleles respectively.
    * 
    * @param aName - Filled in with a value name.
    * @param bName - Filled in with b value name.
    */
    void getAlleleValueNames(std::string &aName, std::string &bName)
    {
        aName = "Signal A";
        bName = "Signal B";
    }

    /** 
     * Get the name of the probeset that these calls are being made for. 
     * @return name of probeset.
     */
    const std::string &getProbeSetName() { return m_ProbesetName; }

    /** 
     * Add an expression reporter to output the summary values for each
     * allele, residuals, etc.
     * 
     * @param reporter - the object that will output all summary values.
     */
    void addExprReporter(QuantMethodReport *reporter) {
        m_Reporters.push_back(reporter);
    }

    /* Set the quantification method to use for summarizing individual alleles */
    void setQuantExprMethod(QuantExprMethod *qMethod) {
        m_QuantMethod = qMethod;
    }
    /* Get the quantification method to use for summarizing individual alleles */
    QuantExprMethod *getQuantExprMethod() {
        return m_QuantMethod;
    }

    double getConfidenceThreshold(){ return(m_ConfidenceThreshold);};
	void setConfidenceThreshold(double v){m_ConfidenceThreshold=v;};

  // We should test for the format of the file here.
  void setPriorFile(const std::string& priorsPath, const string &chipType, 
                    std::set<const char *, Util::ltstr> *probeSetsToLoad, 
                    const std::map<std::string, std::pair< int, int > > &specialSnps);

  void setVerbosity(int v)
  {
    verbosity = v;
  }
   
//   void setClusterOutFile(const std::string& clusterOutPath)
//   {
//     m_ClusterOutFile = clusterOutPath;
//     m_ClusterOstrm.reset(new ofstream(clusterOutPath.c_str()));
//   }
    
  double getMinThresh() { return -DBL_MAX;}
  double getMaxThresh() { return m_ConfidenceThreshold;}

  virtual bool finish() {
      bool result = true;
      if (m_ClusterOstrm.get() != NULL) {
          Fs::carefulClose(*m_ClusterOstrm);
      }

      for(size_t i=0; i<m_Reporters.size(); i++)
          result = result && m_Reporters[i]->finish(*this);
      return result;
  }
};

#endif /* _QUANTBIRDSEEDV2_H_ */

/******************************************************************/
/**************************[END OF QuantBirdseedv2.h]****************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */

