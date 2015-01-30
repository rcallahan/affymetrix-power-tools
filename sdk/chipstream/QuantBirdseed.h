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
 * FILE QuantBirdseed.h
 */

#ifndef _QUANTBIRDSEED_H_
#define _QUANTBIRDSEED_H_


#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethodExprReport.h"
//
#include "birdseed-dev/PriorsReader.h"
#include "util/Util.h"
//

class QuantBirdseed : public QuantGTypeMethod {

  public:
    /** Constructor, currently creates prior estimates from pre-calculated data from R. */
  QuantBirdseed(double confidenceThreshold, double correctionFactor);

  virtual ~QuantBirdseed();
  virtual void setQuantExprMethod(QuantExprMethod *qMethod) = 0;
  virtual QuantExprMethod *getQuantExprMethod() = 0;
  virtual void setVerbosity(int v) = 0;
  virtual void setClusterOutFile(const std::string& clusterOutPath);
  
    virtual std::string getModelFile() { return m_ClusterOutFile; }
    virtual void setGenders(const std::vector<affx::Gender> &genders) {
                m_Genders = genders;
    }
    virtual std::vector<affx::Gender> &getGenders() {
        return m_Genders;
    }
	virtual void setPriorFile(const std::string& priorsPath, const std::string &chipType, 
                    std::set<const char *, Util::ltstr> *probeSetsToLoad, 
                    const std::map<std::string, std::pair< int, int > > &specialSnps) = 0;
	
  /** 
   * Get set up for a run of reporting probesets. Often used to open file
   * streams and print headers to files etc.
   * 
   * @param layout - Where the probesets, probes, etc are on the chip.
   * @param iMart  - intensity holder
   * 
   * @return true if success, false otherwise.
   */
  virtual bool prepare(const IntensityMart &iMart);

    /** 
     * @brief Set up the quantification method given all the data about the probe
     * set, chip layout and data.
     * 
     * @param psGroup - Probes to be used for final estimate.
     * @param iMart - Raw data from chips.
     * @param iTrans - Transformations to be applied to data before use.
     * @param pmAdjust - How to estimate background, or MM probe.
     * @return True if setup sucessful, false otherwise.
     */
    virtual bool setUp(ProbeSetGroup &psGroup, 
                       const IntensityMart &iMart, 
                       std::vector<ChipStream *> &iTrans, 
                       PmAdjuster &pmAdjust);
	
    virtual void blankSelf();


    virtual void setParameters(PsBoard &board);

  protected:
    QuantExprMethod *m_QuantMethod; ///< Quantification method used for summarizing snps.
    std::vector<QuantMethodReport *> m_Reporters;   ///< Report object for outputting results
    double m_ConfidenceThreshold;
    double m_CorrectionFactor;
  
    std::vector<double> m_Confidences; ///< Our resulting confidences in those calls.
    std::vector<std::vector<double> > m_Distances; ///< Our resulting distances for AA, AB, and BB genotypes.
  
    std::vector<affx::Gender> m_Genders;     ///< What gender is each sample?

    std::set<std::string> m_ChrXSNPs;

    std::auto_ptr<birdseed::dev::PriorsReader> m_PriorsReader;
    int verbosity;
    std::auto_ptr<std::ofstream> m_ClusterOstrm;
    
    std::string m_ClusterOutFile;

};

#endif /* _QUANTBIRDSEEDV1_H_ */

/******************************************************************/
/**************************[END OF QuantBirdseed.h]****************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
