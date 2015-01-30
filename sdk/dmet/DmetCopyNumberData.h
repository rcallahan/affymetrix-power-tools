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
#ifndef _COPYNUMBERDATA_H_
#define _COPYNUMBERDATA_H_

#include <cstring>
#include<iostream>
#include<map>
#include<set>
#include<string>
#include<vector>
//

class DmetCopyNumberData{

  public:

//    DmetCopyNumberData();
    ~DmetCopyNumberData(){}

// Region specific set and get functions.
    double getN_0(std::string region);
    double getN_1(std::string region);
    double getN_2(std::string region);
   
    double getM_01_0(std::string region); 
    double getM_01_1(std::string region); 
    double getM_12_1(std::string region); 
    double getM_12_2(std::string region); 

    double getS_01_0(std::string region); 
    double getS_01_1(std::string region); 
    double getS_12_1(std::string region); 
    double getS_12_2(std::string region); 

    double getPrior_0(std::string region);
    double getPrior_1(std::string region);
    double getPrior_2(std::string region);

    double getOutlierProbability(std::string region);
    double getNoCallProbability(std::string region);

    void setN_0(std::string region, double newValue);
    void setN_1(std::string region, double newValue);
    void setN_2(std::string region, double newValue);
 
    void setM_01_0(std::string region, double newValue); 
    void setM_01_1(std::string region, double newValue); 
    void setM_12_1(std::string region, double newValue); 
    void setM_12_2(std::string region, double newValue); 

    void setS_01_0(std::string region, double newValue); 
    void setS_01_1(std::string region, double newValue); 
    void setS_12_1(std::string region, double newValue); 
    void setS_12_2(std::string region, double newValue); 

    void setPrior_0(std::string region, double newValue);
    void setPrior_1(std::string region, double newValue);
    void setPrior_2(std::string region, double newValue);

    void setOutlierProbability(std::string region, double newValue);
    void setNoCallProbability(std::string region, double newValue);

    void addRegion(std::string newRegion);


// Probeset specific set and get functions.
    void addPredictionProbeSet(std::string region, std::string probeSet); 
    void addGenotypeProbeSet(const std::string region, const std::string probeset);
    void setW_01(std::string region, double newValue);
    void setW_12(std::string region, double newValue);
    void setMean_0(std::string region, double newValue);
    void setMean_1(std::string region, double newValue);
    void setMean_2(std::string region, double newValue);

    double getW_01(std::string region) const;
    double getW_12(std::string region) const;
    double getMean_0(std::string region) const;
    double getMean_1(std::string region) const;
    double getMean_2(std::string region) const;

    std::set<std::string> getPredictionProbeSetNamesForRegion (std::string regionName);

    std::vector<std::string> getRegionNames();

    std::vector<std::string> getSampleNames();

    std::vector<std::string> getGenotypeProbeSets(const std::string sampleName);

// Intensity set and get functions.
    double getSampleProbeSetIntensity(std::string sampleName, std::string probeSetName) const ;
    void setSampleProbeSetIntensity(std::string sampleName, std::string probeSetName, double intensity);
    void addSampleName(std::string sampleName);

  protected:

  private:

    std::vector<std::string>	m_Regions; 

//Region specific data, structures below are indexed with region name.
    std::map<std::string, double>	m_N_0; 
    std::map<std::string, double>	m_N_1; 
    std::map<std::string, double>	m_N_2; 

    std::map<std::string, double>	m_M_01_0;
    std::map<std::string, double>	m_M_01_1;
    std::map<std::string, double>	m_M_12_1;
    std::map<std::string, double>	m_M_12_2;
    

    std::map<std::string, double>	m_S_01_0;
    std::map<std::string, double>	m_S_01_1;
    std::map<std::string, double>	m_S_12_1;
    std::map<std::string, double>	m_S_12_2;

    std::map<std::string, double> m_Prior_0;
    std::map<std::string, double> m_Prior_1;
    std::map<std::string, double> m_Prior_2;

    std::map<std::string, double>	m_OutlierProbability;
    std::map<std::string, double> m_NoCallProbability;	

    // Map of region names to genotype probeset ids
    std::map<std::string, std::vector<std::string> > m_GenotypeProbeSetMap;
    // This a map from region names to the set of probeset names used
    // for modelling purposes in that region.
    std::map<std::string, std::set<std::string> >	m_PredictionProbeSetMap;	

// ProbeSet specific data.  Structures are indexed by the probeset name.
    std::map<std::string, double>		m_W_01;
    std::map<std::string, double>		m_W_12;
    std::map<std::string, double>		m_Mean_0;
    std::map<std::string, double>		m_Mean_1;
    std::map<std::string, double>		m_Mean_2;


// Intensity data.  This will be read from a summarization file, or in the future directly from 
// summarization compoutation. 

    std::map<std::string, std::map<std::string, double> >	m_IntensityData;	// This is a map from sample names to a map
								// of probeset names to intensities.  
    std::vector<std::string>		m_SampleNames;		// This is the set of sample names for which
								// intensity data exists. 

};

#endif /* _COPYNUMBERDATA_H_ */
