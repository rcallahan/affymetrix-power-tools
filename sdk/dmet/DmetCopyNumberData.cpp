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
#include "dmet/DmetCopyNumberData.h"
//
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Verbose.h"
//
#include <iostream>
//

using namespace std;

//Region specific set and get functions.

double DmetCopyNumberData::getN_0(string region){
  map<string, double>::iterator iter=m_N_0.find(region);
  if( iter != m_N_0.end()){
    return iter->second;
  }
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getN_1(string region){
  map<string, double>::iterator iter=m_N_1.find(region);
  if( iter != m_N_1.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getN_2(string region){
  map<string, double>::iterator iter=m_N_2.find(region);
  if( iter != m_N_2.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    


double DmetCopyNumberData::getM_01_0(string region){
  map<string, double>::iterator iter=m_M_01_0.find(region);
  if( iter != m_M_01_0.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getM_01_1(string region){
  map<string, double>::iterator iter=m_M_01_1.find(region);
  if( iter != m_M_01_1.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getM_12_1(string region){
  map<string, double>::iterator iter=m_M_12_1.find(region);
  if( iter != m_M_12_1.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getM_12_2(string region){
  map<string, double>::iterator iter=m_M_12_2.find(region);
  if( iter != m_M_12_2.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getS_01_0(string region){
  map<string, double>::iterator iter=m_S_01_0.find(region);
  if( iter != m_S_01_0.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getS_01_1(string region){
  map<string, double>::iterator iter=m_S_01_1.find(region);
  if( iter != m_S_01_1.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getS_12_1(string region){
  map<string, double>::iterator iter=m_S_12_1.find(region);
  if( iter != m_S_12_1.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getS_12_2(string region){
  map<string, double>::iterator iter=m_S_12_2.find(region);
  if( iter != m_S_12_2.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getPrior_0(string region){
  map<string, double>::iterator iter=m_Prior_0.find(region);
  if( iter != m_Prior_0.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}

double DmetCopyNumberData::getPrior_1(string region){
  map<string, double>::iterator iter=m_Prior_1.find(region);
  if( iter != m_Prior_1.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}


double DmetCopyNumberData::getPrior_2(string region){
  map<string, double>::iterator iter=m_Prior_2.find(region);
  if( iter != m_Prior_2.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}



double DmetCopyNumberData::getOutlierProbability(string region){
  map<string, double>::iterator iter=m_OutlierProbability.find(region);
  if( iter != m_OutlierProbability.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    

double DmetCopyNumberData::getNoCallProbability(string region){
  map<string, double>::iterator iter=m_NoCallProbability.find(region);
  if( iter != m_NoCallProbability.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}    


void DmetCopyNumberData::setN_0(string region, double newValue){ m_N_0[region] = newValue;}
void DmetCopyNumberData::setN_1(string region, double newValue){ m_N_1[region] = newValue;}
void DmetCopyNumberData::setN_2(string region, double newValue){ m_N_2[region] = newValue;}

void DmetCopyNumberData::setM_01_0(string region, double newValue){ m_M_01_0[region] = newValue;} 
void DmetCopyNumberData::setM_01_1(string region, double newValue){ m_M_01_1[region] = newValue;} 
void DmetCopyNumberData::setM_12_1(string region, double newValue){ m_M_12_1[region] = newValue;} 
void DmetCopyNumberData::setM_12_2(string region, double newValue){ m_M_12_2[region] = newValue;} 

void DmetCopyNumberData::setS_01_0(string region, double newValue){ m_S_01_0[region] = newValue;} 
void DmetCopyNumberData::setS_01_1(string region, double newValue){ m_S_01_1[region] = newValue;} 
void DmetCopyNumberData::setS_12_1(string region, double newValue){ m_S_12_1[region] = newValue;} 
void DmetCopyNumberData::setS_12_2(string region, double newValue){ m_S_12_2[region] = newValue;} 

void DmetCopyNumberData::setPrior_0(string region, double newValue){ m_Prior_0[region] = newValue;} 
void DmetCopyNumberData::setPrior_1(string region, double newValue){ m_Prior_1[region] = newValue;} 
void DmetCopyNumberData::setPrior_2(string region, double newValue){ m_Prior_2[region] = newValue;} 

void DmetCopyNumberData::setOutlierProbability(string region, double newValue){ m_OutlierProbability[region] = newValue;}
void DmetCopyNumberData::setNoCallProbability(string region, double newValue){ m_NoCallProbability[region] = newValue;}

void DmetCopyNumberData::addRegion(string newRegion){
    m_Regions.push_back(newRegion); 
}


// Probeset specific set and get functions.
void DmetCopyNumberData::addPredictionProbeSet(string region, string probeSet){
  set<string> newSetOfProbes;
  map<string, set<string> >::iterator iter=m_PredictionProbeSetMap.find(region);
  if(iter!=m_PredictionProbeSetMap.end()){
    (iter->second).insert(probeSet);
  } else {
    // We have a new region
    newSetOfProbes.insert(probeSet);
    m_PredictionProbeSetMap.insert(	map<string, set<string> >::
				value_type(region, newSetOfProbes ) );
  } 
}      

void DmetCopyNumberData::addGenotypeProbeSet(string region, string probeset){
  //  cout<< "set:" + region + '\n';
  m_GenotypeProbeSetMap[region].push_back(probeset);
}      


void DmetCopyNumberData::setW_01(string region, double newValue){ m_W_01[region] = newValue;}
void DmetCopyNumberData::setW_12(string region, double newValue){ m_W_12[region] = newValue;}
void DmetCopyNumberData::setMean_0(string region, double newValue){ m_Mean_0[region] = newValue;}
void DmetCopyNumberData::setMean_1(string region, double newValue){ m_Mean_1[region] = newValue;}
void DmetCopyNumberData::setMean_2(string region, double newValue){ m_Mean_2[region] = newValue;}

double DmetCopyNumberData::getW_01( string region) const {
  map<string, double>::const_iterator iter=m_W_01.find(region);
  if( iter != m_W_01.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}

double DmetCopyNumberData::getW_12( string region) const {
  map<string, double>::const_iterator iter=m_W_12.find(region);
  if( iter != m_W_12.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}


double DmetCopyNumberData::getMean_0( string region) const {
  map<string, double>::const_iterator iter=m_Mean_0.find(region);
  if( iter != m_Mean_0.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}


double DmetCopyNumberData::getMean_1( string region) const {
  map<string, double>::const_iterator iter=m_Mean_1.find(region);
  if( iter != m_Mean_1.end())
    return iter->second;
  Err::errAbort("Unable to get region info");
  return 0.0;
}


double DmetCopyNumberData::getMean_2( string region) const {
  map<string, double>::const_iterator iter=m_Mean_2.find(region);
  if( iter != m_Mean_2.end()){
    return iter->second;
  } 
  Err::errAbort("Unable to get region info");
  return 0.0;
}





set<string> DmetCopyNumberData::getPredictionProbeSetNamesForRegion(string regionName){
    map<string, set<string> >::iterator iter=m_PredictionProbeSetMap.find(regionName);
    if( iter != m_PredictionProbeSetMap.end()){
        return iter->second;
    } 
    Err::errAbort("Unable to find probeset list for region " + regionName);
    set<string> nullSet;
    return nullSet;
}

vector<string> DmetCopyNumberData::getRegionNames(){
    return m_Regions;
}

vector<string> DmetCopyNumberData::getSampleNames(){
  return m_SampleNames;

}

vector<string> DmetCopyNumberData::getGenotypeProbeSets(const std::string regionName) {
  //  cout << "get:" + regionName + '\n';
    return m_GenotypeProbeSetMap[regionName];
}


double DmetCopyNumberData::getSampleProbeSetIntensity(string sampleName, string probeSetName) const {
  map<string, map<string, double> >::const_iterator iter1=m_IntensityData.find(sampleName);
  if(iter1!=m_IntensityData.end()){
    map<string,double>::const_iterator iter2=(iter1->second).find(probeSetName);
    if( iter2!=(iter1->second).end() ){
      return iter2->second;
    } else {
      throw "Unable to find summary value for probeset " + probeSetName;
    }
  } else {
    throw "Unable to find summary data for sample " + sampleName;
  }

  throw ("Reached unexpected part of code.");
  return 0.0;
}


void DmetCopyNumberData::setSampleProbeSetIntensity(string sampleName, 
						string probeSetName, 
						double intensity){

  //Verbose::out(1,"Setting data: " + sampleName + ", " + probeSetName + ", " + ToStr(intensity));

  typedef map<string,double>::value_type mapValueType;
  typedef map<string, map<string,double> >::value_type mapOfMapValueType;

  string markerName = probeSetName.substr(0,probeSetName.rfind("-"));
  string alleleName = probeSetName.substr(probeSetName.rfind("-"));

  map<string,double> intensitiesForANewSample; 
  map<string, map<string,double> >::iterator iter=m_IntensityData.find(sampleName);
  if(iter!=m_IntensityData.end()){
    if((iter->second).find(markerName) != (iter->second).end()){
        // we already have data for one or more alleles
        m_IntensityData[sampleName][markerName] += intensity;
    } else {
      // We have a new probeset for a known region.
      (iter->second).insert( mapValueType(mapValueType(markerName, intensity))); 
    }
  } else { 
    // We have a new region 
    intensitiesForANewSample.insert(mapValueType(markerName,intensity)); 
    m_IntensityData.insert(mapOfMapValueType(sampleName, intensitiesForANewSample));
    intensitiesForANewSample.erase(	intensitiesForANewSample.begin(),
  					intensitiesForANewSample.end()); 
  }
}



void DmetCopyNumberData::addSampleName(string sampleName){
  m_SampleNames.push_back(sampleName);
}

