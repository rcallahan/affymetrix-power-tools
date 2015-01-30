////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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
 * @file   DiskIntensityMart.cpp
 * @author Chuck Sugnet
 * @date   Mon Apr  7 07:15:21 PDT 2008
 * 
 * @brief Class for encapsulating microarray data and storing it on disk in a
 * reasonable way to reduce memory impact. Key idea is to reorganize the data on
 * disk in the order which it will be used by the application to avoid lots of
 * seeks. Hopefully with a small number of seeks the performance impact of
 * on-disk stoarage will be minimal.
 */

//
#include "chipstream/DiskIntensityMart.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "file/FileWriter.h"
#include "file5/File5.h"
#include "util/Fs.h"
#include "util/TmpFileFactory.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstdio>

using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;
using namespace std;
using namespace affx;

// int DiskIntensityMart::m_cache_misses = 0;

DiskIntensityMart::DiskIntensityMart(const std::vector<int>& layoutOrder,
                                     const std::vector<std::string>& celNames,
                                     int cacheSize, 
                                     const std::string& tempDir,
                                     const std::string& tempPrefix,
				     bool storeAllCelIntensities)
{
  m_CacheProbeSize = cacheSize;
  m_CacheStart = -1;
  m_CacheEnd = -1;
  m_Cache.resize(celNames.size());  // not necessary to set this here
  m_DeleteFilesWhenFinished = true;
  m_TempDir = tempDir;
  m_CelFiles = celNames;
  m_Flushed = false;
  m_NumChips = 0;
  m_TempPrefix = tempPrefix;
  m_StoreAllCelIntensities = storeAllCelIntensities;

  m_File5 = NULL;
  m_File5Name = "";

  if ( !Fs::dirExists(m_TempDir) ) {
    Fs::mkdirPath(m_TempDir, false);
  }
  GlobalTmpFileFactory()->setTmpdir(m_TempDir);
  m_File5Name = GlobalTmpFileFactory()->genFilename(tempPrefix, ".a5");

  m_AnalysisOrder = layoutOrder;
  m_ExpectedChips = celNames.size();
  m_Size = 0;
  m_UniqueAnalysisOrderSize = 0;
  m_useAuxMemCache = false;
}

/**
 * @brief Basic destructor
 */
DiskIntensityMart::~DiskIntensityMart() {
  for (int i = 0; i < m_TmpVectors.size(); i++) {
    m_TmpVectors[i]->close();
    delete m_TmpVectors[i];
  }  
  deleteTmpfile();
}

/**
 * @brief deep copy constructor -- NOT IN CURRENT USE.  UNTESTED. 
 *
 * @param diskMart - DiskIntensityMart to be copied
 */
DiskIntensityMart::DiskIntensityMart(DiskIntensityMart& diskMart) {
  m_TempDir = diskMart.m_TempDir;
  m_TempPrefix = diskMart.m_TempPrefix;
  m_CacheProbeSize = diskMart.m_CacheProbeSize;
  m_CacheStart = -1;
  m_CacheEnd = -1;
  m_Cache.resize(diskMart.m_CelFiles.size());
  m_DeleteFilesWhenFinished = diskMart.m_DeleteFilesWhenFinished;
  m_CelFiles = diskMart.m_CelFiles;
  m_Flushed = diskMart.m_Flushed;
  m_NumChips = diskMart.m_NumChips;
  m_StoreAllCelIntensities = diskMart.m_StoreAllCelIntensities;
  m_File5 = diskMart.m_File5;
  m_File5Name = diskMart.m_File5Name;
  m_AnalysisOrder = diskMart.m_AnalysisOrder;
  m_Order = diskMart.m_Order;
  m_ExpectedChips = diskMart.m_ExpectedChips;
  m_Map = diskMart.m_Map;
  m_ChipChannelToCacheMap = diskMart.m_ChipChannelToCacheMap;
  m_CacheProbeSize = diskMart.m_CacheProbeSize;
  m_Size = diskMart.m_Size;
  m_UniqueAnalysisOrderSize = 0;
  m_AuxMemCache = diskMart.m_AuxMemCache;

  for (int i = 0; i < diskMart.m_NumChips; i++) {
    std::vector<float> temp = diskMart.getCelData(i);
    setProbeIntensity(i, temp);
  }
}



/**
 * @brief Return number of CEL channels stored in diskmart
 */
int DiskIntensityMart::getCelDataSetCount() const {
  return m_NumChips;
}

/**
 * @brief Return number of multi-CEL files stored in diskmart
 */
int DiskIntensityMart::getCelFileCount() const {
  return m_CelFiles.size();
}

float DiskIntensityMart::fromCache(probeid_t pIx, chipid_t chipIx) const {
  if(m_Map[pIx] < 0) {
    Err::errAbort("DiskIntensityMart::fromCache() - Illegal value in map for probeIx: " + ToStr(pIx));
  }
  if (m_TmpVectors[chipIx] == NULL) {
    Err::errAbort("DiskIntensityMart::fromCache() - Illegal read: No data has been written for chipIx: " +ToStr(chipIx));
  }
  float val = m_Cache[chipIx][m_Map[pIx] - m_CacheStart];
  return val;
}

bool DiskIntensityMart::inCache(probeid_t pIx, chipid_t chipIx) const {
  if(m_Map[pIx] >= m_CacheStart && m_Map[pIx] < m_CacheEnd) {
    return true;
  }
  return false;
}



/**
 * @brief set object that will translate from chip,channel to cache index
 *
 * @param chip_channel_map - First dimension is chipIx, second is
 * channelIx, value is cache index
 */
void DiskIntensityMart::setChipChannelMap(std::vector<std::vector<unsigned int> >& chip_channel_map) {
  m_ChipChannelToCacheMap = chip_channel_map;
}

int DiskIntensityMart::getChannelCount() const {
  if (m_ChipChannelToCacheMap.empty())
    return(1); // 1 channel
  else
    return(m_ChipChannelToCacheMap[0].size()); // number of channels assumed same for all chips
}


const std::vector<std::string> &DiskIntensityMart::getCelFileNames() const {
  return m_CelFiles;
}


/** 
 * @brief Method for setting boolean to indicate if diskMart should store all given intesities or only the ones specified in m_Order.
 * @param flag - boolean indicating desired behavior
 */
void DiskIntensityMart::setStoreAllCelIntensities(bool flag) {
  m_StoreAllCelIntensities = flag;
}


size_t DiskIntensityMart::getProbeCount() const {
  return m_Map.size();
}

void DiskIntensityMart::setChannelMapping(const IdxGroup &idxGroup) {
  m_CelChannels = idxGroup;
  std::vector<std::vector<unsigned int> > TempChannelGroup = m_CelChannels.getGroupingVec("channels");
  setChipChannelMap(TempChannelGroup);
}


/**
 * @brief Method to create empty diskMart with same parameters,
 * except with a new filename.
 *
 * @param tempDir - optional parameter to set new diskMart directory
 * @param tempPrefix - optional parameter to set new diskMart file prefix
 */
// DiskIntensityMart* DiskIntensityMart::copyMetaDataToEmptyMart(const std::string tempPrefix, const std::string tempDir) const {
//   std::string new_tempDir(tempDir);
//   if (new_tempDir.empty()) {
//     new_tempDir = m_TempDir;
//   }
  
//   std::string new_tempPrefix(tempPrefix);
//   if (new_tempPrefix.empty()) {
//     new_tempPrefix = m_TempPrefix;
//   }
DiskIntensityMart* DiskIntensityMart::copyMetaDataToEmptyMart() const {
  DiskIntensityMart* newMart = new DiskIntensityMart(m_AnalysisOrder, m_CelFiles, m_CacheProbeSize, m_TempDir, m_TempPrefix, m_StoreAllCelIntensities);
  
  newMart->m_Order = m_Order;
  newMart->m_Map = m_Map;
  newMart->m_ChipChannelToCacheMap = m_ChipChannelToCacheMap;
  newMart->m_useAuxMemCache = m_useAuxMemCache;
  return newMart;
}


std::string DiskIntensityMart::getFile5Name() const {
  if (m_File5Name=="") {
    m_File5Name = Fs::join(m_TempDir,m_TempPrefix+".a5");
  }
  return m_File5Name;
}

void DiskIntensityMart::openCacheFileToRead() const {
  m_TmpVectors.resize(m_NumChips);
  for(int i = 0; i < m_NumChips; i++) {
    if (m_TmpVectors[i] != NULL) {
      m_TmpVectors[i] = m_File5->openVector(ToStr(i), affx::FILE5_DTYPE_FLOAT, affx::FILE5_RW);
    }
  }
}

void DiskIntensityMart::closeTmpfile() const {
  if (m_File5!=NULL) {
    m_File5->close();
    delete m_File5;
    m_File5=NULL;
  }
}

void DiskIntensityMart::deleteTmpfile() const {
  // should be closed before deleting it.
  closeTmpfile();
  //
  if(m_DeleteFilesWhenFinished) {
    if (m_File5Name!="") {
      if(Fs::rm(m_File5Name, false) != APT_OK)  {
        // Can't throw an exception from destructor...
        Verbose::warn(0, "DiskIntensityMart::deleteTmpfile()"
                      " - Error can't delete file: " + m_File5Name);
      }
    }
  }
}


float DiskIntensityMart::getProbeIntensity(probeid_t pIx, chipid_t chipIx, unsigned int channelIx) const {
//   if(!m_Flushed) {
//     openCacheFileToRead();
//     m_Flushed = true;
//   }

  assert(pIx < m_Map.size() && pIx >= 0);
  assert(m_Map[pIx] >= 0);

  int cacheIx;
  if (m_ChipChannelToCacheMap.empty()) {
    cacheIx = chipIx;
  }
  else {
    cacheIx = m_ChipChannelToCacheMap[chipIx][channelIx];
  }

  if (!(cacheIx < m_NumChips && cacheIx >= 0)) {
    Verbose::out(1, ToStr("cacheIx: ") + ToStr(cacheIx));
    Verbose::out(1, ToStr("m_NumChips: ") + ToStr(m_NumChips));
  }
  assert(cacheIx < m_NumChips && cacheIx >= 0 && m_TmpVectors[cacheIx] != NULL);
  float rtnVal = -1.0;
  if (m_AuxMemCache.size() > cacheIx && 
      m_AuxMemCache[cacheIx].find(pIx) != m_AuxMemCache[cacheIx].end()) {
    std::map<int,float>::const_iterator result = m_AuxMemCache[cacheIx].find(pIx);
    rtnVal = result->second;
  }
  else {
    if (!inCache(pIx,cacheIx)) {
      loadIntoCache(pIx, cacheIx);
    }
    rtnVal = fromCache(pIx, cacheIx);
  }
  return rtnVal;
}


void DiskIntensityMart::setProbeIntensity(const int dataIdx, const std::vector<float> &data) {

  assert(data.size() > 0);

  if (m_Map.size() == 0) {
    // initialize map and prune AnalysisOrder of duplicate ids
    /// @todo - Could be more probes in the layout than the cel file if probes used in more than one probeset. Can we optimize that somehow?
    probeidx_t probeCount = data.size();

    m_Order.reserve(probeCount);
    /// Make the mapping from cel file position to position in the optimized format
    m_Map.resize(probeCount);

    fill(m_Map.begin(), m_Map.end(), -1);
    fill(m_Order.begin(), m_Order.end(), -1);

    int indexCount = 0;
    int i;
    for (i = 0; i < m_AnalysisOrder.size(); i++) {
      int probeIndex = m_AnalysisOrder[i]; 
      // APT-512, 2010/01/10
      // probeIndex comes back as -1 in Windows and has to be checked. 
      // The problem differs between when Cgwin and DOS. 
      // m_AnalysisOrder is returning -1 In DOS on the very last probe set
      // but not for Cygwin. 
      if (probeIndex >= 0) {
        if (m_Map[probeIndex] == -1) {
          m_Order.push_back(probeIndex);
          m_Map[probeIndex] = indexCount++;
        }
        else if (m_useAuxMemCache) {  
          // this probeIndex has already been seen.  store its id for so that
          // later its intensity can be stored in memory
          m_AuxMemCache.resize(1);
          m_AuxMemCache[0].insert(std::pair<int,float>(probeIndex,-1.0));
        }
      }
    }
    // store the size of the analysis order after any duplicate probes
    // ids have been removed.
    m_UniqueAnalysisOrderSize = m_Order.size();
    if (m_StoreAllCelIntensities && probeCount > m_Order.size()) {
      // probe intensities that are not in m_AnalysisOrder are tacked
      // onto the end of m_Order.  in this way, m_AnalysisOrder probe
      // intensities will be cached efficiently, and yet
      // non-m_AnalysisOrder probe intensities will still be available
      // in the diskmart.
      for (i = 0; i < probeCount; i++) {
	if (m_Map[i] == -1) {
	  m_Order.push_back(i);
	  m_Map[i] = indexCount++;
	}
      }
    }
    m_Size = m_Order.size();
  }

  if (m_File5 == NULL) {
    m_File5 = new File5_File();
    m_File5->open(getFile5Name(), affx::FILE5_REPLACE);
  }

  if (dataIdx >= m_NumChips) {
    m_NumChips = dataIdx + 1;
    m_TmpVectors.resize(m_NumChips, NULL);
    m_TmpVectors[dataIdx] = m_File5->openVector(ToStr(dataIdx), affx::FILE5_DTYPE_FLOAT, affx::FILE5_RW);
    // make copies of aux memory cache for each dataset
    if (m_AuxMemCache.size() > 0) {
      m_AuxMemCache.resize(m_NumChips);
      m_AuxMemCache[dataIdx] = m_AuxMemCache[0];
    }
  }
  affx::File5_Vector *tmpVec = m_TmpVectors[dataIdx];
  std::map<int,float>* auxMemCache = NULL;
  if (dataIdx < m_AuxMemCache.size()) {
    auxMemCache = &(m_AuxMemCache[dataIdx]);
  }
  writeReorderedData(m_Order, data, tmpVec, auxMemCache);
}


void DiskIntensityMart::writeReorderedData(const std::vector<int> &order, const std::vector<float>& data, File5_Vector *f5, std::map<int,float>* auxMemCache) {
  std::vector<float> reorderedData;

  int write_data_size = data.size();
  if (!m_StoreAllCelIntensities &&
      write_data_size > m_UniqueAnalysisOrderSize) {
    write_data_size = m_UniqueAnalysisOrderSize;
  }

  // store intensities for duplicated probes in special cache
  for(int i = 0; i < write_data_size; i++) {
    if(order[i] >=0) {
      reorderedData.push_back(data[order[i]]);
      if (auxMemCache != NULL) {
        std::map<int,float>::iterator result = auxMemCache->find(order[i]);
        if (result != auxMemCache->end()) {
          result->second = data[order[i]]; 
        }
      }
    }
  }

  reorderedData.reserve(order.size());
  APT_ERR_ASSERT(reorderedData.size()==write_data_size,"internal error.");

  f5->resize(reorderedData.size());
  int written = f5->write_array(0, reorderedData.size(), &reorderedData[0]);
  assert(written == reorderedData.size() );
}


void DiskIntensityMart::loadIntoCache(probeid_t pIx, chipid_t chipIx) const {
  // cache debug statements
//   cout << "probeId: " << pIx+1 << '\t';
//   cout << "index: " << m_Map[pIx] << '\t';
//   cout << "cache start: " << m_CacheStart << '\t';
//   cout << "cache end: " << m_CacheEnd << '\n';
  int nRows = m_CacheProbeSize / m_TmpVectors.size();
  if (m_CacheProbeSize > static_cast<uint64_t>(m_Size * m_NumChips)) {
    nRows = m_Size;
  }
  int startIx = m_Map[pIx];
  if(startIx < 0) {
    Err::errAbort("DiskIntensityMart::loadIntoCache() - Can't have probe with no position info, did you use the right ChipLayout?");
  }
  int endIx = startIx + nRows;
  // Make sure the cache is valid
  if(m_Cache.size() != m_TmpVectors.size()) {
    m_Cache.resize(m_TmpVectors.size());
  }
  // Read the values from each data file into the cache.
  for(int cIx = 0; cIx < m_TmpVectors.size(); cIx++) {
    if (m_TmpVectors[cIx] != NULL) {
      m_Cache[cIx].resize(nRows);
      // Blank out any existing values
      fill(m_Cache[cIx].begin(), m_Cache[cIx].end(), -1.0f);
      File5_Vector *f5 = m_TmpVectors[cIx];
      int size = min((int)m_Cache[cIx].size(), (int)f5->size() - startIx);
      f5->read_array(startIx, size, &m_Cache[cIx][0]);
    }
  }
  m_CacheStart = startIx;
  m_CacheEnd = endIx;
}


probeidx_t DiskIntensityMart::getProbeCountFromCel(const std::string &celFile) {
  std::string tmp_unc_name=Fs::convertToUncPath(celFile);
  FusionCELData cel;
  cel.SetFileName(tmp_unc_name.c_str());
  probeidx_t result = -1;
  try {
    if(!cel.Read()) 
      Err::errAbort("\nDiskIntensityMart::getProbeCountFromCel() - Can't read cel file: " + cel.GetFileName() + 
                    "\n>>> Error reported: " + StringUtils::ConvertWCSToMBS(cel.GetError()));
    result = cel.GetRows() * cel.GetCols();
    cel.Close();
  }
  catch (...) {
    Err::errAbort("\nDiskIntensityMart::getProbeCountFromCel() - Can't read cel file: " + cel.GetFileName() + "\n");
  }
  return result;
}

std::vector<float> DiskIntensityMart::getCelData(int dataSetIx) {
  assert(dataSetIx < m_NumChips && dataSetIx >= 0);
  File5_Vector *f5 = m_TmpVectors[dataSetIx];
  std::vector<float> data(f5->size());
  f5->read_vector(0, &data);
  std::vector<float> origOrderedData(m_Map.size());
  for (int i = 0; i < m_Map.size(); i++) {
    origOrderedData[i] = data[m_Map[i]];
  }
  return origOrderedData;
}

std::vector<float> DiskIntensityMart::getCelData(chipid_t celIx, unsigned int channelIx) {
  int dataSetIx = 0;
  // convert celname index to cache index, if a mapping exists
  if (m_ChipChannelToCacheMap.empty()) {
    if (channelIx == 0) {
      dataSetIx = celIx;
    }
    else {
      Err::errAbort("\nDiskIntensityMart::getCelData - Accessing CEL data with multi-channel index when no multi-channel data is available.");
    }
  }
  else {
    assert(celIx >= 0 && celIx < m_ChipChannelToCacheMap.size());
    assert(channelIx >= 0 && channelIx < m_ChipChannelToCacheMap[celIx].size());
    dataSetIx = m_ChipChannelToCacheMap[celIx][channelIx];
  }
  std::vector<float> data = getCelData(dataSetIx);
  return data;
}

