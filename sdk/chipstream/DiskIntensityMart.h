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
 * @file   DiskIntensityMart.h
 * @author Chuck Sugnet
 * @date   Mon Apr  7 06:56:22 2008
 * 
 * @brief Class for encapsulating microarray data and storing it on disk in a
 * reasonable way to reduce memory impact. Key idea is to reorganize the data on
 * disk in the order which it will be used by the application to avoid lots of
 * seeks. Hopefully with a small number of seeks the performance impact of
 * on-disk stoarage will be minimal.
 */

#ifndef _DISKINTENSITYMART_H_
#define _DISKINTENSITYMART_H_

//
#include "chipstream/AptTypes.h"
#include "chipstream/CelListener.h"
#include "chipstream/IntensityMart.h"
//
#include "file/CELFileData.h"
#include "file/FileWriter.h"
#include "file5/File5.h"
#include "portability/affy-base-types.h"

#include "util/Util.h"
//
#include <cstring>
#include <string>
#include <vector>
//
#ifdef WIN32
#include <windows.h>
#else
#include  <unistd.h>
#endif /* WIN32 */

/**
 * Very simple version of an IntensityMart based directly on cel files
 * for comparison and troubleshooting.
 * TODO: 
 *  - Constructor with cache files supplied
 *  - Meta information to check if file matches current file.
 *  - Multiple caches in case two parts of chip are being accessed at the same time.
 *  - Ability to specify prefix directory for temp files.
 *  - Delete Temp files when finished.
 *  - Measure performance against Harley's HDF5 layer.
 *  - Ability to store normalized, modified, results to avoid chipstream processing.
 */
class DiskIntensityMart : public IntensityMart {

public :

  /** 
   * Constructor that will build the data files with CEL data reorganized according to the order
   * seen in the ChipLayout object. 
   * 
   * @param layout - Specififes the probesets (ProbeList) and the order that they will be processed in.
   * @param celNames - Cel files that will be processed
   * @param cacheSize - Size of internal cache to be used.
   * @param tempDir - directory to write a5 file
   * @param tempPrefix - beginning of filename of a5 file.  tempPrefix will be appended with an object-specific code and an .a5 extension
   * @param storeAllCelIntensities - if false (default) then only intensities in layoutOrder will be stored.  if true then layoutOrder intensities will be stored, followed by remaining cel intensities.
   */
  
  DiskIntensityMart(const std::vector<int>& layoutOrder,
                    const std::vector<std::string>& celNames,
                    int cacheSize, 
                    const std::string& tempDir,
                    const std::string& tempPrefix="tempcel.",
		    bool storeAllCelIntensities = false);


  /**
   * @brief Basic destructor
   */
  ~DiskIntensityMart();

  /**
   * @brief Copy constructor
   *
   * @param diskMart - DiskIntensityMart to be copied
   */
  DiskIntensityMart(DiskIntensityMart& diskMart);


  /**
   * @brief Method to create empty diskMart with same parameters,
   * except with a new filename.
   *
   * @param tempDir - optional parameter to set new diskMart directory
   * @param tempPrefix - optional parameter to set new diskMart file prefix
   */
/*   DiskIntensityMart* copyMetaDataToEmptyMart(const std::string tempPrefix="", const std::string tempDir="") const; */
  DiskIntensityMart* copyMetaDataToEmptyMart() const;


/*   void newChip(const std::string &celFileName, const std::vector<float>& data); */
  /* void newChip(std::vector<float> data); */
  /** 
   * Read in and reorder the cel files to a more efficient on disk
   * representation based on the order in the ChipLayout probesets.
   * 
   * @param layout - ChipLayout that we are optimizing too.
   * @param celNames - Names of cel files that will be processed.
   */
/*   void reorderCelFilesToDataFiles(ChipLayout &layout, const std::vector<std::string> &celNames); */

  /** 
   * @brief Read in and reorder the cel files to a more efficient on disk
   * representation based on the order in the ChipLayout probesets.
   * 
   * @param layout - ChipLayout that we are optimizing too.
   * @param celNames - Names of cel files that will be processed.
   */
/*   void reorderCelFilesToDataFiles(const std::vector<int> &layoutOrder, const std::vector<std::string> &celNames); */
  
  void writeReorderedData(const std::vector<int> &order, const std::vector<float>& data, affx::File5_Vector *f5, std::map<int,float>* auxMemCache = NULL);

  /**
   * @brief Return number of CEL channels stored in diskmart
   */
  int getCelDataSetCount() const;

  /**
   * @brief Return number of multi-CEL files stored in diskmart
   */
  int getCelFileCount() const;

  float fromCache(probeid_t pIx, chipid_t chipIx) const;

  bool inCache(probeid_t pIx, chipid_t chipIx) const;

  void loadIntoCache(probeid_t pIx, chipid_t chipIx) const;

  // not really const.
  void openCacheFileToRead() const;

  /**
   * @brief set object that will translate from chip,channel to cache index
   *
   * @param chip_channel_map - First dimension is chipIx, second is
   * channelIx, value is cache index
   */
  void setChipChannelMap(std::vector<std::vector<unsigned int> >& chip_channel_map);

  int getChannelCount() const;

  float getProbeIntensity(probeid_t pIx, chipid_t chipIx, unsigned int channelIx = 0) const;
  
  void setProbeIntensity(const int dataIdx, const std::vector<float> &data);
/*     Err::errAbort("DiskIntensityMart::setProbeIntensity() - Not implemented."); */
/*     assert(false); */


  const std::vector<std::string> &getCelFileNames() const;

  /** 
   * @brief Method for getting vector of intensities in original order in CEL file.
   * @param dataSetIx - index of CEL intensity dataset in DiskIntensityMart 
   * @return vector of intensities in CEL file order
   */
  std::vector<float> getCelData(int dataSetIx);

  /** 
   * @brief Method for getting vector of intensities in original order in CEL file.
   * @param celIx - index of input CEL file
   * @param channelIx - index of CEL channel that desired dataSet is on
   * @return vector of intensities in CEL file order
   */
  std::vector<float> getCelData(chipid_t celIx, unsigned int channelIx);

  /** 
   * @brief Method for setting boolean to indicate if diskMart should store all given intesities or only the ones specified in m_Order.
   * @param flag - boolean indicating desired behavior
   */
  void setStoreAllCelIntensities(bool flag);

  /** 
   * Open a CEL file and see how many probes are stored in it.
   * @param celFile - Path to cel file of interest.
   * 
   * @return - Count of probes contained.
   */
  static probeidx_t getProbeCountFromCel(const std::string &celFile);

  size_t getProbeCount() const;

  virtual void setChannelMapping(const IdxGroup &idxGroup);

  void setUseAuxMemCache(bool use) {m_useAuxMemCache = use;}

private: 

  /// Close the tmpfile.
  void closeTmpfile() const;
  /// Close and delete.
  void deleteTmpfile() const;
  /// Get the name of the tmpfile to use.
  std::string getFile5Name() const;

  /// Names of the data files that are being used to read.
  std::vector<std::string> m_CelFiles;

  /// Translate from positions on original disk to those in our data files that
  /// have been reorganized for sequential access. 
  std::vector<probeidx_t> m_Map;
  std::vector<probeidx_t> m_Order;
  /// For what regions of m_Map is the cache valid
  mutable probeidx_t m_CacheStart, m_CacheEnd;

  /// How many probes are stored in the cache at any given time.
  int m_CacheProbeSize;

  /// Cache of data in memory
  mutable std::vector< std::vector<float> > m_Cache;

  /// Directory for temp files
  std::string m_TempDir;

  /// Prefix for temporary files.
  std::string m_TempPrefix;

  /// Delete temporary files when finished?
  bool m_DeleteFilesWhenFinished;

  mutable affx::File5_File *m_File5;
  mutable std::string m_File5Name;
  mutable std::vector<affx::File5_Vector *> m_TmpVectors;

  /// number of intensities stored per dataset
  int m_Size;
  /// number of CELs
  int m_NumChips;
  mutable bool m_Flushed;
  /// probes listed in order of appearance in ChipLayout object
  std::vector<int> m_AnalysisOrder;
  /// number of unique probe ids in m_AnalysisOrder
  int m_UniqueAnalysisOrderSize;
  /// number of CEL files anticipated 
  int m_ExpectedChips;
  /// auxilliary cache of duplicated probes
  std::vector<std::map<int,float> > m_AuxMemCache;
  /// flag to indicate if auxilliary cache should be used
  bool m_useAuxMemCache;

  static int m_cache_misses;
};

#endif /* _DISKINTENSITYMART_H_ */
