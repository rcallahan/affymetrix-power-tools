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
 * @file   DataStore
 * @author Chuck Sugnet
 * @date   Fri Oct  2 11:29:47 PDT 2009
 * 
 * @brief Class for encapsulating microarray data and storing it on disk in a
 * reasonable way to reduce memory impact. Key idea is to reorganize the data on
 * disk in the order which it will be used by the application to avoid lots of
 * seeks. Hopefully with a small number of seeks the performance impact of
 * on-disk stoarage will be minimal. Similar to DiskIntensityMart in spirit but
 * with a couple other features like reorganizing. Ability to serialize out all
 * data including the map and probe ids is also important.
 *
 * Open questions: 
 * - Can be generalized to support formats like a summarized probesets or other
 *   data that is now in chp files?
 * - Support generic sorting for things like chromosome order, or probeset order
 *   or chip order for blemishs
 */

#ifndef _DATASTORE_H_
#define _DATASTORE_H_

class DataStore;

//
#include "chipstream/AptTypes.h"
#include "chipstream/CelListener.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/PsBoard.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/DiskIntensityMart.h"
#include "file/CELFileData.h"
#include "file/FileWriter.h"
#include "file5/File5.h"
#include "portability/affy-base-types.h"

#include "util/Options.h"
#include "util/TmpFileFactory.h"
#include "util/Util.h"
#include "util/Verbose.h"

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

#define CEL_LABEL "cels"
#define PROBE_LABEL "probes"
#define PROBESET_LABEL "probe-sets"

using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;

/**
 * Class that stores data about probes and probesets out of core but loads into
 * RAM as necessary for computation. By ordering data on disk in the same order
 * as processing occurs a large decrease in RAM usage can be obtained with minimal
 * impact in performance times.
 *
 * Currently this class is a 'write once' type affair. It is not recommended to try
 * resorting or reuse an existing object.
 */
class DataStore : public IntensityMart, public CelListener {

private :

  /**
   * Types of probe data that the store can hold
   */
  enum ProbeColumns {
    ProbeGc,
    ProbeGcBgrd,
    ProbePm,
    ProbeMm,
    ProbePmMm,
    ProbeInProbeset,
    ProbePmAlleleMatch,
    // Insert new probe annotations here
    ProbeLastEnum,
  };
                                               
public :

  /** 
   * Constructor.
   */
  DataStore();

  /** 
   * Constructor allowing specification of temporary file to be used
   * 
   * @param fileLocation - Path to file location to be used.
   */
  DataStore(const std::string &fileLocation);

  /** 
   * Destructor
   */
  ~DataStore();

  /** 
   * Initialize class from another DataStore file. Reading in appropriate orders, number of channels, etc.
   * 
   * @param fileLocation - Path to DataStore file to read from.
   */
  void initFromFile(const std::string &fileLocation);

  /** 
   * Initialize values directly from another DataStore object
   * 
   * @param ds - DataStore to read order, celnames, num channels, etc from (but not intensities)
   */
  void initFromDataStore(const DataStore &ds);

  /** 
   * Get the order that the probes are currently stored in.
   * @return - vector with probe ids in the order currently sorted on disk
   */  
  const std::vector<int> &getProbeOrder() const;

  /** 
   * Set the probe order that things will be sorted on disk
   * 
   * @param order - order that things should be stored on disk.
   */
  void setProbeOrder(const std::vector<int> &order);

  /** 
   * Cel file of data.
   * 
   * @param cel - Handle to an open cel file.
   */
  void newChip(affymetrix_fusion_io::FusionCELData *cel);
  
  /** 
   * Utility function to read data from cel files outside of the
   * CelReader/IntensityMart framework.
   * 
   * @param celFiles - Vector of paths to cel files to be read.
   */
  void slurpDataFromCelFiles(const std::vector<std::string> &celFiles);

  /** 
   * Given a chip and a channel calculate the the correct column of our data matrix.
   * Data is stored interleaved so chip 0 channel 0 is column 0 and chip 0 channel 1 is 
   * column 1, etc. 
   * 
   * @param chipIx - What cel file of data is being queried
   * @param channelIx - What channel of that cel file is desired
   * 
   * @return - correcto column  of our data matrix to get the data.
   */
  int getDataSetIx(int chipIx, int channelIx) const;

  /** 
   * Given an data set index return the cel file name that it comes from. For
   * example if the cel files have 2 channels and you ask for the cel name for
   * dataSetIx = 1 or dataSetIx = 2 you'll get the same name.
   * 
   * @param dataSetIx - What cel file of data is being queried
   * 
   * @return - string name of cel file containing that data set index
   */
  std::string getCelName(int dataSetIx) const;
  
  /** 
   * @brief Given the probe index and chip index return the intensity
   * data appropriate for that probe in that chip. Depending on the
   * object that intensity may have been modified from the original
   * found on the chip.
   * 
   * @param probeIx - Probe Index number.
   * @param chipIx - Chip Index number.
   * 
   * @return double - intensity for that position on array.
   */
  virtual float getProbeIntensity(probeid_t probeIx, chipid_t chipIx, unsigned int channelIx = 0) const;

  /** 
   * @brief Given the probe index and chip index return the intensity
   * data appropriate for that probe in that chip. Depending on the
   * object that intensity may have been modified from the original
   * found on the chip.
   * 
   * @param probeIx - Probe Index number.
   * @param chipIx - Chip Index number.
   * 
   * @return double - intensity for that position on array.
   */
  virtual void setProbeIntensity(float value, probeid_t probeIx, chipid_t chipIx);

  /** 
   * @brief Return all of the intensities for given chipIx
   * @param chipIx - Chip Index number.
   * @return data - vector of all intensities in the DataStore for given chipIx.
   */
  virtual std::vector<float> getCelData(int dataSetIx);

  /** 
   * @brief Method for getting vector of intensities in original order in CEL file.
   * @param celIx - index of input CEL file
   * @param channelIx - index of CEL channel that desired dataSet is on
   * @return vector of intensities in CEL file order
   */
  virtual std::vector<float> getCelData(chipid_t celIx,
                                        unsigned int channelIx);

  /** 
   * @brief Return all of the intensities for given chipIx
   * @param chipIx - Chip Index number.
   * @return double - vector of all intensities in the DataStore for given chipIx.
   */
  void fillInCelData(chipid_t chipIx, std::vector<float> &data) const;

  /** 
   * @brief Given a vector of data use it to fill in all of the datapoints
   * that are going to be needed.
   * 
   * @param dataName - Name of the vector of data (usually the cel filename).
   * @param data - cel file intensity data. 
   */
  virtual void setProbeIntensity(const int chipIx, const std::vector<float> &data) {
    writeColumn(chipIx, m_CelNames[chipIx / m_NumChannels], data);
  };

  /** 
   * @brief Get the names (cel files) for the various data that has been seen.
   * @return Reference to all of the filenames.
   */
  virtual const std::vector<std::string> &getCelFileNames() const { return m_CelNames; };

  /** 
   * @brief Get the names (cel files) for the various data that has been seen.
   * @return Reference to all of the filenames.
   */
  virtual void setCelFileNames(const std::vector<std::string> &names);

  /** 
   * @brief Get the names (cel files) for the various data that has been seen.
   * @return Reference to all of the filenames.
   */
  virtual const std::vector<std::string> &getCelFileNames() { return m_CelNames; };

  /** 
   * @brief Get the total number of probes that this mart can supply.
   * @return int - total number of probes.
   */
  virtual size_t getProbeCount() const { return m_OrderProbe.size(); }

  /** 
   * @brief Get the total number of CEL channels that this mart can supply
   * @return int - total number of chips.
   */
  virtual int getCelDataSetCount() const {  
    return m_Chips->getColumnCount(0);
  }

  virtual int getChannelCount() const { return m_NumChannels; }

  /** 
   * @brief Get the total number of CELs/multi-CELs that this mart can supply
   * @return int - total number of chips.
   */
  virtual int getCelFileCount() const { 
    return int(m_CelNames.size());
  }

  /** 
   * Wrapper to check for common setup errors.
   * 
   * @param msg - Where checking is happening.
   */
  void checkSetup(const std::string &msg);

  /** 
   * Load up values from a ChipLayout object and place them on disk in store.
   * 
   * @param layout - ChipLayout with probe data of interest.
   * @param board - Board with other options (eg gc content file) to load data from
   */
  void setValues(ChipLayout &layout, PsBoard &board);


  /**
   * @brief Method to create empty DataStore with same parameters.
   * @todo either replace this method with IntensityMart copy
   * constructor, or flesh it out in a way appropriate to DataStore
   *  usage
   */
  virtual DataStore* copyMetaDataToEmptyMart() const;

  /** 
   * Open up the file5 tsv with the correct number of channels
   * 
   * @param numChannels - Number of channels for this batch of cel files
   */
  void initChannels(int numChannels);

  /** 
   * Initialize the intensity data.
   * 
   * @param numProbes - Number of probes in a particular cel file
   * @param numChannels - Number of channels in the cel files
   */
  void initIntensity(int numProbes, int numChannels);

  /** 
   * Write a column of intensity data to our intensity store.
   * 
   * @param chipIx - Chip number (column index in the file5 tsv)
   * @param colName - Name of the column 
   * @param data - Vector of data to write in cel file order
   */
  void writeColumn(int chipIx, const std::string &colName, const std::vector<float> &data);

  // Accesors for various attributes
  void setGcControlProbes(const std::vector<Probe *> &vec);
  void getGcControlProbes(std::vector<int> &vec);

  bool isSetProbeGc();
  void setProbeGc(const std::vector<char> &vec);
  void getProbeGc(std::vector<char> &vec);
  char getProbeGc(int probeIx);
  void getProbeGcVec(std::vector<int> &vec);

  bool isSetProbeGcBgrd();
  void setProbeGcBgrd(const std::vector<char> &vec);
  void getProbeGcBgrd(std::vector<char> &vec);
  char getProbeGcBgrd(int probeIx);

  bool isSetProbePm();
  void setProbePm(const std::vector<char> &vec);
  void getProbePm(std::vector<char> &vec);
  char getProbePm(int probeIx);
  void getProbePm(std::vector<bool> &pmBool);

  bool isSetProbeMm();
  void setProbeMm(const std::vector<char> &vec);
  void getProbeMm(std::vector<char> &vec);
  char getProbeMm(int probeIx);

  bool isSetProbePmMm();
  void setProbePmMm(const std::vector<int> &vec);
  void getProbePmMm(std::vector<int> &vec);
  int getProbePmMm(int probeIx);

  bool isSetProbeInProbeset();
  void setProbeInProbeset(const std::vector<char> &vec);
  void getProbeInProbeset(std::vector<char> &vec);
  char getProbeInProbeset(int probeIx);

  bool isSetProbePmAlleleMatch();
  void setProbePmAlleleMatch(const std::vector<int> &vec);
  void getProbePmAlleleMatch(std::vector<int> &vec);
  int getProbePmAlleleMatch(int probeIx);

  void setBufferSize(int bytes) const;

protected:

  /** 
   * Initialize class to appropriate values
   */
  void initVariables();

  /** 
   * Read the data from one cel file into the datastore
   * @param celName - Path to file to read
   */
  void readCelFile(const std::string &celName);

  /** 
   * Get the file5 group with the intensity tsv. Opening if necessary
   * @return - Pointer to the appropriate group.
   */
  affx::File5_Group *getIntensityGroup();

  /** 
   * Get the file5 group with the probe tsv. Opening if necessary
   * @return - Pointer to the appropriate group.
   */
  affx::File5_Group *getGroupProbe();

  /** 
   * Get the file5 group with the probeset tsv. Opening if necessary
   * @return - Pointer to the appropriate group.
   */
  affx::File5_Group *getGroupProbeSet();

  /** 
   * Get the file5 tsv  with the probe data. Opening if necessary
   * @return - Pointer to the appropriate tsv.
   */
  affx::File5_Tsv *getTsvProbe();

  /** 
   * Get the file5 tsv  with the probeset data. Opening if necessary
   * @return - Pointer to the appropriate tsv.
   */
  affx::File5_Tsv *getTsvProbeSet();

  /**
   * Cleanup
   */ 
  void closeResources();

  /** 
   * Initialize number of channels and probes from a cel file.
   * 
   * @param celName - Path to location of cel file.
   */
  void initFromCelFile(const std::string &celName);

  /** 
   * Initialize our number of channels and probes from a cel file.
   * 
   * @param - cel file object to get channel and probe number from
   */
  void initFromCelFile(affymetrix_fusion_io::FusionCELData &cel);

protected:
  int m_NumChannels;
  int m_NumChips;
  int m_NumProbe;
  int m_NumProbeSet;

  std::vector<int> m_ColIndexesProbe;
  std::vector<int> m_ColIndexesProbeSet;
  int m_ColCountProbe;
  int m_ColCountProbeSet;
  std::string m_FileStorage;
  affx::File5_File m_File5;	

  affx::File5_Group *m_IntensityGroup;
  mutable affx::File5_Tsv* m_Chips;

  affx::File5_Group *m_GroupProbe;
  affx::File5_Tsv *m_TsvProbe;

  affx::File5_Group *m_GroupProbeSet;
  affx::File5_Tsv *m_TsvProbeSet;
  std::vector<std::string> m_CelNames;
  std::vector<int> m_OrderProbe;
  affx::File5_Vector *m_OrderProbeStorage;
  std::vector<int> m_OrderProbeSet;
  affx::File5_Vector *m_ProbeSetOrderStorage;

  std::vector<probeidx_t> m_MapProbe;
  std::vector<probeidx_t> m_MapProbeSet;
};

#endif /* _DATASTORE_H_ */

