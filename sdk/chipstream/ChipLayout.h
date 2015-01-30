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

#ifndef _CHIPLAYOUT_H_
#define _CHIPLAYOUT_H_

//
#include "chipstream/CumulativeStats.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/ProbeSet.h"
//
#include "calvin_files/fusion/src/FusionCDFData.h"
#include "file/TsvFile/PgfFile.h"
#include "file/TsvFile/SpfFile.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>
//

/**
  ChipLayout - Data structure to represent the probesets encoded in the CDF, SPF or
  PGF files which define the probes that belong which probe sets and if they are
  PM, MM, genotyping etc.  Parse the file and store it in memory for easy access
  by index or probe set name. Note that there is currently no way to reuse a ChipLayout
  object, that is there is no close() function. Delete the ChipLayout and open another
  one if wish to open another cdf/pgf/spf file.

  For the CEL file SDK the order of the probes indexes are:
  \verbatim
  (0,0) (1,0) (2,0) (3,0) ... (numCols,0)
  (0,1) (1,1) (2,1) (3,1) ... (numCols,1)
  (0,2) (1,2) (2,2) (3,2) ... (numCols,2)
  ...
  (0,numRows)       ... (numCols,numRows)
  \endverbatim

  Where in cel file terminology cols are the x dimension and
  rows are the y dimension. So using the standard XYToIndex()
  function in the cel file API you can transform into a single
  vector index i using:
  - i = ((y*numColumns) + x)
  - x = i % numCols
  - y = i / numCols

  So the indexing he numbering when num cols = 10
  \verbatim
  0  1  2  3  4  5  6  7  8  9
  10 11 12 13 14 15 16 17 18 19
  20 21 22 23 24 25 26 27 28 29
  \endverbatim

  Note that these indexes are only for the CEL files as the underlying DAT files
  can be rotated depending on the scanner and application.
 */
class ChipLayout {

public:

  /** Constructor. */
  ChipLayout();

  /*** Destructor. */
  ~ChipLayout();

  /**
   * Are the probesets randomly accessible by name?
   * @return true if can lookup by name, false otherwise.
   */
  inline bool canAccessByName() { 
    //return m_HaveNameMap; 
    return true;
  }

  /**
   * Find a probe set using the name as a key. canAccessByName()
   * must be true for this to work.
   * @param psName - Name of probe set desired.
   * @return Pointer to ProbeList
   */
  ProbeListPacked getProbeListByName(const std::string& name);

  int getProbeSetIndexByName(const std::string& name);

  /**
   * Does the layout have a reference to the probeset with name psName? Must
   * have name map filled to call this function
   *
   * @param psName - identifier of probeset.
   * @return true if in list false otherwise
   */
  bool containsProbeSet(const std::string& psName) {
    // we can alwasy access by name.
    //if(!canAccessByName())
    //Err::errAbort("ChipLayout::containsProbeSet() - Can't access by name unless canAccessByName() == true");
    // use the probelist factory
    // return m_PsNameHash.find(psName.c_str()) != m_PsNameHash.end();
    int idx=m_PlFactory.getProbeListIndexByName(psName);
    return (idx!=-1);
  }

  /**
   * Get a probelist by the index in the layout.
   * @param index - index of probeset in chiplayout list.
   * @return Probelist requested.
   */
  ProbeListPacked getProbeListAtIndex(unsigned int index);

  /**
   * How many probe sets in Pgf file?
   * @return int - number of probe sets defined.
   */
  inline unsigned int getProbeSetCount() const { return m_PlFactory.getProbeSetCount(); }

  /**
   * Open and parse out the data of a file.
   * @param fileName - Name of PGF file to be opened.
   * @param xSize - How many features in x dimension.
   * @param ySize - How many features in y dimension.
   * @param probeSetsToLoad - Vector of probe set ids to be loaded.
   * @param probeSubset - Subset of probes to be loaded.
   * @param chipType - What sort of chip should this pgf file be for?
   * @param killList - Map of individual probes that should not be used.
   * @return bool - true if successful.
   */
  bool openPgf(const std::string& fileName,
            unsigned int xSize, unsigned int ySize,
            const std::set<const char *, Util::ltstr> &probeSetsToLoad,
            std::vector<const char *> *probesetNames,
			std::vector<const char *> *probesetDisplayNames,
            std::vector<bool> &probeSubset,
            const std::string &chipType,
            probeidmap_t &killList,
            bool justStats, bool doPList) ;

  /**
   * Open and parse out the data of a file.
   * @param fileName - Name of file to be opened.
   * @param xSize -
   * @param ySize -
   * @param numProbes - How many probes are on this chip?
   * @param probeSetsToLoad - Vector of probe set ids to be loaded.
   * @param probeSubset - Subset of probes to be loaded.
   * @param chipType - What sort of chip should this pgf file be for?
   * @return bool - true if successful.
   */
  bool openPgf(const std::string& fileName,
               unsigned int xSize,
               unsigned int ySize,
               const std::set<const char *, Util::ltstr> &probeSetsToLoad,
               std::vector<bool> &probeSubset,
               const std::string &chipType) {
    probeidmap_t killList;
    return openPgf(fileName, xSize, ySize, probeSetsToLoad, NULL, NULL, probeSubset,chipType, killList, false, false);
  }

  /**
   * Open a file and parse the probe sets.
   * @param fileName - Name of file to be opened.
   * @return bool - true if success. Errors out on problems.
   */
  bool openPgfAll(const std::string& fileName, unsigned int xSize, unsigned int ySize) {
    std::set<const char *, Util::ltstr> probeSetsToLoad;
    std::vector<bool> probeSubset;
    std::string chipType;
    probeidmap_t killList;
    std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad;
    return openPgf(fileName, xSize, ySize, probeSetsToLoad, probeSubset, chipType);
  }

  /**
   * Open a file and parse the probe sets.
   * @param fileName - Name of file to be opened.
   * @param probeSetsToLoad - Which probe sets should be loaded?
   * @param probeSubset - Subset of probes to be loaded.
   * @param chipType - What sort of chip should this cdf file be for?
   * @param probesetNames - If not null filled in with the
   * probesetname for each probeset even if not actually into
   * chiplayout. Useful for chp files where every single probeset must
   * be included and in order.
   * @param - killList map of probes that are to be excluded from
   * the probesets. A probe must have its id mapped to the probeset name
   * for this to work.
   * @param justStats - just read and generate stats, don't load into memory.
   * @param psTypesToLoad - What types of probesets to load into memory.
   * @return bool - true if success. Errors out on problems.
   */
  bool openCdf(const std::string& fileName,
               const std::set<const char *, Util::ltstr> &probeSetsToLoad,
               std::vector<const char *> *probesetNames,
               std::vector<bool> &probeSubset, const std::string &chipType,
               probeidmap_t &killList,
               bool justStats,
               std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad);
  /**
   * Open a file and parse the probe sets.
   * @param fileName - Name of file to be opened.
   * @param probeSetsToLoad - Which probe sets should be loaded?
   * @param probeSubset - Subset of probes to be loaded.
   * @param chipType - What sort of chip should this cdf file be for?
   * @param probesetNames - If not null filled in with the
   * probesetname for each probeset even if not actually into
   * chiplayout. Useful for chp files where every single probeset must
   * be included and in order.
   * @param - killList map of probes that are to be excluded from
   * the probesets. A probe must have its id mapped to the probeset name
   * for this to work.
   * @return bool - true if success. Errors out on problems.
   */
  bool openCdf(const std::string& fileName,
               const std::set<const char *, Util::ltstr> &probeSetsToLoad,
               std::vector<const char *> *probesetNames,
               std::vector<bool> &probeSubset, const std::string &chipType,
               probeidmap_t &killList, bool justStats=false) {
    std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad;
    return openCdf(fileName, probeSetsToLoad, probesetNames, probeSubset, chipType,
            killList, justStats, psTypesToLoad);
  }

  /**
   * Open a file and parse the probe sets.
   * @param fileName - Name of file to be opened.
   * @param probeSetsToLoad - Which probe sets should be loaded?
   * @param probeSubset - Subset of probes to be loaded.
   * @param chipType - What sort of chip should this cdf file be for?
   * @param probesetNames - If not null filled in with the
   * probesetname for each probeset even if not actually into
   * chiplayout. Useful for chp files where every single probeset must
   * be included and in order.
   * @return bool - true if success. Errors out on problems.
   */
  bool openCdf(const std::string& fileName, const std::set<const char *, Util::ltstr> &probeSetsToLoad,
               std::vector<bool> &probeSubset, const std::string &chipType) {
    probeidmap_t killList;
    std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad;
    return openCdf(fileName, probeSetsToLoad, NULL, probeSubset, chipType, killList, false, psTypesToLoad);
  }

  /**
   * Open a file and parse the probe sets.
   * @param fileName - Name of file to be opened.
   * @return bool - true if success. Errors out on problems.
   */
  bool openCdfAll(const std::string& fileName) {
    std::set<const char *, Util::ltstr> probeSetsToLoad;
    std::vector<bool> probeSubset;
    std::string chipType;
    probeidmap_t killList;
    std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad;
    return openCdf(fileName, probeSetsToLoad, NULL, probeSubset, chipType, killList, false, psTypesToLoad);
  }

  /**
   * Write out the data in class to file in PGF format.
   * @param fileName - Name of file to write data.
   */
  void write(const std::string& fileName);

  /**
   * Write out the contents of probe sets in PGF format.
   * @param fileName - File to write data to.
   * @param vec - Vector of probe sets to be stored.
   */
  static void writePgfFile(const std::string& fileName, const std::vector<ProbeSet *> &vec);

  /**
   * Write out the contents of probe sets in PGF format.
   * @param out - output stream to write data to.
   * @param vec - Vector of probe sets to be stored.
   */
  static void writePgfFile(affx::PgfFile& pgf, const std::vector<ProbeSet *> &vec);

  /**
   * @brief Get the perfect match probes of chip.
   * @return Bitmask representig pm probe on chip.
   */
  const std::vector<bool> &getPmProbeMask() { return m_PmProbes; }

  /**
   * @brief Get the perfect match probes of chip.
   * @return Bitmask representig pm probe on chip.
   */
  const std::vector<bool> &getMmProbeMask() { return m_MmProbes; }

  /**
   * @brief Set the perfect match probes of chip.
   * @param mask with PM probes marked as true, others false.
   */
  void setPmProbeMask(const std::vector<bool> &mask) { m_PmProbes = mask; }

  /**
   * @brief Get the value of a layout key=value header pair.
   * @param key String to lookup value of.
   *
   * @return value of key/value pair.
   */
  const std::vector<std::string> &getHeaderValue(const std::string &key) {
    static std::vector<std::string> tmp;
    std::map<std::string, std::vector<std::string> >::iterator iter;
    iter = m_Header.find(key);
    if(iter == m_Header.end()) {
      Verbose::out(1, "Warning: Don't recognize header key: " + key);
      return tmp;
    }
    return iter->second;
  }

  /**
   * @brief Get the number of x counts
   * @return - int x maximum
   */
  inline int getXCount() { return m_XCount; }

  /**
   * @brief Get the number of y counts
   * @return - int y maximum
   */
  inline int getYCount() { return m_YCount; }

  /**
   * @brief Set the x and y dimensions for the array.
   *
   * @param xCount - Maximum value of x + 1 (corresponds to columns)
   * @param yCount - Maximium value of y + 1 (corresponds to rows)
   */
  inline void setDimensions(unsigned int xCount, unsigned int yCount) {
    m_XCount = xCount;
    m_YCount = yCount;
    // also set the dimensions of the ProbeListFactory.
    m_PlFactory.setDimensions(m_XCount,m_YCount);
  }

  int numChannels() {
    return m_numChannels;
  }

  /**
   * Get the total number of probes on this chip.
   * @return - Number of probes on array.
   */
  inline int getProbeCount() { return m_YCount * m_XCount; }

  /**
   * Get the total number of probes actually used in the cdf or pgf
   * @return - Number of probes in a probeset.
   */
  int getProbeUsedCount();

  /**
   * @brief Given probe id supply the index.
   * @param probeId - Id of probe.
   * @return Index in chip.
   */
  inline unsigned int probeIdToIndex(unsigned int probeId) { return probeId -1; }

  /**
   * @brief Given probe index supply the id
   * @param probeIndex - Index of probe
   * @return Probe id.
   */
  inline unsigned int indexToProbeId(unsigned int probeIndex) { return probeIndex + 1; }

  /**
   * @brief Translate an index (id) to the x and y coordinates on a chip.
   *
   * @param probeId - Id of probe of interest.
   * @param numProbes - Number of total probes on the chip
   * @param x - X coordinate of interest.
   * @param y - Y coordinate of interest.
   */
  static void indexToXY(unsigned int probeIndex, int numCol, unsigned int &x, unsigned int &y) {
    x = probeIndex % (numCol);
    y = probeIndex / (numCol);
  }

  /**
   * @brief Translate an index (id) to the x and y coordinates on a chip.
   *
   * @param probeId - Id of probe of interest.
   * @param x - X coordinate of interest.
   * @param y - Y coordinate of interest.
   */
  inline void indexToXY(unsigned int probeIndex, unsigned int &x, unsigned int &y) {
    indexToXY(probeIndex, m_XCount, x, y);
    y = probeIndex / (m_XCount);
  }


  /**
   * @brief Translate an index (id) to the x and y coordinates on a chip. (just like the unsigned version.)
   *
   * @param probeId - Id of probe of interest.
   * @param x - X coordinate of interest.
   * @param y - Y coordinate of interest.
   */
  inline void indexToXY(unsigned int probeIndex, int &x, int &y) {
    x = probeIndex % (m_XCount);
    y = probeIndex / (m_XCount);
  }

  /**
   * Turn the x and y cordinates of a probe on the array into an index. Should
   * be exactly the same index as provided by cel files.
   *
   * @param x - x coordinate on array.
   * @param y - y coordinate on array.
   * @return - index (id) of feature in cel file and other arrays.
   */
  int xyToIndex(int x, int y);


  /// Typedef for utility.
  typedef std::vector<ProbeSet *>::const_iterator psVecIter;

  /**
   * Convert from numeric representation to string.
   * @param type - numeric representation.
   */
  const char *typeFromCdfType(unsigned short type);

  /**
   * Is this probe a perfect match probe?
   * @param cel - probe information from cdf.
   * @return true if probe is a perfect match, false otherwise.
   */
  static bool isPm(const affymetrix_fusion_io::FusionCDFProbeInformation &cel);

  /**
   * Is a particular probe in the kill list (i.e. not to be used?)
   * @param probeIx - Id of probe.
   * @param psName - Name of probeset that contains probe.
   * @param killList - Map of probes not to be used.
   * @return true if probe id is in kill list, false otherwise.
   */
  bool inKillList(probeid_t, const std::string& psName, probeidmap_t &killList);

  /**
   * Read the contents of rowFile into a vector of
   * probe sets.
   * @param pgf - Opened Pgf file to read from
   * @param numProbes - How many probes are on this chip?
   * @param psVec - vector to be filled in.
   * @param probeSetsToLoad - Vector of probe set ids to be loaded.
   * @param probeSubset - Subset of probes to be loaded.
   * @param killList - List of probes not to be used.
   */
  void readProbeListPgfFileKillList(affx::PgfFile &pgf,
                                    unsigned int numProbes,
                                    //std::vector<ProbeListPacked> &plVec,
                                    ProbeListFactory &plFactory,
                                    const std::set<const char *, Util::ltstr> &probeSetsToLoad,
                                    std::vector<const char *> *probesetNames,
                                    std::vector<const char *> *probesetDisplayNames,
                                    std::vector<bool> &probeSubset,
                                    probeidmap_t &killList,
                                    bool justStats);

  /**
   * Read a vector of ProbeSets from a cdf file representation.
   *
   * @param cdf - cdf file object with probe set information.
   * @param psVec - vector to put newly created probe sets in.
   * @param probeSetsToLoad - which probe sets should be converted.
   * @param probeSubSet - which individual probes were loaded.
   * @param probesetNames - If not null filled in with the
   * probesetname for each probeset even if not actually into
   * chiplayout. Useful for chp files where every single probeset must
   * be included and in order.
   * @param killList - map of probes ids that shouldn't be used.
   * @param justStats - just read and generate stats, don't load into memory.
   * @param psTypesToLoad - What types of probesets to load into memory.
   */
  void readProbeListCdfFileKillList(affymetrix_fusion_io::FusionCDFData &cdf,
                                    //std::vector<ProbeListPacked> &plVec,
                                    ProbeListFactory &plFactory,
                                    const std::set<const char *, Util::ltstr> &probeSetsToLoad,
                                    std::vector<const char *> *probesetNames,
                                    std::vector<bool> &probeSubSet,
                                    probeidmap_t &killList,
                                    bool justStats,
                                    std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad);

  /**
   * Read cdf probe information into a Probe object.
   *
   * @param p - probe object to be filled in.
   * @param parent - Atom that is parent of this probe.
   * @param numCol - column on chip.
   * @param numRow - row on chip.
   * @param direction - synthesis direction.
   * @param group - group information from cdf.
   * @param cel -probe information from cdf.
   */
  void readCdfProbe(Probe *p, Atom *parent, unsigned int numCol, unsigned int numRow,
                    unsigned int direction,
                    const affymetrix_fusion_io::FusionCDFProbeGroupInformation &group,
                    const affymetrix_fusion_io::FusionCDFProbeInformation &cel);

  /**
   * Check to see if a charater is a valid DNA base of 'a','t','g', or 'c'
   * @param c - base to check.
   * @return - return true if c is in 'atgcATGC', false otherwise.
   */
  bool validBase(int c);

  /**
   * Build map between pm and mm probes.
   */
  void buildPmMmMap();

  /**
   * @brief Function to return the mm id associated with a given pm id.
   * @param probeIx - Index of perfect match probe.
   *
   * @return - The index of the associated mismatch probe, if any.
   */
  const int mmId(const int probeIx) const;

  static ProbeSet::ProbeSetType getProbeSetType(ProbeSet *ps) {
    ProbeSet::ProbeSetType type = (ProbeSet::ProbeSetType) ps->psType;
    if ((type == ProbeSet::GenoType) && (ps->atoms.size() != 2) && (ps->atoms.size() != 4)) {
      return ProbeSet::Expression;
    }
    return type;
  }

  /**
   * Get the probeset type for a particular probeset based on its
   * index in layout file. Note that probesets marked as "GenoType"
   * but only having one block are changed here to
   * "Expression". Currently only supported for cdf files.
   *
   * @param psIndex - Index of probeset.
   * @return - Type of probeset (i.e. ProbeSet::Genotype)
   */
  ProbeSet::ProbeSetType getProbeSetType(int psIndex) {
    ProbeListPacked plp=m_PlFactory.getProbeListAtIndex(psIndex);
    if (plp.isNull()) {
      APT_ERR_ABORT("psIndex out of range.");
    }
    ProbeSet::ProbeSetType ps_type=(ProbeSet::ProbeSetType)plp.get_type();
    int ps_block_cnt=plp.block_cnt();
    // do the conversion as described above.
    // if ((ps_type==ProbeSet::GenoType)&&(plp.block_cnt()==1)) {
    //   ps_type=ProbeSet::Expression;
    // }
    // This is rays version.
    if ((ps_type == ProbeSet::GenoType) && (ps_block_cnt != 2) && (ps_block_cnt != 4)) {
      return ProbeSet::Expression;
    }
    //
    return ps_type;
  };

  /**
   * @brief Determine if the probeset at the given index is a genotype probeset that has more than one annotated allele
   * @param psIndex - Index of probeset.
   * @return - boolean
   */
  bool isGenotypable(int psIndex) {
    bool genotypable = false;
    ProbeListPacked plp=m_PlFactory.getProbeListAtIndex(psIndex);
    if (plp.isNull()) {
      APT_ERR_ABORT("psIndex out of range.");
    }
    ProbeSet::ProbeSetType ps_type=(ProbeSet::ProbeSetType)plp.get_type();
    int block_count = plp.block_cnt();
    if (block_count > 1) {
      int first_allele = plp.get_blockAllele(0);
      bool more_than_one_allele = false;
      for (int i = 1; i < block_count; i++) {
        if (plp.get_blockAllele(i) != first_allele) {
          more_than_one_allele = true;
          break;
        }
      }
      if (more_than_one_allele &&
          ((ps_type == ProbeSet::GenoType && 
            (block_count == 2 || block_count == 4)
            )
           ||
           ps_type == ProbeSet::Marker
           ||
           ps_type == ProbeSet::MultichannelMarker
           )
          ) {
        genotypable = true;
      }    
    }
    return genotypable;
  }

  void readCdfProbeSet(const std::string& name,
                       affymetrix_fusion_io::FusionCDFProbeSetInformation &set,
                       unsigned int psIx, int numCol, int numRow,
                       std::vector<bool> &probeSubSet,
                       bool &keepProbeSet,
                       probeidmap_t &killList,
                       ProbeSet &ps, int &psNameSize,
                       std::vector<Atom> **locAtoms,
                       std::vector<Probe> **locProbes,
                       int& errorCount);

  std::vector<int> getPmMmVec() {
    return m_PmMmVec;
  }

  std::vector<int> getPmAlleleMatchVec() {
    return m_PmAlleleMatchVec;
  }

  const std::vector<char> & getGcProbeVec() {
    return m_ProbeGcVec;
  }

  std::vector<int> getProbesetProbeCounts() {
    return m_ProbeCounts;
  }

  std::vector<bool> getProbesetProbes() {
    return m_ProbesetProbes;
  }

  double meanProbesPerProbeset() {
    return m_ProbesPerProbeset.getMean();
  }

  double meanMemPerProbeset() {
    return m_ProbesetMemSizes.getMean();
  }

  void setNeedMismatch(bool needMisMatchMap) {
    m_NeedMismatch = needMisMatchMap;
  }

  void setNeedPmAlleleMatch(bool needPmAlleleMatch) {
    m_NeedPmAlleleMatch = needPmAlleleMatch;
  }

  void setNeedGc(bool needGcMap) {
    m_NeedGc = needGcMap;
  }

  /// @brief     Set the name of the spf output file
  /// @param     fileName   Filename to write.
  void setSpfFileName(const std::string &fileName) {
    m_SpfFileName = fileName;
  }

  //
  void openSpfForWrite(affx::SpfFile& spfFile,const std::string& spfFilename,int spfVersion=2);
  //
  void writeSpfProbeList(const std::string& spfFilename,int spfVersion=2);
  //
  static void writeSpfProbeSetRecord(affx::SpfFile& tsv,const ProbeSet& ps);
  static void writeSpfProbeListRecord(affx::SpfFile& tsv,const ProbeListPacked& pL);

  /**
   * Read in probelist data from a simple probe format file.
   *
   * @param fileName - Name of file to be opened.
   * @param probeSetsToLoad - Which probe sets should be loaded?
   * @param probesetNames - If not null filled in with the
   * probesetname for each probeset even if not actually into
   * chiplayout. Useful for chp files where every single probeset must
   * be included and in order.
   * @param probeSubset - Subset of probes to be loaded.
   * @param chipTypeExpected - What sort of chip should this cdf file be for?
   * @param justStats - just read and generate stats, don't load into memory.
   * @param psTypesToLoad - What types of probesets to load into memory.
   */
  /// @todo openSpf should allow for a kill list
  void openSpf(const std::string& fileName,
               const std::set<const char *, Util::ltstr>& probeSetsToLoad,
               std::vector<const char *> *probesetNames,
               std::vector<bool> &probeSubSet,
               const std::string &chipType,
               bool justStats,
               std::set<affxcdf::GeneChipProbeSetType>& psTypesToLoad);

  void openSpf_v2(const std::string& fileName,
                  const std::set<const char *, Util::ltstr>& probeSetsToLoad,
                  std::vector<const char *> *probesetNames,
                  std::vector<bool> &probeSubSet,
                  const std::string &chipType,
                  bool justStats,
                  std::set<affxcdf::GeneChipProbeSetType>& psTypesToLoad);
  void openSpf_v3(const std::string& fileName,
                  const std::set<const char *, Util::ltstr>& probeSetsToLoad,
                  std::vector<const char *> *probesetNames,
                  std::vector<bool> &probeSubSet,
                  const std::string &chipType,
                  bool justStats,
                  std::set<affxcdf::GeneChipProbeSetType>& psTypesToLoad);

  void openSpf_v4(const std::string& fileName,
                  const std::set<const char *, Util::ltstr>& probeSetsToLoad,
                  std::vector<const char *> *probesetNames,
                  std::vector<bool> &probeSubSet,
                  const std::string &chipType,
                  bool justStats,
                  std::set<affxcdf::GeneChipProbeSetType>& psTypesToLoad);

  void openSpfAll(const std::string& fileName) {
    std::set<const char *, Util::ltstr> probeSetsToLoad;
    std::vector<bool> probeSubset;
    std::string chipType;
    std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad;

    openSpf(fileName, probeSetsToLoad, NULL, probeSubset, chipType, false, psTypesToLoad);
  }

  /**
   * Reads in probes to kill.
   *
   * @param fileName - file containing BadProbes
   * @param killList - things to kill, determined from file
   * @param rows - number of rows in cel file
   * @param cols - number of cols in cel file
   */
  static void fillInKillList(const std::string& fileName, probeidmap_t &killList, int rows, int cols);

  const std::vector<int> & getProbeLayoutOrder() {
    return m_ProbeLayoutOrder;
  }

  /**
   * Parse the results from one record into a ProbeList if desired,
   * othewise just keep track of statistics
   *
   * @param name - Probeset identifier.
   * @param type - Class of probeset (i.e. genotyping, expression, copynumber) as in ProbeSet::Type enumeration
   * @param numBlocks - How many blocks are there in probeset.
   * @param blockSizes - Comma separated size in probes of each block.
   * @param blockAnns - Comma separated annotation for each block as in ProbeSet::BlockAnnotation enumeration
   * @param numMatch - How many matches (1 for pm only, 2 for pm-mm).
   * @param numProbes - Total number of probes.
   * @param probesStr - Comma separated list of probe ids (indexes).
   * @param lineNumber - Line number in file.
   * @param keep - Should this probelist be kept after getting statistics.
   * @param probeSubSet - Bool array with probes to be loaded marked as true.
   * @param keepOrder - Should the order of these probes be saved for layout order
   */
/*   void readProbeListRecord(std::string &name, int type, int numBlocks, std::string &blockSizesStr, */
/*                            std::string &blockAnnsStr, int numMatch, int numProbes, std::string &probesStr, */
/*                            int lineNumber, bool keep, std::vector<bool> &probeSubSet, bool keepOrder); */
  void parseProbeListRecord(const std::string& name,
                            int type,
                            int numBlocks,
                            const std::string& blockSizesStr,
                            const std::string& blockAnnsStr,
                            bool have_blockAlleleStr,
                            const std::string& blockAlleleStr,
                            bool have_blockContextStr,
                            const std::string& blockContextStr,
                            bool have_blockChannelStr,
                            const std::string& blockChannelStr,
                            bool have_blockRepTypeStr,
                            const std::string& blockRepTypeStr,
                            int numMatch,
                            int numProbes,
                            const std::string& probesStr,
                            int lineNumber,
                            bool keep,
                            std::vector<bool> &probeSubSet, bool keepOrder);

  /**
   * Create a probelist allocating from the factory with the attributes
   * passed in.
   *
   * @param name - Probeset identifier.
   * @param type - Class of probeset (i.e. genotyping, expression, copynumber) as in ProbeSet::Type enumeration
   * @param numBlocks - How many blocks are there in probeset.
   * @param blockSizes - Size in probes of each block.
   * @param blockAnns - Annotation for each block as in ProbeSet::BlockAnnotation enumeration
   * @param numMatch - How many matches (1 for pm only, 2 for pm-mm).
   * @param numProbes - Total number of probes.
   * @param probes - Vector of probe ids (indexes).
   *
   * @return Created ProbeList.
   */
  ProbeListPacked makeProbeListPacked(const std::string &name,
                                      int type, int numBlocks,
                                      const std::vector<int> &blockSizes,
                                      const std::vector<int> &blockAnns,
                                      const std::vector<int> &blockAlleles,
                                      const std::vector<int> &blockContexts,
                                      const std::vector<int> &blockChannels,
                                      const std::vector<int> &blockRepTypes,
                                      int numMatch, int numProbes,
                                      const std::vector<int> &probes);

  /**
   * Make maps of the probe sets in m_PsVec by id and name if possible.
   */
  void makePsMaps();

  void resizeStats(uint32_t numProbes) {
    m_PmProbes.resize(numProbes, false);
    if(m_NeedMismatch) {
      m_MmProbes.resize(numProbes, false);
      m_PmMmVec.resize(numProbes, -1);
    }
    if(m_NeedGc) {
      m_ProbeGcVec.resize(numProbes, (char)NULLPROBEGC);
    }
    if(m_NeedPmAlleleMatch) {
      m_PmAlleleMatchVec.resize(numProbes, -1);
    }
    m_ProbesetProbes.resize(numProbes, false);
  }

  void resizeAtoms(std::vector<Atom> **old, int32_t newSize, ProbeSet *ps);

  void resizeProbes(std::vector<Probe> **old, int32_t newSize,
                    std::vector<Atom> *currentAtoms, uint32_t currentAtomCount);

  void clearProbe(Probe &p);

  void clearAtom(Atom &a);

  void clearProbeSet(ProbeSet &ps);

  void trackGlobalStats(ProbeSet &ps, int numBytes);

  /// pass these along to the ProbeListFactory for it to handle.
  /// @brief     query the Apid of a a probe (unique for this run only.)
  ///            The qPid, qAllele, qContext and qChannel must match.
  /// @param     qName     name of probeset
  /// @param     qPid      the query probeid (0-based)
  /// @param     qAllele   the query allele
  /// @param     qContext  the query context 
  /// @param     qChannel  the query channel
  /// @return    "-1" if not found.  The Apid if a probe matching is found.
  int findProbeApidByName(const std::string& qName,int qPid, int qAllele, int qContext, int qChannel) const {
    return m_PlFactory.findProbeApidByName(qName,qPid,qAllele,qContext,qChannel);
  }
  /// @brief     query the Apid of a a probe (unique for this run only.)
  ///            The qPid, qAllele, qContext and qChannel must match.
  ///            Just like findProbeApidByName, but the first arg is the probelist index.
  /// @param     qPlIdx    
  /// @param     qPid      
  /// @param     qAllele   
  /// @param     qContext  
  /// @param     qChannel  
  /// @return    "-1" if not found.  The Apid if a probe matching is found.
  int findProbeApidByIdx (int qPlIdx, int qPid, int qAllele, int qContext, int qChannel) const {
    return m_PlFactory.findProbeApidByIdx(qPlIdx,qPid,qAllele,qContext,qChannel);
  }


  /// @brief     Find the Apid for a probe given its ProbeId.  
  ///            This only works in the single channel case.
  ///            (Otherwise it is not a one-to-one map.)
  /// @param     qPid
  /// @return    -1 if not found.
  int findProbeApidByProbeId(int qPid) const  {
    return m_PlFactory.findProbeApidByProbeId(qPid);
  }


//  int findProbeApid2(const std::string& qName,int qPid,
//                     int& qAllele, int& qContext, int& qChannel) const {
//    return m_PlFactory.findProbe2(qName,qPid,qAllele,qContext,qChannel);
//  }

  // If we only have PlIdx and Pid, we could look it up with a method
  // which looked like this.
  // int findProbeApid2(int qPlIdx,int qPid,
  //                    int& oAllele, int& oContext, int& oChannel) const {
  //   return m_PlFactory.findProbe2(qName,qPid,oAllele,oContext,oChannel);
  // }

  /**
   * Keep track of global stats as we parse through file.
   *
   * @param name - Name of probeset
   * @param type - Class of probeset (as in ProbeSet::Type)
   * @param numBlocks - How many blocks are in the probeset.
   * @param numMatch - How many matches are there (1 for pm only, 2 for pm-mm)
   * @param probes - Vector of the probe ids/indexes in the array.
   */
  void trackGlobalStats(const std::string &name,
                        int type, int numBlocks,
                        const std::vector<int> &blockSizes,
                        int numMatch,
                        const std::vector<int> &probes,
                        bool multAlleles = false);

  static bool sameValues(std::vector<int> &v) {
    bool same = true;
    if(v.empty())
      return same;
    int first = v[0];
    for(size_t i = 1; i < v.size(); i++) {
      if(v[i] != first) {
        return false;
      }
    }
    return same;
  }
  
  int getMaxNameLength() {
    return m_MaxNameLength;
  }

  void incrementTypes(int type, int numBlocks, bool multipleAlleles);

  /// Number of expression probesets
  int getNumExpressionPSets() { return m_NumExpression; }
  /// Number of genotyping expression
  int getNumGenotypeingPSets() { return m_NumGenotyping; }

  /**
   * @brief compare probeset names according to respective probeset ids
   * @param a, b - probeset names to be ordered
   */
/*   static bool compareProbesetNameByProbesetIdx(std::string a, std::string b); */

  void sortProbesetNameByProbesetIdx(std::vector<std::string>::iterator begin, std::vector<std::string>::iterator end);

  ///< Do we need mismatch probes?
  bool m_NeedMismatch;
  ///< Do we need gc counts?
  bool m_NeedGc;
  ///< Do we need pm matching?
  bool m_NeedPmAlleleMatch;
  ///< Name of file.
  std::string m_FileName;
  ///< Header parameters from %thisKey=thisParam
  std::map<std::string, std::vector<std::string> > m_Header;
  ///< Verbatim lines from header.
  std::vector<std::string> m_HeaderLines;
  /// @todo: harley remove this
  ///< Vector of actual probe sets that got read information
  //std::vector<ProbeListPacked> m_PlVec;

  ProbeListFactory m_PlFactory;
  /// @todo change to colrow_t
  /// Number of values in X and Y axis of chip.
  unsigned int m_XCount, m_YCount;
  // this is a duplicate of the one in ProbeListFactory
  int m_numChannels;
  //int m_numCols; // X
  //int m_numRows; // Y
  /// Iterator to walk through the map.
  //typedef std::map<const char *, int, Util::ltstr>::iterator psNameIter;
  /// @todo: harley remove this
  /// map of probe sets by names.
  // std::map<const char *, int, Util::ltstr> m_PsNameHash;
  /// map of probe sets by ids
  // bool m_HaveNames;
  /// Are probe sets accesible by name?
  // bool m_HaveNameMap;

  /// Number of expression probesets
  int m_NumExpression;
  /// Number of genotyping expression
  int m_NumGenotyping;
  /// Pm probes in this layout (useful for some methods like RMA)
  std::vector<bool> m_PmProbes;
  /// Mm probes in this layout
  std::vector<bool> m_MmProbes;
  /// Map between PM and MM index on chip.
  std::vector<int> m_PmMmVec;
  /// Map between PM and PM for other allele
  std::vector<int> m_PmAlleleMatchVec;
  /// GC count of each probe as indexed by position.
  std::vector<char> m_ProbeGcVec;
  /// Sizes of various probesets (in probes)
  std::vector<int> m_ProbeCounts;
  /// Bitmask of probes that are actually in probesets.
  std::vector<bool> m_ProbesetProbes;
  CumulativeStats<double> m_ProbesetMemSizes;
  CumulativeStats<double> m_ProbesPerProbeset;
  int m_MaxNameLength;
  /// Vector of the types for each probeset, used by calvin chp files to
  /// determine size of file. Note that probesets marked as "GenoType" but only
  /// having one block are changed here to "Expression"
  // std::vector<char> m_PsTypes;

  /// Was the chip layout read in with a kill list
  bool m_haveKillList;
  std::vector<int> m_ProbeLayoutOrder;
  //
  std::string   m_SpfFileName;
  // @todo Do we need "m_SpfFile"? Cant we just allocate it when needed?
  affx::SpfFile m_SpfFile;
  // How many warnings about complete atom removal
  int m_warningCount;
};

#endif /* _CHIPLAYOUT_H_ */
