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

/// @file   MidasEngine.h
/// @brief  Headers for MidasEngine.
///         Iterator which reads input data files, calls
///         the MidasSpliceDetector for each probeset_id (gene) and
///         probeset_list_id (exon) data pair, writes output to file.

#ifndef MIDASENGINE_H
#define MIDASENGINE_H

//
#include "midas/MidasSpliceDetector.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/PgOptions.h"
//
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
//

// hard coded column names
#define PROBESET_ID "probeset_id"
#define PROBESET_LIST "probeset_list"
#define PROBESET_LIST_ID "probeset_list_id"

/// @brief define the midas options.
void define_midas_opts(PgOptions* opts);

/**
 *  midasEngine
 *  \brief Object for midas engine.
 *
 */
class midasEngine
{
public:
  /** Constructor.
   * @param celFiles List of cel_file names.
   * @param groups List of group indices.
   * @param metaFileName Meta probeset file name.
   * @param geneDataFileName Gene data file name.
   * @param exonDataFileName Exon data file name.
   * @param pvaluesFileName Output file name for p values.
   * @param fstatsFileName Output file name for F statistics.
   * @param normalizedFileName Output file name for normalized exon signal.
   * @param wantPvalues Caller requests p values.
   * @param wantFstats Caller requests F statistics.
   * @param wantNormalized Caller requests normalized exon signals.
   * @param logStabilize Amount to add to data before taking logarithm.
   * @param noLogTransform Do not log transform data.
   */

  /* Note: the parameters wantPvalues, wantFstats, wantNormalized, and noLogTransform,
     user configurable options, are functionally either true or false and could be
     treated as bools.  Here they are defined as ints so that this code could be
     swig wrapped; in Visual Studio builds of swig wrapped code, using bools for these
     parameters produced unresolved external symbol errors at link time. */

  midasEngine (const std::vector<std::string>& celFiles, const std::vector<int>& groups,
      std::string metaFileName, std::string geneDataFileName, std::string exonDataFileName,
      std::string pvaluesFileName, std::string fstatsFileName, std::string normalizedFileName,
      int wantPvalues, int wantFstats, int wantNormalized, const float& logStabilize,
      int noLogTransform);

  /** Destructor.
   */
  ~midasEngine ();

  /** Read, process data, write output.
   */
  void analyze();

  // a place to park metatags
  affx::TsvFile metatagBucketTsv;

private:

  /** Open output file to append
   * @param fileName File name.
   * @param Ofstream Ofstream.
   */
  //void openOutputAppend (const char* fileName, ofstream& Ofstream);
  void openOutput_Pvalues(std::string filename);
  void openOutput_Fstats(std::string filename);
  void openOutput_Normalized(std::string filename);

  /** Read meta probeset file
   */
  void readMeta ();

  /** Read, process gene data
   */
  void readGeneData ();

  /** Process exon data
   */
  void processExonData ();

  /// private data
  /// list of cel_file names used
  const std::vector<std::string> celFiles;
  /// experimental groups
  std::vector<int> groups;
  /// meta probeset FileData object
  affx::TsvFile metaTsv;
  /// meta probeset file name
  std::string metaFileName;
  /// gene FileData object
  affx::TsvFile geneDataTsv;
  /// gene data file name
  std::string geneDataFileName;
  /// exon FileData object
  affx::TsvFile exonDataTsv;
  /// exon data file name
  std::string exonDataFileName;
  /// pvalues output file name
  std::string pvaluesFileName;
  /// fstats output file name
  std::string fstatsFileName;
  /// normalized output file name
  std::string normalizedFileName;

  /// pvalues output fstream
  affx::TsvFile pvaluesTsv;
  /// fstats output fstream
  affx::TsvFile fstatsTsv;
  /// normalized output fstream
  affx::TsvFile normalizedTsv;

  /// caller requests p values
  const int wantPvalues;
  /// caller requests F statistics
  const int wantFstats;
  /// caller requests normalized exon signal
  const int wantNormalized;
  /// Amount to add to data before taking logarithm
  const float logStabilize;
  /// caller requests no log transform of data
  const int noLogTransform;
  /// SpliceDetector object
  midasSpliceDetector spliceDetector;
  /// number of experimental groups
  const unsigned int numGroups;
  /// map from probeset_list id to vector of probeset ids
  std::map<int, std::vector<int>* > probesetListIdToProbesetIdVectorMap;
  /// map from gene id to row number
  std::map<int, int> geneIdToRowMap;
  /// vector of vectors containing gene data
  std::vector<std::vector<float> > geneData;
  /// fstream for output
  std::ofstream out;
};

#endif /* MIDASENGINE_H */
