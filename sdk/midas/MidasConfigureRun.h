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

/// @file   MidasConfigureRun.h
/// @brief  Headers for configuring, running the
///         MidasEngine iterator of the
///         MidasSpliceDetector statistical analysis
///         of exon data.

#ifndef MIDASCONFIGURE_H
#define MIDASCONFIGURE_H

#include "file/TsvFile/TsvFile.h"

#define PVALUES_OUTPUT    "midas.pvalues.txt"
#define FSTATS_OUTPUT     "midas.fstats.txt"
#define NORMALIZED_OUTPUT "midas.normalized.txt"

std::string MidasFileRoot(const std::string& str);

/**
 *  midasConfigureRun
 *  @brief Object for midasConfigureRun.
 *
 */
class midasConfigureRun
{
public:
  /** Constructor.
   * @param celsFileName Cels file name.
   * @param geneDataFileName Gene data file name.
   * @param exonDataFileName Exon data file name.
   * @param metaFileName Meta probeset file name.
   * @param outputDirectory Output directory.
   * @param wantPvalues Caller requests p values.
   * @param wantFstats Caller requests F statistics.
   * @param wantNormalized Caller requests normalized exon signal.
   * @param logStabilize Amount to add to data before taking logarithm.
   * @param commandLine User command line.
   * @param execVersion Version, cvs id string.
   * @param noLogTransform Do not log transform data.
   * @param keepPath Keep cel file path.
   */
  midasConfigureRun (std::string celsFileName, std::string geneDataFileName,
                     std::string exonDataFileName, std::string metaFileName, 
                     std::string outputDirectory, const bool wantPvalues,
                     const bool wantFstats, const bool wantNormalized, 
                     const float& logStabilize, const std::string& commandLine, 
                     const std::string& execVersion, const bool noLogTransform = false,
                     const bool keepPath = false);

  /** Set up for run.
   * @return pointer to non-fatal warning message
   */
  std::string* configure();

  /** Read, process data, write output.
   */
  void run();

  /** Delete output files.
   */
  void deleteOutputs_1(bool want,std::string filename);
  void deleteOutputs(); // unused

  void throwIfFileExists(std::string fileName);

private:

  /// @brief     Add a set of metaTags to the tsvfile
  /// @param     tsvfile   file to add the tags to
  void addMetaTags(affx::TsvFile& tsvfile);

  /** checkMetaFile
   * Validate meta probeset file - must contain probeset_id, probeset_list columns.
   */
  void checkMetaFile();

  /** readCelsFile
   * Read cels.txt file, generate map between cel_file name and group_id string.
   */
  void readCelsFile();

  /** readGeneDataFile
   * Read genedata file - for each cel_file name found in both the genedata.txt and
   * cels.txt files, assign an (int) group_id index, write the cel_file name to
   * the geneDataCelFiles vector and the (int) group_id index to the groupId Ints
   * vector.  celFilesUsed is the size of both of these vectors.
   */
 void readGeneDataFile();

  /** readExonDataFile
   * read exon data file header - here we're only validating that the list
   * of cel_file names matches those found in the gene data file.
   */
  void readExonDataFile();

  /** openOutput
   * Attempt to open requested output file, write headers.
   * @param fileName string containing name of file.
   * @param Ofstream output fstream to open, write to.
   */
  //void openOutput (const std::string& fileName, std::ofstream& Ofstream);

  /** addWarning
   * If a warning string is already present, concatenate the message to it,
   * else create a new warning string containing the message
   */
  void addWarning (const std::string& message);

  /** propagateMetaTags
   * Propagate meta tags from a given input file, prepending the supplied tagPrefix.
   * If a tagSubset vector is supplied, propagate only those tags.
   * @param tagPrefix string to prepend before existing tags
   * @param tagSubset optional subset of tags to propagate
   * @param inputFileTSV TSV object of file to propagate tags from
   * @param optionalString optional tag string, e.g., data summary file name
   */
  void propagateMetaTags (affx::TsvFile& destTsv,std::string prefix, affx::TsvFile& srcTsv);
  //affxtsv::CTSVFileData& intputFileTSV, const std::string& optionalString);

  // map between group_id string and group index (int)
  std::map <std::string,int> groupIdStringToIntMap;
  /// The end of the group ids (as used above)
  int groupIdMax;

  int allocGroupId(std::string groupIdStr);

  /// private data
  /// Cels file name
  std::string celsFileName;
  affx::TsvFile celsTsv;
  /// Gene data file name
  std::string geneDataFileName;
  /// Exon data file name
  std::string exonDataFileName;
  /// meta probeset file name
  std::string metaFileName;
  /// Output directory
  std::string outputDirectory;
  /// caller requests p values
  const bool wantPvalues;
  /// caller requests F statistics
  const bool wantFstats;
  /// caller requests normalized signal
  const bool wantNormalized;
  /// Amount to add to data before taking logarithm
  const float logStabilize;
  /// User command line
  const std::string commandLine;
  /// caller requests no log transform of data
  const bool noLogTransform;
  /// meta FileData object
  affx::TsvFile metaTsv;
  /// gene data FileData object
  affx::TsvFile geneDataTsv;
  /// exon data FileData object
  affx::TsvFile exonDataTsv;
  /// p values output stream
  std::ofstream pvaluesOfstream;
  /// fstats output stream
  std::ofstream fstatsOfstream;
  /// normalized output stream
  std::ofstream normalizedOfstream;
  /// map between cel_file and group_id string
  std::map <std::string, std::string> celFileGroupIdStringMap;
  /// map between cel_file name and group index (int)
  /// this map contains only the cel_file names found in both the cels.txt file
  /// and the genedata files, mapped to the (int) group index assigned by
  /// readGeneDataFile()
  std::map <std::string, int> celFileGroupIdIntMap;
  /// vector of cel_file names found in the geneData file
  /// this is passed to the midasEngine
  std::vector <std::string> geneDataCelFiles;
  /// vector of group_id ints in the same order as geneData cel_file names
  /// this is passed to the midasEngine
  std::vector <int> groupIdInts;
  /// warning message to return to caller
  std::string* pWarningMessage;
  /// ascii time string
  std::string timeString;
  /// globally unique identifier for this run
  std::string execGuid;
  /// count of cel_files used (data columns) in gene data
  /// used to validate cel_files used in exon data
  int celFilesUsedInGeneData;
  /// count of distinct group ids used (distinct groups)
  int groupIdsUsed;
  /// pvalues output file name
  std::string pvaluesFileName;
  /// fstats output file name
  std::string fstatsFileName;
  /// normalized output file name
  std::string normalizedFileName;
  /// meta tags to propagate
  affx::TsvFile metatagBucketTsv;
  
  /// comment string for headers
  const std::string comment;
  /// Tags to propagate from gene, exon data files
  std::vector<std::string> dataFileTags;
  /// version, cvs_id string for exec_version
  const std::string execVersion;
  /// keep cel file path
  const bool keepPath;
};

#endif /* MIDASCONFIGURE_H */
