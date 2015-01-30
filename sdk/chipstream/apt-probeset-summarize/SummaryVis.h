////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

#ifndef SUMMARY_VIS_H
#define SUMMARY_VIS_H

//
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h"
#include "util/PgOptions.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cassert>
#include <cstring>
#include <ostream>
#include <string>
#include <vector>
//

/**
 * @file  SummaryVis.h
 * @brief Headers for SummaryVis.cpp.
 */
class summaryVis
{
public:
  /** Constructor.
   * @param argc Number of command line arguments.
   * @param argv Command line arguments.
   * @param version Version string.
   */

  summaryVis (int argc,
              const char* argv[],
              const std::string& version);

  /** Destructor.
   */
  ~summaryVis();

  /** Read, process data, write output.
   */
  void run();

private:

  /** Clear data.
   */
  void clear();

  /** Read optional transcript cluster id and probeset id files.
   */
  void readIdFiles (void);

  /** Begin output.
   */
  void beginOutput (void);

  /** Open input files.
   */
  void openInputFiles (void);

  /** Begin output file header.
   */
  void beginOutputHeader (void);

  /** Write header lines copied from input file.
   * @param tsv Tsv object for file.
   * @param idx File index.
   */
  void writeHeadersFromFile (affx::TsvFile& tsv, const unsigned int idx);

  /** Write column names.
   */
  void writeColumnNames (void);

  /** Write matching lines.
   */
  void writeMatches (void);

  /** Write first part of file related header line.
   * @param i File index.
   */
  void headerLine (const unsigned int i);

  /** fileRoot
   * Strip path (return root) from input file.
   * Prefixes using either windows or unix format path separators will be stripped,
   * but not a combination of the two.
   * @param str input string.
   * @return string root.
   */
//  std::string fileRoot (const std::string& str)
//  {
//    const char* unixSeparator = "/";
//    const char* windowsSeparator = "\\";
//    std::string::size_type pos = 0;
//    pos = str.rfind (unixSeparator);
//    if (pos != std::string::npos)
//      // Return root based on unix separators if any were found.
//      return str.substr (pos+1);
//
//    pos = str.rfind (windowsSeparator);
//    if (pos != std::string::npos)
//      // Else return root based on windows separators if any were found.
//      return str.substr (pos+1);
//
//    // Return full name if no separators were found.
//    return str;
//  }

  /**
   * @brief Class to store output data for sorting.
   */
  class summaryVisOutput
  {
  public:
    /** Constructor.
     * @param seqName Chromosome name.
     * @param start Start position.
     * @param stop Stop position.
     * @param data Summary data value.
     */
    summaryVisOutput (const std::string& seqName, const int start, const int stop, const double& data);
    const std::string seqName;
    const int start;
    const int stop;
    const double data;
  };

  /**
   * @brief Class for function predicate, to sort summaryVisOutput by the seqName,
   * then by start position.
   */
  class outputSortCriterion
  {
  public:
   bool operator() (const summaryVisOutput* o1, const summaryVisOutput* o2) const;
  };
  
  /// private data
  /// Version string.
  const std::string& m_Version;
  /// Command line options.
  PgOptions* m_Opts;
  /// Command line as a string.
  std::string m_CommandLine;
  /// Genome position file name.
  std::string m_GenomePosFileName;
  /// Transcript cluster id file name(s).
  std::vector<std::string> m_TranscriptClusterIdFileNames;
  /// Probeset id file name(s).
  std::vector<std::string> m_ProbesetIdFileNames;
  /// Output file name.
  std::string m_Outfile;
  /// Should filenames be prepended?
  bool m_Prepend;
  /// Should wiggle format output be generated?
  bool m_Wiggle;
  /// Column name to use for wiggle format output.
  std::string m_WiggleColName;
  /// Column index to use for wiggle format output.
  unsigned int m_WiggleColIndex;
  /// Column name used for wiggle format output.
  std::string m_WiggleUseColName;
  /// Wiggle column data.
  double m_WiggleData;
  /// Number of input files.
  unsigned int m_FileCount;
  /// Egr format version
  int m_EgrVersion;
  /// Summary file names.
  std::vector<std::string> m_FileNames;
  /// Output file stream.
  std::ofstream m_Out;
  /// Genome coordinates file tsv object.
  affx::TsvFile m_GenomePosTsv;
  /// Base summary file tsv object.
  affx::TsvFile m_BaseTsv;
  /// Base file column names.
  std::vector<std::string> m_BaseColNames;
  /// Base file probeset id.
  int m_BaseProbesetId;
  /// Base file data values for row.
  std::vector<double> m_BaseData;
  /// Output order for base file data.
  std::vector<int> m_BaseOutputOrder;
  /// Number of data columns in base file.
  unsigned int m_BaseDataColCount;
  /// Extra summary file tsv objects.
  std::vector<affx::TsvFile*> m_ExtraTsvs;
  /// Extra file column names.
  std::vector<std::vector<std::string> > m_ExtraColNames;
  /// Extra file probeset ids.
  std::vector<int> m_ExtraProbesetIds;
  /// Extra file data values for row.
  std::vector<std::vector<double> > m_ExtraData;
  /// Output order for extra file data.
  std::vector<std::vector<int> > m_ExtraOutputOrder;
  /// Number of data columns in extra files.
  std::vector<unsigned int> m_ExtraDataColCount;
  /// Execution guid.
  std::string m_ExecGuid;
  /// Output section separator.
  const std::string m_CommentLine;
  /// Map for transcript cluster ids.
  std::map<int, bool> m_TranscriptClusterIDs;
  /// Map for probeset ids.
  std::map<int, bool> m_ProbesetIDs;
  /// Chromosome name.
  std::string m_Seqname;
  /// Chromosome strand.
  std::string m_Strand;
  /// Start position.
  int m_Start;
  /// Stop position.
  int m_Stop;
  /// Genome coordinates file probeset_id.
  int m_GenomePosProbesetId;
  /// Genome coordinates file transcript_cluster_id.
  int m_GenomePosTranscriptClusterId;
  /// Genome version.
  std::string m_GenomeVersion;

};

#endif /* SUMMARY_VIS_H */
