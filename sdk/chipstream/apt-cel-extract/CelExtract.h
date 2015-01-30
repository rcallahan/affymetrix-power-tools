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

/// @file  CelExtract.h
/// @brief Headers for CelExtract.cpp.

#ifndef CEL_EXTRACT_H
#define CEL_EXTRACT_H

#include "chipstream/AnalysisStream.h"
#include "chipstream/AnalysisStreamExpression.h"
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/CelReader.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/IntensityMart.h"
//
#include "file/TsvFile/TsvFile.h"
//
#include <cstring>
#include <string>
#include <vector>
//

class CelExtractOptions
{
public:

  // Input files
  std::string cdfFile;                   ///< Name of cdf-file.
  std::string spfFile;                   ///< Name of spf file
  std::string pgfFile;                   ///< Name of pgf file
  std::string clfFile;                   ///< Name of clf file
  std::string bgpFile;                   ///< Name of bgp file
  std::vector<std::string> celFiles;     ///< Cel files to analyze.

  // Basic options control runtime
  std::string help;                      ///< Are we doing help? What kind?
  int verbosity;                         ///< Our verbosity level.
  bool force;                            ///< Don't check chip type.
  
  int m_output_precision; ///< precision of output

  // Set output content
  std::vector<std::string> probesetIdsFiles;           ///< File of probeset IDs/Names to load
  std::vector<std::string> probeIdsFiles;              ///< File of probe IDs to load

  bool pmOnly;                           ///< Do we only dump PM probes
  bool ignoreProbesWithoutMM;            ///< Do we ignore probes that do not have MM
  bool pairWithBackground;               ///< Do we show the background values (ie mm, gcbg, ...)
  bool subtractBackground;               ///< subtract background from value reported
  std::string outFile;                   ///< Output file

  // Analysis Stream Options
  std::string analysisString;            ///< An analysis string (no quant method) to use to extract intensities
  std::string sketchInFile;              ///< File specifying

  // Things filled in
  std::string progName;                  ///< Name of this program (generally what would appear in argv[0])
  std::string version;                   ///< Version of program.
  std::string cvsId;                     ///< Cvs id string.
  std::string commandLine;               ///< String we were called with.
  std::string timeStr;                   ///< Date and time of execution.
  std::string execGuid;                  ///< Execution id.

  // Disk Mode
  bool useDisk;
  std::string diskDir;
  int diskCache;

  /**
   * Constructor
   */
  CelExtractOptions();

  /**
   * Destructor
   */
  ~CelExtractOptions();

};

class CelExtract
{
public:
  /**
   * Constructor.
   */
  CelExtract(CelExtractOptions &o);

  /**
   * Destructor.
   */
  ~CelExtract();

  /**
   * Do the extraction
   */
  void extract();

  /**
   * Dump intensities iterating over probes
   * @param iMart - the intensity mart to use -- already allocated, but not loaded
   * @param tsv - the tsv object to write to
   */
  void doProbe(IntensityMart *iMart, std::vector<affx::TsvFile*>& tsv, bool bDoSaturationReport=false);

  /**
   * Dump intensities iterating over probe sets
   * @param iMart - the intensity mart to use -- already allocated, but not loaded
   * @param tsv - the tsv object to write to
   */
  void doProbeSet(IntensityMart *iMart, std::vector<affx::TsvFile*>& tsv);

  /**
   * Setup the TSV file for output, write headers, define columns, open ostream
   * @param tsv - TsvFile object
   */
  void defineExtractTsv(std::vector<affx::TsvFile*>& tsv, bool bDoSaturationReport=false);

  /**
   * Set the various probe level annotations for a specific probeset
   * @param tsv - the TsvFile object
   * @param columnCount - the column to start with. will be updated.
   * @param ps - pointer to the probeset
   * @param aIx - the atom index (relative to probeset)
   * @param plIx - the probe index (relative to probeset)
   */
  void tsvSetProbeAnnotations(affx::TsvFile* tsv, int &columnCount, ProbeSet *ps, int aIx, int plIx);

private:
  /**
   * Check the options for correctness
   */
  void checkOptions(CelExtractOptions &o);

  /**
   * Create an analysis stream object
   */
  AnalysisStreamExpression *newAnalysisStreamObject();

  /// order of probesets if probeset list given
  std::vector<int> m_PsOrder_jhg;
  /// our copy of the options
  CelExtractOptions m_Options;
  /// the number of rows on the chip
  int m_NumRows;
  /// the number of columns on the chip
  int m_NumCols;
  /// the number of features on the chip
  int m_NumFeatures;
  /// number of data channels in CEL file
  int m_NumChannels;
  /// column index of first CEL file
  int m_FirstCelColIndex;
  /// vector of output file names for each channel
  std::vector<std::string> m_OutFiles;
  /// mask of which probes to dump
  std::vector<bool> m_ProbesToDump;
  /// our chip layout info. NULL if no layout.
  ChipLayout *m_Layout;
  /// pointer to analysis stream
  AnalysisStreamExpression *m_as;
  /// pointer to analysis stream
  AnalysisStreamFactory *m_asFactory;
  /// list of gc probe IDs
  std::vector<int> m_gcProbes;
  /// Order of probes for disk mart
  std::vector<int> m_ProbeOrder;

};

#endif //CEL_EXTRACT_H
