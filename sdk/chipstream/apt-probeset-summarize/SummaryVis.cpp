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

/// @file   SummaryVis.cpp
/// @brief  Class for merging multiple tab separated text files.

//
#include "chipstream/apt-probeset-summarize/SummaryVis.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/Util.h"
//
#include <algorithm>
#include <cstring>
#include <set>
#include <string>

using namespace std;
using namespace affx;

summaryVis::summaryVisOutput::summaryVisOutput (const std::string& seqName, const int start, const int stop, const double& data)
      : seqName (seqName), start (start), stop (stop), data (data)
    {}


bool summaryVis::outputSortCriterion::operator() (const summaryVis::summaryVisOutput* o1, const summaryVis::summaryVisOutput* o2) const
{
  if (o1->seqName == o2->seqName)
    return o1->start < o2->start;
  // Want upper case, e.g. X or Y, last.
  const char* upperCase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  if (o1->seqName.find_first_of (upperCase) != std::string::npos)
    {
      if (o2->seqName.find_first_of (upperCase) != std::string::npos)
        // If both contain upper case, do lexicographical comparison.
        return o1->seqName < o2->seqName;
      else
        // The first seqName contains upper case, not the second.
        return false;
    }
  else if (o2->seqName.find_first_of (upperCase) != std::string::npos)
    // The second seqName contains upper case, not the first.
    return true;
  // No upper case - want shorter first, e.g. chr9 before chr10.
  else
    {
      const unsigned int o1Size = o1->seqName.size();
      const unsigned int o2Size = o2->seqName.size();
      if (o1Size != o2Size)
        return o1Size < o2Size;
      else
        return o1->seqName < o2->seqName;
    }
}

/** Destructor.
 */
summaryVis::~summaryVis()
{
    clear();
}

/** Clear data.
 */
void summaryVis::clear()
{
  delete m_Opts;
  for (unsigned int i = 0; i < m_ExtraTsvs.size(); ++i)
    {
      m_ExtraTsvs[i]->close();
      delete m_ExtraTsvs[i];
    }
}

/** Write first part of file related header line.
 * @param i File index.
 */
void summaryVis::headerLine (const unsigned int i)
{
  m_Out << "# file" << i << "_";
}

void define_summaryvis_optsion(PgOptions* opts)
{
  opts->setUsage("apt-summary-vis - Create an IGB egr file from summary and annotation files.\n"
                 "Usage:\n"
                 "   apt-summary-vis -g genome_position_file -o output_file summary_file [...]");
  
  opts->defineOption("g", "genome-position", PgOpt::STRING_OPT,
                     "Name of the file containing genome position information. "
                     "At a minimum this file should contain columns probeset_id, "
                     "seqname, strand, start, and stop.",
                     "");
  opts->defineOption("o", "out-file", PgOpt::STRING_OPT,
                     "Output file to send the merged results to.",
                     "");
  opts->defOptMult("", "transcript-cluster-ids", PgOpt::STRING_OPT,
                   "Optional name of a file containing transcript_cluster_ids "
                   "to include in the egr file.  Can be specified multiple times.",
                   "");
  opts->defOptMult("", "probeset-ids", PgOpt::STRING_OPT,
                   "Optional name of a file containing probeset_ids to include "
                   "in the egr file.  Can be specified multiple times.",
                   "");
  opts->defineOption("", "egr-version", PgOpt::INT_OPT,
                     "Optional type of egr file to generate. For version 2 "
                     "both the probeset id and genome coordinates are included."
                     "Version 3 omits genome coordinates "
                     "but requires that the design be loaded into IGB "
                     "before the egr file. Version 1 is deprecated; it "
                     "has genome coordinates but no probeset id.",
                     "2");
  opts->defineOption("p", "prepend-filename", PgOpt::BOOL_OPT,
                     "Prepend the filename to the column header.",
                     "false");
  opts->defineOption("", "wiggle-col-name", PgOpt::STRING_OPT,
                     "Generate wiggle format output, compatible with "
                     "the UCSC genome browser, using the summary file "
                     "column with the specified name. The generated "
                     "wiggle output will not work in the UCSC genome "
                     "browser if the extracted probesets overlap.",
                     "");
  opts->defineOption("", "wiggle-col-index", PgOpt::INT_OPT,
                     "Generate wiggle format output using the summary file "
                     "column with the specified index. Numbering is "
                     "one-based, ignoring the probeset_id column. "
                     "The generated wiggle output will not work in "
                     "the UCSC genome browser if the extracted "
                     "probesets overlap.",
                     "0");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "Print help message.",
                     "false");
}

/**
 *  @brief Constructor.
 *
 *  @param argc Number of command line arguments.
 *  @param argv Command line arguments.
 *  @param version Version string.
 *
 *  Errors: throw exception to display help messages, if too few input files.
*/
summaryVis::summaryVis (int argc,
                        const char* argv[],
                        const std::string& version)
  : m_Version (version),
  m_CommentLine ("############################################################\n")
{

  // Prefer throw() to exit().
  Err::setThrowStatus (true);

  m_WiggleData = 0;
  m_BaseProbesetId = -1;
  m_BaseDataColCount = 0;
  m_Start = 0;
  m_Stop = 0;
  m_GenomePosProbesetId = -1;
  m_GenomePosTranscriptClusterId = -1;

  //
  m_Opts = new PgOptions();
  define_summaryvis_optsion(m_Opts);
  m_Opts->parseArgv(argv);

  // Optionally display usage message.
  if (m_Opts->getBool("help") || argc == 1)
  {
    m_Opts->usage();
    string msg = "version: " + version + "\n";
    cout << msg;
    exit(0);
  }
  // Optionally display version.
  if (m_Opts->getBool("version"))
  {
    string msg = "version: " + version + "\n";
    cout << msg;
    exit(0);
  }

  // Require genome position file.
  m_GenomePosFileName = m_Opts->get("genome-position");
  if (m_GenomePosFileName.empty())
  {
    string msg = "FATAL: Must provide --genome-position file.";
    Err::errAbort (msg);
  }
  // Require at least one summary file.
  m_FileCount = (unsigned int)m_Opts->getArgCount();
  if (m_FileCount < 1)
  {
    string msg = "FATAL: Must provide summary files.";
    Err::errAbort (msg);
  }
  for (unsigned int i = 0; i < m_FileCount; ++i)
    m_FileNames.push_back (string (m_Opts->getArg(i)));

  // Egr format version.
  m_EgrVersion = m_Opts->getInt("egr-version");
  if (m_EgrVersion == 0 || m_EgrVersion > 3)
  {
    string msg = "FATAL: Unknown egr format version " + ToStr (m_EgrVersion) + ".";
    Err::errAbort (msg);
  }

  // Require a writeable output file.
  m_Outfile = m_Opts->get("out-file");
  if (m_Outfile.empty())
  {
    string msg = "FATAL: Must provide an output file.";
    Err::errAbort (msg);
  }
  else
    Fs::mustOpenToWrite(m_Out, m_Outfile);

  // Save optional prepend-filename flag.
  m_Prepend = m_Opts->getBool("prepend-filename");

  // Set up for optional wiggle format output.
  m_WiggleColName = m_Opts->get("wiggle-col-name");
  m_WiggleColIndex = m_Opts->getInt("wiggle-col-index");
  if (! m_WiggleColName.empty() || m_WiggleColIndex > 0)
    m_Wiggle = true;
  else
    m_Wiggle = false;
  if (! m_WiggleColName.empty() && m_WiggleColIndex > 0)
  {
    string msg = "FATAL: For wiggle format output, provide either a wiggle-col-name ";
    msg += "or a wiggle-col-index, but not both.";
    Err::errAbort (msg);
  }

  // Save optional id file names.
  //for (PgOpt* fileOpt = m_Opts->findOpt ("transcript-cluster-ids"); fileOpt != NULL; fileOpt = fileOpt->next) {
  //  printf("transcript-cluster-ids=%p\n",fileOpt);
  //  // Discard empty strings.
  //  if ((fileOpt->getValue()!=NULL) && (strlen(fileOpt->getValue())!=0))
  //    m_TranscriptClusterIdFileNames.push_back (fileOpt->getValue());
  //}
  PgOpt* opt;
  opt=m_Opts->mustFindOpt("transcript-cluster-ids");
  for (int i=0;i<opt->getValueCount();i++) {
    printf("transcript-cluster-ids=%p\n",opt->getValue(i).c_str());
    m_TranscriptClusterIdFileNames.push_back(opt->getValue(i));
  }
  //  for (PgOpt* fileOpt = m_Opts->findOpt ("probeset-ids"); fileOpt != NULL; fileOpt = fileOpt->next)
  //    if (strlen (fileOpt->getValue()) != 0)
  //      m_ProbesetIdFileNames.push_back (fileOpt->getValue());
  opt=m_Opts->mustFindOpt("probeset-ids");
  opt->push_user_values_into(m_ProbesetIdFileNames);

  //
  m_CommandLine=m_Opts->commandLine();
}

/**
 *  @brief Read, process data, write output.
 *
 *  Errors: abort if unable to read any input file.
*/
void summaryVis::run()
{
  // Write initial output.
  beginOutput();

  // Read optional transcript cluster id and probeset id files.
  readIdFiles();

  // Open input files.
  openInputFiles();

  // Generate initial output header lines.
  beginOutputHeader();

  // Write file specific header lines.
  unsigned int idx = 0;
  writeHeadersFromFile (m_BaseTsv, idx);
  for (idx = 0; idx < m_FileCount - 1; ++idx)
    writeHeadersFromFile (*m_ExtraTsvs[idx], idx + 1);
  m_Out << m_CommentLine;

  // Write genome version if available.
  if (! m_GenomeVersion.empty())
  {
    Verbose::out (1, "Genome position file '" + m_GenomePosFileName
      + "' has a genome version tag of '" + m_GenomeVersion + "'.");
    m_Out << "# genome_version = " << m_GenomeVersion << "\n";
  }
  else
    Verbose::out (1, "WARNING: No genome_version tag found in genome position file '"
      + m_GenomePosFileName + "'.");

  // Write matching lines.
  writeMatches();

  // Close output file.
  m_Out.close();
}

/**
 *  @brief Read optional transcript cluster id and probeset id files.
*/
void summaryVis::readIdFiles()
{
  // Transcript cluster id file(s).
  for (unsigned int i = 0; i < m_TranscriptClusterIdFileNames.size(); ++i)
  {
    string& clusterFileName = m_TranscriptClusterIdFileNames[i];
    Verbose::out (1, "Processing " + clusterFileName + " for transcript clusters.");
    TsvFile tsv;
    int transcriptClusterId;
    tsv.bind (0, "transcript_cluster_id", &transcriptClusterId, TSV_BIND_REQUIRED);
    if (tsv.open (clusterFileName) != TSV_OK)
      Err::errAbort ("Problem opening transcript cluster file " + clusterFileName);
    while (tsv.nextLevel (0) == TSV_OK)
      m_TranscriptClusterIDs [transcriptClusterId] = true;
    tsv.close();
  }
  // Probeset id file(s).
  for (unsigned int i = 0; i < m_ProbesetIdFileNames.size(); ++i)
  {
    string& probesetIDFileName = m_ProbesetIdFileNames[i];
    Verbose::out (1, "Processing " + probesetIDFileName + " for probesets.");
    TsvFile tsv;
    int probesetId;
    tsv.bind (0, "probeset_id", &probesetId, TSV_BIND_REQUIRED);
    if (tsv.open (probesetIDFileName) != TSV_OK)
      Err::errAbort ("Problem opening probeset id file " + probesetIDFileName);
    while (tsv.nextLevel (0) == TSV_OK)
      m_ProbesetIDs [probesetId] = true;
    tsv.close();
  }
}

/**
 *  @brief Begin output.
*/
void summaryVis::beginOutput()
{
  Verbose::out (1, "MODULE: " + m_Version);
  Verbose::out (1, "CMD: " + m_CommandLine);
  m_ExecGuid = affxutil::Guid::GenerateNewGuid();
  Verbose::out (1, "exec_guid " + m_ExecGuid);
  if (m_Wiggle)
    Verbose::out (1, "Generating wiggle format output.");
}

/**
 *  @brief Open input files.
 *
 *  Errors: abort if a file could not be opened or if
 *  a required column name was not found.
*/
void summaryVis::openInputFiles()
{
  // Bind columns of interest in genome position file, open it.
  Verbose::out (1, "Opening genome coordinates file " + m_GenomePosFileName + ".");
  m_GenomePosTsv.bind (0, "probeset_id", &m_GenomePosProbesetId, TSV_BIND_REQUIRED);
  m_GenomePosTsv.bind (0, "seqname", &m_Seqname, TSV_BIND_REQUIRED);
  m_GenomePosTsv.bind (0, "strand", &m_Strand, TSV_BIND_REQUIRED);
  m_GenomePosTsv.bind (0, "start", &m_Start, TSV_BIND_REQUIRED);
  m_GenomePosTsv.bind (0, "stop", &m_Stop, TSV_BIND_REQUIRED);
  // Bind transcript_cluster_id if the user has provided a transcript cluster id list.
  if (! m_TranscriptClusterIDs.empty())
    m_GenomePosTsv.bind (0, "transcript_cluster_id", &m_GenomePosTranscriptClusterId, TSV_BIND_REQUIRED);
  if (m_GenomePosTsv.open (m_GenomePosFileName) != TSV_OK)
    Err::errAbort ("Problem opening file " + m_GenomePosFileName + ".");
  // Will be indexing over probeset_id.
  m_GenomePosTsv.defineIndex (0, "probeset_id", TSV_INDEX_INT, 0);

  // Bind probeset_id, open, get column names of base summary file.
  string& baseFileName = m_FileNames[0];
  Verbose::out (1, "Opening summary file " + baseFileName + ".");
  m_BaseTsv.bind (0, "probeset_id", &m_BaseProbesetId, TSV_BIND_REQUIRED);
  if (m_BaseTsv.open (baseFileName) != TSV_OK)
    Err::errAbort ("Problem opening file " + baseFileName);

  unsigned int colCount = m_BaseTsv.getColumnCount(0);
  // All summary columns except the probeset_id are data.
  m_BaseColNames.resize (colCount - 1);
  if (! m_Wiggle)
  {
    m_BaseData.resize (colCount - 1);
    m_BaseOutputOrder.resize (colCount - 1);
  }

  bool foundWiggleCol = false;
  // Do not require that probeset_id be the first column.
  unsigned int k = 0;
  for (unsigned int i = 0; i < colCount; ++i)
  {
    string colName;
    m_BaseTsv.cidx2cname (0, i, colName);
    if (colName != "probeset_id")
    {
      m_BaseColNames[k] = colName;
      if (m_Wiggle)
      {
	++k;
	if ( (m_WiggleColName == colName) || (m_WiggleColIndex == k) )
	{
	  // Found match.
	  foundWiggleCol = true;
	  m_WiggleUseColName = colName;
	  // Bind only this column.
	  m_BaseTsv.bind (0, i, &m_WiggleData);
	}
      }
      else
	m_BaseTsv.bind (0, i, &m_BaseData[k++]);
    }
  }
  // Remember the number of data columns.
  m_BaseDataColCount = colCount - 1;
  // Offset into extra files for wiggle column.
  int wiggleColOffset = m_WiggleColIndex - m_BaseDataColCount;

  // Store column names sorted - make copy to build
  // index between data read from file, values sorted by
  // column name.
  vector<string> colNames (m_BaseColNames);
  // Not checking that column names are unique.
  sort (m_BaseColNames.begin(), m_BaseColNames.end());

  typedef vector<string>::const_iterator VSIter_t;
  // Don't set up output order for wiggle output.
  if (! m_Wiggle)
  {
    const VSIter_t unsortedBegin = colNames.begin();
    const VSIter_t unsortedEnd = colNames.end();
    for (unsigned int i = 0; i < m_BaseDataColCount; ++i)
    {
      const string& colName = m_BaseColNames[i];
      m_BaseOutputOrder[i] = distance (unsortedBegin, find (unsortedBegin, unsortedEnd, colName));
    }
  }

  assert (m_FileCount > 0);
  m_ExtraProbesetIds.resize (m_FileCount - 1);
  // Open extra summary files if present.
  for (unsigned int i = 0; i < m_FileCount - 1; ++i)
  {
    string& extraFileName = m_FileNames[i + 1];
    Verbose::out (1, "Opening summary file " + extraFileName + ".");
    TsvFile* extraTsv = new TsvFile;
    m_ExtraTsvs.push_back (extraTsv);
    extraTsv->bind (0, "probeset_id", &m_ExtraProbesetIds[i], TSV_BIND_REQUIRED);
    if (extraTsv->open (extraFileName) != TSV_OK)
      Err::errAbort ("Problem opening file " + extraFileName);

    colCount = extraTsv->getColumnCount(0);
    m_ExtraColNames.push_back (vector<string> (colCount - 1));
    if (! m_Wiggle)
    {
      m_ExtraData.push_back (vector<double> (colCount - 1));
      m_ExtraOutputOrder.push_back (vector<int> (colCount - 1));
    }
    // Set up index.
    extraTsv->defineIndex (0, "probeset_id", TSV_INDEX_INT, 0);
    // Remember the number of data columns.
    m_ExtraDataColCount.push_back (colCount - 1);

  }  // end for (unsigned int i = 0; i < m_FileCount - 1; ++i)

  // Bind to the data vectors after setting up them up: push_back(), also resize(),
  // may cause reallocation, which will change the addresses.
  for (unsigned int i = 0; i < m_FileCount - 1; ++i)
  {
    TsvFile* extraTsv = m_ExtraTsvs[i];
    k = 0;
    // Add one since had stored number of data columns.
    unsigned int colCount = m_ExtraDataColCount[i] + 1;

    for (unsigned int j = 0; j < colCount; ++j)
    {
      string colName;
      extraTsv->cidx2cname (0, j, colName);
      if (colName != "probeset_id")
      {
	m_ExtraColNames[i][k] = colName;
	if (m_Wiggle)
	{
	  ++k;
	  // Not checking that extra file column names are different from
	  // base column names or each other.
	  if ( (m_WiggleColName == colName) || (wiggleColOffset == (int)k) )
	  {
	    // Found match.
	    foundWiggleCol = true;
	    m_WiggleUseColName = colName;
	    // Bind only this column.
	    extraTsv->bind (0, j, &m_WiggleData);
	  }
	}
	else
	  // These binds must be done after m_ExtraData allocation is complete.
          extraTsv->bind (0, j, &m_ExtraData[i][k++]);
      }
    }
    // Sort column names, build output order index.
    vector<string> extraColNames (m_ExtraColNames[i]);
    sort (m_ExtraColNames[i].begin(), m_ExtraColNames[i].end());
    if (! m_Wiggle)
    {
      const VSIter_t extraUnsortedBegin = extraColNames.begin();
      const VSIter_t extraUnsortedEnd = extraColNames.end();
      const unsigned int extraColCount = m_ExtraDataColCount[i];
      for (unsigned int l = 0; l < extraColCount; ++l)
      {
        const string& colName = m_ExtraColNames[i][l];
        m_ExtraOutputOrder[i][l]
  	= distance (extraUnsortedBegin, find (extraUnsortedBegin, extraUnsortedEnd, colName));
      }
    } // end for (unsigned int j = 0; j < colCount; ++j)
    wiggleColOffset -= m_ExtraDataColCount[i];
  }  // end for (unsigned int i = 0; i < m_FileCount - 1; ++i)

  // Quit if wiggle output requested and the requested column was not found.
  if (m_Wiggle && ! foundWiggleCol)
    Err::errAbort ("No match found for the column requested for wiggle output");
}

/**
 *  @brief Begin output file header.
*/
void summaryVis::beginOutputHeader()
{
  // Generic apt meta tags.
  const string guid = affxutil::Guid::GenerateNewGuid();
  m_Out << "# guid = " << guid << "\n";
  m_Out << "# exec_guid = " << m_ExecGuid << "\n";
  m_Out << "# exec_version = " << m_Version << "\n";
  m_Out << "# create_date = " << Util::getTimeStamp() << "\n";
  m_Out << "# cmd = " << m_CommandLine << "\n";
  m_Out << m_CommentLine;

  // Genome position file meta tags.
  m_Out << "# genome_pos_file_name = " << m_GenomePosFileName << "\n";
  // Copy genome position file header meta tags to output.
  string key, value;
  m_GenomePosTsv.headersBegin();
  while (m_GenomePosTsv.headersNext (key, value) == TSV_OK)
  {
    // Comment lines among the headers give rise to a pair with a null key.
    if (! key.empty())
    {
      // Save genome version.
      if ( (key == "genome-version") && m_GenomeVersion.empty())
	m_GenomeVersion = value;
      // In key, replace '-' with '_'.
      string::size_type pos = 0;
      for (;;)
      {
	pos = key.find ('-', pos);
	if (pos == string::npos)
	  break;
	key[pos] = '_';
      }
      m_Out << "# genome_pos_file_" << key << " = " << value << "\n";
    }
  }
}

/**
 *  @brief Write header lines copied from input file.
*/
void summaryVis::writeHeadersFromFile (TsvFile& tsv, const unsigned int idx)
{
  m_Out << m_CommentLine;
  headerLine (idx);
  m_Out << "name = " << m_FileNames[idx] << "\n";

  // Copy file header meta tags to output.
  string key, value;
  tsv.headersBegin();
  while (tsv.headersNext (key, value) == TSV_OK)
    // Comment lines among the headers give rise to a pair with a null key.
    if (! key.empty())
    {
      // In key, replace '-' with '_'.
      string::size_type pos = 0;
      for (;;)
      {
	pos = key.find ('-', pos);
	if (pos == string::npos)
	  break;
	key[pos] = '_';
      }
      headerLine (idx);
      m_Out << key << " = " << value << "\n";
    }
}

/**
 *  @brief Write summary file column names.
*/
void summaryVis::writeColumnNames()
{
  unsigned int i = 0;

  // Columns of the base file.
  for (i = 0; i < m_BaseDataColCount; ++i)
  {
    m_Out << "# score" << i << " = ";
    // Prepend file name if requested, without path.
    if (m_Prepend)
      m_Out << ": " << Fs::basename(m_FileNames[0]) << ": ";
    // Strip path.
    m_Out << Fs::basename (m_BaseColNames[i]) << "\n";
  }

  // Columns of extra files.
  for (unsigned int idx = 0; idx < m_FileCount - 1; ++idx)
  {
    const unsigned int extraColCount = m_ExtraDataColCount[idx];
    for (unsigned int j = 0; j < extraColCount; ++j)
    {
      // Continue count from base file.
      m_Out << "# score" << i++ << " = ";
      if (m_Prepend) {
        m_Out << ": " << Fs::basename(m_FileNames[idx + 1]) << ": ";
      }
      m_Out << Fs::basename(m_ExtraColNames[idx][j]) << "\n";
    }
  }
}

/**
 *  @brief Write matching lines.
*/
void summaryVis::writeMatches()
{
  // Write column names only if find at least one match.
  bool columnNamesWritten = false;

  const map<int, bool>::const_iterator probesetIDsEnd = m_ProbesetIDs.end();
  const map<int, bool>::const_iterator transcriptClusterIDsEnd = m_TranscriptClusterIDs.end();

  // Wiggle output must be sorted; use a set with functor as sort criterion.
  typedef set<summaryVisOutput*, outputSortCriterion> outputSet;
  outputSet wiggleOutput;

  // Read base file.
  while (m_BaseTsv.nextLevel (0) == TSV_OK)
  {
    // Skip if we're subsetting probeset_ids and this one is not
    // in the list.
    if (! m_ProbesetIDs.empty() && m_ProbesetIDs.find (m_BaseProbesetId) == probesetIDsEnd)
      continue;

    // Proceed if no extra files are present.
    bool foundMatch = true;

    // Look for matches in each extra file.
    for (unsigned int idx = 0; idx < m_FileCount - 1; ++idx)
    {
      TsvFile* extraTsv = m_ExtraTsvs[idx];
      if (extraTsv->findBegin (0, "probeset_id", TSV_OP_EQ, m_BaseProbesetId) != TSV_OK)
	Err::errAbort ("Problem reading extra file " + m_FileNames[idx + 1]);

      const int resultCount = extraTsv->findResultsCount();
      // Require a match in all summary files, if present.
      if (resultCount == 0)
      {
        foundMatch = false;
	break;
      }
      // Fatal error if more than one match was found.
      else if (resultCount > 1)
      {
	string msg = "FATAL: probeset_id '" + ToStr (m_BaseProbesetId);
	msg += "' is not a unique index. Duplicate probeset_id found, [";
	msg += ToStr (m_BaseProbesetId) + "] for file " + m_FileNames[idx + 1] + ".";
	Err::errAbort (msg);
      }
      // Found one match - read data
      if (extraTsv->findNext() != TSV_OK)
	Err::errAbort ("Problem reading extra file " + m_FileNames[idx + 1]);
      foundMatch = true;
    } // end for (unsigned int idx = 0; idx < m_FileCount - 1; ++idx)

    if (foundMatch)
    {
      // Found the probeset_id in all extra files, if any - get entry from the genome position file.
      if (m_GenomePosTsv.findBegin (0, "probeset_id", TSV_OP_EQ, m_BaseProbesetId) != TSV_OK)
	Err::errAbort ("Problem reading genome position file " + m_GenomePosFileName);
      const int resultCount = m_GenomePosTsv.findResultsCount();
      // Skip if no entry in the genome position file.
      if (resultCount == 0)
	continue;
      if (m_GenomePosTsv.findNext() != TSV_OK)
	Err::errAbort ("Problem reading genome position file " + m_GenomePosFileName);

      // Skip if we're subsetting by transcript cluster id and this one not in the list.
      if (! m_TranscriptClusterIDs.empty()
	  && m_TranscriptClusterIDs.find (m_GenomePosTranscriptClusterId) == transcriptClusterIDsEnd)
        continue;

      // If found at least one match, write column names and optional wiggle track definition line once.
      if (! columnNamesWritten)
      {
	writeColumnNames();
	if (m_Wiggle)
	  m_Out << "track type=wiggle_0 name=\"" << m_WiggleUseColName << "\"\n";
	columnNamesWritten = true;
      }
      // Write version dependent output.
      if (m_Wiggle)
      {
	// Wiggle output must be sorted: save in a set, output later.
	summaryVisOutput* wiggleLine = new summaryVisOutput (m_Seqname, m_Start, m_Stop, m_WiggleData);
	wiggleOutput.insert (wiggleLine);
	continue;
      }
      if (m_EgrVersion == 2)
      {
	m_Out << m_BaseProbesetId << "\t" << m_Seqname << "\t" << m_Start;
	m_Out << "\t" << m_Stop << "\t" << m_Strand;
      }
      else if (m_EgrVersion == 3)
	m_Out << m_BaseProbesetId;
      else if (m_EgrVersion == 1)
      {
	m_Out << m_Seqname << "\t" << m_Start;
	m_Out << "\t" << m_Stop << "\t" << m_Strand;
      }

      // Write columns from base file, in column sorted order.
      for (unsigned int i = 0; i < m_BaseDataColCount; ++i)
	m_Out << "\t" << m_BaseData [m_BaseOutputOrder [i]];

      // Write columns from extra files if any.
      for (unsigned int idx = 0; idx < m_FileCount - 1; ++idx)
      {
	const unsigned int extraCols = m_ExtraDataColCount[idx];
	for (unsigned int i = 0; i < extraCols; ++i)
	  m_Out << "\t" << m_ExtraData[idx] [m_ExtraOutputOrder [idx][i]];
      }
      m_Out << "\n";
    } // end if (foundMatch)
  } // end while (m_BaseTsv.nextLevel (0) == TSV_OK)

  // Write sorted wiggle output.
  if (m_Wiggle)
  {
    const outputSet::const_iterator wiggleOutputEnd = wiggleOutput.end();
    for (outputSet::const_iterator iter = wiggleOutput.begin(); iter != wiggleOutputEnd; ++iter)
    {
      m_Out << (*iter)->seqName << " " << (*iter)->start << " " << (*iter)->stop;
      m_Out << " " << (*iter)->data << "\n";
    }
    // Clean up.
    for (outputSet::iterator iter = wiggleOutput.begin(); iter != wiggleOutput.end(); ++iter)
      delete *iter;
  }
}
