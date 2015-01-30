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

/**
/// @file   CelExtract.cpp
/// @brief  Class for extracting probe level intensities from cel files.
*/

//
#include "chipstream/apt-cel-extract/CelExtract.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/DiskIntensityMart.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/SparseMart.h"
#include "file/TsvFile/ClfFile.h"
#include "util/Fs.h"

using namespace std;
using namespace affx;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;

/**
 * Constructor
 */
CelExtractOptions::CelExtractOptions()
{
  /* Some defaults... */
  verbosity = 1;
  force = false;
  pmOnly = false;
  pairWithBackground = false;
  ignoreProbesWithoutMM = false;
  subtractBackground = false;
  useDisk = true;
  diskCache = 50;
  m_output_precision=2;
}

/**
 * Destructor
 */
CelExtractOptions::~CelExtractOptions()
{
}

/**
 * Destructor.
 */
CelExtract::~CelExtract() {}

/**
 * Check options.
 */
void CelExtract::checkOptions(CelExtractOptions &o)
{
  // make sure we have cel file
  if (o.celFiles.size() < 1)
    Err::errAbort("Must provide at least 1 cel file.");

  // do we have chip layout info
  bool haveLayout = !(o.pgfFile == "" && o.cdfFile == "" && o.spfFile == "");

  // need library file if using probe id or probeset id file
  if (!haveLayout && (o.probesetIdsFiles.size() != 0 || o.probeIdsFiles.size() != 0))
    Err::errAbort("Must provide chip layout information (ie cdf, spf or pgf/clf file) when using probe or probeset list.");

  // user must specify outFile
  if (o.outFile == "")
    Err::errAbort("Must provide an output file.");

  // pgf needs clf
  if (o.pgfFile != "" && o.clfFile == "")
    Err::errAbort("Must provide a clf file when using pgf file.");

  // only pgf or cdf not both
  if (o.pgfFile != "" && o.cdfFile != "")
    Err::errAbort("Cannot provide both a CDF and a PGF file.");

  // cannot request pmOnly if no layout file
  if (o.pmOnly && !haveLayout)
    Err::errAbort("Must provide chip layout information (ie cdf, spf or pgf/clf file) when using PM only option.");

  // cannot request itnoreProbesWithoutMM if no layout file
  if (o.ignoreProbesWithoutMM && !haveLayout)
    Err::errAbort("Must provide chip layout information (ie cdf, spf or pgf/clf file) when using PM with MM only option.");

  // must have layout to do analysis
  if (o.analysisString != "" && !haveLayout)
    Err::errAbort("Must provide layout information (ie cdf, spf or pgf/clf file) when requesting an analysis.");

  // cannot request pairWithBackground if no layout and no analysis stream
  if ((o.pairWithBackground || o.subtractBackground) && !haveLayout)
    Err::errAbort("Must provide chip layout information (ie cdf, spf or pgf/clf file) when using one of the background options.");

  // cannot request pairWithBackground if no analysis
  if ((o.pairWithBackground || o.subtractBackground) && o.analysisString == "")
    Err::errAbort("Must provide an analysis string when using one of the background options.");

  // cannot request both options
  if (o.pairWithBackground && o.subtractBackground)
    Err::errAbort("You cannot request background values to both be reported and subtracted.");

  // have to have an analysis string to use with target-sketch-file
  if (o.sketchInFile.size() && ! o.analysisString.size()) {
    Err::errAbort("The --analysis option must be used when specifying an input target-sketch-file");
  }

  // chip type cel files all match
  if (o.force) {
    Verbose::out(1, "Force option given, skipping chip type checks.");
  } else {
    Verbose::out(1, "Checking chip type on all cel files.");
    vector<string> chipTypes;
    colrow_t numRows, numCols;
    int probeCount;
    int probeSetCount;

    if (o.cdfFile != "") {
      EngineUtil::getCdfChipType(chipTypes, numRows, numCols, probeCount, probeSetCount, o.cdfFile);
    } else if (o.pgfFile != "") {
      EngineUtil::getPgfChipType(chipTypes, numRows, numCols, probeCount, o.pgfFile, o.clfFile);
    } else if (o.spfFile != "") {
      EngineUtil::getSpfChipType(chipTypes, numRows, numCols, probeCount, probeSetCount, o.spfFile);
    } else {
      FusionCELData cel;
      cel.SetFileName(o.celFiles[0].c_str());
      if (!cel.ReadHeader())
        Err::errAbort("Cannot read cel file: " + o.celFiles[0]);
      chipTypes.push_back(StringUtils::ConvertWCSToMBS(cel.GetChipType()));
      if (chipTypes[0] == "")
        Err::errAbort("Unable to get chip type from cel file " + o.celFiles[0]);
      numCols = cel.GetCols();
      numRows = cel.GetRows();
      probeCount = numCols * numRows;
      cel.Close();
    }
    EngineUtil::checkCelChipTypes(chipTypes, probeCount, o.celFiles, numRows, numCols);
  }
}

/**
 * Constructor.
 */
CelExtract::CelExtract(CelExtractOptions &o)
{
  checkOptions(o);
  m_as = NULL;
  m_asFactory = NULL;
  m_FirstCelColIndex = 0;

  // make a copy of our options
  m_Options = o;

  // figure out basic stats from CEL file
  FusionCELData cel;
  cel.SetFileName(o.celFiles[0].c_str());
  if (!cel.ReadHeader()) {
    Err::errAbort("Cannot read cel file: " + o.celFiles[0]);
  }
  m_NumCols = cel.GetCols();
  m_NumRows = cel.GetRows();
  m_NumFeatures = m_NumCols * m_NumRows;
  m_NumChannels = cel.GetChannels().size();
  if (m_NumChannels == 0) {
    // not a multi-channel CEL file; treat as a single-channel CEL file
    m_NumChannels = 1;
  }
  cel.Close();

  m_OutFiles.resize(m_NumChannels, o.outFile);
  if (m_NumChannels > 1) {
    for (int channelIx = 0; channelIx < m_NumChannels; channelIx++) {
      m_OutFiles[channelIx].append(".channel" + ToStr(channelIx));
    }
  }

  Verbose::out(1, "Expecting " + ToStr(m_NumCols) + " columns and " +
               ToStr(m_NumRows) + " rows in each CEL file.");

  /// what probesets to load, we leave this empty == load all
  /// or fill in to restrict if probeset list file given
  set<const char *, Util::ltstr> probeSetsToLoad;

  for (int i = 0; i < m_NumFeatures; i++)
    m_ProbesToDump.push_back(false);

  // set from probe file
  if (o.probeIdsFiles.size() != 0) {
    int count = 0;
    vector<probeid_t> probesToExtract;
    for (int i = 0;i < o.probeIdsFiles.size();i++) {
      Verbose::out(1, "Loading probes to extract from: " + o.probeIdsFiles[i]);
      EngineUtil::readProbeFile(probesToExtract, o.probeIdsFiles[i]);
    }
    for (int i = 0; i < probesToExtract.size(); i++) {
      if (probesToExtract[i] >= m_NumFeatures)
        Err::errAbort("Invalid probe ID requested, " + ToStr(probesToExtract[i]));
      if (probesToExtract[i] < 0)
        Err::errAbort("Invalid probe ID requested, " + ToStr(probesToExtract[i]));
      m_ProbesToDump[probesToExtract[i]] = true;
      count++;
    }
    Verbose::out(1, "Loaded " + ToStr(count) + " probes.");
  }
  // set from probeset file
  vector<string> psVec;
  if (o.probesetIdsFiles.size() != 0) {
    for (int i = 0;i < o.probesetIdsFiles.size();i++) {
      Verbose::out(1, "Loading probesets to extract from: " + o.probesetIdsFiles[i]);
      EngineUtil::readProbeSetFile(psVec, o.probesetIdsFiles[i]);
    }
    for (int i = 0; i < psVec.size(); i++)
      probeSetsToLoad.insert(psVec[i].c_str());
    Verbose::out(1, "Will load " + ToStr(probeSetsToLoad.size()) + " probesets.");
  }
  // set default list of probe IDs from cel file
  if ((o.probeIdsFiles.size() == 0) && (o.probesetIdsFiles.size()) == 0) {
    Verbose::out(1, "Defaulting to extract all probes.");
    for (int i = 0; i < m_NumFeatures; i++)
      m_ProbesToDump[i] = true;
  }

  // load up layout
  m_Layout = NULL;
  /// what probes to ignore, we leave this empty == none to kill
  probeidmap_t killList;
  /// what probeset types (ie genotype) to load, we leave this empty == load all
  set<affxcdf::GeneChipProbeSetType> psTypesToLoad;

  // make sure gc control probes get loaded, if requested
  std::vector<bool> probesToLoad = m_ProbesToDump;
  if (o.bgpFile != "") {
    int count = 0;
    Verbose::out(1, "Loading probes to extract from " + o.bgpFile);
    vector<probeid_t> bgpProbes;
    EngineUtil::readProbeFile(bgpProbes, o.bgpFile.c_str());
    for (int i = 0; i < bgpProbes.size(); i++) {
      if (bgpProbes[i] >= m_NumFeatures)
        Err::errAbort("Invalid probe ID requested, " + ToStr(bgpProbes[i]));
      if (bgpProbes[i] < 0)
        Err::errAbort("Invalid probe ID requested, " + ToStr(bgpProbes[i]));
      count++;
      probesToLoad[bgpProbes[i]] = true;
      m_ProbeOrder.push_back(bgpProbes[i]);
    }
    Verbose::out(1, "Loaded " + ToStr(count) + " bgp probes.");
  }

  if (o.cdfFile != "" || o.pgfFile != "" || o.spfFile != "") {
    m_Layout = new ChipLayout;
    m_Layout->setNeedGc(true);
    m_Layout->setNeedMismatch(true);
  }

  if (o.cdfFile != "") {
    Verbose::out(1, "Opening layout file: " + o.cdfFile);
    if (!m_Layout->openCdf(o.cdfFile, probeSetsToLoad, NULL, probesToLoad, "", killList, false, psTypesToLoad))
      Err::errAbort("Could not open layout file " + o.cdfFile);
    Verbose::out(1, "Done loading chip layout info.");
  } else if (o.pgfFile != "") {
    Verbose::out(1, "Opening layout file: " + o.clfFile);
    affx::ClfFile clf;
    if (!clf.open(o.clfFile))
      Err::errAbort("Could not open layout file " + o.clfFile);
    /* Make sure our assumptions about CLF file are true */
    Err::check(clf.getSequential() == 1,
               "ProbesetSummarizeEngine::loadPgfLayout() - unable to handle clf file without sequential set to 1.");
    Err::check(clf.getOrder().compare("col_major") == 0 || clf.getOrder().compare("row_major") == 0,
               "Unable to handle clf file without order set to row_major (old mislabeled 'col_major' accepted due to earlier bug.)");

    Verbose::out(1, "Opening layout file: " + o.pgfFile);
    if (!m_Layout->openPgf(o.pgfFile.c_str(), clf.getXMax() + 1, clf.getYMax() + 1, probeSetsToLoad, NULL, NULL, probesToLoad, "", killList, false, false))
      Err::errAbort("Could not open layout file " + o.pgfFile);
    Verbose::out(1, "Done loading chip layout info.");
  } else if (o.spfFile != "") {
    Verbose::out(1, "Opening layout file: " + o.spfFile);
    m_Layout->openSpf(o.spfFile, probeSetsToLoad, NULL, probesToLoad, "", false, psTypesToLoad);
    Verbose::out(1, "Done loading chip layout info.");
  } else {
    Verbose::out(1, "No chip layout information available. Doing basic extraction.");
  }

  // now we need to flag the probes in the requested probesets
  // to be dumped
  if (o.probesetIdsFiles.size() != 0) {
    Verbose::out(1, "Flagging probes from requested probesets to be extracted.");
    assert(m_Layout);
    set<const char *, Util::ltstr>::iterator iter;
    for (iter = probeSetsToLoad.begin(); iter != probeSetsToLoad.end(); iter++) {
      ProbeListPacked pList = m_Layout->getProbeListByName(*iter);
      if (pList.isNull()) {
        Verbose::out(1, "Unable to find probeset " + ToStr(*iter) + " in layout file. Skipping.");
      } else {
        ProbeSet *ps = ProbeListFactory::asProbeSet(pList);
        //printf("m_PsOrder_jhg:1: %d\n",m_Layout->getProbeSetIndexByName(*iter));
        m_PsOrder_jhg.push_back(m_Layout->getProbeSetIndexByName(*iter));
        for (int aIx = 0; aIx < ps->atoms.size(); aIx++)
          for (int plIx = 0; plIx < ps->atoms[aIx]->probes.size(); plIx++) {
            m_ProbesToDump[ps->atoms[aIx]->probes[plIx]->id] = true;
            m_ProbeOrder.push_back(ps->atoms[aIx]->probes[plIx]->id);
          }
        Freez(ps);
      }
    }
  } else if (m_Layout != NULL) {
    // iterate over everything in the layout
    for (int psIx = 0; psIx < m_Layout->getProbeSetCount(); psIx++) {
      ProbeListPacked pList = m_Layout->getProbeListAtIndex(psIx);
      ProbeSet *ps = ProbeListFactory::asProbeSet(pList);
      //printf("m_PsOrder_jhg:2: %d\n",m_Layout->getProbeSetIndexByName(ps->name));
      m_PsOrder_jhg.push_back(m_Layout->getProbeSetIndexByName(ps->name));
      for (int aIx = 0; aIx < ps->atoms.size(); aIx++)
        for (int plIx = 0; plIx < ps->atoms[aIx]->probes.size(); plIx++)
          m_ProbeOrder.push_back(ps->atoms[aIx]->probes[plIx]->id);
      Freez(ps);
    }
  }

  // does the user only want a certain type of probe
  if (o.pmOnly || o.ignoreProbesWithoutMM) {
    if (o.ignoreProbesWithoutMM)
      Verbose::out(1, "Figuring out which probes are PM with MM to restrict extraction.");
    else if (o.pmOnly)
      Verbose::out(1, "Figuring out which probes are PM to restrict extraction.");
    assert(m_Layout);
    // swap probes to dump - we need to know what is true so far
    vector<bool> originalProbesToDump(m_NumFeatures, false);
    m_ProbesToDump.swap(originalProbesToDump);
    // iterate over everything in the layout
    for (int psIx = 0; psIx < m_Layout->getProbeSetCount(); psIx++) {
      ProbeListPacked pList = m_Layout->getProbeListAtIndex(psIx);
      ProbeSet *ps = ProbeListFactory::asProbeSet(pList);
      for (int aIx = 0; aIx < ps->atoms.size(); aIx++) {
        for (int plIx = 0; plIx < ps->atoms[aIx]->probes.size(); plIx++) {
          if (o.ignoreProbesWithoutMM) {
            // if has MM and is PM and was true before, set true
            if (ps->hasMM() && Probe::isPm(*(ps->atoms[aIx]->probes[plIx]))) {
              if (originalProbesToDump[ps->atoms[aIx]->probes[plIx]->id] == true)
                m_ProbesToDump[ps->atoms[aIx]->probes[plIx]->id] = true;
            }
          } else if (o.pmOnly) {
            // if is PM and was true before, set true
            if (Probe::isPm(*(ps->atoms[aIx]->probes[plIx]))) {
              if (originalProbesToDump[ps->atoms[aIx]->probes[plIx]->id] == true)
                m_ProbesToDump[ps->atoms[aIx]->probes[plIx]->id] = true;
            }
          }
        }
      }
      Freez(ps);
    }
  }
  Verbose::out(1, "Done with extraction setup.");
}

void CelExtract::defineExtractTsv(std::vector<affx::TsvFile*>& tsv, bool bDoSaturationReport/*=false*/)
{

  for (int f = 0; f < tsv.size(); f++) {
    tsv[f] = new affx::TsvFile();
    int columnCount = 0;

    // add some nice header info
    tsv[f]->addHeader("program", "apt-cel-extract");
    tsv[f]->addHeader("command-line",m_Options.commandLine);
    tsv[f]->addHeader("exec-guid",m_Options.execGuid);

    //
    if (m_Layout != NULL) {
      tsv[f]->addHeader("num-cols", int(m_Layout->m_XCount));
      tsv[f]->addHeader("num-rows", int(m_Layout->m_YCount));
    }
    else {
      tsv[f]->addHeader("num-features", m_NumFeatures);
      tsv[f]->addHeader("num-cols", m_NumCols);
      tsv[f]->addHeader("num-rows", m_NumRows);
    }

    // full input path.
    for (int i = 0; i < m_Options.celFiles.size(); i++) {
      tsv[f]->addHeader("input-" + ToStr(i), m_Options.celFiles[i]);
    }

    // primary key
	if (bDoSaturationReport == false)
	{
		tsv[f]->defineColumn(0, columnCount++, "probe_id");
		tsv[f]->defineColumn(0, columnCount++, "x");
		tsv[f]->defineColumn(0, columnCount++, "y");
	}
	else
	{
		tsv[f]->defineColumn (0, columnCount++, "fraction_of_features_saturated", TSV_TYPE_FLOAT);
	}
    // set up probe annotations column labels
    if (m_Layout != NULL) {
      tsv[f]->defineColumn(0, columnCount++, "probe_type");
      tsv[f]->defineColumn(0, columnCount++, "probeset_id");
      tsv[f]->defineColumn(0, columnCount++, "probeset_type");
      tsv[f]->defineColumn(0, columnCount++, "block");
    }

    m_FirstCelColIndex = columnCount;

    // column for each cel file
    for (int i = 0; i < m_Options.celFiles.size(); i++) {
      int cidx=columnCount++;
      tsv[f]->defineColumn(0, cidx, Fs::basename(m_Options.celFiles[i]));
      tsv[f]->setPrecision(0, cidx, m_Options.m_output_precision);
    }
    tsv[f]->writeTsv_v1(m_OutFiles[f]);
  }
}

/**
 * Do the extraction
 */
void CelExtract::extract()
{
	bool bDoSaturationReport = false;
	if (strstr (m_Options.outFile.c_str (),"SaturationReport"))
		bDoSaturationReport = true;
	
  // create analysis stream if requested
  if (m_Options.analysisString != "")
    m_as = newAnalysisStreamObject();

  // cache for probe intensity data
  IntensityMart *iMart = NULL;
  string diskDir = m_Options.diskDir;
  if (m_Options.diskDir == "") {
    diskDir = Fs::join(".","temp");
  }
  if (m_Layout != NULL) {
    if (m_Options.useDisk == true) {
      DiskIntensityMart* diskMart = 
        new DiskIntensityMart(m_ProbeOrder, 
                              m_Options.celFiles, 
                              m_Options.diskCache * 1048576,
                              diskDir, 
                              "apt-extract.tmp",
                              true);
      iMart = diskMart;
    } else {
      SparseMart* sparseMart = new SparseMart(m_ProbeOrder, 
                                              m_Options.celFiles, true);
      iMart = sparseMart;
    }
  } else {
    vector<int> order;
    for (int pIx = 0; pIx < m_ProbesToDump.size(); pIx++)
      order.push_back(pIx);
    iMart = new DiskIntensityMart(order, m_Options.celFiles, m_Options.diskCache * 1048576, diskDir, "apt.extract.tmp", true);
  }

  // Setup the output file
  std::vector<affx::TsvFile*> tsv(m_NumChannels);
  defineExtractTsv(tsv,bDoSaturationReport);

  // Iterate over probesets
  if (m_Layout != NULL) {
    if (m_PsOrder_jhg.size() == 0) {
      for (int psIx = 0; psIx < m_Layout->getProbeSetCount(); psIx++) {
        //printf("m_PsOrder_jhg:3: %d\n",psIx);
        m_PsOrder_jhg.push_back(psIx);
      }
    }
    doProbeSet(iMart, tsv);
    // Iterate over probes
  } else {
    doProbe(iMart, tsv, bDoSaturationReport);
  }

  // close analysis stream
  if (m_as != NULL)
    m_as->finish();

  // close and free output tsv file
  for (int f = 0; f < tsv.size(); f++) {
    tsv[f]->close();
    Freez(tsv[f]);
  }

  // free up memory
  Freez(iMart);
  Freez(m_as);
  Freez(m_asFactory);
  Freez(m_Layout);

  /* Remove the temp dir if it is empty */
  if (Fs::isWriteableDir(diskDir))
    Fs::rmdir(diskDir, false);

}

void setupAnalysisStreamFactory(AnalysisStreamFactory &aFactory, CelExtractOptions &opts)
{
  if (opts.sketchInFile != "") {
    Verbose::out(1, "Opening target normalization file: " + Fs::basename(opts.sketchInFile.c_str()));
    aFactory.readTargetSketchFromFile(opts.sketchInFile.c_str());
  }
  if (opts.bgpFile != "") {
    Verbose::out(1, "Opening bgp file: " + Fs::basename(opts.bgpFile.c_str()));
    aFactory.readControlProbes(opts.bgpFile.c_str());
  }
}

AnalysisStreamExpression *CelExtract::newAnalysisStreamObject()
{

  // Construct the analysis stream factory
  AnalysisStreamExpression *as = NULL;
  map<string, string> stdMethods; // there are none for this app
  string emptyString;
  AnalysisStreamFactory *asFactory = new AnalysisStreamFactory(QuantMethodFactory::Expression);

  setupAnalysisStreamFactory(*asFactory, m_Options);

  // Construct the analysis stream
  string analysisString = m_Options.analysisString + ",med-polish"; // add quant method, although we won't use it
  as = asFactory->constructExpressionAnalysisStream(analysisString.c_str(), *m_Layout, stdMethods, emptyString);

  // get list control probes
  vector<Probe *> gcProbes = asFactory->getGcControlProbes();
  for (int i = 0; i < gcProbes.size(); i++) {
    m_gcProbes.push_back(gcProbes[i]->id);
  }

  // We cannot blow away the factory yet as it owns some of the
  // memory used in the analysis stream. At a minimum the probes
  // instances for the GC background probes are owned by the
  // factory.
  m_asFactory = asFactory;

  return as;
}

void CelExtract::doProbeSet(IntensityMart *iMart, std::vector<affx::TsvFile*>& tsv)
{
  assert(m_Layout);

  Verbose::out(2, "Reading cel files.");
  vector<string> fileCopy = m_Options.celFiles;
  CelReader cReader;
  cReader.registerIntensityMart(iMart);
  if (m_as != NULL)
    m_as->registerChipStreamObjs(cReader);
  cReader.setFiles(fileCopy);
  cReader.readFiles();

  // dump out intensities
  unsigned int dotMod = max((int)m_PsOrder_jhg.size() / 40, 1);
  Verbose::progressBegin(1, "Dumping intensities", 40, dotMod, m_PsOrder_jhg.size());
  for (int index = 0; index < m_PsOrder_jhg.size(); index++) {
    Verbose::progressStep(1);

    int psIx = m_PsOrder_jhg[index];
    ProbeListPacked pList = m_Layout->getProbeListAtIndex(psIx);
    ProbeSet *ps = ProbeListFactory::asProbeSet(pList);

    for (int aIx = 0; aIx < ps->atoms.size(); aIx++) {
      unsigned int channelIx = ps->atoms[aIx]->getChannelCode();
      for (int plIx = 0; plIx < ps->atoms[aIx]->probes.size(); plIx++) {
        if (m_ProbesToDump[ps->atoms[aIx]->probes[plIx]->id]) {
          // set probe ID
          int columnCount = 0;
          tsv[channelIx]->set(0, columnCount++, ps->atoms[aIx]->probes[plIx]->id + 1);

          // set probe x/y
          int x = -1, y = -1;
          m_Layout->indexToXY(ps->atoms[aIx]->probes[plIx]->id, x, y);
          tsv[channelIx]->set(0, columnCount++, x);
          tsv[channelIx]->set(0, columnCount++, y);

          tsvSetProbeAnnotations(tsv[channelIx], columnCount, ps, aIx, plIx);

          // set probe intensities
          for (int cIx = 0; cIx < m_Options.celFiles.size(); cIx++) {
            if (m_as != NULL) {
              float pm = - 1;
              pm = QuantMethod::transformPrimaryData(ps->atoms[aIx]->probes[plIx]->id, cIx, *iMart, *(m_as->getChipStream()), channelIx);
              if (m_Options.pairWithBackground || m_Options.subtractBackground) {
                float mm = -1;
                m_as->getPmAdjuster()->pmAdjustment(ps->atoms[aIx]->probes[plIx]->id, cIx, *iMart, *(m_as->getChipStream()), pm, mm);
                if (m_Options.subtractBackground) {
                  tsv[channelIx]->set(0, columnCount++, pm - mm);
                } else {
                  tsv[channelIx]->set(0, columnCount++, ToStr((int)pm) + "," + ToStr((int)mm));
                }
              } else {
                tsv[channelIx]->set(0, columnCount++, pm);
              }
            } else {
              tsv[channelIx]->set(0, columnCount++,
                                  iMart->getProbeIntensity(ps->atoms[aIx]->probes[plIx]->id, cIx, channelIx));
            }
          }
          // write
          tsv[channelIx]->writeLevel(0);
        } // end if probe to dump
      } // end for each probe
    } // end for each atom
    Freez(ps);
  } // end for each probeset
  Verbose::progressEnd(1, "Done.");
}

void CelExtract::doProbe(IntensityMart *iMart, std::vector<affx::TsvFile*>& tsv, bool bDoSaturationReport/*=false*/)
{
	Verbose::out(2, "Reading cel files.");
	vector<string> fileCopy = m_Options.celFiles;
	CelReader cReader;
	cReader.registerIntensityMart(iMart);

	if (m_as != NULL)
		m_as->registerChipStreamObjs(cReader);

	cReader.setFiles(fileCopy);
	cReader.readFiles();

	if (bDoSaturationReport == false)
	{
		// dump out intensities
		unsigned int dotMod = max((int)m_ProbesToDump.size() / 40, 1);
		Verbose::progressBegin(1, "Dumping intensities", 40, dotMod, m_ProbesToDump.size());
		for (int channelIx = 0; channelIx < m_NumChannels; channelIx++) 
		{
			for (int pIx = 0; pIx < m_ProbesToDump.size(); pIx++) 
			{
				Verbose::progressStep(1);

				if (!m_ProbesToDump[pIx])
				break;

				// set probe ID
				int columnCount = 0;
				tsv[channelIx]->set(0, columnCount++, pIx + 1);

				// set probe x/y
				unsigned int x = pIx % (m_NumCols);
				unsigned int y = pIx / (m_NumCols);
				tsv[channelIx]->set(0, columnCount++, x);
				tsv[channelIx]->set(0, columnCount++, y);

				// set probe intensities
				for (int cIx = 0; cIx < m_Options.celFiles.size(); cIx++) 
				{
					tsv[channelIx]->set(0,
								columnCount++,
								iMart->getProbeIntensity(pIx, cIx, channelIx));
				}

				if (bDoSaturationReport == false)
					tsv[channelIx]->writeLevel(0);
			}
		}
	} 
	else
	{  // bDoSaturationReport == true

		unsigned int dotMod = max((int)m_ProbesToDump.size() / 40, 1);
		Verbose::progressBegin(1, "Computing saturation", 40, dotMod, m_ProbesToDump.size());

		for (int ifile = 0; ifile < m_Options.celFiles.size(); ifile++)
		{
			for (int ichannel = 0; ichannel < m_NumChannels; ichannel++) 
			{
				int ntotal = 0;
				int nsat = 0;

				for (int ifeature = 0; ifeature < m_ProbesToDump.size(); ifeature++) 
				{
					if (iMart->getProbeIntensity (ifeature, ifile, ichannel) > 3800)
						nsat++;

					ntotal++;
				}

				tsv[ichannel]->set(0,ifile+1,(float)nsat/(float)ntotal);
			}
		}

		for (int ichannel = 0; ichannel < m_NumChannels; ichannel++)
		{
			tsv[ichannel]->writeLevel (0);
		}
	}
 
	Verbose::progressEnd(1, "Done.");
}

void CelExtract::tsvSetProbeAnnotations(affx::TsvFile* tsv, int &columnCount, ProbeSet *ps, int aIx, int plIx)
{
  if (ps != NULL) {
    Probe *p = ps->atoms[aIx]->probes[plIx];

    int block = -1;
    int atomSum = 0;
    for (int gIx = 0; gIx < ps->numGroups; gIx++) {
      int atomsInGroup = ps->atomsPerGroup[gIx];
      if (aIx >= atomSum && aIx < atomSum + atomsInGroup) {
        block = gIx;
        break;
      }
      atomSum += atomsInGroup;
    }

    // set values -- we do not use Probe::stringForType
    // as the types are screwed up due to the ProbeList refactor.
    // ie every probe is ST (ie PMST) even if it is an AT probe.
    string type = "unknown";
    if (Probe::isPm(*p))
      type = "pm";
    else if (Probe::isMm(*p))
      type = "mm";
    tsv->set(0, columnCount++, type);
    tsv->set(0, columnCount++, ps->name);
    // detect bogus genotype probesets and change them to expression
    if (ps->psType == ProbeSet::GenoType && ps->numGroups != 2 && ps->numGroups != 4) {
      // there is a bug in some cdf files where controls are marked as genotyping...
      Verbose::out(4, "Labelling probeset " + ToStr(ps->name) + " as expression rather than genotyping because it does not have 2 or 4 groups.");
      tsv->set(0, columnCount++, ProbeSet::typeEnumToString(ProbeSet::Expression));
    } else {
      tsv->set(0, columnCount++, ProbeSet::typeEnumToString(ps->psType));
    }
    tsv->set(0, columnCount++, block);
  } else {
    // set null values
    for (int i = 0; i < m_FirstCelColIndex; i++)
      tsv->set(0, columnCount++, "NULL");
  }
}
