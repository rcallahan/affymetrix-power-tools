////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License
// (version 2.1) as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

//
#include "calvin_files/fusion/src/FusionCHPLegacyData.h"
#include "chipstream/EngineUtil.h"
#include "file/TsvFile/TsvFile.h"
#include "util/AptVersionInfo.h"
#include "util/AffxByteArray.h"
#include "util/CalvinToText.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/Verbose.h"


struct opt {
  std::string help;
  int verbose;
  std::string progName;
  std::string version;
  std::string execGuid;
  std::string timeStr;
  std::string commandLine;
  std::string outDir;
  bool no_header;
  bool no_body;
  std::string cdfFile;
  vector<std::string> chp_files;
};

string getOutFileName(string chpName, struct opt *o);
void run(opt& o);
bool doLegacy(opt& o, int chp_index);
void doLegacyExpression(opt& o, affx::TsvFile& chpTsv, affymetrix_fusion_io::FusionCHPLegacyData* data);
void doLegacyGenotyping(opt& o, affx::TsvFile& chpTsv, affymetrix_fusion_io::FusionCHPLegacyData* data);
void doLegacyUniversal(opt& o, affx::TsvFile& chpTsv, affymetrix_fusion_io::FusionCHPLegacyData* data);
void doLegacyReseq(opt& o, affx::TsvFile& chpTsv, affymetrix_fusion_io::FusionCHPLegacyData* data);

void define_apt_chp_to_txt_options(PgOptions* opts)
{
  ///@todo add force option
  opts->setUsage(
    "apt-chp-to-txt - Writes copynumber CHP files in tab-delimited text format.\n\n"
    "usage:\n"
    "	apt-chp-to-txt  \\\n"
    "       --o results --chp-files chpfiles.txt");
  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "This message.",
                     "false");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
  opts->defineOption("v", "verbose", PgOpt::INT_OPT,
                     "How verbose to be with status messages 0 = quiet, 1 - usual messages, 2 - more messages.",
                     "1");
  opts->defineOption("o", "out-dir", PgOpt::STRING_OPT,
                     "Directory to write result files into.",
                     ".");
  opts->defineOption("", "cdf-file", PgOpt::STRING_OPT,
                     "CDF file to pull probe set names from.",
                     "");
  opts->defineOption("", "no-header", PgOpt::BOOL_OPT, "If set to true, then the file header will not be output.", "false");
  opts->defineOption("", "no-body", PgOpt::BOOL_OPT, "If set to true, then the file body will not be output.", "false");
  opts->defineOption("", "chp-files", PgOpt::STRING_OPT,
                     "Text file specifying chp files to process, one per line with the first line being 'chp_files'.",
                     "");
  opts->defineOption("", "log-file", PgOpt::STRING_OPT,
                     "The name of the log file. Generally defaults to the program name in the out-dir folder.",
                     "");
}

/**
 * Print out a little ditty for the given type of self documenting object.
 *
 * @param query - type of analysis/adjustment/summarization information is desired for.
 */

void fillInOptions(PgOptions* pgopts, opt& opt)
{
  // Standard stuff...
  opt.progName = pgopts->getProgName();
  opt.progName = Fs::basename(opt.progName);
  opt.timeStr = Util::getTimeStamp();
  if (pgopts->getBool("help") || pgopts->argc() == 1) {
    pgopts->usage();
    cout << endl << "version: " << opt.version << endl;
    exit(0);
  }
  if (pgopts->getBool("version")) {
    Verbose::out(0, "version: " + opt.version); exit(0);
  }
  opt.verbose = pgopts->getInt("verbose");
  Verbose::setLevel(opt.verbose);
  opt.outDir = pgopts->get("out-dir");

  opt.cdfFile = pgopts->get("cdf-file");
  opt.no_header = pgopts->getBool("no-header");
  opt.no_body = pgopts->getBool("no-body");
  if ((opt.no_header) && (opt.no_body)) {
    Err::errAbort("If you specify both --no-header and --no-body, then there will be no output.");
  }

  // Read in chp file list from other file if specified.
  std::string chpFiles_filename = pgopts->get("chp-files");
  if (chpFiles_filename != "") {
    affx::TsvFile::extractColToVec(chpFiles_filename, "chp_files", &opt.chp_files, 1);
    //
    Verbose::out(1, "Read " + ToStr(opt.chp_files.size()) + " chp files from: "
                 "'" + Fs::basename(chpFiles_filename) + "'");
  }
  // The rest of the args are celfiles.
  else {
    opt.chp_files = pgopts->getArgVector();
  }
}

int main(int argc, char* argv[])
{
  try {
    struct opt opt;
    PgOptions* pgopts = new PgOptions;

    define_apt_chp_to_txt_options(pgopts);
    pgopts->parseArgv(argv);

    opt.version = AptVersionInfo::versionToReport();
    opt.execGuid = affxutil::Guid::GenerateNewGuid();
    opt.commandLine = pgopts->commandLine();

    fillInOptions(pgopts, opt);

    time_t startTime = time(NULL);
    ofstream ostrmLog;

    Fs::ensureWriteableDirPath(opt.outDir);

    string strLogFileName;
    if (pgopts->get("log-file") != "") {
      strLogFileName = pgopts->get("log-file");
    } else {
      strLogFileName = Fs::join(opt.outDir, "apt-chp-to-txt.log");
    }

    Fs::mustOpenToWrite(ostrmLog, strLogFileName);
    LogStream log(3, &ostrmLog);
    Verbose::pushMsgHandler(&log);
    Verbose::pushProgressHandler(&log);
    Verbose::pushWarnHandler(&log);
    Verbose::out(3, "version=" + opt.version);

#ifndef WIN32
    char name[8192];
    if (gethostname(name, ArraySize(name)) == 0) {
      Verbose::out(3, "host=" + ToStr(name));
    } else {
      Verbose::out(3, "host=unknown");
    }
    if (getcwd(name, ArraySize(name)) != NULL) {
      Verbose::out(3, "cwd=" + ToStr(name));
    } else {
      Verbose::out(3, "cwd=unknown");
    }
#endif /* WIN32 */

    Verbose::out(3, "command-line=" + pgopts->commandLine());
    delete pgopts; pgopts = NULL;

    run(opt);

    time_t endTime = time(NULL);
    int t = int((float)(endTime - startTime) / 60.0 * 100);   // convert to minutes
    Verbose::out(1, ToStr("Run took approximately: ") + ToStr((float)t / 100) + ToStr(" minutes."));
    Verbose::out(1, "Closing log.");
    ostrmLog.close();
    return 0;
  } catch (...) {
    Verbose::out(1, "Unexpected Error: uncaught exception.");
    return 1;
  }
  return 1;
}

bool isCalvinFile(const AffxString& strFileName)
{
  AffxBinaryFile file;
  if (!file.open(strFileName, AffxBinaryFile::LOAD)) {
    return false;
  }
  unsigned char ucFileType = file.readUnsignedChar();
  file.close();
  if (ucFileType == 59) {
    return true;
  } else {
    return false;
  }
}

void run(opt& o)
{
  for (unsigned int chp_index = 0; chp_index < o.chp_files.size(); chp_index++) {
    Verbose::out(1, "Processing file: " + o.chp_files[chp_index]);
    try {
      bool bSuccess = false;
      if (isCalvinFile(o.chp_files[chp_index])) {
        AffxByteArray chpFileName;
        chpFileName.assign(o.chp_files[chp_index]);
        CalvinToText converter;
        std::string outFile=Fs::join(o.outDir,Fs::basename(chpFileName.toString())+".txt");
        converter.run(o.chp_files[chp_index], 
                      outFile,
                      !o.no_header, 
                      !o.no_body, 
                      false);
      } else {
        bSuccess = doLegacy(o, chp_index);
      }
    } catch (...) {
      // Ios exception messages aren't very informative - just write
      Err::errAbort("A problem occurred while processing file (Exception thrown) " + o.chp_files[chp_index]);
    }
  }
}

bool readCdfFile(const AffxString& strFileName, bool bMapping, std::vector<std::string>& vProbeSetNames)
{
  vProbeSetNames.clear();
  affxcdf::CCDFFileData cdf;
  cdf.SetFileName(strFileName.c_str());
  if (cdf.Read()) {
    int iUnitCount = cdf.GetHeader().GetNumProbeSets();
    vProbeSetNames.resize(iUnitCount);
    for (int iUnitIndex = 0; (iUnitIndex < iUnitCount); iUnitIndex++) {
      vProbeSetNames[iUnitIndex] = cdf.GetProbeSetName(iUnitIndex);
    }
  } else {
    Err::errAbort("Cannot read CDF file: " + strFileName);
  }
  return true;
}


bool doLegacy(opt& o, int chp_index)
{
  bool bSuccessful = false;
  try {
    affymetrix_fusion_io::FusionCHPData* chp = NULL;
    affymetrix_fusion_io::FusionCHPLegacyData* data = NULL;
    chp = affymetrix_fusion_io::FusionCHPDataReg::Read(o.chp_files[chp_index]);
    if (chp == NULL)
      Err::errAbort("Problem reading file '" + o.chp_files[chp_index] + "'");
    else {
      data = affymetrix_fusion_io::FusionCHPLegacyData::FromBase(chp);
      if (data != NULL) {
        if (o.cdfFile != "") {
          colrow_t rowCount, colCount;
          int probeCount;
          int psCount;
          vector<string> chipTypes;
          EngineUtil::getCdfChipType(chipTypes, rowCount, colCount, probeCount, psCount, o.cdfFile);
          string observedChipType = StringUtils::ConvertWCSToMBS(data->GetHeader().GetChipType());
          if (!EngineUtil::checkChipTypeOk(observedChipType, chipTypes)) {
            string chipTypeString = "";
            for (int i = 0; i < chipTypes.size(); i++) {
              chipTypeString += chipTypes[i];
              if (i + 1 != chipTypes.size())
                chipTypeString += ", ";
            }
            Verbose::out(1, "CHP File: " + o.chp_files[chp_index] + " (" + observedChipType + ")");
            Verbose::out(1, "CDF File: " + o.cdfFile + " (" + chipTypeString + ")");
            Err::errAbort("CDF chip type does not match CHP file.");
          }
        }

        AffxString chpFileName = (o.outDir != "" ? o.outDir + "/" : std::string("")) + Fs::basename(o.chp_files[chp_index]) + std::string(".txt");
        affx::TsvFile chpTsv;
        Verbose::out(1, "Writing file: " + chpFileName);

        if (!o.no_header) {
          ParameterNameValueType param;
          chpTsv.addHeader("Magic", ::getInt(data->GetHeader().GetMagic()));
          chpTsv.addHeader("ChipType", StringUtils::ConvertWCSToMBS(data->GetHeader().GetChipType()));
          chpTsv.addHeader("AlgName", StringUtils::ConvertWCSToMBS(data->GetHeader().GetAlgName()));
          chpTsv.addHeader("AlgVersion", StringUtils::ConvertWCSToMBS(data->GetHeader().GetAlgVersion()));
          chpTsv.addHeader("ProgID" , StringUtils::ConvertWCSToMBS(data->GetHeader().GetProgID()));
          chpTsv.addHeader("ParentCellFile", StringUtils::ConvertWCSToMBS(data->GetHeader().GetParentCellFile()));
          chpTsv.addHeader("Cols" , ::getInt(data->GetHeader().GetCols()));
          chpTsv.addHeader("Cols" , ::getInt(data->GetHeader().GetCols()));
          chpTsv.addHeader("NumProbeSets" , ::getInt(data->GetHeader().GetNumProbeSets()));
          affymetrix_fusion_io::FusionTagValuePairTypeList list;
          data->GetHeader().AlgorithmParameters(list);
          affymetrix_fusion_io::FusionTagValuePairTypeList::iterator iter;
          for (iter = list.begin(); iter != list.end(); ++iter) {
            chpTsv.addHeader("affymetrix-algorithm-param-" + StringUtils::ConvertWCSToMBS(iter->Tag) , StringUtils::ConvertWCSToMBS(iter->Value));
          }
          data->GetHeader().SummaryParameters(list);
          for (iter = list.begin(); iter != list.end(); ++iter) {
            chpTsv.addHeader("affymetrix-chipsummary-" + StringUtils::ConvertWCSToMBS(iter->Tag), StringUtils::ConvertWCSToMBS(iter->Value));
          }
          affxchp::BackgroundZoneInfo info;
          data->GetHeader().GetBackgroundZoneInfo(info);
          chpTsv.addHeader("NumberBackgroundZones" , ::getInt(info.number_zones));
          chpTsv.addHeader("SmoothFactor" , ::getDouble(info.smooth_factor, 6));
          affxchp::BackgroundZoneTypeList::iterator iterBackground;
          int iIndex = 1;
          for (iterBackground = info.zones.begin(); iterBackground != info.zones.end(); ++iterBackground) {
            chpTsv.addHeader("ZoneType" + ::getInt(iIndex) , ::getDouble(iterBackground->centerx, 6) + "," + ::getDouble(iterBackground->centery, 6) + "," + ::getDouble(iterBackground->background, 6));
            iIndex++;
          }
          chpTsv.writeTsv_v2(chpFileName);
        } else {
          chpTsv.writeTsv(chpFileName);
        }

        if (!o.no_body) {
          doLegacyExpression(o, chpTsv, data);
          doLegacyGenotyping(o, chpTsv, data);
          doLegacyGenotyping(o, chpTsv, data);
          doLegacyUniversal(o, chpTsv, data);
          //        doLegacyReseq(o, chpTsv, data);
        }
        chpTsv.clear();
        delete chp;
      } else {
        Err::errAbort("A problem occurred while processing file (chp = NULL) " + o.chp_files[chp_index]);
      }
      bSuccessful = (data != NULL);
    }
  } catch (...) {
    Err::errAbort("A problem occurred while processing file (Exception thrown) " + o.chp_files[chp_index]);
  }
  return (bSuccessful);
}


void doLegacyExpression(opt& o, affx::TsvFile& chpTsv, affymetrix_fusion_io::FusionCHPLegacyData* data)
{
  std::vector<std::string> vProbeSetNames;
  if (o.cdfFile != "") {
    //Verbose::out(1, "Reading CDF file: " + o.cdfFile);
    readCdfFile(o.cdfFile, false, vProbeSetNames);
  }
  affymetrix_fusion_io::FusionExpressionProbeSetResults entry;
  bool bFound = true;
  int iIndex = 0;
  std::vector< std::string > arHeaders;
  arHeaders.push_back("Index");
  arHeaders.push_back("ProbeSetName");
  arHeaders.push_back("Signal");
  arHeaders.push_back("DetectionPValue");
  arHeaders.push_back("Detection");
  arHeaders.push_back("DetectionString");
  arHeaders.push_back("NumPairs");
  arHeaders.push_back("NumUsedPairs");
  arHeaders.push_back("NumCommonPairs");
  arHeaders.push_back("NumChangePValue");
  arHeaders.push_back("Change");
  arHeaders.push_back("ChangeString");
  arHeaders.push_back("SignalLogRatioLow");
  arHeaders.push_back("SignalLogRatioHigh");

  while ((bFound = data->GetExpressionResults(iIndex, entry))) {
    int iCidx = 0;
    if (iIndex == 0) {
      if (!entry.HasCompResults()) {
        arHeaders.resize(8);
      }
      if (vProbeSetNames.size() == 0) {
        arHeaders.erase(arHeaders.begin() + 1);
      }
      for (int i = 0; i < arHeaders.size(); i++) {
        chpTsv.defineColumn(0, i, arHeaders[i]);
      }
      chpTsv.writeColumnHeaders_clvl(0);
      chpTsv.flush();
    }

    chpTsv.set(0, iCidx++, ::getInt(iIndex).c_str());
    if (vProbeSetNames.size() > 0) {
      if (iIndex >= (int)vProbeSetNames.size()) {
        Err::errAbort("A problem occurred while processing file (Wrong CDF file?)");
      }
      chpTsv.set(0, iCidx++, vProbeSetNames[iIndex].c_str());
    }
    chpTsv.set(0, iCidx++, ::getDouble(entry.GetSignal(), 6).c_str());
    chpTsv.set(0, iCidx++, ::getDouble(entry.GetDetectionPValue(), 6).c_str());
    chpTsv.set(0, iCidx++, ::getInt(entry.GetDetection()).c_str());
    chpTsv.set(0, iCidx++, entry.GetDetectionString().c_str());
    chpTsv.set(0, iCidx++, ::getInt(entry.GetNumPairs()).c_str());
    chpTsv.set(0, iCidx++, ::getInt(entry.GetNumUsedPairs()).c_str());
    if (entry.HasCompResults()) {
      chpTsv.set(0, iCidx++, ::getInt(entry.GetNumCommonPairs()).c_str());
      chpTsv.set(0, iCidx++, ::getDouble(entry.GetChangePValue(), 6).c_str());
      chpTsv.set(0, iCidx++, ::getInt(entry.GetChange()).c_str());
      chpTsv.set(0, iCidx++, entry.GetChangeString().c_str());
      chpTsv.set(0, iCidx++, ::getDouble(entry.GetSignalLogRatio(), 6).c_str());
      chpTsv.set(0, iCidx++, ::getDouble(entry.GetSignalLogRatioLow(), 6).c_str());
      chpTsv.set(0, iCidx++, ::getDouble(entry.GetSignalLogRatioHigh(), 6).c_str());
    }
    chpTsv.writeLevel(0);
    iIndex++;
  }
  if ((iIndex > 0) && (vProbeSetNames.size() > 0)) {
    if (iIndex != (int)vProbeSetNames.size()) {
      chpTsv.clear();
      Err::errAbort("A problem occurred while processing file (Wrong CDF file?)");
    }
  }
}

void doLegacyGenotyping(opt& o, affx::TsvFile& chpTsv, affymetrix_fusion_io::FusionCHPLegacyData* data)
{
  std::vector<std::string> vProbeSetNames;
  if (o.cdfFile != "") {
    //Verbose::out(1, "Reading CDF file: " + o.cdfFile);
    readCdfFile(o.cdfFile, true, vProbeSetNames);
  }
  affymetrix_fusion_io::FusionGenotypeProbeSetResults g;
  bool bFound = true;
  int iIndex = 0;
  std::vector< std::string > arHeaders;
  arHeaders.push_back(std::string("Index"));
  arHeaders.push_back(std::string("ProbeSetName"));
  arHeaders.push_back(std::string("AlleleCall"));
  arHeaders.push_back(std::string("AlleleCallString"));
  arHeaders.push_back(std::string("Confidence"));
  arHeaders.push_back(std::string("PValueAA"));
  arHeaders.push_back(std::string("PValueAB"));
  arHeaders.push_back(std::string("PValueBB"));
  arHeaders.push_back(std::string("RAS1"));
  arHeaders.push_back(std::string("RAS2"));

  while ((bFound = data->GetGenotypingResults(iIndex, g))) {
    int iCidx = 0;
    if (iIndex == 0) {
      if (vProbeSetNames.size() == 0) {
        arHeaders.erase(arHeaders.begin() + 1);
      }
      for (int i = 0; i < arHeaders.size(); i++) {
        chpTsv.defineColumn(0, i, arHeaders[i]);
      }
      chpTsv.writeColumnHeaders_clvl(0);
      chpTsv.flush();
    }
    chpTsv.set(0, iCidx++, ::getInt(iIndex).c_str());
    if (vProbeSetNames.size() > 0) {
      if (iIndex >= (int)vProbeSetNames.size()) {
        Err::errAbort("A problem occurred while processing file (Wrong CDF file?)");
      }
      chpTsv.set(0, iCidx++, vProbeSetNames[iIndex].c_str());
    }
    chpTsv.set(0, iCidx++, ::getInt(g.GetAlleleCall()).c_str());
    chpTsv.set(0, iCidx++, g.GetAlleleCallString().c_str());
    chpTsv.set(0, iCidx++, ::getDouble(g.GetConfidence(), 6).c_str());
    chpTsv.set(0, iCidx++, ::getDouble(g.GetPValueAA(), 6).c_str());
    chpTsv.set(0, iCidx++, ::getDouble(g.GetPValueAB(), 6).c_str());
    chpTsv.set(0, iCidx++, ::getDouble(g.GetPValueBB(), 6).c_str());
    //      chpTsv.set(0,iCidx++,::getDouble(g.GetPValueNoCall(), 6).c_str());
    chpTsv.set(0, iCidx++, ::getDouble(g.GetRAS1(), 6).c_str());
    chpTsv.set(0, iCidx++, ::getDouble(g.GetRAS2(), 6).c_str());
    chpTsv.writeLevel(0);
    iIndex++;

  }
  if ((iIndex > 0) && (vProbeSetNames.size() > 0)) {
    if (iIndex != (int)vProbeSetNames.size()) {
      chpTsv.clear();
      Err::errAbort("A problem occurred while processing file (Wrong CDF file?)");
    }
  }
}

void doLegacyUniversal(opt& o, affx::TsvFile& chpTsv, affymetrix_fusion_io::FusionCHPLegacyData* data)
{
  std::vector<std::string> vProbeSetNames;
  if (o.cdfFile != "") {
    //Verbose::out(1, "Reading CDF file: " + o.cdfFile);
    readCdfFile(o.cdfFile, false, vProbeSetNames);
  }
  affymetrix_fusion_io::FusionUniversalProbeSetResults entry;
  std::vector< std::string > arHeaders;
  arHeaders.push_back(std::string("Index"));
  arHeaders.push_back(std::string("ProbeSetName"));
  arHeaders.push_back(std::string("BackGround"));

  bool bFound = true;
  int iIndex = 0;
  while ((bFound = data->GetUniversalResults(iIndex, entry))) {
    int iCidx = 0;
    if (iIndex == 0) {
      if (vProbeSetNames.size() == 0) {
        arHeaders.erase(arHeaders.begin() + 1);
      }
      for (int i = 0; i < arHeaders.size(); i++) {
        chpTsv.defineColumn(0, i, arHeaders[i]);
      }
      chpTsv.writeColumnHeaders_clvl(0);
      chpTsv.flush();
    }
    chpTsv.set(0, iCidx++, ::getInt(iIndex).c_str());
    if (vProbeSetNames.size() > 0) {
      if (iIndex >= (int)vProbeSetNames.size()) {
        Err::errAbort("A problem occurred while processing file (Wrong CDF file?)");
      }
      chpTsv.set(0, iCidx++, vProbeSetNames[iIndex].c_str());
    }
    chpTsv.set(0, iCidx++, ::getDouble(entry.GetBackground(), 6).c_str());
    chpTsv.writeLevel(0);
    iIndex++;
  }
  if ((iIndex > 0) && (vProbeSetNames.size() > 0)) {
    if (iIndex != (int)vProbeSetNames.size()) {
      chpTsv.clear();
      Err::errAbort("A problem occurred while processing file (Wrong CDF file?)");
    }
  }
}

void doLegacyReseq(opt& o, affx::TsvFile& chpTsv, affymetrix_fusion_io::FusionCHPLegacyData* data)
{
  std::vector<std::string> vProbeSetNames;
  if (o.cdfFile != "") {
    //Verbose::out(1, "Reading CDF file: " + o.cdfFile);
    readCdfFile(o.cdfFile, false, vProbeSetNames);
  }
  affymetrix_fusion_io::FusionResequencingResults entry;
  bool bFound = true;
  int iIndex = 0;
  std::vector< std::string > arHeaders;
  arHeaders.push_back(std::string("Index"));
  arHeaders.push_back(std::string("ProbeSetName"));

  while ((bFound = data->GetReseqResults(entry))) {
    int iCidx = 0;
    if (iIndex == 0) {
      if (vProbeSetNames.size() == 0) {
        arHeaders.erase(arHeaders.begin() + 1);
      }
      for (int i = 0; i < arHeaders.size(); i++) {
        chpTsv.set(0, i, arHeaders[i]);
      }
      chpTsv.writeLevel(0);
    }
    chpTsv.set(0, iCidx++, ::getInt(iIndex).c_str());
    if (vProbeSetNames.size() > 0) {
      if (iIndex >= (int)vProbeSetNames.size()) {
        Err::errAbort("A problem occurred while processing file (Wrong CDF file?)");
      }
      chpTsv.set(0, iCidx++, vProbeSetNames[iIndex].c_str());
    }

    chpTsv.writeLevel(0);
    iIndex++;
  }

  if ((iIndex > 0) && (vProbeSetNames.size() > 0)) {
    if (iIndex != (int)vProbeSetNames.size()) {
      chpTsv.clear();
      Err::errAbort("A problem occurred while processing file (Wrong CDF file?)");
    }
  }
}
