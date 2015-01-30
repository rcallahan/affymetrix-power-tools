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
 * @file   EngineUtil.cpp
 * @author Chuck Sugnet
 * @date   Tue Jul 11 12:47:10 2006
 * 
 * @brief  Some common functions that application executables make use of.
 */

//
#include "chipstream/EngineUtil.h"
//
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/ChipStreamFactory.h"
#include "chipstream/GcAdjust.h"
#include "chipstream/PmAdjusterFactory.h"
#include "chipstream/QuantDabg.h"
#include "chipstream/QuantMethodFactory.h"
//
#include "calvin_files/fusion/src/FusionCDFData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "file/CELFileData.h"
#include "file/TsvFile/ClfFile.h"
#include "file/TsvFile/PgfFile.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Fs.h"
#include "util/PgOptions.h"
#include "util/Util.h"
#include "util/Verbose.h"

using namespace std;
using namespace affx; 
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_io;

/**
 * Read in a file of probe IDs
 */
void EngineUtil::readProbeFile( std::vector< std::vector< probeid_t > > &probeIds, 
				const std::string& probeFileName) {

    affx::TsvFile tsv;
    probeid_t probeId;
    int channelIx = 0;
    /* sanity checks */
    if (tsv.open(probeFileName) != TSV_OK) {
        Err::errAbort("Couldn't open file: " + string(probeFileName) + " to read.");
    }

    if (tsv.cname2cidx(0, "probe_id") != TSV_ERR_NOTFOUND) {
        tsv.bind(0, "probe_id", &probeId, TSV_BIND_REQUIRED);
    }
    else {
        Err::errAbort("Unable to find probe_id column in file '"+ToStr(probeFileName)+"'");
    }

    if (tsv.cname2cidx(0, "channel") != TSV_ERR_NOTFOUND) {
        tsv.bind(0, "channel", &channelIx, TSV_BIND_REQUIRED);
    }
//     else
//         Err::errAbort("Unable to find channel column in file '"+ToStr(probeFileName)+"'");

    while (tsv.nextLevel(0) == TSV_OK) {
        if (channelIx >= probeIds.size()) {
            probeIds.resize(channelIx + 1);
        }
        probeIds[channelIx].push_back(TO_ZEROBASED(probeId));
    }
    tsv.close();
}

void EngineUtil::readProbeFile(std::vector<probeid_t> &probeIds, 
                               const std::string& probeFileName) {

    affx::TsvFile tsv;
    probeid_t probeId;
    /* sanity checks */
    if (tsv.open(probeFileName) != TSV_OK) 
        Err::errAbort("Couldn't open file: " + string(probeFileName) + " to read.");
    if (tsv.cname2cidx(0, "probe_id") != TSV_ERR_NOTFOUND)
        tsv.bind(0, "probe_id", &probeId, TSV_BIND_REQUIRED);
    else
        Err::errAbort("Unable to find probe_id column in file '"+ToStr(probeFileName)+"'");

    while (tsv.nextLevel(0) == TSV_OK) {
        probeIds.push_back(TO_ZEROBASED(probeId));
    }
    tsv.close();
}

void EngineUtil::readProbeFiles(std::vector<probeid_t> &probeIds,
                                const vector<std::string>& probeFileNames) {
    for (int i=0;i<probeFileNames.size();i++) {
        readProbeFile(probeIds,probeFileNames[i]);
    }
}


/**
 * Read in a file of probeset names or IDs. IDs trumps names.
 */
void EngineUtil::readProbeSetFile(std::vector<std::string> &probeSets, 
                                  const std::string& probeSetFileName) {
    affx::TsvFile tsv;
    string probeset;

    //
    if (tsv.open(probeSetFileName) != TSV_OK) 
        Err::errAbort("Couldn't open file: " + probeSetFileName + " to read.");

    //
    int cidx;
    // we try all of these in a case insenstive manner as other people like to use spaces.
    // we (APT) dont like spaces in column names as that causes problems
    // (Some programs treat spaces as tabs.)
    if ((cidx=tsv.cname2cidx(0, "probeset_id",affx::TSV_OPT_CASEINSENSTIVE)) != TSV_ERR_NOTFOUND) {
      tsv.bind(0, cidx, &probeset,TSV_BIND_REQUIRED);
    }
    else if ((cidx=tsv.cname2cidx(0, "probeset id",affx::TSV_OPT_CASEINSENSTIVE)) != TSV_ERR_NOTFOUND) {
      tsv.bind(0, cidx, &probeset, TSV_BIND_REQUIRED);
    }
    else if ((cidx=tsv.cname2cidx(0, "probe set id",affx::TSV_OPT_CASEINSENSTIVE)) != TSV_ERR_NOTFOUND) {
      tsv.bind(0, cidx, &probeset, TSV_BIND_REQUIRED);
    }
    else if ((cidx=tsv.cname2cidx(0, "probeset_name",affx::TSV_OPT_CASEINSENSTIVE)) != TSV_ERR_NOTFOUND) {
      tsv.bind(0, cidx, &probeset, TSV_BIND_REQUIRED);
    }
    else
        Err::errAbort("Unable to find probeset_id or probeset_name column in file '"+probeSetFileName+"'");

    // slurp the list.
    while (tsv.nextLevel(0) == TSV_OK)
        probeSets.push_back(probeset);
    
    tsv.close();
}

void EngineUtil::readProbeSetFiles(std::vector<std::string> &probeSets, 
                                   const std::vector<std::string>& probeSetFileNames) {
    for (int i=0;i<probeSetFileNames.size();i++) {
        readProbeSetFile(probeSets,probeSetFileNames[i]);
    }
}

/**
 * Read in the probe class file
 *
 * @param file   - The file to read
 * @param masks  - Vector of probe masks to fill in
 */
void EngineUtil::readProbeClassFile(const string &file, const int arraySize, 
                                    map<string, vector<bool> > &classMap) {
    affx::TsvFile tsv;
    int probeId = 0;
    string probeClass;
    map<string, vector<bool> >::iterator classMapIter;


    tsv.bind(0, "probe_id", &probeId, TSV_BIND_REQUIRED);
    tsv.bind(0, "class", &probeClass, TSV_BIND_REQUIRED);
    if (tsv.open(file) != TSV_OK)
        Err::errAbort("Couldn't open file: '" + file + "' to read.");


    while (tsv.nextLevel(0) == TSV_OK) {
        classMapIter = classMap.find(probeClass);
        if (classMapIter == classMap.end()) {
            classMap[probeClass].resize(arraySize,false);
        }
        classMap[probeClass].at(probeId-1) = true;
    }
    tsv.close();
}


/**
 * Return true if the qType of chip is found in the list of possible
 * types that are ok. Note that this is a case sensitive operation and
 * is meant to be so.
 * We want to be restrictive to avoid problems later on.
 * Filesystems we care about are at least case-preserving (HFS+)
 */
bool EngineUtil::checkChipTypeOk(const std::string &qType, const std::vector<std::string> &possibleTypes) {
    for (unsigned int i = 0; i < possibleTypes.size(); i++) {
        if (qType == possibleTypes[i]) {
            return true;
        }
    }
    // maybe they have the wrong case?
    // we always return false, but we can say what is wrong.
    std::string qType_lower=Util::downcaseString(qType);
    for (unsigned int i = 0; i < possibleTypes.size(); i++) {
        std::string possibleType_lower=Util::downcaseString(possibleTypes[i]);
        if (qType_lower == possibleType_lower) {
            Verbose::out(1,"Case mismatch in chiptypes. ('"+qType+"'!='"+possibleTypes[i]+"') Case must match.");
            return false;
        }
    }

    return false;
}

/**
 * Check chip types in TSV file against those in chipTypes vector
 */
void EngineUtil::checkTsvFileChipType(const string &path, 
                                      const std::vector<std::string> &chipTypes,
                                      bool mustSee) {
    affx::TsvFile tsv;
    string chip;
    vector<string> chipSeen;
    if (tsv.open(path) != TSV_OK)
        Err::errAbort("Couldn't open file: '" + path + "' to read.");
    bool chipMatch = false;
    if (chipTypes[0] == "")
        chipMatch = true; // treat empty as wildcard.
    while (tsv.headersFindNext("chip-type", chip) == TSV_OK) {
        chipSeen.push_back(chip);
        for (int i = 0; i<chipTypes.size(); i++) {
            if (chip == chipTypes[i])
                chipMatch = true;
        }
    }
    if (!chipMatch && (chipSeen.size() > 0 || mustSee)) {
        string foundTypes, expectedTypes;
        for (int i=0; i<chipSeen.size(); i++)
            foundTypes += chipSeen[i] + ", ";
        for (int i=0; i<chipTypes.size(); i++)
            expectedTypes += chipTypes[i] + ", ";

        Err::errAbort("Chiptype(s) [" + foundTypes + "] in file '" + path + "' does not match the set of valid chip types [" + expectedTypes + "]");
    }
}

/**
 * Check to make sure that there is at least one match between the two
 * vectors of chip types. Call Err::errAbort() if not.
 */
void EngineUtil::checkChipTypeVectors(std::vector<std::string> &chipTypesOk,
                                      std::vector<std::string> &chipTypesCheck) {
    for (unsigned int i = 0; i < chipTypesCheck.size(); i++)
        if (checkChipTypeOk(chipTypesCheck[i],chipTypesOk))
            return;
    std::string chipTypeString = "";
    for (unsigned int i = 0; i < chipTypesOk.size(); i++)
        chipTypeString += chipTypesOk[i] + ", ";
    Verbose::out(1,"Expected Chip Types: " + chipTypeString);
    chipTypeString = "";
    for (unsigned int i = 0; i < chipTypesCheck.size(); i++)
        chipTypeString += chipTypesCheck[i] + ", ";
    Verbose::out(1,"Found Chip Types: " + chipTypeString);
    Err::errAbort("No matching chip types found.");
}

/** 
 * Check to make sure all the cel files are compatible with the
 * chiptypes supported and number or probes on them. Treats the empty
 * string "" as a wildcard, so if chipTypes contains "" then any chip
 * type is accepted, but only if number of probes is correct. Calls
 * Err::errAbort() if something is wrong.
 * 
 * @param chipTypes - Names of chiptypes that are compatible.
 * @param probeCount - How many probes (aka features, ada cels) should be on array.
 * @param celFiles - Vector of cel file names to check.
 */
void EngineUtil::checkCelChipTypes(std::vector<std::string> &chipTypes, uint32_t probeCount,
                                   std::vector<std::string> &celFiles, uint32_t rows, uint32_t cols) {
    map<string,string> seenFiles;
    string root;
    FusionCELData cel;
  
    vector<string>::iterator celIx;
    bool allOk = true;
    int maxWarnings = 0;
    int vLevel = 1;
    for (celIx = celFiles.begin(); celIx != celFiles.end(); celIx++) {
        root = Fs::basename(*celIx);
        /* check to see we haven't seen this file before. */
        if (seenFiles.find(root) != seenFiles.end()) {
            Err::errAbort("Filename: '" + root + "' seen multiple times. Filenames must be unique.");
        }
        else {
            seenFiles[root] = "seen";
        }
  
        try {
          std::string tmp_unc_name=Fs::convertToUncPath(*celIx);
          cel.SetFileName(tmp_unc_name.c_str());
            if (!cel.ReadHeader())
              Err::errAbort("Can't read file: "+FS_QUOTE_PATH(tmp_unc_name));

            bool hack = false;
            std::vector<std::wstring> data_channels = cel.GetChannels();
            if (data_channels.empty()) {
                std::wstring temp = L"Default Group";
                data_channels.push_back(temp);
            }
            else if (data_channels.size() > 1) {
                affymetrix_calvin_io::GenericData *gd = cel.GetGenericData();
                if (gd != NULL) {
                    if (gd->DataGroupCnt() != data_channels.size()) {
                        // Falcon Cel files with QC info have 3 data groups and two channels
                        // Internal single channel cel files from MC system at one point reported two channels, but only had one data group
                        if (gd->DataGroupCnt() == 1) {
                            Verbose::out(1,"WARNING: Bogus multi channel cel file hack triggered for cel file '"+*celIx+"'.");
                            hack = true;
                            std::wstring temp = L"Default Group";
                            data_channels.clear();
                            data_channels.push_back(temp);
                        }
                    }
                }
            }
            for (int chanIx = 0; chanIx < data_channels.size(); chanIx++) {
                if (!hack) {
                    cel.SetActiveDataGroup(data_channels[chanIx]);
                }

                if (cel.GetNumCells() != probeCount) {
                    Verbose::out(vLevel, "Wrong number of cells in channel " + ToStr(chanIx) + ": Expecting: " + ToStr(probeCount) + " and " +
                                 *celIx + " has: " + ToStr(cel.GetNumCells()));
                    allOk = false;
                }
                else if (cel.GetCols() != cols || cel.GetRows() != rows) {
                    Verbose::out(vLevel, "Row/Col mismatch in channel " + ToStr(chanIx) + ". Expecting " + ToStr(rows) + " rows and " + 
                                 ToStr(cols) + " cols. " + *celIx + " has " + ToStr(cel.GetRows()) +
                                 " rows and " + ToStr(cel.GetCols()) + " cols.");
                    allOk = false;
                }
            }
            if (!checkChipTypeOk("", chipTypes) && 
                !checkChipTypeOk(StringUtils::ConvertWCSToMBS(cel.GetChipType()), chipTypes)) {
                Verbose::out(vLevel, "Wrong CEL ChipType: expecting: '" + chipTypes[0]  + "' and " +
                             *celIx + " is: '" + StringUtils::ConvertWCSToMBS(cel.GetChipType()) + "'");
                allOk = false;
            }
    
            if (allOk == false) { 
                maxWarnings++; 
            }
            if (maxWarnings > 10) {
                vLevel = 2;
            }
            cel.Close();
        }
        catch (...)
        {
            Err::errAbort("Can't read file: " + ToStr(*celIx));
        }
    }
    if (!allOk) {
        Err::errAbort("Different chiptypes or probe counts found. Only first 10 errors printed, see log for rest.");
    }
}

/** 
 * @brief Get a list from a text file to select a subset of probe sets for
 * analysis.
 * 
 * @param fileName - Path to file for reading probe subset from.
 * @param psGroups - Groups to be filled in, must be deleted later
 * @param psNameLoadMap - Map of probesets names to be loaded which is
 * filled in.
 */
void EngineUtil::makePSetNameSubset(const std::string &fileName, 
                                    std::vector<ProbeSetGroup *> &psGroups, 
                                    std::map<std::string, bool> &psNameLoadMap) {
    affx::TsvFile tsv;
    std::string probeset;
    tsv.bind(0, "probeset_id", &probeset, TSV_BIND_REQUIRED);
    if (tsv.open(fileName.c_str()) != TSV_OK) 
        Err::errAbort("Couldn't open file: " + fileName + " to read.");

    while (tsv.nextLevel(0) == TSV_OK) {
        psNameLoadMap[probeset] = true;
        ProbeSetGroup *ps = new ProbeSetGroup();
        ps->probeSetNames.push_back(Util::cloneString(probeset.c_str()));
        psGroups.push_back(ps);
    }
    tsv.close();
}

/** 
 * Open up a pgf and clf file pair and read the valid chiptypes
 * (i.e. HuEx-1_0-st-v2, HuEx-1_0-st-v1, HuEx-1_0-st-ta1) and number
 * of probes expected from them.
 * 
 * @param chipTypes - Names of chiptypes that this pgf file is compatible with.
 * @param probeCount - Number of probes that are on chip according to clf file.
 * @param pgfFile - Name of pgf file to be opened and read from.
 * @param clfFile - Name of clf file to be opened and read from.
 */
void EngineUtil::getPgfChipType(std::vector<std::string> &chipTypes, 
                                colrow_t& rowCount,
                                colrow_t& colCount,
                                int& probeCount,
                                const std::string &pgfFile, 
                                const std::string &clfFile) {
    PgfFile pgf;
    ClfFile clf;
    /* Fill in the chiptypes from the pgf file. */
    if (pgf.open(pgfFile) != TSV_OK) {
        Err::errAbort("EngineUtil::getChipType() - Can't open Pgf file: '" + pgfFile + "' to read.");
    }
    chipTypes.clear();
    string chip;
    while (pgf.m_tsv.headersFindNext("chip_type", chip) == TSV_OK) {
        chipTypes.push_back(chip);
    }
    pgf.close();

    /* Fill in the probecount from the clf file. */
    if (clf.open(clfFile) != TSV_OK) {
        Err::errAbort("EngineUtil::getChipType() - Can't open clf file: '" + clfFile + "' to read.");
    }
    int rows = 0, cols = 0;
    if (clf.m_tsv.getHeader("rows", rows) != TSV_OK) 
        Err::errAbort("EngineUtil::getPgfChipType() - Clf file format not valid, no 'rows' field.");
    if (clf.m_tsv.getHeader("cols", cols) != TSV_OK) 
        Err::errAbort("EngineUtil::getPgfChipType() - Clf file format not valid, no 'cols' field.");
    rowCount = rows;
    colCount = cols;
    probeCount = (uint32_t)rows * cols;
    clf.close();
}

/** 
 * Open up a spf (simple probe format) file and read the valid chiptypes 
 * (i.e. HuEx-1_0-st-v2, HuEx-1_0-st-v1, HuEx-1_0-st-ta1) and number
 * of probes expected.
 * 
 * @param chipTypes - Names of chiptypes that this pgf file is compatible with.
 * @param probeCount - Number of probes that are on chip according to spf file.
 * @param spfFile - Name of spf file to be opened and read from.
 */
void EngineUtil::getSpfChipType(std::vector<std::string>& chipTypes, 
                                colrow_t& rowCount, 
                                colrow_t& colCount,
                                int& probeCount,
                                int& probeSetCount,
                                const std::string& spfFile) {
    TsvFile tsv;
    /* Fill in the chiptypes from the pgf file. */
    if (tsv.open(spfFile) != TSV_OK) {
        Err::errAbort("EngineUtil::getSpfChipType() - Can't open file: '" + spfFile + "' to read.");
    }
    chipTypes.clear();
    string chip;
    while (tsv.headersFindNext("chip_type", chip) == TSV_OK) {
        chipTypes.push_back(chip);
    }
    int rows = 0, cols = 0, probesets = 0;
    if (tsv.getHeader("num-rows", rows) != TSV_OK) 
        Err::errAbort("EngineUtil::getSpfChipType() - Spf file format not valid, no 'num-rows' field.");
    rowCount = rows;
    if (tsv.getHeader("num-cols", cols) != TSV_OK) 
        Err::errAbort("EngineUtil::getSpfChipType() - Spf file format not valid, no 'num-cols' field.");
    colCount = cols;
    probeCount = (uint32_t)rows * cols;
    if (tsv.getHeader("num-probesets", probesets) != TSV_OK) 
        Err::errAbort("EngineUtil::getSpfChipType() - Spf file format not valid, no 'num-probesets' field.");
    probeSetCount = probesets;

    tsv.close();

  
}
  
/** 
 * Open up a cdf file and get the chip type (i.e. 'Mapping250K_Sty') and number
 * of probes on the chip.
 * 
 * @param chipTypes - Names of chiptypes that this pgf file is compatible with.
 * @param rowCount - Number of rows on the chip.
 * @param colCount - Number of columns on chip.
 * @param probeCount - Number of probes that are on chip according to clf file.
 * @param cdfFile - Name of cdf file to be opened and read from.
 */
void EngineUtil::getCdfChipType(std::vector<std::string> &chipTypes, 
                                colrow_t& rowCount, 
                                colrow_t& colCount,
                                int& probeCount,
                                int& probeSetCount,
                                const std::string &cdfFile) {
    FusionCDFData cdf;
    std::string tmp_unc_name=Fs::convertToUncPath(cdfFile);
    cdf.SetFileName(tmp_unc_name.c_str());
    try {
        if (!cdf.ReadHeader()) {
            Err::errAbort("EngineUtil::getCdfChipType() - Can't read header for " +
                          FS_QUOTE_PATH(tmp_unc_name) +
                          ". Description: " + cdf.GetError());
        }
    }
    catch(...) {
      Err::errAbort("Can't read file: " + FS_QUOTE_PATH(tmp_unc_name));
    }
    chipTypes.clear();
    chipTypes = cdf.GetChipTypes();
    FusionCDFFileHeader &cdfHeader = cdf.GetHeader();
    rowCount = cdfHeader.GetRows();
    colCount = cdfHeader.GetCols();
    probeCount = (uint32_t)rowCount * colCount;
    probeSetCount = cdfHeader.GetNumProbeSets();
}

/**
 * @brief Print out a header and then a series of self documenting
 * objects.
 * @param header - Descriptive header for this collection of self documenting objects.
 * @param doc - Collection of self documenting objects.
 */
void EngineUtil::printSelfDocs(const char *header, vector<SelfDoc> doc) {
    unsigned int i = 0;

    int maxLength = 0;
    int extraChars = 3;
    cout << endl << header << endl;
    for (i=0; i < doc.size(); i++) {
        string s = doc[i].getDocName();
        maxLength = max((int)maxLength, (int)s.length());
    }
    maxLength += 3 + extraChars;
    maxLength = min(maxLength, 25);
    for (i=0; i < doc.size(); i++) {
        int currentLength = 0;
        string s = doc[i].getDocName();
        cout << "   " << s;
        currentLength = s.length() + 3;
        while (currentLength < maxLength) {
            cout.put(' ');
            currentLength++;
        }
        PgOptions::printStringWidth(doc[i].getDocDescription().c_str(), maxLength, currentLength);
        cout << endl;
    }
}

/**
 * Print out a little ditty for the given type of self documenting object.
 *
 * @param query - type of analysis/adjustment/summarization information is desired for.
 */
void EngineUtil::explainParameter(const std::string& query) {
    ChipStreamFactory cFactory;
    PmAdjusterFactory aFactory;
    QuantMethodFactory qFactory(QuantMethodFactory::Expression);
    QuantMethodFactory gFactory(QuantMethodFactory::GenoType);
    AnalysisStreamFactory asFactory;
    vector<vector<SelfDoc> > docs;

  
    docs.push_back(cFactory.getDocs());
    docs.push_back(aFactory.getDocs());
    docs.push_back(qFactory.getDocs());
    docs.push_back(gFactory.getDocs());
    docs.push_back(asFactory.getDocs());

    for (unsigned int i = 0; i < docs.size(); i++) {
        for (unsigned int j = 0; j < docs[i].size(); j++) {
            if (docs[i][j].getDocName() == query) {
                SelfDoc::printExplanation(docs[i][j],cout);
                exit(0);
            }
        }
    }
    Err::errAbort("Didn't find any documentation for: '" + ToStr(query) + "'");
}


void EngineUtil::getCelFiles(std::vector<std::string> &celFiles, BaseEngine *engine) {
    if (engine->getOpt("cel-files")!="") {
        affx::TsvFile tsv;
#ifdef WIN32
        tsv.m_optEscapeOk = false;
#endif
        std::string celFilesFile = engine->getOpt("cel-files");
        string file;
        tsv.bind(0, "cel_files", &file, TSV_BIND_REQUIRED);
        if (tsv.open(celFilesFile) != TSV_OK) {
            Err::errAbort("Couldn't open cel-files file: " + celFilesFile);
        }
        tsv.rewind();
        while (tsv.nextLevel(0) == TSV_OK) {
            celFiles.push_back(file);
        }
        tsv.close();
        Verbose::out(1, "Read " + ToStr(celFiles.size()) + " cel files from: " + Fs::basename(celFilesFile));
    } 
    else {
        celFiles = engine->getOptVector("cels");
        if (celFiles.size() == 0) {
            for (vector<const char *>::size_type i = 0; i < engine->getArgCount(); i++)
                celFiles.push_back(engine->getArg(i));
        }
#ifdef WIN32
        // Windows doesn't do wildcard expansion in windows, let the user know.
        if (celFiles.size() > 0 && celFiles[0].find("\*") != string::npos) {
            Err::errAbort("Wildcard ('*') expansion not available in windows, please use --cel-files option.");
        }
#endif
    }
}

void EngineUtil::getChpFiles(std::vector<std::string> &chpFiles, BaseEngine *engine) {
    if (engine->getOpt("chp-files")!="") {
        affx::TsvFile tsv;
#ifdef WIN32
        tsv.m_optEscapeOk = false;
#endif
        std::string chpFilesFile = engine->getOpt("chp-files");
        string file;
        tsv.bind(0, "chp_files", &file, TSV_BIND_REQUIRED);
        if (tsv.open(chpFilesFile) != TSV_OK) {
            Err::errAbort("Couldn't open chp-files file: " + chpFilesFile);
        }
        tsv.rewind();
        while (tsv.nextLevel(0) == TSV_OK) {
            chpFiles.push_back(file);
        }
        tsv.close();
        Verbose::out(1, "Read " + ToStr(chpFiles.size()) + " chp files from: " + Fs::basename(chpFilesFile));
    } 
    else {
        chpFiles = engine->getOptVector("chps");
        if (chpFiles.size() == 0) {
            for (vector<const char *>::size_type i = 0; i < engine->getArgCount(); i++)
                chpFiles.push_back(engine->getArg(i));
        }
#ifdef WIN32
        // Windows doesn't do wildcard expansion in windows, let the user know.
        if (chpFiles.size() > 0 && chpFiles[0].find("\*") != string::npos) {
            Err::errAbort("Wildcard ('*') expansion not available in windows, please use --chp-files option.");
        }
#endif
    }
}
