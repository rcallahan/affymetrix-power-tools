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

/// @file  apt-data-subset.cpp
/// @brief Program make a smaller sample data set from a large one

//
#include "calvin_files/data/src/CELData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCelFileWriter.h"
#include "calvin_files/writers/src/CombinedDatConstants.h"
#include "calvin_files/writers/src/MultiChannelCelFileCollater.h"
#include "chipstream/ChipLayout.h"
#include "file/TsvFile/BgpFile.h"
#include "file/TsvFile/ClfFile.h"
#include "file/TsvFile/SpfFile.h"
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "stats/affy_random_sample.h"
#include "util/AptVersionInfo.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/MixedFileCheck.h"
#include "util/PgOptions.h"
#include "util/RowFile.h"
#include "util/TableFile.h"
#include "util/TextFileCheck.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

using namespace std;
using namespace affx;
using namespace affymetrix_calvin_io;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;

#define NO_WAVELENGTH L"wavelength unavailable"
#define MULTI_CHANNEL_FILE_ID "affymetrix-calvin-multi-intensity"

/**
   \page apt-data-subset MANUAL: apt-data-subset (\aptversion)

  \par This application is supported on Linux and OS X platforms only. It has not been tested on Windows platforms.
   

*/


/**
 * Helper class for meta probesets in apt-data-subset.cpp. Holds a 
 * meta probeset name and a list of sub-probesets making up the meta set.
 */
class Meta {
public :
  string id;  ///< Name/Identifier of the meta probset
  vector<string> probeSets; ///< Ids/names of the sub probesets that make up this meta-set

  /// Standard sorting version
  bool operator< (const Meta &r) const {
    return id < r.id;
  }

};


/**
 * Helper class containing information about the subset transformation
 * like the set of probesets to be kept and the mapping of original
 * probe ids to new probe ids.
 */
class FileTranslation {

public :

  /// Default constructor
  FileTranslation() {
    totalProbes = 0;
    nRow = 0;
    nCol = 0;
    pOffset = 1;
  }
 
  vector<int> probeMask;  ///< Mapping of original probe_ids (zero based) to new ids. -1 means not in subset
  set<string> keepProbeSets; ///< Which probesets will be kept in the new subset
  string probeField; ///< Identifier for probe id column in tsv files, usually "probe_id"
  string probeSetField; ///< Identifier for probeset id column in tsv files, usually "probeset_id"
  string chipType; ///< New chiptype
  int totalProbes; ///< How many probes are currently in subset  
  int nRow; ///< Number of rows on our spoofed cel file
  int nCol; ///< Number of columns on our spoofed cel file
  int pOffset; ///< Probe offset, most files this is 1 as they are 1 based, but some can be zero
  vector<string> magicOrder;
  vector<int> magicOrderTranslation;
  map<string,int> magicOrderMap;

  /// Add a probe to our transformation with default offset
  bool addProbe(int pId) {
    return addProbe(pId, pOffset);
  }

  /// Add a probe id to our mapping and transformation array with specified offset
  bool addProbe(int pId, int offset) {
    pId = pId - offset;
    if(probeMask[pId] < 0) {
      probeMask[pId] = totalProbes++;
      return true;
    }
    return false;
  }

  /**
   * Loop through a chip layout and make a mapping with probesets that have been
   * chosen to be kept making a map from current probe ids to new probe ids.
   */
  void fillInFromProbeSets(ChipLayout &layout) {
    set<string> keepersSeen;
    for(int i = 0; i < layout.getProbeSetCount(); i++) {
      ProbeListPacked pl = layout.getProbeListAtIndex(i);

      // Keep the probesets with probes identified
      for(int j = 0; j < pl.probe_cnt(); j++) {
        if(probeMask[pl.get_probeId(j)] >= 0) {
          keepProbeSets.insert(pl.get_name_string());
        }
      }

      // For probesets identified mark the probes
      if( keepProbeSets.find(pl.get_name_string()) != keepProbeSets.end()) {
        keepersSeen.insert(pl.get_name_string());
        for(int j = 0; j < pl.probe_cnt(); j++) {
          int pId = pl.get_probeId(j);
          addProbe(pId, 0);
        } 
      }
    }
    Verbose::out(2, "Expeting to get: " + ToStr(keepProbeSets.size()) + " but only saw: " + ToStr(keepersSeen.size()) + " probesets in the layout.");
    keepProbeSets = keepersSeen;
  }

  /// Make a new square chip based on current number of probes.
  void calcDimensions(int totalProbes) {
      double sq = sqrt(totalProbes);
      int width = Util::round(ceil(sq));
      nRow = width;
      nCol = width;
  }

};

/**
 * Hold all the options for the program
 */
class DsOptions {

public:
  string libOutputPrefix;
  string celOutputPrefix;
  string spfFileStr;
  string pgfFileStr;
  string clfFileStr;
  string probesetFileStr;
  string filePrefix;
  string celFileFile;
  string tsvConvertFileFile;
  string tsvSelectFileFile;
  string file5ConvertFileFile;
  string file5SelectFileFile;
  string bgpFile;
  string qccFileStr;
  string metaFile;
  string killFileStr;
  string magicOrderColStr;
  string magicOrderFileStr;
  string magicOrderGroupStr;
  string chipRoot;
  
  int numKeep;

  DsOptions(PgOptions *opts) {
   
    libOutputPrefix = opts->get("out-lib-dir");
    celOutputPrefix = opts->get("out-cel-dir");
    spfFileStr = opts->get("spf-file");
    pgfFileStr = opts->get("pgf-file");
    clfFileStr = opts->get("clf-file");
    filePrefix = opts->get("prefix");
    probesetFileStr = opts->get("probeset-file");
    celFileFile = opts->get("cel-files");
    tsvConvertFileFile = opts->get("tsv-convert-files");
    tsvSelectFileFile = opts->get("tsv-select-files");
    file5ConvertFileFile = opts->get("file5-convert-files");
    file5SelectFileFile = opts->get("file5-select-files");
    bgpFile = opts->get("bgp-file");
    qccFileStr = opts->get("qcc-file");
    metaFile = opts->get("meta-file");
    killFileStr = opts->get("kill-list-file");
    magicOrderColStr = opts->get("magic-order-column");
    magicOrderFileStr = opts->get("magic-order-file");
    magicOrderGroupStr = opts->get("magic-order-group");
    chipRoot = spfFileStr;
    if(spfFileStr=="") 
      chipRoot = pgfFileStr;
    numKeep = opts->getInt("num");
  }
};

/**
 * File5 options
 */ 
class File5Options {

public :
  string fileName;
  string groupName;
  string probeId;
  string probeSetId;
  int idOffset;

};

/**
 * User interface options
 */
void defineDataSubsetOptions(PgOptions* opts) {
  opts->setUsage("apt-data-subset - Create a sample data subset from a larger set of files.\n"
                 "usage:\n"
                 "   apt-data-subset --spf file.spf --num 10000 --cel-file celfiles.txt --out-dir mydatadir --prefix 10k\n");

  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "Display program options and extra documentation about usage.",
                     "false");
  opts->defineOption("v", "verbose", PgOpt::INT_OPT,
                     "How verbose to be with status messages 0 - quiet, 1 - usual messages, 2 - more messages.",
                     "1");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");
  opts->defineOption("", "cel-files", PgOpt::STRING_OPT,
                     "Text file specifying cel files to process, one per line with the first line being 'cel_files'.",
                     "");
  opts->defineOption("", "out-lib-dir", PgOpt::STRING_OPT,
                     "Directory to write result library files into.",
                     "");
  opts->defineOption("", "out-cel-dir", PgOpt::STRING_OPT,
                     "Directory to write result celf files into.",
                     "");
  opts->defineOption("", "spf-file", PgOpt::STRING_OPT,
                     "Spf file to read probesets data from.",
                     "");
  opts->defineOption("", "probeset-file", PgOpt::STRING_OPT,
                     "List of probesets to use as sample.",
                     "");
  opts->defineOption("", "prefix", PgOpt::STRING_OPT,
                     "Prefix to tag on before file type suffix",
                     "");

  opts->defineOption("", "probeset-field", PgOpt::STRING_OPT,
                     "Column name for probeset fields, usually probeset_id or probeset_name",
                     "probeset_id");

  opts->defineOption("", "bgp-file", PgOpt::STRING_OPT,
                     "File with GC background probes",
                     "");

  opts->defineOption("", "meta-file", PgOpt::STRING_OPT,
                     "Meta-probset file to be sampled from.",
                     "");

  opts->defineOption("", "pgf-file", PgOpt::STRING_OPT,
                     "pgf-probset file to be sampled from.",
                     "");

  opts->defineOption("", "clf-file", PgOpt::STRING_OPT,
                     "clf-probset file to be sampled from.",
                     "");

  opts->defineOption("", "kill-list-file", PgOpt::STRING_OPT,
                     "file with kill-list probes to be sampled from.",
                     "");
  opts->defineOption("", "magic-order-column", PgOpt::STRING_OPT,
                     "What column has the magic order for copynumber.",
                     "");
  opts->defineOption("", "magic-order-group", PgOpt::STRING_OPT,
                     "What group in the file5 has the magic order for copynumber.",
                     "");
  opts->defineOption("", "magic-order-file", PgOpt::STRING_OPT,
                     "What group in the file5 has the magic order for copynumber.",
                     "");

  opts->defineOption("", "qcc-file", PgOpt::STRING_OPT,
                     "file with QC probesets to be sampled from.",
                     "");

  opts->defineOption("", "tsv-convert-files", PgOpt::STRING_OPT,
                     "List of tsv based files to read from and transform.",
                     "");

  opts->defineOption("", "tsv-select-files", PgOpt::STRING_OPT,
                     "List of tsv based files to set probesets from.",
                     "");

  opts->defineOption("", "file5-convert-files", PgOpt::STRING_OPT,
                     "List of file5 based files to read from and transform.",
                     "");

  opts->defineOption("", "file5-select-files", PgOpt::STRING_OPT,
                     "List of file5 based files to set probesets from.",
                     "");

  opts->defineOption("", "num", PgOpt::INT_OPT,
                     "Number of probesets to sample for creation.",
                     "-1");

  opts->defineOption("", "spf-version", PgOpt::INT_OPT,
                     "Which spf version to write.",
                     "1");

  opts->defineOption("", "num-channels", PgOpt::INT_OPT,
                     "How many channels for the spf-file.",
                     "1");

}


/// Open an Spf file and write out the headers
void openSpfForWrite(affx::SpfFile& spffile,const std::string& fileName, 
                     const std::string &chipType, int version, int nCol, int nRow, int psCount,
                     int nChannels) {
  // the standard set of columns
  spffile.define_file(version);
  // stuff what will appear in the headers into the vars...
  spffile.m_header_chipTypes.clear();
  spffile.m_header_chipTypes.push_back(chipType);
  spffile.m_header_numProbesets=psCount;
  spffile.m_header_numCols=nCol;
  spffile.m_header_numRows=nRow;
  spffile.m_header_numChannels=nChannels;
  if(version == 4) {
    spffile.m_has_allele_info = 1;
    spffile.m_has_context_info = 1;
    spffile.m_has_channel_info = 1;
    spffile.m_has_rep_type_info = 1;
  }
  // ...and this will add the headers in the standard order.
  spffile.addStandardHeaders();
  // ready for data.
  spffile.writeSpf(fileName, version);
}

/// Open a tsv file and read the cel files we want to process
void fillInCelFiles(vector<string> &celFiles, const string &fileName) {
  celFiles.clear();
  affx::TsvFile tsv;
#ifdef WIN32
  tsv.m_optEscapeOk = false;
#endif
  string file;
  tsv.bind(0, "cel_files", &file, TSV_BIND_REQUIRED);
  if(tsv.open(fileName) != TSV_OK) {
    Err::errAbort("Couldn't open cell-files file: " + fileName);
  }
  tsv.rewind();
  while(tsv.nextLevel(0) == TSV_OK) {
    celFiles.push_back(file);
  }
  tsv.close();
  Verbose::out(1, "Read " + ToStr(celFiles.size()) + " cel filenames from: " + Fs::basename(fileName));
}

/// Open up a tsv file and read the collection of tsv files that we are going to want to process
void fillInTsvFiles(vector<string> &tsvFiles, vector<string> &ids, const string &fileName) {
  tsvFiles.clear();
  affx::TsvFile tsv;
#ifdef WIN32
  tsv.m_optEscapeOk = false;
#endif
  string file;
  string id;
  tsv.bind(0, "tsv_files", &file, TSV_BIND_REQUIRED);
  tsv.bind(0, "probeset_field", &id, TSV_BIND_OPTIONAL);
  if(tsv.open(fileName) != TSV_OK) {
    Err::errAbort("Couldn't open tsv-files file: " + fileName);
  }
  tsv.rewind();
  while(tsv.nextLevel(0) == TSV_OK) {
    tsvFiles.push_back(file);
    ids.push_back(id);
  }
  tsv.close();
  Verbose::out(1, "Read " + ToStr(tsvFiles.size()) + " tsv filenames from: " + Fs::basename(fileName));
}

/// Open up a tsv file and read the collection of tsv files that we are going to want to process
void fillInFile5Files(vector<File5Options> &file5Files,  const string &inputName) {
  file5Files.clear();
  affx::TsvFile tsv;
#ifdef WIN32
  tsv.m_optEscapeOk = false;
#endif
  string file;
  string group;
  string probeId;
  string probesetId;
  int idOffset;
  tsv.bind(0, "file5_group", &group, TSV_BIND_REQUIRED);
  tsv.bind(0, "probe_id_field", &probeId, TSV_BIND_REQUIRED);
  tsv.bind(0, "probe_id_offset", &idOffset, TSV_BIND_REQUIRED);
  tsv.bind(0, "probeset_field", &probesetId, TSV_BIND_REQUIRED);
  tsv.bind(0, "file5_file", &file, TSV_BIND_REQUIRED);
  if(tsv.open(inputName) != TSV_OK) {
    Err::errAbort("Couldn't open tsv-files file: " + inputName);
  }
  tsv.rewind();
  while(tsv.nextLevel(0) == TSV_OK) {
    File5Options f5Opts;
    f5Opts.fileName = file;
    f5Opts.groupName = group;
    f5Opts.probeId = probeId;
    f5Opts.idOffset = idOffset;
    f5Opts.probeSetId = probesetId;

    file5Files.push_back(f5Opts);
  }
  tsv.close();
  Verbose::out(1, "Read " + ToStr(file5Files.size()) + " file5 filenames from: " + Fs::basename(inputName));
}



/// Write a subset of cel file to a new file
void makeCelFileSubset(const string &orig, const string &outCel, 
                        const FileTranslation &trans) {
  FusionCELData cel;
  cel.SetFileName(orig.c_str());
  cel.Read();
  if(!cel.Read()) 
    Err::errAbort("\nCan't read cel file: " + cel.GetFileName() + 
                  "\n>>> Error reported: " + StringUtils::ConvertWCSToMBS(cel.GetError()));
  int numCells = cel.GetNumCells();
  int newNumCells = trans.nRow * trans.nCol;

  vector<float> intensities(newNumCells);
  fill(intensities.begin(), intensities.end(), 0.0f);
  vector<float> stdevs(newNumCells);
  fill(stdevs.begin(), stdevs.end(), 0.0f);
  vector<int16_t>  nPixels(newNumCells);
  fill(nPixels.begin(), nPixels.end(), 0);

  CelFileData out(outCel.c_str());
  out.SetCols(trans.nCol);
  out.SetRows(trans.nRow);
  out.SetArrayType(StringUtils::ConvertMBSToWCS(trans.chipType));
  out.SetIntensityCount(newNumCells);
  out.SetStdDevCount(newNumCells);
  out.SetPixelCount(newNumCells);
  CelFileWriter writer(out);

  for(int featIx = 0; featIx < numCells; featIx++) {
    if(trans.probeMask[featIx] >= 0) {
      intensities[trans.probeMask[featIx]] = cel.GetIntensity(featIx);
      stdevs[trans.probeMask[featIx]] = cel.GetStdv(featIx);
      nPixels[trans.probeMask[featIx]] = cel.GetPixels(featIx);
    }
  }
  writer.WriteStdDevs(stdevs);
  writer.WriteIntensities(intensities);
  writer.WritePixels(nPixels);
}

// void CopyHeaders(const std::string& celFile, const std::wstring& channel, FileHeader& outFileHdr) {
//   GenericData gReader;
//   std::wstring wavelength;
//   //   OpenReader(celFile, gReader);
//   gReader.Clear();
//   GenericFileReader reader;
//   reader.SetFilename(celFile);
//   reader.Open(gReader);
//   reader.Close();


//   // CopyFileHdr(gReader.Header(), outFileHdr);
//   GenDataHdrVectorIt beginVec, endVec;
//   gReader.GetGenericDataHdr()->GetParentIterators(beginVec, endVec);
//   while(beginVec != endVec) {
//     outFileHdr.GetGenericDataHdr()->AddParent(*beginVec);
//     beginVec++;
//   }
//   // CopyGenericHeader(gReader.GetGenericDataHdr(), outFileHdr.GetGenericDataHdr());
//   outFileHdr->SetFileTypeId(MULTI_CHANNEL_FILE_ID);
//   ParameterNameValueTypeIt beginNv, endNv;
//   gReader.GetGenericDataHdr()->GetNameValIterators(beginNv, endNv);
//   while(beginNv != endNv) {
//     if(beginNv->GetName() == AFFY_FILTER_WAVELENGTH)   {
//       wavelength = begin->GetValueText();
//     }
//     else     {
//       outFileHdr.GetGenericDataHdr()->AddNameValParam(*beginNv);
//     }
//     beginNv++;
//   }

//   DataGroupHeader dstGrp;
//   DataGroupHeader* srcGrp = gReader.Header().FindDataGroupHeader(channel);
//   //  CopyDataGroupHdr(*srcGrp, dstGrp);
//   dstGroup.SetName(wavelength);
//   DataSetHdrIt beginDs, endDs;
//   srcGrp->GetDataSetIterators(beginDs, endDs);
//   while(begin != end) {
//     DataSetHeader dstSet;
//     CopyDataSetHdr(*begin, dstSet);

// 	dstSet.SetName(srcGrp->GetName());
// 	dstSet.SetRowCnt(srcGrp->GetRowCnt());
// 	for(int i = 0; i < srcGrp->GetColumnCnt(); i++)
// 	{
// 		dstSet.AddColumn(srcGrp->GetColumnInfo(i));
// 	}
// 	ParameterNameValueTypeIt beginNvT, endNvT;
// 	srcGrp->GetNameValIterators(beginNvT, endNvT);
// 	while(beginNvT != endNvT)
// 	{
// 		dstSet.AddNameValParam(*begin);
// 		begin++;
// 	}
//     dstGroup.AddDataSetHdr(dstSet);
//     beginNvT++;
//   }
//   outFileHdr.AddDataGroupHdr(dstGrp);
//   gReader.Clear();
// }

// int32_t GetGrpIndex(FileHeader& hdr, const std::wstring& channel) {
//   for(int i = 0; i < hdr.GetDataGroupCnt(); i++)
//     {
//       if(hdr.GetDataGroup(i).GetName() == channel)
//         {
//           return i;
//         }
//     }
//   return -1;
// }

// void ExtractCel(const std::string& inFile, const std::wstring& channel, const std::string& outFile,
//                 const FileTranslation &trans) {
//   FileHeader outFileHdr;
//   outFileHdr.SetFilename(outFile);
//   CopyHeaders(inFile, channel, outFileHdr);
//   //  UpdateFileHeader(outFileHdr, channel);
//   outFileHeader.GetDataGroup(0).SetName(CelDataGroupName);
//   //  AddWavelengthParam(hdr.GetGenericDataHdr(), channel);
//   ParameterNameValueType p;
//   p.SetName(FILTER_PARAM_NAME);
//   p.SetValueText(channel);
//   outFileHeader.GetGenericDataHdr()->AddNameValParam(p);

//   outFileHeader.GetGenericDataHdr()->SetFileTypeId(INTENSITY_DATA_TYPE);

//   GenericFileWriter writer(&outFileHdr);
//   writer.WriteHeader();
//   GenericFileReader reader;
//   GenericData data;
//   reader.SetFilename(inFile);
//   reader.Open(data);
//   DataGroupWriter grpWriter = writer.GetDataGroupWriter(0);
//   int32_t grpIndex = GetGrpIndex(data.Header(), channel);
//   WriteDataGroup(reader.GetDataGroupReader(grpIndex), data.Header().GetDataGroup(grpIndex), grpWriter, true);
// }

/// Write a subset of cel file to a new file
void makeMultiCelFileSubset(const string &orig, const string &outCel, 
                        const FileTranslation &trans) {
  FusionCELData cel;
  cel.SetFileName(orig.c_str());
  cel.Read();
  
  if(!cel.Read()) 
    Err::errAbort("\nCan't read cel file: " + cel.GetFileName() + 
                  "\n>>> Error reported: " + StringUtils::ConvertWCSToMBS(cel.GetError()));
  WStringVector channelVec = cel.GetChannels();
  if(channelVec.size() != 2) {
    Err::errAbort("Expecting 2 channels but got: " + ToStr(channelVec.size()) + " in file: " + orig);
  }
  vector<string> tmpFileNamesToDelete;
  for(int i = 0; i < channelVec.size(); i++) {
    cel.SetActiveDataGroup(channelVec[i]);
    int numCells = cel.GetNumCells();
    int newNumCells = trans.nRow * trans.nCol;
    
    vector<float> intensities(newNumCells);
    fill(intensities.begin(), intensities.end(), 0.0f);
    vector<float> stdevs(newNumCells);
    fill(stdevs.begin(), stdevs.end(), 0.0f);
    vector<int16_t>  nPixels(newNumCells);
    fill(nPixels.begin(), nPixels.end(), 0);
    string tempFile = "celFile." + StringUtils::ConvertWCSToMBS(channelVec[i]) + "." + ToStr(i) + ".cel";
    tmpFileNamesToDelete.push_back(tempFile);
    CelFileData out(tempFile);

    /* Have to copy over all the parent headers as the CELData.cpp code won't recognize
       this as a multi-channel cel file unless there is a parent GenericDataHeader named
       MULTI_SCAN_ACQUISITION_DATA_TYPE ("affymetrix-calvin-multi-scan-acquisition"). 
       Apparently this is from the original dat file used for this cel file. Implying that
       Cel files only come from DAT files... */
    GenericDataHeader *inSrc =  cel.GetGenericData()->Header().GetGenericDataHdr();
    GenericDataHeader *outDst = out.GetGenericData().Header().GetGenericDataHdr();
    GenDataHdrVectorIt begin, end;
    inSrc->GetParentIterators(begin, end);
    while(begin != end)      {
      outDst->AddParent(*begin);
      begin++;
    }    
    //    out.SetActiveChannel(channelVec[i]);
    
    GenericDataHeader *hdr = out.GetGenericData().Header().GetGenericDataHdr();
    ParameterNameValueType pFilt;
    pFilt.SetName(FILTER_PARAM_NAME);
    pFilt.SetValueText(channelVec[i]);
    hdr->AddNameValParam(pFilt);
    ParameterNameValueType p;
    //    hdr->SetFileTypeId(INTENSITY_DATA_TYPE);
    p.SetName(AFFY_FILTER_WAVELENGTH);
    p.SetValueText(channelVec[i]);
    hdr->AddNameValParam(p);
    out.SetCols(trans.nCol);
    out.SetRows(trans.nRow);
    out.SetArrayType(StringUtils::ConvertMBSToWCS(trans.chipType));
    out.SetIntensityCount(newNumCells);
    out.SetStdDevCount(newNumCells);
    out.SetPixelCount(newNumCells);
    CelFileWriter writer(out);

    for(int featIx = 0; featIx < numCells; featIx++) {
      if(trans.probeMask[featIx] >= 0) {
        intensities[trans.probeMask[featIx]] = cel.GetIntensity(featIx);
        stdevs[trans.probeMask[featIx]] = cel.GetStdv(featIx);
        nPixels[trans.probeMask[featIx]] = cel.GetPixels(featIx);
      }
    }
    writer.WriteStdDevs(stdevs);
    writer.WriteIntensities(intensities);
    writer.WritePixels(nPixels);
  }
  MultiChannelCelFileCollater collator;
  collator.Collate(tmpFileNamesToDelete, outCel);
   { for( int i=0; i < tmpFileNamesToDelete.size();  i++ ) { Fs::rm(tmpFileNamesToDelete[i], false);} } ;
}

void subsetCelFile(const string &orig, const string &outCel, 
                        const FileTranslation &trans) {
  FusionCELData cel;
  cel.SetFileName(orig.c_str());
  cel.Read();
  if(!cel.Read()) 
    Err::errAbort("\nCan't read cel file: " + cel.GetFileName() + 
                  "\n>>> Error reported: " + StringUtils::ConvertWCSToMBS(cel.GetError()));
  WStringVector channelVec = cel.GetChannels();
  cel.Close();
  if(channelVec.size() <= 1) {
    makeCelFileSubset(orig, outCel, trans);
  }
  else {
    makeMultiCelFileSubset(orig, outCel,trans);
  }
}

/// Update the bgp files with new probe ids for appropriate bgp entries
void convertBgpFile(vector<int> &probeMask, const string &chipType,
                    const string &bgpFileIn, const string &bgpFileOut ) {
  affx::BgpFile bgp;
  affx::BgpFile bgpOut;
  bgpOut.m_tsv.defineFile("probe_id\ttype\tgc_count\n");
  bgpOut.m_tsv.addHeader("chip_type",chipType);
  bgpOut.m_tsv.addHeader("lib_set_name",chipType);

  bgpOut.m_tsv.writeTsv_v1(bgpFileOut);
  if (bgp.open(bgpFileIn)!=TSV_OK) {
    Err::errAbort("Unable to open "+ToStr(bgpFileIn));
  }

  while (bgp.next_bgprobe()==TSV_OK) {
    int pId = probeMask[bgp.probe_id-1];  // NOTE: -1 because this probe_id is 1-based.
    if(pId < 0) {
      Err::errAbort("Can't have pId < 0");
    }
    bgpOut.m_tsv.set(0,"probe_id",pId);
    bgpOut.m_tsv.set(0, "type", bgp.type);
    bgpOut.m_tsv.set(0,"gc_count",bgp.gc_count);
    bgpOut.m_tsv.writeLevel(0);
  }
  bgpOut.close();
  bgp.close();
}

/// Set the bgp probes and probesets to be kept
void setBgpFile(FileTranslation &trans, string &bgpFile) {
  affx::BgpFile bgp;
  if (bgp.open(bgpFile)!=TSV_OK) {
    Err::errAbort("Unable to open "+ToStr(bgpFile));
  }

  while (bgp.next_bgprobe()==TSV_OK) {
    int pId = bgp.probe_id;
    trans.addProbe(pId);
    trans.keepProbeSets.insert(ToStr(bgp.probeset_id_1));
  }
  bgp.close();
}

/// Copy a tsv record from one tsv-file to another
void copyRow(TsvFile &in, TsvFile &out, int level, const FileTranslation &trans) {
  string s;
  int pId = -1;
  for(int cidx = 0; cidx < in.getColumnCount(level); cidx++) {
    if(in.get(level,cidx,s) == TSV_OK) {
      if(in.getColumnName(level, cidx) == trans.probeField) {
        pId = Convert::toInt(s) - trans.pOffset;
        pId = trans.probeMask[pId];
        if(pId < 0) {
          Err::errAbort("Shouldn't get -1 in the probeMask for pId: " + s + " at column: " + in.getColumnName(level, cidx) + " at " + ToStr(level) + " in: " + in.getFileName() + " with offset: " + ToStr(trans.pOffset) + " at line: " + ToStr(in.lineNumber()));
        }
        s = ToStr(pId + trans.pOffset);
      }
      else if((in.getColumnName(level, cidx) == "y" || in.getColumnName(level,cidx) == "Probe Y") && pId >= 0) {
        int y = pId / trans.nCol;
        s = ToStr(y);
      }
      else if((in.getColumnName(level, cidx) == "x" || in.getColumnName(level,cidx) == "Probe X") && pId >= 0) {
        int x = pId % trans.nCol;
        s = ToStr(x);
      }
      out.set(level, cidx, s);
    }
    else {
      Err::errAbort("Error - Bad read in: " + in.getFileName() + " at line number " + ToStr(in.lineNumber()) +  " at level: " + ToStr(level));
    }
  }
  out.writeLevel(level);
}

/// Recursively copy rows from one tsv-file to another
void recursiveWrite(TsvFile &in, TsvFile &out, int level, const FileTranslation &trans) {
  string s;
  while(in.nextLevel(level) == TSV_OK) {
    copyRow(in, out, level, trans);
    if(level < in.getLevelCount() -1) {
      recursiveWrite(in, out, level+1, trans);
    }
  }
}

/**
 * Copy all the headers with the exception of chip_type which is supplied
 */
void copyTsvHeaders(TsvFile &in, TsvFile &out, const map<string,string> &replace) {
  string key;
  string value;
  while(in.headersNext(key, value) == TSV_OK) {
    if(replace.find(key) != replace.end()) {
      value = replace.find(key)->second;
    }
    out.addHeader(key,value);
  }
}

/**
 * Copy over the "format" of one tsvfile to another
 */ 
void copyTsvDefinition(TsvFile &in, TsvFile &out) {
  out.copyFormat(in);
}

/// Based on the probesets and probes in the translation object determine if we're keeping this tsv row
bool keepRow(TsvFile &in, int level,
             const FileTranslation &trans) {
  string s;
  string pString;
  if(trans.probeSetField == "*") 
    return true;
  for(int cidx = 0; cidx < in.getColumnCount(level); cidx++) {
    if(in.get(level,cidx,s) == TSV_OK) {
      string colName = in.getColumnName(level, cidx);
      if(colName == trans.probeField) {
        pString = s;
        int pId = Convert::toInt(s);
        pId = trans.probeMask[pId-trans.pOffset];
        if(pId >= 0) {
          return true;
        }
      }
      if(colName == trans.probeSetField) {
        string s2;
        size_t x = s.rfind(":");
        if(x != string::npos) {
          s2 = s.substr(0,x);
        }
        if(trans.keepProbeSets.find(s) != trans.keepProbeSets.end() || trans.keepProbeSets.find(s2) != trans.keepProbeSets.end()) {
          return true;
        }
      }
    }
  }
  return false;
}

/// Based on the probesets and probes in the translation object determine if we're keeping this tsv row
bool keepRow(File5_Tsv &in, int level, File5Options &f5Opts, const FileTranslation &trans, int lineCount) {
  string s;
  string pString;
  string colName;
  int pId = 0;
  if(f5Opts.probeSetId == "*") 
    return true;
  for(int cidx = 0; cidx < in.getColumnCount(level); cidx++) {
    affx::File5_dtype_t type = in.getColumnDtype(level,cidx);
    if(type == FILE5_DTYPE_INT) {
      in.getColumnName(level, cidx, &colName);
      if(colName == f5Opts.probeId) {
        in.get(level, cidx, &pId);
        pId = trans.probeMask[pId - f5Opts.idOffset];
        if(pId >= 0) {
          return true;
        }
      }
    }
    else if(type == FILE5_DTYPE_STRING) {
      in.getColumnName(level, cidx, &colName);
      if(colName == f5Opts.probeSetId) {
        in.get(level,cidx,&s);
        string s2;
        size_t x = s.rfind(":");
        if(x != string::npos) {
          s2 = s.substr(0,x);
        }
        if(trans.keepProbeSets.find(s) != trans.keepProbeSets.end() || trans.keepProbeSets.find(s2) != trans.keepProbeSets.end()) {
          return true;
        }
      }
    }
    else if(type == FILE5_DTYPE_DOUBLE) {
      in.getColumnName(level, cidx, &colName);
      if(colName.find('X') == 0) {
        return trans.magicOrderTranslation[lineCount] >= 0;
      }
    }
  }
  return false;
}

int adjustColumn(File5_Tsv &in, string &colName, int value, File5Options &f5Opts, const FileTranslation &trans) {
  int adjusted = value;
  if(colName == f5Opts.probeId) {
    adjusted = value - f5Opts.idOffset;
    adjusted = trans.probeMask[adjusted];
    if(adjusted < 0) {
      Err::errAbort("Shouldn't get -1 in the probeMask for pId: " + ToStr(value) + " at column: " + colName + " with offset: " + ToStr(trans.pOffset));
    }
    adjusted = adjusted + f5Opts.idOffset;
  }
  else if((colName == "y" || colName == "Probe Y") && adjusted >= 0) {
    adjusted = adjusted / trans.nCol;
    
  }
  else if((colName == "x" || colName == "Probe X") && adjusted >= 0) {
    adjusted = adjusted % trans.nCol;
  }
  return adjusted;
}

void checkRet(int retVal, int level, int col, File5_Tsv &tsv5) {
  if(retVal != FILE5_OK) {
    string colName;
    tsv5.getColumnName(level, col, &colName);
    Err::errAbort("Error: bad return value at column: " + colName + " level " + ToStr(level) + " column: " + ToStr(col));
  }
}

void copyColumn(File5_Tsv &in, File5_Tsv &out, File5Options &f5Opts,
                int level, int colIx, const FileTranslation &trans) {

  affx::File5_dtype_t type = in.getColumnDtype(level,colIx);
  int ival = 0;
  float fval = 0.0f;
  double dval = 0.0f;
  string sval;
  string colName;
  int ret = 0;
  switch (type) {
    
  case FILE5_DTYPE_INT :
    colName = "";
    ret = in.getColumnName(level, colIx, &colName);
    checkRet(ret, level, colIx, in);
    ret = in.get(level,colIx,&ival);
    checkRet(ret, level, colIx, in);
    ival = adjustColumn(in, colName, ival, f5Opts, trans);
    if(colName == "ProbeSetIndex") {
      int psIndex = trans.magicOrderTranslation[ival];
      ret = out.set_i(level,colIx,psIndex);
    }
    else 
      ret = out.set_i(level,colIx,ival);
    checkRet(ret, level, colIx, out);
    break;
  case FILE5_DTYPE_FLOAT :
    ret = in.get(level,colIx,&fval);
    checkRet(ret, level, colIx, in);
    ret = out.set_f(level,colIx,fval);
    checkRet(ret, level, colIx, out);
    break;
  case FILE5_DTYPE_DOUBLE :
    ret = in.get(level,colIx,&dval);
    checkRet(ret, level, colIx, in);
    ret = out.set_d(level,colIx,dval);
    checkRet(ret, level, colIx, out);
    break;
  case FILE5_DTYPE_STRING : 
    ret = in.get(level,colIx,&sval);
    checkRet(ret, level, colIx, in);
    ret = out.set_string(level,colIx,sval);
    checkRet(ret, level, colIx, out);
    break;
  case FILE5_DTYPE_UNKNOWN:
  case FILE5_DTYPE_ERR:
  case FILE5_DTYPE_ANY:
  default:
    Err::errAbort("Unknown Type - Unhandled case.");

  }
}

void copyRow(File5_Tsv &in, File5_Tsv &out, File5Options &f5Opts, int level, const FileTranslation &trans) {
  for(int cidx = 0; cidx < in.getColumnCount(level); cidx++) {
    copyColumn(in, out, f5Opts, level, cidx, trans);
  }
  out.writeLevel(level);
}

/// Recursively copy rows from one tsv-file to another
void recursiveWrite(File5_Tsv &in, File5_Tsv &out, File5Options &f5Opts, int level, const FileTranslation &trans) {
  string s;
  while(in.nextLevel(level) == FILE5_OK) {
    copyRow(in, out, f5Opts, level, trans);
    if(level < in.getLevelCount() -1) {
      recursiveWrite(in, out, f5Opts, level+1, trans);
    }
  }
}

/// Copy the appropriate entries from the input tsv file into the output tsv file.
void copyTsvFile(const string &in, const string &out, const map<string,string> &replace,
                 const FileTranslation &trans) {
  TsvFile iTsv;
  TsvFile oTsv;
  iTsv.open(in);
  copyTsvDefinition(iTsv, oTsv);
  copyTsvHeaders(iTsv, oTsv, replace);
  if(iTsv.getLevelCount() == 1) {
    oTsv.writeTsv_v1(out);
  }
  else {
    oTsv.writeTsv(out);
  }
  int level = 0;
  int seen = 0, kept = 0;
  while(iTsv.nextLevel(level) == TSV_OK) {
    seen++;
    if(keepRow(iTsv, level, trans)) {
      kept++;
      copyRow(iTsv, oTsv, level, trans);
      recursiveWrite(iTsv, oTsv,level+1, trans);
    }
  }
  iTsv.close();
  oTsv.close();
  Verbose::out(2, "Copy tsv: Saw " + ToStr(seen) + " and kept " + ToStr(kept) + " in file " + in);
}

/// Copy the appropriate entries from the input tsv file into the output tsv file.
void copyFile5File(File5Options f5Opts, const string &out,
                   const FileTranslation &trans) {
  int ret;
  affx::File5_File iFile5;	
  Verbose::out(2,"Opening " + f5Opts.fileName + " with group " + f5Opts.groupName + " to read.");
  ret = iFile5.open(f5Opts.fileName, affx::FILE5_OPEN_RO);
  if(ret != FILE5_OK) {
    Err::errAbort("Couldn't open file: " + f5Opts.fileName);
  }
  Verbose::out(2,"Opening " + out + " with group " + f5Opts.groupName + " to write.");
  affx::File5_File oFile5;	
  ret = oFile5.open(out, affx::FILE5_OPEN_CREATE);
  if(ret != FILE5_OK) {
    Err::errAbort("Couldn't open file: " + out);
  }
  
  affx::File5_Group* iGroup5 = iFile5.openGroup(f5Opts.groupName, affx::FILE5_OPEN);
  affx::File5_Group* oGroup5 = oFile5.openGroup(f5Opts.groupName, affx::FILE5_OPEN_CREATE);
  affx::File5_Tsv* iTsv5 = iGroup5->openTsv(f5Opts.groupName, affx::FILE5_OPEN);
  affx::File5_Tsv* oTsv5 = oGroup5->openTsv(f5Opts.groupName, affx::FILE5_OPEN_CREATE);

  oTsv5->copyFormat(*iTsv5);
  int level = 0;
  int seenCount = 0;
  int keepCount = 0;
  while(iTsv5->nextLevel(level) == affx::FILE5_OK) {

    // @todo - fix this hack to figure out what ProbeSetIndex the values for reference file are...
    if(keepRow(*iTsv5, level, f5Opts, trans, seenCount++)) {
      keepCount++;
      copyRow(*iTsv5, *oTsv5, f5Opts, level, trans);
    }
  }
  Verbose::out(2, "Saw " + ToStr(seenCount) + " and kept " + ToStr(keepCount) + " records.");
  iTsv5->close();
  oTsv5->close();
  delete iTsv5;
  delete oTsv5;

  iGroup5->close();
  oGroup5->close();
  delete iGroup5;
  delete oGroup5;

  iFile5.close();
  oFile5.close();
}

/// Read a tsv row and see if it has a probe or probeset we want to keep
bool readTsvProbesetRow(TsvFile &in, int level, FileTranslation &trans, set<string> &probeSets) {
  string s;
  bool ret = false;
  for(int cidx = 0; cidx < in.getColumnCount(level); cidx++) {
    if(in.get(level,cidx,s) == TSV_OK) {
      string colName = in.getColumnName(level, cidx);
      if(colName == trans.probeSetField) {
        probeSets.insert(s);
        ret = true;
      }
      else if(colName == trans.probeField) {
        int pId = Convert::toInt(s);
        trans.addProbe(pId);
        ret = true;
      }
    }
  }
  return ret;
}

/// Recursively loop through tsv entries looking for probes and probesets to keep
void recursiveRead(TsvFile &in, int level, FileTranslation &trans, set<string> &probeSets) {
  string s;
  while(in.nextLevel(level) == TSV_OK) {
    readTsvProbesetRow(in, level, trans, probeSets);
    if(level < in.getLevelCount() -1) {
      recursiveRead(in, level+1, trans, probeSets);
    }
  }
}

//
void readProbesetTsvFile(const string &fileName, FileTranslation &trans, int numKeep) {
  TsvFile in;
  in.open(fileName);
  int level = 0;
  set<string> probeSets;
  int count = 0;
  int kept = 0;
  while(in.nextLevel(level) == TSV_OK) {
    if(readTsvProbesetRow(in, level, trans, probeSets)) {
      kept++;
    }
    recursiveRead(in, level+1, trans, probeSets);
    count++;
  }
  in.close();
  Verbose::out(2, "Saw " + ToStr(count) + " and kept " + ToStr(kept) + " in file " + fileName);
  vector<string> toAdd(min((size_t)numKeep,probeSets.size()));
  srand(1);
  random_sample(probeSets.begin(), probeSets.end(), toAdd.begin(), toAdd.end());
  sort(toAdd.begin(), toAdd.end());
  for(vector<string>::iterator iter = toAdd.begin(); iter != toAdd.end(); iter++) {
    trans.keepProbeSets.insert(*iter);
  }
}

void readFile5ProbesetRow(File5_Tsv &iTsv5, File5Options &f5Opts, int level, FileTranslation &trans,
                          int &probeId, int &probesetId) {
  int nothingCount;
  for(int cidx = 0; cidx < iTsv5.getColumnCount(level); cidx++) {
    affx::File5_dtype_t type = iTsv5.getColumnDtype(level,cidx);
    int ival = 0;
    string sval;
    string colName;
    int ret = 0;
    switch (type) {
    case FILE5_DTYPE_INT :
      colName = "";
      ret = iTsv5.getColumnName(level, cidx, &colName);
      checkRet(ret, level, cidx, iTsv5);
      ret = iTsv5.get(level,cidx,&ival);
      checkRet(ret, level, cidx, iTsv5);
      if(colName == f5Opts.probeId) {
        if(trans.addProbe(ival, f5Opts.idOffset))
          probeId++;
      }
      break;
    case FILE5_DTYPE_STRING : 
      colName = "";

      ret = iTsv5.getColumnName(level, cidx, &colName);
      checkRet(ret, level, cidx, iTsv5);

      ret = iTsv5.get(level,cidx,&sval);
      checkRet(ret, level, cidx, iTsv5);
      
      if(colName == f5Opts.probeSetId) {
        pair<set<string>::iterator,bool> rVal = trans.keepProbeSets.insert(sval);
        if(rVal.second) {
          probesetId++;
        }
      }

      break;
    case FILE5_DTYPE_UNKNOWN:
    case FILE5_DTYPE_ERR:
    case FILE5_DTYPE_ANY:
    default:
      // do nothing
      nothingCount++;
    }
  }
}

//
void readProbesetFile5File(File5Options &f5Opts, FileTranslation &trans) {
  int ret;
  affx::File5_File iFile5;	
  Verbose::out(2,"Opening " + f5Opts.fileName + " with group " + f5Opts.groupName + " to read.");
  ret = iFile5.open(f5Opts.fileName, affx::FILE5_OPEN_RO);
  if(ret != FILE5_OK) {
    Err::errAbort("Couldn't open file: " + f5Opts.fileName);
  }
  affx::File5_Group* iGroup5 = iFile5.openGroup(f5Opts.groupName, affx::FILE5_OPEN);
  affx::File5_Tsv* iTsv5 = iGroup5->openTsv(f5Opts.groupName, affx::FILE5_OPEN);
  int level = 0;
  set<string> probeSets;
  int probesetCount = 0, probeCount = 0, seenCount = 0;
  while(iTsv5->nextLevel(level) == affx::FILE5_OK) {
    readFile5ProbesetRow(*iTsv5, f5Opts, level, trans, probeCount, probesetCount);
    seenCount++;
  }
  Verbose::out(2, "Saw " + ToStr(seenCount) + " and kept " + ToStr(probeCount) + " probe records and " + ToStr(probesetCount) + " probeset records");
  iTsv5->close();
  delete iTsv5;
  iGroup5->close();
  delete iGroup5;
  
  iFile5.close();
}

void readMetaProbesetList(vector<Meta> &toLoad, const string &metaFile, int numLoad){
  string probeset, probesetList;
  vector<Meta> sets;
  vector<string> words;
  TsvFile tsv;
  // open file and make sure that it has at least a probeset_id column
  if(tsv.open(metaFile) != TSV_OK) 
    Err::errAbort("Couldn't open file: '" + ToStr(metaFile) + "'");
  tsv.bind(0, "probeset_id", &probeset, TSV_BIND_REQUIRED);
  // list of probesets can be either called probeset_list or 
  if(tsv.cname2cidx(0, "probeset_list") != TSV_ERR_NOTFOUND)
    tsv.bind(0, "probeset_list", &probesetList, TSV_BIND_REQUIRED);
  else if(tsv.cname2cidx(0, "probeset_ids") != TSV_ERR_NOTFOUND)
    tsv.bind(0, "probeset_ids", &probesetList, TSV_BIND_REQUIRED);
  else
    Err::errAbort("Meta probeset list must contain 'probeset_list'");

  set<string> seen;
  while(tsv.nextLevel(0) == TSV_OK) {
    Meta m;
    m.id = probeset;
    Util::chopString(probesetList, ' ', m.probeSets);
    sets.push_back(m);
  }
    
  tsv.close();
  toLoad.clear();
  toLoad.resize(numLoad);
  srand(1);
  random_sample(sets.begin(), sets.end(), toLoad.begin(), toLoad.end());
  sort(toLoad.begin(), toLoad.end());
}

void fillPsNamesFromLayout(set<string> &keepProbeSets, ChipLayout &layout, int numKeep) {
  // - For each probsets to be kept mark the probes to be kept
  vector<string> psNames(layout.getProbeSetCount());
  vector<string> psSample(numKeep);
  fill(psSample.begin(), psSample.end(), -1);
  for(int i = 0; i < layout.getProbeSetCount(); i++) 
    psNames[i] = layout.getProbeListAtIndex(i).get_name_string();
  srand(1);
  random_sample(psNames.begin(), psNames.end(), psSample.begin(), psSample.end());
  //  sort(psSample.begin(), psSample.end());
  for(int i = 0; i < psSample.size(); i++) {
    keepProbeSets.insert(psSample[i]);
  }

}

void writeMetaProbesets(vector<Meta> &metaSets, const string &fileName, const string &chipType) {
  TsvFile tsv;
  tsv.defineFile("probeset_id\tprobeset_list\n");
  tsv.addHeader("chip_type", chipType);
  tsv.addHeader("lib_set_name", chipType);
  tsv.writeTsv_v1(fileName);
  for(int i = 0; i < metaSets.size(); i++) {
    tsv.set(0,"probeset_id",metaSets[i].id);
    string list = metaSets[i].probeSets[0];
    for(int j = 1; j < metaSets[i].probeSets.size(); j++) {
      list += " " + metaSets[i].probeSets[j];
    }
    tsv.set(0,"probeset_list", list);
    tsv.writeLevel(0);
  }
  tsv.close();
}

void writeSpfFile(ChipLayout &layout, 
                  const string &outputSpf, 
                  const FileTranslation &trans,
                  int version, int nChannels) {
  SpfFile outSpf;
  openSpfForWrite(outSpf, outputSpf, trans.chipType, version, trans.nCol, trans.nRow, trans.keepProbeSets.size(), nChannels);
  int written = 0;
  set<string> working;
  for(int i = 0; i < layout.getProbeSetCount(); i++) {
    ProbeListPacked pl = layout.getProbeListAtIndex(i);
    if( trans.keepProbeSets.find(pl.get_name_string()) != trans.keepProbeSets.end()) {
      for(int j = 0; j < pl.probe_cnt(); j++) {
        int pId = pl.get_probeId(j);
        pl.set_probeId(j,trans.probeMask[pId]);
      }
      if(version == 4) {
        ProbeListFactory::writeProbeListToSpfFile_v4(outSpf, pl);
      }
      else {
        ChipLayout::writeSpfProbeListRecord(outSpf, pl);
      }
      working.insert(pl.get_name_string());
      written++;
    }
  }
  Verbose::out(2, "Wrote " + ToStr(written) + " probesets");
  outSpf.close();
  for(set<string>::iterator iter = trans.keepProbeSets.begin(); iter != trans.keepProbeSets.end(); iter++) {
    if(working.find(*iter) == working.end()) {
      Verbose::out(2,"Couldn't find: " + *iter);
    }
  }
}

void openChipLayout(ChipLayout &layout, string &pgfFileStr, string &clfFileStr, string &spfFileStr) {
  if(spfFileStr != "") {
    Verbose::out(1, "Reading probesets from spf file: " + spfFileStr);
    layout.openSpfAll(spfFileStr);
  }
  else if(pgfFileStr != "") {
    if(clfFileStr == "") {
      Err::errAbort("If specifying --pgf-file must also specifying --clf-file");
    }
    ClfFile clfFile;
    clfFile.open(clfFileStr);
    int xSize = clfFile.getXMax()+1;
    int ySize = clfFile.getYMax()+1;
    Verbose::out(1, "Reading probesets from pgf file: " + pgfFileStr);
    layout.openPgfAll(pgfFileStr, xSize, ySize);
  }
  else {
    Err::errAbort("Must specify specify pgf file or clf file");
  }
}

string makeOutputName(const string &dir, const string &file, const string &prefix, const string &suffix) {
  string s;
  if(prefix == "") {
    s = Fs::join(dir,Util::chopSuffix(Fs::basename(file)) + "." + suffix);
  }
  else {
    s = Fs::join(dir,Util::chopSuffix(Fs::basename(file)) + "." + prefix + "." + suffix);
  }
  return s;
}

void readProbesetList(set<string> &toLoad, const string &psFile) {
  string probeset;
  toLoad.clear();
  TsvFile tsv;
  // open file and make sure that it has at least a probeset_id column
  if(tsv.open(psFile) != TSV_OK) 
    Err::errAbort("Couldn't open file: '" + ToStr(psFile) + "'");
  tsv.bind(0, "probeset_name", &probeset, TSV_BIND_REQUIRED);

  set<string> seen;
  while(tsv.nextLevel(0) == TSV_OK) {
    toLoad.insert(probeset);
  }
    
  tsv.close();
}

void fillInProbesets(const string &metaFile, vector<Meta> &metaSets,
                     ChipLayout &layout, FileTranslation &trans, DsOptions &o) {
  if(o.probesetFileStr != "") {
    readProbesetList(trans.keepProbeSets, o.probesetFileStr);
  }
  else if(metaFile.size() > 0) {
    Verbose::out(2,"Loading meta probset file: " + o.metaFile);
    readMetaProbesetList(metaSets, o.metaFile, o.numKeep);
    for(int i = 0; i < metaSets.size(); i++) {
      for(int j = 0; j < metaSets[i].probeSets.size(); j++) {
        trans.keepProbeSets.insert(metaSets[i].probeSets[j]);
      }
    }
  }
  else {
    Verbose::out(2, "Selecting " + ToStr(o.numKeep) + " probesets from layout file");
    fillPsNamesFromLayout(trans.keepProbeSets, layout, o.numKeep);
  }
}

void loadMagicOrder(FileTranslation &trans, string &col, string &file, string &group) {
  
  affx::File5_File iFile5;
  Verbose::out(2,"Opening " + file + " with group " + group + " to read.");
  int ret = 0;
  ret = iFile5.open(file, affx::FILE5_OPEN_RO);
  if(ret != FILE5_OK) {
    Err::errAbort("Couldn't open file: " + file);
  }
  affx::File5_Group* iGroup5 = iFile5.openGroup(group, affx::FILE5_OPEN);
  affx::File5_Tsv* iTsv5 = iGroup5->openTsv(group, affx::FILE5_OPEN);
  trans.magicOrder.clear();
  trans.magicOrderTranslation.clear();
  trans.magicOrderMap.clear();
  int cIdx = -1;
  int level = 0;
  int seenCount = 0;
  int convertedCount = 0;
  while(iTsv5->nextLevel(level) == FILE5_OK) {
    if(cIdx < 0) {
      for(int i = 0; i < iTsv5->getColumnCount(level); i++) {
        string colName;
        iTsv5->getColumnName(level, i, &colName);
        if(colName == col) {
          cIdx = i;
          break;
        }
      }
      if(cIdx < 0) {
        Err::errAbort("Index column not found.");
      }
      else {
        Verbose::out(2, "Magic is index: " + ToStr(cIdx));
      }
    }
    string val;
    iTsv5->get(level, cIdx, &val);
    if(trans.keepProbeSets.find(val) != trans.keepProbeSets.end()) {
      trans.magicOrderTranslation.push_back(convertedCount);
      trans.magicOrderMap[val] = convertedCount;
      trans.magicOrder.push_back(val);
      convertedCount++;
    }
    else {
      trans.magicOrderTranslation.push_back(-1);
    }
    seenCount++;
  }
  Verbose::out(2, "Magic - Saw: " + ToStr(seenCount) + " converted: " + ToStr(convertedCount));
  iTsv5->close();
  delete iTsv5;
  iGroup5->close();
  delete iGroup5;
  
  iFile5.close();
}

void dataSubset(PgOptions *opts) {

  DsOptions o(opts);

  if ( !Fs::dirExists(o.libOutputPrefix) ) {
    Fs::mkdirPath(o.libOutputPrefix, false);
  }
  if ( !Fs::dirExists(o.celOutputPrefix) ) {
    Fs::mkdirPath(o.celOutputPrefix, false);
  }
  int numChannels = opts->getInt("num-channels");
  int spfVersion = opts->getInt("spf-version");
  vector<Meta> metaSets;
  FileTranslation trans;
  trans.probeSetField = opts->get("probeset-field");
  trans.probeField = "probe_id";
  trans.chipType = Util::chopSuffix(Fs::basename(o.chipRoot));
  if(o.filePrefix != "") 
    trans.chipType += "." + o.filePrefix;

  // - Make vector for size of chip, init all to -1
  ChipLayout layout;
  openChipLayout(layout,o.pgfFileStr, o.clfFileStr, o.spfFileStr);

  trans.probeMask.resize(layout.getXCount() * layout.getYCount());
  fill(trans.probeMask.begin(),trans.probeMask.end(),-1);

  /* --- Determine probesets to be kept. --- */

  // First probesets via layout or meta-probesets
  if(o.probesetFileStr != "") {
    readProbesetList(trans.keepProbeSets, o.probesetFileStr);
    o.numKeep = trans.keepProbeSets.size();
  }
  else {
    if(o.numKeep > layout.getProbeSetCount()) 
      Err::errAbort("Asked for " + ToStr(o.numKeep) + " but only " + ToStr(layout.getProbeSetCount()) + " probesets.");
    fillInProbesets(o.metaFile, metaSets, layout, trans, o);
  }



  // Add in any qcc probesets
  if(o.qccFileStr != "") {
    Verbose::out(2, "Selecting " + ToStr(o.numKeep) + " probesets from qcc file");
    readProbesetTsvFile(o.qccFileStr, trans, o.numKeep);
  }

  if(o.file5SelectFileFile != "") {
    vector<File5Options> file5Conv;
    fillInFile5Files(file5Conv, o.file5SelectFileFile);
    for(int i = 0; i < file5Conv.size(); i++) {
      readProbesetFile5File(file5Conv[i], trans);
    }
  }
  
  if(o.tsvSelectFileFile != "") {
    vector<string> tsvFiles;
    vector<string> ids;
    fillInTsvFiles(tsvFiles,ids, o.tsvSelectFileFile);
    for(int i = 0; i < tsvFiles.size(); i++) {
      trans.pOffset = 1;
      string id = trans.probeSetField;
      if(ids[i] != "")
        trans.probeSetField = ids[i];
      readProbesetTsvFile(tsvFiles[i], trans, o.numKeep);
      trans.probeSetField = id;
    }
  }
  

  // - Mark the bgp probes as needing to be kept as well
  if(o.bgpFile.size() > 0) {
    Verbose::out(1, "Reading bgp file: " + o.bgpFile);
    setBgpFile(trans, o.bgpFile);
  }

  // - Make a mask mapping old probes to new ones from list of probesets
  trans.fillInFromProbeSets(layout);
  
  // - Determine new size of cel file
  trans.calcDimensions(trans.totalProbes);

  // - Figure out magic orderings for copynumber if desired
  if(o.magicOrderColStr != "") {
    loadMagicOrder(trans, o.magicOrderColStr, o.magicOrderFileStr, o.magicOrderGroupStr);
  }

  /*  ---- For each library file write a smaller one ---- */

  // Chip Layout file, either pgf or spf
  if(o.spfFileStr != "") {
    string outputSpf = makeOutputName(o.libOutputPrefix, o.spfFileStr, o.filePrefix, "spf");
    Verbose::out(1,"Writing new spf file: " + outputSpf);
    writeSpfFile(layout, outputSpf, trans, spfVersion, numChannels);
  }
  if(o.pgfFileStr != "") {
    string pgfFileOut = makeOutputName(o.libOutputPrefix, o.pgfFileStr, o.filePrefix, "pgf");
    string clfFileOut = makeOutputName(o.libOutputPrefix, o.clfFileStr, o.filePrefix, "clf");
    map<string,string> pgfHeaders;
    pgfHeaders.insert(pair<string,string>("chip_type", trans.chipType));
    pgfHeaders.insert(pair<string,string>("rows", ToStr(trans.nRow)));
    pgfHeaders.insert(pair<string,string>("cols", ToStr(trans.nCol)));
    pgfHeaders.insert(pair<string,string>("probesets", ToStr(trans.keepProbeSets.size())));
    trans.pOffset = 1;
    Verbose::out(1,"Writing new pgf file: " + pgfFileOut);
    copyTsvFile(o.pgfFileStr, pgfFileOut, pgfHeaders, trans);
    Verbose::out(1,"Writing new clf file: " + clfFileOut);
    copyTsvFile(o.clfFileStr, clfFileOut, pgfHeaders, trans);

  }

  // Write out meta sets if set
  if(metaSets.size() > 0) {
    string metaOutput =  makeOutputName(o.libOutputPrefix, o.metaFile, o.filePrefix, "mps");
    Verbose::out(1,"Writing new meta-probeset file: " + metaOutput);
    writeMetaProbesets(metaSets, metaOutput, trans.chipType);
  }

  // Write a new qcc file if supplied
  if(o.qccFileStr != "") {
    string qccOutput = makeOutputName(o.libOutputPrefix, o.qccFileStr,  o.filePrefix, "qcc");
    map<string,string> qccHeaders;
    qccHeaders.insert(pair<string,string>("chip_type", trans.chipType));
    Verbose::out(1, "Writing new qcc file: " + qccOutput);
    trans.pOffset = 0;
    copyTsvFile(o.qccFileStr, qccOutput, qccHeaders, trans);
  }

  // - Write new bgp file if supplied
  if(o.bgpFile.size() > 0) {
    string bgpOutput = makeOutputName(o.libOutputPrefix, o.bgpFile, o.filePrefix, "bgp");
    map<string,string> bgpHeaders;
    bgpHeaders.insert(pair<string,string>("chip_type", trans.chipType));
    Verbose::out(1, "Writing new bgp file: " + bgpOutput);
    trans.pOffset = 1;
    copyTsvFile(o.bgpFile, bgpOutput, bgpHeaders, trans);
  }
  // - Write new kill list file if supplied
  if(o.killFileStr.size() > 0) {
    string killListOutput =  makeOutputName(o.libOutputPrefix, o.killFileStr, o.filePrefix, "txt");
    map<string,string> killHeaders;
    Verbose::out(1, "Writing new kill-list file: " + killListOutput);
    trans.pOffset = 1;
    copyTsvFile(o.killFileStr, killListOutput, killHeaders, trans);
  }

  if(o.file5ConvertFileFile != "") {
    vector<File5Options> file5Conv;
    fillInFile5Files(file5Conv, o.file5ConvertFileFile);
    for(int i = 0; i < file5Conv.size(); i++) {
      string suffix = file5Conv[i].fileName.substr(file5Conv[i].fileName.rfind(".")+1,file5Conv[i].fileName.size());
      string outputName = makeOutputName(o.libOutputPrefix, file5Conv[i].fileName, o.filePrefix, suffix);
      Verbose::out(2, "Starting processing of file: " + file5Conv[i].fileName + " group: " + file5Conv[i].groupName);
      copyFile5File(file5Conv[i], outputName, trans);
      Verbose::out(2, "Finished processing file: " + file5Conv[i].fileName + " group: " + file5Conv[i].groupName);
    }
  }

  if(o.tsvConvertFileFile != "") {
    vector<string> tsvFiles;
    vector<string> ids;
    fillInTsvFiles(tsvFiles, ids, o.tsvConvertFileFile);
    for(int i = 0; i < tsvFiles.size(); i++) {
      string suffix = tsvFiles[i].substr(tsvFiles[i].rfind(".")+1,tsvFiles[i].size());
      string outputName = makeOutputName(o.libOutputPrefix, tsvFiles[i], o.filePrefix, suffix);
      Verbose::out(2, "Writing from: " + Fs::basename(tsvFiles[i]) + " into: " + outputName);
      trans.pOffset = 1;
      map<string,string> pgfHeaders;
      pgfHeaders.insert(pair<string,string>("chip_type", trans.chipType));
      pgfHeaders.insert(pair<string,string>("rows", ToStr(trans.nRow)));
      pgfHeaders.insert(pair<string,string>("cols", ToStr(trans.nCol)));
      pgfHeaders.insert(pair<string,string>("probesets", ToStr(trans.keepProbeSets.size())));
      string id = trans.probeSetField;
      if(ids[i] != "")
        trans.probeSetField = ids[i];
      copyTsvFile(tsvFiles[i], outputName, pgfHeaders, trans);
      trans.probeSetField = id;
    }
  }

  // - For each cel file write a smaller one
  vector<string> celFiles;
  fillInCelFiles(celFiles, o.celFileFile);
  for(int i = 0; i < celFiles.size(); i++) {
    // string rootName = Util::chopSuffix(Fs::basename(celFiles[i]));
    string suffix = celFiles[i].substr(celFiles[i].rfind(".")+1,celFiles[i].size());
    string outputName = makeOutputName(o.celOutputPrefix, celFiles[i], o.filePrefix, suffix);
    Verbose::out(2, "Writing from: " + Fs::basename(celFiles[i]) + " into: " + outputName);
    //makeCelFileSubset(celFiles[i], outputName, trans);
    subsetCelFile(celFiles[i], outputName, trans);
  }
  Verbose::out(1, "Done.");
}

/** Everybody's favorite function. */
int main(int argc, const char *argv[]) {
  const string version = AptVersionInfo::versionToReport();

  PgOptions *opts = NULL;
  opts = new PgOptions();
  defineDataSubsetOptions(opts);
  opts->parseArgv(argv);
  Verbose::setLevel(opts->getInt("verbose"));
  // Do we need help?
  if(opts->getBool("help") || argc == 1) {
    opts->usage();
    cout << "version: " << version << endl;
  }
  else if(opts->getBool("version"))
    cout << "version: " << version << endl;
  else
    dataSubset(opts);
  delete opts;
  return 0;
}

