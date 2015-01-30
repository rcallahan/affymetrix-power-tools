////////////////////////////////////////////////////////////////
//
// Copyright (C) 2013 Affymetrix, Inc.
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

#include "chipstream/apt-cel-transformer/CELTransformerEngine.h"
#include "calvin_files/converters/cel/src/CELFileVersion.h"
#include "calvin_files/data/src/CELData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/parsers/src/CelFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCelFileWriter.h"
#include "calvin_files/writers/src/MultiChannelCelFileCollater.h"
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/CelReader.h"
#include "chipstream/DiskIntensityMart.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/SparseMart.h"
#include "file/CELFileWriter.h"
#include "file/TsvFile/ClfFile.h"
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/PgOptions.h"
#include "util/TmpFileFactory.h"
//
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace affx;
using namespace affymetrix_calvin_exceptions;

CELTransformerEngine::Reg CELTransformerEngine::reg;

CELTransformerEngine * CELTransformerEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == CELTransformerEngine::EngineName())
		return (CELTransformerEngine *)engine;
	return NULL;
}

CELTransformerEngine::CELTransformerEngine()
{
	defineOptions();
}

CELTransformerEngine::~CELTransformerEngine()
{
}

void CELTransformerEngine::defineOptions()
{
	setUsage("apt-cel-transformer - Utility program for taking cel files and converting them\n"
		"using a ChipStream path and saving the results to a modified set of cel files.\n"
		"An example would include rma background subtracting and normalizing a set\n"
		"of cel files that could then be used in another application without those\n"
		"capabilities.\n"
		"\n"
		" *Warning* - If the output directory is the same as the input directory the\n"
		"the original cel files will be overwritten with the transformed ones.\n"
		"\n"
		"A chipstream pathway is specified by a comma separated list of transformations\n"
		"with specific parameters passed as key value pairs. For example:\n"
		"rma-bg,quant-norm.usepm=true.sketch=50000 -> rma background subtract and then\n"
		"                                             quantile normalize using only PM\n"
		"                                             probes and a sketch size of 50K\n"
		"rma-bg,med-norm -> rma background subtraction and then median normalization.\n"
		"\n"
		"usage:\n"
		"   apt-cel-transformer -d chip.cdf -c rma-bg,quant-norm.sketch=50000 -o out-dir \\ \n"
		"      celFile1.cel celFile2.cel ... celFileN.cel\n");


	defineOption("", "target-sketch", PgOpt::STRING_OPT,
		"File specifying a target distribution to use for quantile normalization.",
		"");
	defineOption("", "a5-target-sketch", PgOpt::STRING_OPT,
		"File specifying an target distribution to use for quantile normalization in A5 format.",
		"");
	defineOption("c", "chipstream", PgOpt::STRING_OPT,
		"String representing tranformation desired. "
		"For example: 'quant-norm.sketch=50000' does a quantile normalization using 50000 points.",
		"");
	defineOption("", "probe-norm-file", PgOpt::STRING_OPT,
		"File specifying index of probes (1 based) to use for normalization.",
		"");
	defineOption("", "use-disk", PgOpt::BOOL_OPT,
		"Store CEL intensities to be analyzed on disk.", "true");
	defineOption("", "disk-cache", PgOpt::INT_OPT,
		"Size of intensity memory cache in millions of intensities (when --use-disk=true).",
		"50");

	defineOption("", "cel-files", PgOpt::STRING_OPT,
		"Text file specifying cel files to process, one per line with the first line being 'cel_files'.",
		"");
	defOptMult("", "cels", PgOpt::STRING_OPT,
        "Cel files to process.",
        "");
	defineOption("d", "cdf-file", PgOpt::STRING_OPT,
		"File defining probe sets",
		"");
	defineOption("", "clf-file", PgOpt::STRING_OPT,
		"File defining probes",
		"");
	defineOption("", "kill-list", PgOpt::STRING_OPT,
		"File defining probes to remove",
		"");
	defineOption("", "pgf-file", PgOpt::STRING_OPT,
		"File defining probe sets",
		"");
	defineOption("", "spf-file", PgOpt::STRING_OPT,
		"File defining probe sets",
		"");
	defineOption("", "write-sketch", PgOpt::BOOL_OPT,
		"Write the quantile normalization distribution (or sketch) "
		"to a file in output-dir for reuse with target-sketch option.",
		"false");
	defineOption("", "a5-write-sketch", PgOpt::BOOL_OPT,
		"Write the quantile normalization distribution (or sketch) "
		"to a file in output-dir for reuse with target-sketch option.",
		"false");
	defineOption("", "update-header", PgOpt::BOOL_OPT, 
		"Change guid and add apt-cel-transform meta info.", "true");
}

/** 
* Figure out chip dimensions from the cel file specified.
* 
* @param fileName - Name of cel file to read header from.
* @param numCols - Number of x columns.
* @param numRows - Number of y rows.
*/
static void getDimensions(string fileName, colrow_t& numCols, colrow_t& numRows) {
	FusionCELData cel;
	cel.SetFileName(fileName.c_str());
	if (!cel.ReadHeader()) 
		Err::errAbort("Can't read cel file: " + cel.GetFileName());
	numCols = cel.GetCols();
	numRows = cel.GetRows();
	cel.Close();
}

/**
* Transform AGCC CEL Files
*/
void CELTransformerEngine::transformNew(int celIx, vector<ChipStream *> &cStreamVec) {

	string outDir = getOpt("out-dir");
	string outCelFileName = Fs::join(outDir , Fs::basename(celFiles[celIx]));
	vector<string> tmpFileNames;
	GlobalTmpFileFactory()->setTmpdir(outDir);

	try {
		// read agcc cel file
		affymetrix_calvin_io::CelFileReader reader;
		reader.SetFilename(celFiles[celIx]);
		affymetrix_calvin_io::CelFileData cel;
		reader.Read(cel);
		int numCells;
		WStringVector channelVec = cel.GetChannels();
		for (int channelIx = 0; channelIx < channelVec.size(); ++channelIx) {
			cel.SetActiveChannel(channelVec[channelIx]);
			numCells = cel.GetNumCells();
			vector<float> intensities;
			if (!cel.GetIntensities(0,numCells,intensities))
				Err::errAbort("Unable to retrieve intensities for " + celFiles[celIx]);
			vector<float> stdev;
			if (!cel.GetStdev(0,numCells,stdev))
				Err::errAbort("Unable to retrieve stdev for " + celFiles[celIx]);
			vector<int16_t>  npixels;
			if (!cel.GetNumPixels(0,numCells,npixels))
				Err::errAbort("Unable to retrieve num pixels for " + celFiles[celIx]);
			vector<XYCoord> outliers;
			cel.GetOutlierCoords(outliers);

			vector<XYCoord> masked;
			cel.GetMaskedCoords(masked);

			string workingCelFileName = outCelFileName;
			// if multi-channel Cel, then temp single-channel Cels have to be
			// written and then collated later
			if (channelVec.size() > 1) {
				// make temporary output CEL file
				workingCelFileName = GlobalTmpFileFactory()->genFilename(StringUtils::ConvertWCSToMBS(channelVec[channelIx]) + ".", ".cel");
				tmpFileNames.push_back(workingCelFileName);
			}

			affymetrix_calvin_io::CelFileData outCel(workingCelFileName);

			affymetrix_calvin_io::GenericDataHeader *inHdr = cel.GetFileHeader()->GetGenericDataHdr();
			affymetrix_calvin_io::GenericDataHeader *outHdr = outCel.GetFileHeader()->GetGenericDataHdr();
			ParameterNameValueType param;
			int headerCount = inHdr->GetNameValParamCnt();
			for (int hdrIx = 0; hdrIx < headerCount; ++hdrIx) {
				param = inHdr->GetNameValParam(hdrIx);
				outHdr->AddNameValParam(param);
			}
			affymetrix_calvin_io::GenDataHdrVectorIt begin, end;
			inHdr->GetParentIterators(begin, end);
			while (begin != end) {
				// this snippet of commented code was useful for generating
				// test data.  it writes the temporary single-channel CELs as
				// proper CEL files that do not require the CelReader HACK.
				// (hard-coded adjustment of the number of channels in a
				// single-channel CEL that was derived from a multi-channel
				// CEL).  as the code stands, if this snippet were uncommented
				// this tool would break.  the code is left here for now as a
				// reference to how I made the test data.
				//                 if (begin->FindNameValParam(AFFY_WAVELENGTH, param)) {
				//                     param.SetValueText(channelVec[channelIx]);
				//                     begin->AddNameValParam(param);
				//                 }
				outHdr->AddParent(*begin);
				begin++;
			}    

			// add wavelength info if this is multichannel -- the
			// collator uses this field to name the data groups
			if (channelVec.size() > 1) {
				param.SetName(AFFY_FILTER_WAVELENGTH);
				param.SetValueText(channelVec[channelIx]);
				outHdr->AddNameValParam(param);
			}

			if (getOptBool("update-header") == true) {
				// change guid and store old one in header
				affymetrix_calvin_utilities::AffymetrixGuidType oldGuid = inHdr->GetFileId();
				param.SetName(L"affymetrix-previous-guid");
				param.SetValueText(StringUtils::ConvertMBSToWCS(oldGuid));
				outHdr->AddNameValParam(param);
				// store command-line in new CEL header
				param.SetName(L"affymetrix-algorithm-param-command-line");
				param.SetValueText(StringUtils::ConvertMBSToWCS(commandLine()));
				outHdr->AddNameValParam(param);
			}

			outCel.SetIntensityCount(intensities.size());
			outCel.SetStdDevCount(stdev.size());
			outCel.SetPixelCount(npixels.size());
			outCel.SetOutlierCount(outliers.size());
			outCel.SetMaskCount(masked.size());

			//  writer
			affymetrix_calvin_io::CelFileWriter writer(outCel);
			// get and transform intensities.
			for (int featIx = 0; featIx < numCells; featIx++) {
				intensities[featIx] = (*cStreamVec[cStreamVec.size() - 1]).getTransformedIntensity(featIx, celIx, channelIx);
			}

			// write out new intensity values
			writer.WriteIntensities(intensities);

			// copy over other parts
			writer.WriteStdDevs(stdev);
			writer.WritePixels(npixels);
			if (outliers.size() > 0) {
				writer.WriteOutlierCoords(outliers);
			}
			if (masked.size() > 0) {
				writer.WriteMaskCoords(masked);
			}
		}

		if (channelVec.size() > 1) {
			affymetrix_calvin_io::MultiChannelCelFileCollater collator;
			collator.Collate(tmpFileNames, outCelFileName);
			{ for( int i=0; i < tmpFileNames.size();  i++ ) { Fs::rm(tmpFileNames[i]);} } ;
		}
	}
	catch(Except &e) {
		Err::errAbort(e.what());
	}
	catch(const std::bad_alloc &e) {
		{ for( int i=0; i < tmpFileNames.size();  i++ ) { Fs::rm(tmpFileNames[i], false);} } ;
		Err::errAbort("Ran out of memory. "
			"Try quitting other applications.");
	}
	catch(CalvinException &ce) {
		{ for( int i=0; i < tmpFileNames.size();  i++ ) { Fs::rm(tmpFileNames[i], false);} } ;
		Err::errAbort("Affymetrix GeneChip Command Console library has thrown an exception. "
			"Description: '" +StringUtils::ConvertWCSToMBS(ce.Description()) + "'");
	}
	catch(const std::exception &e) {
		Err::errAbort("Exception caught. "
			"Most problems with this program are memory related. "
			"Message is: " + ToStr(e.what()));
	}
	catch(const BaseException &e) {
		Err::errAbort("newmat issue: " + ToStr(e.what()));
	}
	catch(...) {
		Err::errAbort("Unknown exception caught "
			"(most problems with this program are memory related).");
	}
}

/**
* Transform old-style CEL files (text, xda)
*/
void CELTransformerEngine::transformOld(int celIx, vector<ChipStream *> &cStreamVec) {
	affxcel::CCELFileWriter celWriterOld;
	celWriterOld.SetFileName(celFiles[celIx].c_str());
	if (!celWriterOld.Read()) 
		Err::errAbort("Can't read cel file: " + celWriterOld.GetFileName());

	string format = "unknown";
	if (celWriterOld.IsXDACompatibleFile())
		format = "xda";
	else if (celWriterOld.IsTranscriptomeBcelFile())
		format = "bcel";
	else if (celWriterOld.IsCompactCelFile())
		format = "ccel";
	else if (celWriterOld.IsVersion3CompatibleFile())
		format = "text";

	// the celfile might be mapped read-only -- make it writeable
	celWriterOld.EnsureNotMmapped();
	int numCells = celWriterOld.GetNumCells();
	for (int featIx = 0; featIx < numCells; featIx++) {
		float intensity = (*cStreamVec[cStreamVec.size() - 1]).getTransformedIntensity(featIx, celIx);
		celWriterOld.SetIntensity(featIx, intensity);
	}
	string outDir = getOpt("out-dir");
	string outCel = Fs::join(outDir , Fs::basename(celWriterOld.GetFileName()));
	celWriterOld.SetFileName(outCel.c_str());

	bool success = false;
	if (format == "text") 
		success = celWriterOld.WriteTextCel();
	else if (format == "xda")
		success = celWriterOld.WriteXDABCel();
	else if (format == "bcel")
		success = celWriterOld.WriteTranscriptomeBCel();
	else if (format == "ccel") 
		success = celWriterOld.WriteCompactBCel();
	else 
		Err::errAbort("Don't recognize type: " + string(format));

	if (!success)
		Err::errAbort("could not save file: " + outCel);
	celWriterOld.Close();
}

/** 
* Transform all of the cel files specified using the chipstream modification stream
* specified by the chipStreamStr and save the new modified cel files in the output directory.
*/
void CELTransformerEngine::runImp() {
	/* get basic layout info */
	colrow_t numCols, numRows;
	int probeCount, psCount;
	vector<string> chipTypes;

	bool force = getOptBool("force");
	string cdfFile = getOpt("cdf-file");
    string clfFile = getOpt("clf-file");
    string pgfFile = getOpt("pgf-file");
    string spfFile = getOpt("spf-file");

	if (cdfFile != "") {
		EngineUtil::getCdfChipType(chipTypes, numRows, numCols, probeCount, psCount, cdfFile);
		if (!force)
			EngineUtil::checkCelChipTypes(chipTypes, probeCount, celFiles, numRows, numCols);
	}
	else if (pgfFile != "") {
		EngineUtil::getPgfChipType(chipTypes, numRows, numCols, probeCount, pgfFile, clfFile);
		if (!force)
			EngineUtil::checkCelChipTypes(chipTypes, probeCount, celFiles, numRows, numCols);
	}
	else if (spfFile != "") {
		EngineUtil::getSpfChipType(chipTypes, numRows, numCols, probeCount, psCount, spfFile);
		if (!force)
			EngineUtil::checkCelChipTypes(chipTypes, probeCount, celFiles, numRows, numCols);
	}
	else {
		getDimensions(celFiles[0], numCols, numRows);
		probeCount = numCols * numRows;
	}

	/* Get the layout for the chip. */
	// Specifies probesets, locations of features on chip.
	ChipLayout layout; 
	// empty list of probesets, ie we want to load all
	std::set<const char *, Util::ltstr> probeSetsToLoad; 
	// empty list, ie load everything
	std::vector<bool> probeSubset; 
	// we have already check chip type, so use empty string;
	string chipType; 
	//
	probeidmap_t killList;

	if (cdfFile != "") {
		Verbose::out(1, "Opening cdf file: " + Fs::basename(cdfFile));
		if (!layout.openCdf(cdfFile, probeSetsToLoad, NULL, probeSubset, chipType, killList, false)) 
			Err::errAbort("Couldn't open cdf file: " + cdfFile);
	}
	else if (pgfFile != "") {
		affx::ClfFile clf;

		/* open clf file. */
		if (clfFile=="") 
			Err::errAbort("Must specify a clf file for chip.");
		Verbose::out(1, "Opening clf file: " +  Fs::basename(clfFile));
		if (!clf.open(clfFile))
			Err::errAbort("Couldn't open clf file: " + clfFile);
		layout.setDimensions(clf.getXMax() + 1, clf.getYMax() + 1);

		/* Make sure our assumptions about CLF file are true */
		Err::check(clf.getSequential() == 1,
			"ProbesetSummarizeEngine::loadPgfLayout() - unable to handle clf file without sequential set to 1.");
		Err::check(clf.getOrder().compare("col_major") == 0 || clf.getOrder().compare("row_major") == 0,
			"Unable to handle clf file without order set to row_major (old mislabeled 'col_major' accepted due to earlier bug.)");

		/* open pgf */
		Verbose::out(1, "Opening pgf file: " + Fs::basename(pgfFile));
		if (!layout.openPgf(pgfFile, clf.getXMax() + 1, clf.getYMax() + 1, probeSetsToLoad, NULL, NULL, probeSubset, chipType, killList, false, false))
			Err::errAbort("Couldn't open PGF file: " + pgfFile);
	}
	else if (spfFile != "") {
		std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad;
		Verbose::out(1, "Opening spf file: " + Fs::basename(cdfFile));
		layout.openSpf(spfFile, probeSetsToLoad, NULL, probeSubset, chipType, false, psTypesToLoad);
	}
	else {
		Verbose::out(1,"Setting layout dimensions from cel file. Assuming all probes are PM.");
		layout.setDimensions(numCols, numRows);
		vector<bool> mask(layout.getProbeCount(), true);
		layout.setPmProbeMask(mask);
	}

	/* get list of probes to kill */
	/* we do this after the chip layout load as we don't want these to drop on the floor */
	string killListFile = getOpt("kill-list");
	if (killListFile != "")
		ChipLayout::fillInKillList(killListFile, killList, numRows, numCols);

	/* Setup chipstream factory */
	ChipStreamFactory cFactory;

	string targetSketch = getOpt("target-sketch");
	string probeNormFile = getOpt("probe-norm-file");
	string outDir = getOpt("out-dir");
	if (targetSketch != "") 
		cFactory.readTargetSketchFromFile(targetSketch);
	if (getOptBool("write-sketch") == true) 
		cFactory.setWriteSketchDir(outDir + Fs::osPathSep());
	if (getOptBool("a5-write-sketch") == true) 
		cFactory.setWriteSketchDir_a5(outDir + Fs::osPathSep());
	if (probeNormFile != "")
		cFactory.readNormProbesFromFile(probeNormFile);
	if (killListFile != "") {
		/// @todo allow provision of kill list to Chipstream Factor
		///cFactory.setKillList(killList);
		Err::errAbort("Kill list not yet implemented. Try probe-norm-file instead.");
	}

	/* Create our chipstream objects. */
	std::string chipStreamStr = getOpt("chipstream");
	ChipStream *cStream = NULL;      // Our chipstream for transforming data.
	vector<ChipStream *> cStreamVec;
	vector<string> words;

	Util::chopString(chipStreamStr, ',', words);
	for (unsigned int i = 0; i < words.size(); i++) {
		string dummy;
		cStream = cFactory.chipStreamForString(words[i], layout, dummy );
		if (!cStreamVec.empty()) {
			cStreamVec[cStreamVec.size() - 1]->registerStream(cStream);
			cStream->registerParent(cStreamVec[cStreamVec.size() - 1]);
		}
		cStreamVec.push_back(cStream);
	}

	/* Create intensity mart */
	std::vector<int> dummy_order(probeCount);
	for (int i = 0; i < probeCount; i++) {
		dummy_order[i] = i;
	}

	IntensityMart* iMart = NULL;
	bool useDisk = getOptBool("use-disk");
	if (useDisk) {
		int cacheSize = getOptInt("disk-cache");
		string tempDir = getOpt("temp-dir");
		DiskIntensityMart* diskMart = new DiskIntensityMart(dummy_order,
			celFiles, 
			cacheSize * 1048576, 
			tempDir, 
			"apt-transform.tmp.", 
			true);
		iMart = diskMart;
	}
	else {
		SparseMart* sparseMart = new SparseMart(dummy_order,celFiles,true);
		iMart = sparseMart;
	}

	/* set cel files */
	CelReader cReader;               // Object for getting data from cel files.

	cReader.setFiles(celFiles);

	/* Let reader know to pass data to our chipstream object. */
	cReader.registerStream(cStreamVec[0]);
	cReader.registerIntensityMart(iMart);
	Verbose::out(1,"Reading cel files.");

        /* Open, read, and send cel file data one by one to chipstream. */
        cReader.readFiles();

        /* Transform the data and write it out. */
        unsigned int dotMod = max((int)celFiles.size()/20, 1);
        Verbose::progressBegin(1, "Adjusting and writing files", 20, (int)dotMod, (int)celFiles.size());

        for (int i = 0; i < celFiles.size(); i++) {
                Verbose::progressStep(1);
                if (affymetrix_cel_converter::CELFileVersion::DetermineCELFileVersion(celFiles[i].c_str())==affymetrix_cel_converter::Calvin_Version1)
                        transformNew(i,cStreamVec);
                else
                        transformOld(i,cStreamVec);
        }

        delete cStream;

        Freez(iMart);

        /* Remove the temp dir if it is empty */
        string tempDir = getOpt("temp-dir");
        if (Fs::isWriteableDir(tempDir))
                Fs::rmdir(tempDir,false);

        Verbose::progressEnd(1, "Done.");
}

void CELTransformerEngine::checkOptionsImp()
{
        EngineUtil::getCelFiles(celFiles, this);
        
        /* Sanity checks... */
        if (getOpt("chipstream") == "") 
                Err::errAbort("Must specify a chipstream transformer.");
        if (getOpt("out-dir") == "") 
                Err::errAbort("Must specify an output directory.");
        if (getOpt("probe-norm-file") != "" && getOpt("kill-list") != "") 
                Err::errAbort("Cannot provide a kill list and a probe norm list.");
        if (celFiles.empty()) 
                Err::errAbort("Must specify at least one cel file to open.");

        /* Make our directory. */
        string outDir = getOpt("out-dir");
        Util::chompLastIfSep(outDir);
        setOpt("out-dir", outDir);
        if (!Fs::isWriteableDir(outDir)) {
                if (Fs::mkdirPath(outDir, false) != APT_OK) {
                        APT_ERR_ABORT("Can't make or write to directory: " + outDir);
                }
        }

        string tempDir = getOpt("temp-dir");
        if (tempDir == "") {
                tempDir = Fs::join(outDir,"temp"); 
                setOpt("temp-dir", tempDir);
        }
}
