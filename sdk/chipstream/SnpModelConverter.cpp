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

//
#include "chipstream/SnpModelConverter.h"
//
#include "chipstream/QuantLabelZ.h"
#include "chipstream/QuantLabelZIO.h"
#include "chipstream/SnpModelDb.h"
#include "file5/File5.h"
#include "util/Convert.h"
#include "util/Fs.h"
#include "util/RowFile.h"
#include "util/md5sum.h"
#include "util/SQLite.h"
//
#include <sstream>

using namespace std;
using namespace affx;

const std::string FILE_INPUT = "from-model-file";
const std::string TOTAL_MODEL_COUNT = "total-model-count";
const std::string USED_FOR_PRIOR = "used-for-prior";
const std::string OFFICIAL_AFFY = "official-affy-model";
const std::string CHIP_ANALYSIS_TYPE = "chip-analysis-type";
const std::string CHIP_TYPE = "chip-type";
const std::string ARRAY_TYPE = "array-type";
const std::string MODEL_TYPE = "model-type";
const std::string ALGORITHM_NAME = "algorithm-name";
const std::string FILE_TYPE = "file-type";
const std::string REAGENT_KEY_VERSION = "reagent-key-version";
const std::string FILE_TYPE_MODELS = "models";
const std::string UNKNOWN_CHIP_TYPE = "UNKNOWN_CHIP";

const char *SnpModelFileType_str[]={
  "_UNKNOWN_FILE",
  "_BIRDSEED_PRIOR",
  "_BIRDSEED_POSTERIOR",
  "_BRLMMP_MODEL",
  "_BRLMMP_MODEL_HDF5",
  "_SQLITE_MODEL_FILE"
};

string nullStr = "null";

void SnpModelConverter::headersToMeta(TsvFile &tsv, SQLiteDatabase &db) {
    string key;
    string value;
    tsv.headersBegin();
    while (tsv.headersNext(key, value) == TSV_OK) {
        string sql = "insert into snp_model_meta values ('" + key + "','" + value + "');";
        db.execute(sql);
    }
}

SnpModelConverter::SnpModelFileType SnpModelConverter::guessFileType(const std::string &file) {
    if (File5_File::isHdf5file(file)) {
        return BRLMMP_MODEL_HDF5;
    }
    else if (SnpModelDb::isSqlModelFile(file)) { 
        return SQLITE_MODEL_FILE;
    }
    else {
        return guessTsvFileType(file);
    }
    return UNKNOWN_FILE;
}

void SnpModelConverter::convertToDbModel(const std::string &fileIn, const std::string &fileOut, const std::string &file5Path) {
    if (m_FileType == UNKNOWN_FILE) {
        m_FileType = guessFileType(fileIn);
    }
    Verbose::out(2, "File " + fileIn + " is type is: " + ToStr(m_FileType));
    if (m_FileType == SnpModelConverter::BIRDSEED_PRIOR) {
        convertBirdseedPriorToDbModel(fileIn, fileOut);
    }
    else if (m_FileType == SnpModelConverter::BIRDSEED_POSTERIOR) {
        convertBirdseedPosteriorToDbModel(fileIn, fileOut);
    }
    else if (m_FileType == SnpModelConverter::BRLMMP_MODEL) {
        convertBrlmmPToDbModel(fileIn, fileOut);
    }
    else if (m_FileType == SnpModelConverter::BRLMMP_MODEL_HDF5) {
        convertBrlmmPFile5ToDbModel(fileIn, file5Path, fileOut);
    }
    else if (m_FileType == SnpModelConverter::SQLITE_MODEL_FILE) {
        Verbose::out(2, "File " + fileIn + " is already a db.");
    }
    else {
        APT_ERR_ABORT("File " + fileIn + " is an unknown model type: " + ToStr(m_FileType));
    }
}

SnpModelConverter::SnpModelFileType SnpModelConverter::guessTsvFileType(const std::string &file) {
    RowFile rf;
    rf.open(file);
    vector<string> words;
    rf.nextRow(words);
    if (words.size() >= 5 && words[0] == "id" && words[1] == "BB" && words[2] == "AB" && words[3] == "AA" && words[4] == "CV") {
        return BRLMMP_MODEL;
    }
    else if (words.size() >= 5 && words [0] == "probeset_id" && words[1] == "copy_number" && words[2] == "BB" && words[3] == "AB" && words[4] =="AA") {
        return BIRDSEED_PRIOR;
    }
    else if (words.size() == 1 && words[0].find("SNP_") != string::npos && words[0].find(';') != string::npos) {
        return BIRDSEED_POSTERIOR;
    }
    return UNKNOWN_FILE;
}

void SnpModelConverter::convertBrlmmPToDbModel(const std::string &fileIn, const std::string &fileOut) {
    SnpModelDb snpDb(fileOut);
    void setupTables();  
    
    affx::TsvFile tsv;
    tsv.open(fileIn);

    snpDb.setupTables();
    string key, value;
    tsv.headersBegin();
    while (tsv.headersNext(key, value) == TSV_OK) {
        snpDb.addMeta(key, value);
    }
    vector<string> chipTypes;
    tsv.getHeader("chip_type", chipTypes);
    if (chipTypes.size() == 0) {
        tsv.getHeader(CHIP_TYPE, chipTypes);
    }
    int numProbesets;
    if (TSV_OK != tsv.getHeader(TOTAL_MODEL_COUNT, numProbesets)) {
        enum SnpChipType type = m_ChipType;
        if (type == UNKNOWN_CHIP) {
            type = guessChipType(m_FileType, chipTypes, "");
        }
        numProbesets = guessNumProbesets(type);
    }
    string aa, ab, bb, id,cv;
    // This is the "v1" format with comma seperated values...
    tsv.bind(0,"id", &id, affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"BB", &bb, affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"AB", &ab, affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"AA", &aa, affx::TSV_BIND_REQUIRED);

    int fc_cidx=tsv.cname2cidx(0,"CV");
    if (fc_cidx!=affx::TSV_ERR_NOTFOUND) {
        tsv.bind(0,"CV", &cv, affx::TSV_BIND_REQUIRED);
        ;  }
    else {
        cv = "0,0,0";
    }
    int iCount = 0;
    string first;
    try {
        unsigned int dotMod = max(numProbesets/20,1);
        Verbose::progressBegin(1, "Indexing " + ToStr(numProbesets) + " Probesets", 20, dotMod, numProbesets);
        snpDb.prepareToWriteModels(numProbesets);
        while (tsv.nextLevel(0) == affx::TSV_OK) {
            Verbose::progressStep(1);
            if (first.empty()) {
                first = id;
            }
            snp_labeled_distribution p = QuantLabelZ__CreatePriorFromStrings(bb,ab,aa,cv);
            p.probeset_id = id;
            snpDb.writeSnpDistribution(p);
            iCount++;
        }
        snpDb.doneWriteModels();
        snpDb.addMeta(TOTAL_MODEL_COUNT, ToStr(iCount));
        string file = Fs::basename(fileIn);
        snpDb.addMeta(FILE_INPUT, file);
        snpDb.addMeta(MODEL_TYPE, ToStr(m_FileType));
        if (m_ChipType == UNKNOWN_CHIP) {
            m_ChipType = guessChipType(m_FileType, chipTypes, first);
        }
        snpDb.addMeta(CHIP_ANALYSIS_TYPE, ToStr(m_ChipType));
		string arrayType = UNKNOWN_CHIP_TYPE;
		if (chipTypes.size() > 0)
			arrayType = chipTypes[0];
        snpDb.addMeta(ARRAY_TYPE, arrayType);
        Verbose::progressEnd(1, "Done.");
    }
    catch (SQLiteException &e) {
        Err::errAbort("Got Exception at record: " + ToStr(iCount) + " " + e.getMessage());
    }

    tsv.close();  
}

bool SnpModelConverter::vecContains(const std::string &needle, const std::vector<std::string> &haystack) {
    for (int i = 0; i < haystack.size(); i++) {
        if (haystack[i].find(needle) != string::npos) {
            return true;
        }
    }
    return false;
}

int SnpModelConverter::guessNumProbesets(enum SnpModelConverter::SnpChipType type) {
    int num = 500000;
    if (type == AXIOM) {
        num = 582406;
    }
    else if (type == SNP6) {
        num = 972549;
    }
    else if (type == SNP5) {
        num = 513865;
    }
    Verbose::out(2, "Num probesets guess is: " + ToStr(num));
    return num;
}

enum SnpModelConverter::SnpChipType SnpModelConverter::guessChipType(enum SnpModelConverter::SnpModelFileType modelType,
                                                                     const std::vector<std::string> &chipTypes, 
                                                                     const std::string &firstPsName) {
    enum SnpModelConverter::SnpChipType type = UNKNOWN_CHIP;
    // Believe the chip-type
    if (vecContains("Axiom", chipTypes)) {
        type = AXIOM;
    }
    else if(vecContains("GenomeWideSNP_6", chipTypes) || vecContains("GenomeWideEx_6", chipTypes)) {
        type = SNP6;
    }
    else if(vecContains("GenomeWideSNP_5", chipTypes)) {
        type = SNP5;
    }

    // Guess based on type, birdseed only really for snp6
    else if (modelType == BIRDSEED_PRIOR || modelType == BIRDSEED_POSTERIOR) {
        type = SNP6;
    }
    
    // Guess based on first probeset name
    if (firstPsName == "SNP_A-1780520") {
        type = SNP5;
    }
    else if (firstPsName == "AFFX-SNP_10000979") {
        type = SNP6;
    }
    else if (firstPsName == "AX-11086525" || firstPsName.find("AX-") == 0) {
        type = AXIOM;
    }
    Verbose::out(2, "Chip Type is: " + ToStr(type));
    return type;
}


void SnpModelConverter::convertBrlmmPFile5ToDbModel(const std::string &fileIn, const std::string &file5Path, const std::string &fileOut) {
    SnpModelDb snpDb(fileOut);
    void setupTables();  

    affx::File5_File* file5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;
    string path = file5Path;
    file5 = new affx::File5_File();
    file5->open(fileIn,affx::FILE5_OPEN_RO);
    if (path.length() == 0) {
        affx::File5_Group* group5 = NULL;
        group5 = file5->openGroup("/",affx::FILE5_OPEN_RO);
        vector<string> names = group5->listNames();
        for (int i = 0; i < names.size(); i++) {
            if (names[i].find(".snp-posteriors") != string::npos) {
                path = names[i];
                break;
            }
            
        }
        if (path.length() == 0) {
            APT_ERR_ABORT("Don't recognize this type of posterior file.");
        }
        group5->close();
        Freez(group5);
    }

    tsv5 = file5->openTsv(path);

    int numProbesets = tsv5->getLineCount();
    int iCount = 0;
    try {
        snpDb.setupTables();
        string key, value;
        for (int i = 0; i < tsv5->getHeaderCount(); i++) {
            tsv5->getHeaderByIdx(i, &key, &value);
            snpDb.addMeta(key, value);
        }
        
        vector<string> chipTypes;
        string type;
        tsv5->getHeader("chip_type", &type);
        chipTypes.push_back(type);

        string first;
        snp_labeled_distribution p;
        p.Dist.Clear();
        snpDb.prepareToWriteModels(numProbesets);
        unsigned int dotMod = max(numProbesets/20,1);
        Verbose::progressBegin(1, "Indexing " + ToStr(numProbesets) + " Probesets", 20, dotMod, numProbesets);
        while (tsv5->nextLevel(0) == affx::FILE5_OK) {
            Verbose::progressStep(1);
            tsv5->get(0,0, &p.probeset_id);
            if (first.empty()) {
                first = p.probeset_id;
            }
            QuantLabelZ__readSnpPrior_tsv5_v2(tsv5, p.Dist);
            snpDb.writeSnpDistribution(p);
            iCount++;
        }
        snpDb.doneWriteModels();
        snpDb.addMeta(TOTAL_MODEL_COUNT, ToStr(iCount));
        string file = Fs::basename(fileIn);
        snpDb.addMeta(FILE_INPUT, file);
        snpDb.addMeta(MODEL_TYPE, ToStr(m_FileType));
        if (m_ChipType == UNKNOWN_CHIP) {
            m_ChipType = guessChipType(m_FileType, chipTypes, first);
        }
        snpDb.addMeta(CHIP_ANALYSIS_TYPE, ToStr(m_ChipType));
        snpDb.addMeta(ARRAY_TYPE, type);
        Verbose::progressEnd(1, "Done.");
    }
    catch (SQLiteException &e) {
        Err::errAbort("Got Exception at record: " + ToStr(iCount) + " " + e.getMessage());
    }
    tsv5->close();
    delete tsv5;
    file5->close();
    delete file5;
}


void SnpModelConverter::birdseedStringToCluster(const std::string &s, char sep, cluster_data &d) {
    vector<string> words;
    if (s != "null") {
        Util::chopString(s, sep, words);
        APT_ERR_ASSERT(words.size() == 6, "Birdseed cluster malformed: " + s);
        d.m = Convert::toDouble(words[0]);
        d.ym = Convert::toDouble(words[1]);
        d.ss = Convert::toDouble(words[2]);
        d.xyss = Convert::toDouble(words[3]);
        d.yss = Convert::toDouble(words[4]);
        d.k = Convert::toDouble(words[5]);
        d.v = Convert::toDouble(words[5]);
    }
    else {
        d.m = 0;
        d.ym = 0;
        d.ss = 0;
        d.xyss = 0;
        d.yss = 0;
        d.k = 0;
        d.v = 0;
    }
}

void SnpModelConverter::birdseedModelFromStrings(const std::string &name,
                                                 int copynumber,
                                                 const std::string &bb, 
                                                 const std::string &ab, 
                                                 const std::string &aa,
                                                 char sep,
                                                 snp_labeled_distribution &p) {
    p.Dist.Clear();
    if (copynumber == 1) {
        p.probeset_id = name + ":" + ToStr(copynumber);
    }
    else {
        p.probeset_id = name;
    }
    birdseedStringToCluster(bb, sep, p.Dist.bb);
    birdseedStringToCluster(ab, sep, p.Dist.ab);
    birdseedStringToCluster(aa, sep, p.Dist.aa);
}

void SnpModelConverter::convertBirdseedPriorToDbModel(const std::string &fileIn, const std::string &fileOut) {
    SnpModelDb snpDb(fileOut);
    void setupTables();  
    
    affx::TsvFile tsv;
    tsv.open(fileIn);


    snpDb.setupTables();
    string key, value;
    tsv.headersBegin();
    while (tsv.headersNext(key, value) == TSV_OK) {
        snpDb.addMeta(key, value);
    }
    vector<string> chipTypes;
    tsv.getHeader("chip_type", chipTypes);    
    if (chipTypes.size() == 0) {
        tsv.getHeader(CHIP_TYPE, chipTypes);
    }
    int numProbesets;
    if (TSV_OK != tsv.getHeader(TOTAL_MODEL_COUNT, numProbesets)) {
        enum SnpChipType type = m_ChipType;
        if (type == UNKNOWN_CHIP) {
            type = guessChipType(m_FileType, chipTypes, "");
        }
        numProbesets = guessNumProbesets(type);
    }

    string aa, ab, bb, id;
    int copynumber;
    // This is the format with comma seperated values...
    tsv.bind(0,"probeset_id", &id, affx::TSV_BIND_REQUIRED);
    tsv.bind(0, "copy_number", &copynumber, affx::TSV_BIND_REQUIRED);
    // swap AA and BB when writing to database.  This is compensating for the
    // fact that the AA and BB clusters have swapped labels in the priors file.
    tsv.bind(0,"BB", &aa, affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"AB", &ab, affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"AA", &bb, affx::TSV_BIND_REQUIRED);
    int iCount = 0;
    snp_labeled_distribution p;
    p.Dist.Clear();
    string first;
    try {
        unsigned int dotMod = max(numProbesets/20,1);
        Verbose::progressBegin(1, "Indexing " + ToStr(numProbesets) + " Probesets", 20, dotMod, numProbesets);
        snpDb.prepareToWriteModels(numProbesets);
        while (tsv.nextLevel(0) == affx::TSV_OK) {
            Verbose::progressStep(1);
            if (first.empty()) {
                first = id;
            }
            birdseedModelFromStrings(id, copynumber, bb, ab, aa, ',',  p);
            snpDb.writeSnpDistribution(p);
            iCount++;
        }
        snpDb.doneWriteModels();
        snpDb.addMeta(TOTAL_MODEL_COUNT, ToStr(iCount));
        string file = Fs::basename(fileIn);
        snpDb.addMeta(FILE_INPUT, file);
        snpDb.addMeta(MODEL_TYPE, ToStr(m_FileType));
        if (m_ChipType == UNKNOWN_CHIP) {
            m_ChipType = guessChipType(m_FileType, chipTypes, first);
        }
        snpDb.addMeta(CHIP_ANALYSIS_TYPE, ToStr(m_ChipType));
		string arrayType = UNKNOWN_CHIP_TYPE;
		if (chipTypes.size() > 0)
			arrayType = chipTypes[0];
        snpDb.addMeta(ARRAY_TYPE, arrayType);
        Verbose::progressEnd(1, "Done.");
    }
    catch (SQLiteException &e) {
        APT_ERR_ABORT("Got Exception at record: " + ToStr(iCount) + " " + e.getMessage());
    }
    tsv.close();  
}

void SnpModelConverter::addBirdseedPosteriorFromStrings(std::string &id, const std::string &aa, 
                                                        const std::string &ab, const std::string &bb,
                                                        snp_labeled_distribution &p) {
    size_t offset = id.rfind('-');
    APT_ERR_ASSERT(offset != string::npos, "Couldn't find copynumber in: " + id );
    string cpnumber = id.substr(offset+1,1);
    int copynumber = Convert::toInt(cpnumber);
    id = id.substr(0, offset);
    if (id == "SNP_A-1984828") {
        Verbose::out(1, "Here we come...");
    }
    if (copynumber == 1) {
        birdseedModelFromStrings(id, copynumber, ab, nullStr, aa, ' ', p);
    }
    else {
        birdseedModelFromStrings(id, copynumber, bb, ab, aa, ' ', p);
    }
    
}

void SnpModelConverter::convertBirdseedPosteriorToDbModel(const std::string &fileIn, const std::string &fileOut) {
    SnpModelDb snpDb(fileOut);
    void setupTables();  
    
    affx::TsvFile tsv;
    tsv.m_optFieldSep = ';';
    tsv.m_optAutoSenseSep = false;
    tsv.open(fileIn);
    int numProbesets;

    //    tsv.getHeader("rows", numProbesets);

    snpDb.setupTables();
    string key, value;
    tsv.headersBegin();
    while (tsv.headersNext(key, value) == TSV_OK) {
        snpDb.addMeta(key, value);
    }
    vector<string> chipTypes;
    tsv.getHeader("chip_type", chipTypes);
    if (chipTypes.size() == 0) {
        tsv.getHeader(CHIP_TYPE, chipTypes);
    }
    if (TSV_OK != tsv.getHeader(TOTAL_MODEL_COUNT, numProbesets)) {
        enum SnpChipType type = m_ChipType;
        if (type == UNKNOWN_CHIP) {
            type = guessChipType(m_FileType, chipTypes, "");
        }
        numProbesets = guessNumProbesets(type);
    }
    string aa, ab, bb, id;
    // This is the format with comma seperated values...
    tsv.bind(0, 0, &id, affx::TSV_BIND_REQUIRED);
    tsv.bind(0, 1, &aa, affx::TSV_BIND_REQUIRED);
    tsv.bind(0, 2, &ab, affx::TSV_BIND_REQUIRED);
    tsv.bind(0, 3, &bb, affx::TSV_BIND_REQUIRED);

    int iCount = 0;
    snp_labeled_distribution p;
    p.Dist.Clear();

    string first;
    try {
        unsigned int dotMod = max(numProbesets/20,1);
        Verbose::progressBegin(1, "Indexing " + ToStr(numProbesets) + " Probesets", 20, dotMod, numProbesets);
        snpDb.prepareToWriteModels(numProbesets);
        // This is an odd little bit as TsvFile thinks that the first row is the header...
        id = tsv.getColumnName(0,0);
        aa = tsv.getColumnName(0,1);
        ab = tsv.getColumnName(0,2);
        bb = tsv.getColumnName(0,3);
        first = id;
        addBirdseedPosteriorFromStrings(id, aa, ab, bb, p);
        snpDb.writeSnpDistribution(p);
        while (tsv.nextLevel(0) == affx::TSV_OK) {
            Verbose::progressStep(1);
            addBirdseedPosteriorFromStrings(id, aa, ab, bb, p);
            snpDb.writeSnpDistribution(p);
            iCount++;
            aa.clear();
            ab.clear();
            bb.clear();
        }
        snpDb.doneWriteModels();
        snpDb.addMeta(TOTAL_MODEL_COUNT, ToStr(iCount));
        string file = Fs::basename(fileIn);
        snpDb.addMeta(FILE_INPUT, file);
        snpDb.addMeta(MODEL_TYPE, ToStr(m_FileType));
        if (m_ChipType == UNKNOWN_CHIP) {
            m_ChipType = guessChipType(m_FileType, chipTypes, first);
        }
        snpDb.addMeta(CHIP_ANALYSIS_TYPE, ToStr(m_ChipType));
		string arrayType = UNKNOWN_CHIP_TYPE;
		if (chipTypes.size() > 0)
			arrayType = chipTypes[0];
        snpDb.addMeta(ARRAY_TYPE, arrayType);
		Verbose::progressEnd(1, "Done.");
    }
    catch (SQLiteException &e) {
        APT_ERR_ABORT("Got Exception at record: " + ToStr(iCount) + " " + e.getMessage());
    }
    tsv.close();  
}

std::string SnpModelConverter::getHd5ArrayTypeAndAlgorithm(const std::string &fileIn, enum SnpModelFileType ft, std::string &alg)
{
	string arrayType = UNKNOWN_CHIP_TYPE;
	affx::File5_File* file5 = NULL;
	affx::File5_Tsv* tsv5 = NULL;
	file5 = new affx::File5_File();
	file5->open(fileIn,affx::FILE5_OPEN_RO);
	string path = "";
	affx::File5_Group* group5 = file5->openGroup("/",affx::FILE5_OPEN_RO);
	vector<string> names = group5->listNames();
	for (int i = 0; i < names.size(); i++) {
		if (names[i].find(".snp-posteriors") != string::npos) {
			path = names[i];
			break;
		}
	}
	if (path.length() > 0)		// valid .A5 (posterior) file
	{
		tsv5 = file5->openTsv(path);
        string type;
		tsv5->getHeader("chip_type", &type);
		if (!type.empty())
			arrayType = type;
	    tsv5->close();
	    delete tsv5;
	}
	string tempAlg;
	tsv5->getHeader(ALGORITHM_NAME, &tempAlg);
	if (!tempAlg.empty())
		alg = tempAlg;
	group5->close();
	Freez(group5);
    file5->close();
    delete file5;
	return arrayType;
}

std::string SnpModelConverter::getTsvArrayTypeAndAlgorithm(const std::string &fileIn, enum SnpModelFileType ft, std::string &alg)
{
	string arrayType = UNKNOWN_CHIP_TYPE;
	affx::TsvFile tsv;
	tsv.open(fileIn);
    tsv.headersBegin();
	vector<string> chipTypes;
	tsv.getHeaderMatchingKeySubstr("chip_type", chipTypes);
	if (chipTypes.size() == 0)
		tsv.getHeaderMatchingKeySubstrAppend(CHIP_TYPE, chipTypes);
	if (chipTypes.size() > 0)
		arrayType = chipTypes[0];
	string tempAlg;
	if (tsv.getHeader(ALGORITHM_NAME, tempAlg) == TSV_OK)
		alg = tempAlg;
	return arrayType;
}

void SnpModelConverter::getArrayTypeAndAlgorithm(const std::string &fileIn, std::vector<std::string> &result) {
	string alg = "";
	string arrayType = "";
	try {
		result.clear();
		enum SnpModelFileType ft = guessFileType(fileIn);
		if (ft == SQLITE_MODEL_FILE)
		{
			SnpModelDb snpDb;
			snpDb.open(fileIn);
			vector<string> values;
			snpDb.getMetaValue(ALGORITHM_NAME, values);
			if (values.size() == 0) {
				snpDb.getMetaValue(MODEL_TYPE, values);
				if (values.size() == 0) {
					alg = SnpModelFileType_str[UNKNOWN_FILE];
				}
				else {
					alg = SnpModelFileType_str[atoi(values[0].c_str())];
				}
			}
			else {
				alg = values[0];
			}
			snpDb.getMetaValue(ARRAY_TYPE, values);
			snpDb.close();
			if (values.size() == 0) {
				arrayType = UNKNOWN_CHIP;
			}
			else {
				arrayType = values[0];
			}

		}
		else
		{
			alg = SnpModelFileType_str[ft];
			if (ft == BRLMMP_MODEL_HDF5) {
				alg = SnpModelFileType_str[BRLMMP_MODEL];
				arrayType = getHd5ArrayTypeAndAlgorithm(fileIn, ft, alg);
			}
			else {
				arrayType = getTsvArrayTypeAndAlgorithm(fileIn, ft, alg);
			}
		}
	}
    catch (SQLiteException &e) {
      Verbose::out(2,"SnpModelConverter::getArrayTypeAndAlgorithm got exception: " + e.getMessage());
    }
    catch (Exception &e){
      string s = string(e.what());
      Verbose::out(2,"SnpModelConverter::getArrayTypeAndAlgorithm got exception: " + s);
    }
    catch (...) {
      Verbose::out(2,"SnpModelConverter::getArrayTypeAndAlgorithm got exception");
    }
	result.push_back(arrayType);
	result.push_back(alg);
}

/////////////////////////////////////////////////////////////////////////////////
// SnpModelConverter::GetModelHeaders CODE START
/////////////////////////////////////////////////////////////////////////////////

enum SMC_GH_ModelType {
  SMC_GH_MODEL_UNKNOWN,
  SMC_GH_MODEL_PRIOR,
  SMC_GH_MODEL_POSTERIOR,
  SMC_GH_MODEL_LAST_MODEL,
};

enum SMC_GH_AlgorithmType  {
  SMC_GH_ALGORITHM_UNKNOWN,
  SMC_GH_ALGORITHM_BIRDSEED,
  SMC_GH_ALGORITHM_BRLMMP,
  SMC_GH_ALGORITHM_LAST_ALGORITHM,
};

const std::string SNP_MODEL_ALGORITHM_TYPE_GTC[] = {
  "unknown",
  "birdseed",
  QUANTBRLMMP,
};

const std::string SNP_MODEL_MODEL_TYPE_GTC[] = {
  "unknown",
  "prior",
  "posterior",
};

class SnpModelConverterValidation
{
public:
  std::string  m_name;
  bool         m_emptyOk;
  std::string  m_description;
  pcrecpp::RE  m_validHeaderRE;
  pcrecpp::RE  m_validFileNameRE;
  std::string  m_value; // run-time only.
  //pcrecpp::RE  m_valueRE; ?? do we validate the value as well?
};

const SnpModelConverterValidation SMC_VALIDATION[] = {
  {
    std::string(FILE_TYPE),
    false,
    std::string("Set when valid model file detected."),
    pcrecpp::RE("file-type|affymetrix-algorithm-param-apt-opt-model-file-(brlmm|brlmmp|birdseed)|SnpPosteriorFormatVer"),
    pcrecpp::RE("^$"),
    std::string(""),
  },
  {
    std::string(MODEL_TYPE),
    false,
    std::string("One of either 'prior' or 'posterior' type model file."),
    pcrecpp::RE("model-type|SnpPosteriorFormatVer"),
    pcrecpp::RE("$^"),
    std::string(""),
  },
  {
    std::string(CHIP_TYPE),
    false,
    std::string("Chip type [GenomeWideSNP_5 | GenomeWideSNP_6 | Axiom_GW_Hu_SNP ..]."),
    pcrecpp::RE(".*chip[-_]type"),
    pcrecpp::RE("^$"),
    std::string(""),
  },
  {
    std::string(ALGORITHM_NAME),
    false,
    std::string("Algorithm or analysis name [birdseed, brlmm, brlmmp, ...]."),
    pcrecpp::RE("algorithm-name|affymetrix-algorithm-param-apt-opt-analysis-name"),
    pcrecpp::RE("^$"),
    std::string(""),
  },
  {
    std::string(REAGENT_KEY_VERSION),
    true,
    std::string("Model file type headers detected."),
    pcrecpp::RE("reagent-key-version|affymetrix-algorithm-param-apt-opt-reagent-key-version"),
    pcrecpp::RE("^$"),
    std::string(""),
  },

};

const unsigned int SMC_VALIDATION_SIZE = sizeof(SMC_VALIDATION) /
    sizeof(SMC_VALIDATION[0]);

struct SnpModelConverterAddDerivedHeadersCallback {
  int m_type;
  void (*m_callback)(std::map< std::string, std::string *> & headerMap,
                     const std::vector< SnpModelConverterValidation > & validations,
                     std::map< std::string, std::string > & derivedHeaders,
                     const char * headerMapValueDelimiter);
};

struct SnpModelConverterAlgorithmColumns {
  std::string m_algorithm;
  std::string m_type;
  int         m_numColumnsToCheck;
  pcrecpp::RE m_col1;
  pcrecpp::RE m_col2;
  pcrecpp::RE m_col3;
  pcrecpp::RE m_col4;
  pcrecpp::RE m_col5;

};

const SnpModelConverterAlgorithmColumns SMC_ALGORITHM_COLUMNS[] = {
  {
    std::string(SNP_MODEL_ALGORITHM_TYPE_GTC[SMC_GH_ALGORITHM_BRLMMP]),
    std::string(SNP_MODEL_MODEL_TYPE_GTC[SMC_GH_MODEL_UNKNOWN]),
    5,
    pcrecpp::RE("id"),
    pcrecpp::RE("BB"),
    pcrecpp::RE("AB"),
    pcrecpp::RE("AA"),
    pcrecpp::RE("CV"),
  },
  {
    std::string(SNP_MODEL_ALGORITHM_TYPE_GTC[SMC_GH_ALGORITHM_BIRDSEED]),
    std::string(SNP_MODEL_MODEL_TYPE_GTC[SMC_GH_MODEL_PRIOR]),
    5,
    pcrecpp::RE("probeset_id"),
    pcrecpp::RE("copy_number"),
    pcrecpp::RE("BB"),
    pcrecpp::RE("AB"),
    pcrecpp::RE("AA"),
  },
  //SNP_A-1984828-2;951.04 381.871 7471.51 -2071.62 4962.37 0.65;709.288 965.904 5277.17 2515.01 10355 0.3;424.846 1464.56 4263.25 3602.01 26292.1 0.05
  {
    std::string(SNP_MODEL_ALGORITHM_TYPE_GTC[SMC_GH_ALGORITHM_BIRDSEED]),
    std::string(SNP_MODEL_MODEL_TYPE_GTC[SMC_GH_MODEL_POSTERIOR]),
    1,
    pcrecpp::RE("(AFFX-)?SNP_[^;]+;(-?\\d+(\\.\\d+)?[\\s;])+-?\\d+(\\.\\d+)?"),
    pcrecpp::RE(""),
    pcrecpp::RE(""),
    pcrecpp::RE(""),
    pcrecpp::RE(""),
  },
};

const int SMC_ALGORITHM_COLUMNS_SIZE =
  sizeof(SMC_ALGORITHM_COLUMNS)  / sizeof(SMC_ALGORITHM_COLUMNS[0]);

static void _SnpModelConverter_AddHeaderToHeaderMap(std::map< std::string, std::string *> & headerMap, const char * headerMapValueDelimiter, const std::string & key, const std::string & value);
// GTC 4.1 Mapped header callback.
static void addDerivedHeadersCallback_GTC4_1(std::map< std::string, std::string *> & headerMap,
    const std::vector< SnpModelConverterValidation > & validations,
                                             std::map< std::string, std::string > & derivedMap,
                                             const char * headerMapValueDelimiter);

// Register new Header Type callbacks here.
const SnpModelConverterAddDerivedHeadersCallback SMC_ADD_DERIVED_HEADERS_CALLBACK[] = {
  { SnpModelConverter::NONE, NULL },
  { SnpModelConverter::GTC4_1, addDerivedHeadersCallback_GTC4_1 },
};
/*****************************************************************************
 * SnpModelGetHeaders: Header Type Call Back Functions.
 *****************************************************************************/
/*****************************************************************************/
/**
 * addDerivedHeadersCallback_GTC4_1:
 * This api assumes the headers and data have been validated.
 * Further the derived data values of algorithm and model type
 * may also be provided for possible use.
 *
 * @param headerMap     - file headers.
 * @param validations   - normalized headers
 * @param derivedMap    - data headers.
 *
 * @return headerMap    - Final mapping including the derived values.
 *
 */
/*****************************************************************************/
static void addDerivedHeadersCallback_GTC4_1(std::map< std::string, std::string *> & headerMap,
    const std::vector< SnpModelConverterValidation > & validations,
                                             std::map< std::string, std::string > & derivedMap,
                                             const char * headerMapValueDelimiter)
{


  std::vector< SnpModelConverterValidation >::const_iterator it;
  std::map< std::string, std::string > newHeaderMap;
  
  // FIX UP DERIVED HEADERS
  if(derivedMap.count(ALGORITHM_NAME)) {
    std::string algorithmName(derivedMap[ALGORITHM_NAME]);
    for(int i = 0; i < SnpModelConverter::LAST_MODEL_TYPE; i++) {
      if(algorithmName == SnpModelFileType_str[i]) {
        switch(i) {
        case SnpModelConverter::BIRDSEED_PRIOR:
          derivedMap[ALGORITHM_NAME] = SNP_MODEL_ALGORITHM_TYPE_GTC[SMC_GH_ALGORITHM_BIRDSEED];
          derivedMap[MODEL_TYPE] = SNP_MODEL_MODEL_TYPE_GTC[SMC_GH_MODEL_PRIOR];
          break;
        case SnpModelConverter::BIRDSEED_POSTERIOR:
          derivedMap[ALGORITHM_NAME] = SNP_MODEL_ALGORITHM_TYPE_GTC[SMC_GH_ALGORITHM_BIRDSEED];
          derivedMap[MODEL_TYPE] = SNP_MODEL_MODEL_TYPE_GTC[SMC_GH_MODEL_POSTERIOR];
          break;
        case SnpModelConverter::BRLMMP_MODEL:
        case SnpModelConverter::BRLMMP_MODEL_HDF5:
          derivedMap[ALGORITHM_NAME] = SNP_MODEL_ALGORITHM_TYPE_GTC[SMC_GH_ALGORITHM_BRLMMP];
          derivedMap[MODEL_TYPE] = SNP_MODEL_MODEL_TYPE_GTC[SMC_GH_MODEL_UNKNOWN];
          break;
        case SnpModelConverter::SQLITE_MODEL_FILE:
        case SnpModelConverter::UNKNOWN_FILE:
        default:
          break;
        }
        break;
      }
    }
  }
  // FIX UP HEADER MAP VALUES
  if ( !headerMap.count(FILE_TYPE) ) {
    headerMap[FILE_TYPE] = new std::string(FILE_TYPE_MODELS);
  }
  else if(derivedMap.count(MODEL_TYPE) && !headerMap.count(MODEL_TYPE)) {
      headerMap[MODEL_TYPE] = new std::string(derivedMap[MODEL_TYPE]);
  }
  else {
    headerMap[MODEL_TYPE] = new std::string(SNP_MODEL_MODEL_TYPE_GTC[SMC_GH_MODEL_UNKNOWN]);
  }

  if ( !headerMap.count(REAGENT_KEY_VERSION) ) {
    headerMap[REAGENT_KEY_VERSION] = new std::string("");
  }
  // Always defer to the data version
  if(derivedMap.count(ALGORITHM_NAME) ) {
    headerMap[ALGORITHM_NAME] = new std::string(derivedMap[ALGORITHM_NAME]);
  }

  if(derivedMap.count(CHIP_TYPE) && !headerMap.count(CHIP_TYPE) ) {
    headerMap[CHIP_TYPE] = new std::string(derivedMap[CHIP_TYPE]);
  }

  //All validated headers.
  //The validation routine matched all ambiguous headers, not
  //just the duplicates.
  //File headers always take precedence. 
  for(it = validations.begin(); it != validations.end(); it++) {
    if ( (it->m_name == ALGORITHM_NAME) && headerMap.count(ALGORITHM_NAME) ) {
      continue;
    }
    if ( !it->m_value.empty() && (it->m_name != FILE_TYPE) ) {
      if ( headerMap[it->m_name] != NULL ) {
        delete headerMap[it->m_name];
      }
      headerMap[it->m_name] = new std::string(it->m_value);
      if((it->m_name == MODEL_TYPE)  && pcrecpp::RE("\\d+").FullMatch(it->m_value) && headerMap.count("SnpPosteriorFormatVer") ){
        *(headerMap[it->m_name]) = SNP_MODEL_MODEL_TYPE_GTC[SMC_GH_MODEL_POSTERIOR];
      }
    }
  }

}
// end addDerivedHeadersCallback_GTC4_1
/*****************************************************************************/
/*****************************************************************************/
/**
 * SnpModelConverter::GetModelHeaders, Add Multiple Values Delimited HELPER
 * _SnpModelConverter_AddHeaderToHeaderMap
 *
 * @param modelFilePath - the models file, posterior or prior
 * @param headerMap     - returned with updated values for key.
 * @param key           - header key to update.
 * @param value         - possible multi-value to append to the key's
 *                        value.
 * @param headerMapValueDelimiter
 *                      - multiple value delimiter
 *
 * @return headerMap    - key/value pair set..
 *
 */
/*****************************************************************************/
static void _SnpModelConverter_AddHeaderToHeaderMap(std::map< std::string, std::string *> & headerMap, const char * headerMapValueDelimiter, const std::string & key, const std::string & value)
{
  if(headerMap.count(key) == 0) {
    headerMap[key] = new std::string(value);
  } else if(headerMapValueDelimiter != NULL) {
    std::stringstream newVal;
    newVal << *(headerMap[key]) << headerMapValueDelimiter << value;
    *(headerMap[key]) = newVal.str();
  }
}

static void _SnpModelConverter_AddHeaderToHeaderMap(std::vector< std::string > & headerMatch, const char * headerMapValueDelimiter, const std::string & key, const std::string & value)
{

  if (  headerMatch.size() != 2) {
    return;
  }
  
  std::map< std::string, std::string *> dummy;

  dummy[headerMatch[0]] = &headerMatch[1];
  _SnpModelConverter_AddHeaderToHeaderMap(dummy, headerMapValueDelimiter, key, value);
  
}
// end _SnpModelConverter_AddHeaderToHeaderMap
/*****************************************************************************/
/*****************************************************************************/
/**
 * SnpModelConverter::GetModelHeaders, Hdf5 header and data parser HELPER.
 * _SnpModelConverter_Hdf5GetModelHeaders
 *
 * @param modelFilePath - the models file, posterior or prior
 * @param derivedMap     - returned derived headers.
 * @param headerMapValueDelimiter
 *                      - returned header values multiple value delimiter
 *
 * @return              - true if the file headers qualify as a models file,
 *                        false otherwise.
 * @return isModelData  - true if first record of data validates as
 *                        the correct schema.
 * @return headerMap    - filled in with file headers.
 * @return derivedMap   - possibley filled in with algorthem and model type
 *                        derived headers.
 *
 */
/*****************************************************************************/
static bool _SnpModelConverter_Hdf5GetModelHeaders(bool & isModelData, const std::string & modelFilePath, std::map< std::string, std::string *> & headerMap, const char * headerMapValueDelimiter, std::map< std::string, std::string > & derivedMap)
{

  affx::File5_File* file5 = NULL;
  affx::File5_Tsv* tsv5 = NULL;
  bool isModel = true;
  std::string hdf5path;

  file5 = new affx::File5_File();
  file5->open(modelFilePath, affx::FILE5_OPEN_RO);

  affx::File5_Group* group5 = NULL;
  group5 = file5->openGroup("/", affx::FILE5_OPEN_RO);

  vector<string> names = group5->listNames();
  for(int i = 0; i < names.size(); i++) {
    if(names[i].find(".snp-posteriors") != string::npos) {
      hdf5path = names[i];
      break;
    }

  }
  group5->close();
  Freez(group5);
  if(hdf5path.empty()) {
    isModel = false;
  } else {

    tsv5 = file5->openTsv(hdf5path);

    string key, value;
    for(int i = 0; i < tsv5->getHeaderCount(); i++) {
      tsv5->getHeaderByIdx(i, &key, &value);
      _SnpModelConverter_AddHeaderToHeaderMap(headerMap, headerMapValueDelimiter, key, value);
    }

    // DATA CHECK
    snp_labeled_distribution p;
    p.Dist.Clear();
    if((tsv5->nextLevel(0) == affx::FILE5_OK) &&
        (tsv5->get(0, 0, &p.probeset_id) > -1)) {

      QuantLabelZ__readSnpPrior_tsv5_v2(tsv5, p.Dist);
      isModelData = (p.Dist.aa.m   != 0.0 ||
                     p.Dist.aa.ss  != 0.0 ||
                     p.Dist.aa.k   != 0.0 ||
                     p.Dist.aa.v   != 0.0 ||
                     p.Dist.aa.ym  != 0.0 ||
                     p.Dist.aa.yss != 0.0 ||
                     p.Dist.aa.xyss != 0.0);
    }

    if(!isModelData) {
      Verbose::out(4, "TSV5 snp data row does not match a valid format.");
    }

    // Currently only BLRMMP is suppported with A5 files.
    derivedMap[ALGORITHM_NAME] = SNP_MODEL_ALGORITHM_TYPE_GTC[SMC_GH_ALGORITHM_BRLMMP];
    derivedMap[MODEL_TYPE] = SNP_MODEL_MODEL_TYPE_GTC[SMC_GH_MODEL_UNKNOWN];

    tsv5->close();
    delete tsv5;
  }
  file5->close();
  delete file5;

  return isModel;
}
// end _SnpModelConverter_Hdf5GetModelHeaders
/*****************************************************************************/
/*****************************************************************************/
/**
 * SnpModelConverter::GetModelHeaders, SQL header and data parser HELPER.
 * _SnpModelConverter_SqlGetModelHeaders
 *
 * @param modelFilePath - the models file, posterior or prior
 * @param headerMapValueDelimiter
 *                      - returned header values multiple value delimiter
 *
 * @return              - true if the file headers qualify as a models file,
 *                        false otherwise.
 * @return isModelData  - true if first record of data validates as
 *                        the correct schema.
 * @return headerMap    - filled in with file and SQL specific headers.
 *
 */
/*****************************************************************************/
static bool _SnpModelConverter_SqlGetModelHeaders(bool & isModelData, const std::string & modelFilePath, std::map< std::string, std::string *> & headerMap, const char * headerMapValueDelimiter)
{


  // DATA CHECK
  // The SnpModelDb::isSqlModelFile call validates schema.
  isModelData = SnpModelDb::isSqlModelFile(modelFilePath);

  if(!isModelData) {
    return isModelData;
  }

  bool isModelHeader = true;

  try {
    SQLiteDatabase sqliteConn;

    sqliteConn.open(modelFilePath, true);
    std::string key;
    std::string value;
    SQLiteRecordset rset(sqliteConn);
    string sql = "select key, value from snp_model_meta;";
    rset.open(sql);
    while(rset.fetch()) {
      key = rset.getString(0);
      value = rset.getString(1);
      _SnpModelConverter_AddHeaderToHeaderMap(headerMap, headerMapValueDelimiter, key, value);
    }
    rset.close();
    sqliteConn.close();
  } catch(...) {
    isModelHeader = false;
  }
  return isModelHeader;
}
// end _SnpModelConverter_SqlGetModelHeaders
/*****************************************************************************/
/*****************************************************************************/
/**
 * SnpModelConverter::GetModelHeaders, TSV header and data parser HELPER.
 * _SnpModelConverter_TsvGetModelHeaders
 *
 * @param modelFilePath - the models file, posterior or prior
 * @param derivedMap     - returned derived headers.
 * @param headerMapValueDelimiter
 *                      - returned header values multiple value delimiter
 *
 * @return              - true if the file headers qualify as a models file,
 *                        false otherwise.
 * @return isModelData  - true if first record of data validates as
 *                        the correct schema.
 * @return headerMap    - filled in with file headers.
 * @return derivedMap   - possibley filled in with algorthem and model type
 *                        derived headers.
 *
 */
/*****************************************************************************/
static bool _SnpModelConverter_TsvGetModelHeaders(bool & isModelData, const std::string & modelFilePath, std::map< std::string, std::string *> & headerMap, const char * headerMapValueDelimiter, std::map< std::string, std::string > & derivedMap)
{

  affx::TsvFile tsv;

  // OPEN FILE
  tsv.m_optAutoTrim = true;
  tsv.m_optAbortOnError = false;

  try {
    int rv = tsv.openTable(modelFilePath);
    if(rv != affx::TSV_OK) {
      return false;
    }
    // DATA ACCESS: Make sure we can access data.
    std::string col;
    rv = tsv.nextLevel(0);
    if(rv != affx::TSV_OK) {
      return false;
    }

    // FIRST ROW SCHEMA CHECK
    std::vector< std::string > cols;
    for(int i = 0; i < tsv.getColumnCount(0); i++) {
      std::string col;
      tsv.get(0, i, col);
      cols.push_back(col);
    }

    SnpModelConverterAlgorithmColumns *algorithmColumns = NULL;
    for(int i = 0 ; !algorithmColumns && (i < SMC_ALGORITHM_COLUMNS_SIZE); i++) {
      const SnpModelConverterAlgorithmColumns &acr = SMC_ALGORITHM_COLUMNS[i];
      if(tsv.getColumnCount(0) < acr.m_numColumnsToCheck) {
        continue;
      }
      bool match = true;
      for(int j = 0; match && (j < acr.m_numColumnsToCheck); j++) {
        switch(j) {
        case 0:
          match = acr.m_col1.FullMatch(cols[j]);
          break;
        case 1:
          match = acr.m_col2.FullMatch(cols[j]);
          break;
        case 2:
          match = acr.m_col3.FullMatch(cols[j]);
          break;
        case 3:
          match = acr.m_col4.FullMatch(cols[j]);
          break;
        case 4:
          match = acr.m_col5.FullMatch(cols[j]);
          break;
        default:
          break;
        }
      }
      if(match) {
        algorithmColumns = (SnpModelConverterAlgorithmColumns *)&acr;
      }
    }
    if(algorithmColumns != NULL) {
      Verbose::out(4, "TSV file first row of data validated as a models file.");
      derivedMap[ALGORITHM_NAME] = algorithmColumns->m_algorithm;
      derivedMap[MODEL_TYPE] = algorithmColumns->m_type;
      isModelData = true;
    } else {
      isModelData = false;
      Verbose::out(4, "TSV Column headers or first row do not match valid formats.");
    }

    // GET HEADERS: and comments.
    headerMap.clear();
    tsv.headersBegin();
    std::string key, val;
    TsvFileHeaderLine *next = NULL;
    while((next = tsv.nextHeaderPtr()) != NULL) {
      key.clear();
      val.clear();
      if(next->m_key.empty()) {
        key = next->m_value;
        val = "comment only.";
      } else {
        key = next->m_key;
        val = next->m_value;
      }
      _SnpModelConverter_AddHeaderToHeaderMap(headerMap, headerMapValueDelimiter, key, val);

    }

    tsv.clear();
  }
  catch (...) {
    return false;
  }

  return true;
}
// end _SnpModelConverter_TsvGetModelHeaders
/*****************************************************************************/
/*****************************************************************************/
/**
 * SnpModelConverter::GetModelHeaders, VALIDATION HELPER
 * _SnpModelConverter_ValidateGetModelHeaders
 *
 * @param modelFilePath - the models file, posterior or prior
 * @param headerMapType - 1.) empty for just the raw headers
 *                        2.) "GTC4.1" - GTC 4.1 map
 * @param headerMapValueDelimiter
 *                      - returned header values multiple value delimiter
 * @param derivedMap     - Derived headers.
 *
 * @return              -  @validations filled in with normalized header values.
 * @return              - true if the file headers qualifyas a models file,
 *                        false otherwise.
 *
 */
/*****************************************************************************/
static bool _SnpModelConverter_ValidateGetModelHeaders(const std::string modelFilePath, const std::map< std::string, std::string *> & headerMap, std::map< std::string, std::string > & derivedMap,  const char * headerMapValueDelimiter, std::vector< SnpModelConverterValidation > & validations)
{

  bool isModel = true;

  std::map < std::string, std::string >::const_iterator its;
  for ( its = derivedMap.begin(); its != derivedMap.end(); its++ ) {
    Verbose::out(5, ToStr("derivedMap ") +  its->first +  "=" + its->second );
  }

  // Ok, validate only the first occurence
  // In the case of multiple headers matching we ignore
  // everything but the first.
  for(int i = 0; i < SMC_VALIDATION_SIZE ; i ++) {
    std::map< std::string, std::string * >::const_iterator it;
    std::vector < std::string > derivedMatch;
    std::vector < std::string > fileMatch;
    std::vector < std::string > headerMatch;

    
    // DERIVED FROM DATA
    if(derivedMap.count(SMC_VALIDATION[i].m_name) && !derivedMap[SMC_VALIDATION[i].m_name].empty()) {
      derivedMatch.push_back(SMC_VALIDATION[i].m_name);
      derivedMatch.push_back(derivedMap[SMC_VALIDATION[i].m_name]);
    }

    // HEADER MATCH:  no choice but to walk through the list.
    for(it = headerMap.begin(); it != headerMap.end(); it++) {
      if(SMC_VALIDATION[i].m_validHeaderRE.FullMatch(it->first)) {
        if (headerMatch.empty() ) {
          headerMatch.push_back(SMC_VALIDATION[i].m_name);
          headerMatch.push_back(*(it->second));
        }
        else if ( headerMapValueDelimiter != NULL ) {
          _SnpModelConverter_AddHeaderToHeaderMap(headerMatch, headerMapValueDelimiter, SMC_VALIDATION[i].m_name, *(it->second));
        }
        else {
          break;
        }
      }
    }

    // FILE NAME MATCH
    std::string fileToMatch = modelFilePath;
    if(headerMap.count(FILE_INPUT)) {
      it = headerMap.find(FILE_INPUT);
      fileToMatch = *(it->second);
    }

    if(SMC_VALIDATION[i].m_validFileNameRE.FullMatch(Fs::basename(fileToMatch))) {
      fileMatch.push_back(SMC_VALIDATION[i].m_name);
      fileMatch.push_back(SMC_VALIDATION[i].m_value);
    }

    // VALIDATATION RESULTS CHECK
    if(!headerMatch.size() && !derivedMatch.size() && !fileMatch.size()) {
      if(! SMC_VALIDATION[i].m_emptyOk) {
        Verbose::out(4, SMC_VALIDATION[i].m_name +
                     ToStr(" header missing and is required. ") +
                     SMC_VALIDATION[i].m_description);
        isModel = false;
      }
    } else if(!headerMatch.size() || (headerMatch[1].empty())) {
      if(derivedMatch.size()) {
        headerMatch = derivedMatch;
      } else if(fileMatch.size()) {
        headerMatch = fileMatch;
      }
    }
    // SAVE VALIDATION FOR OUTPUT IF OK
    if(! headerMatch.empty()) {
      validations.push_back(SMC_VALIDATION[i]);
      validations.back().m_value = headerMatch[1];
    }
  }

  return isModel;
}
// end _SnpModelConverter_ValidateGetModelHeaders
/*****************************************************************************/
/*****************************************************************************/
/**
 * SnpModelConverter::GetModelHeaders
 * Given a potentional model file return the headers.
 * Model files can be of type db, A5, or TSV.
 *
 * @param modelFilePath - the models file, posterior or prior
 * @param headerMapType - 1.) empty for just the raw headers
 *                        2.) "GTC4.1" - GTC 4.1 map
 * @param headerMap     - return header key/value pairs
 * @param headerMapValueDelimiter
 *                      - return multiple header values using delimiter
 * If this value is passed in then all ambiguous matching
 * header values are returned delimited by the passed in value.
 * If this value is not passed (NULL) then only the first value
 * of the first qualifying header is returned.
 *
 * @return              - true if the file qualifies as a models file,
 *                        false otherwise.
 */
/*****************************************************************************/
bool SnpModelConverter::GetModelHeaders(const std::string & modelFilePath, const std::string & headerMapType, std::map< std::string, std::string *> & headerMap, const char * headerMapValueDelimiter)
{

  bool isModelHeader = false;
  bool isModelData = false;
  int derivedType = SnpModelConverter::NONE;
  std::map< std::string, std::string > derivedMap;
  std::vector< std::string > arrayResults;
  std::string uncModelFilePath = Fs::Unc(modelFilePath);

  Verbose::out(4, "SnpModelConverter::GetModelHeaders BEGIN");
  // GET HEADERS
  if(!Fs::fileExists(uncModelFilePath)) {
    APT_ERR_ABORT(ToStr(" SnpModelConverter::GetModelHeaders : ") + uncModelFilePath + " file not found.");
  }

  if(!headerMapType.empty()) {
    bool okType = false;
    for(int i = 0; i < SnpModelConverter::LAST_DERIVED_TYPE; i++) {
      if(headerMapType == SNP_MODEL_DERIVED_HEADER_TYPE_STRINGS[i]) {
        derivedType = i;
        okType = true;
        break;
      }
    }
    if(!okType) {
      Err::errAbort(ToStr("SnpModelConverter::GetModelHeaders invalid header type passed: ") + headerMapType);
    }
  }
    
  if(File5_File::isHdf5file(uncModelFilePath)) {
    isModelHeader = _SnpModelConverter_Hdf5GetModelHeaders(isModelData, uncModelFilePath, headerMap, headerMapValueDelimiter, derivedMap);
  } else if(SnpModelDb::isSqlModelFile(uncModelFilePath)) {
    isModelHeader = _SnpModelConverter_SqlGetModelHeaders(isModelData, uncModelFilePath,  headerMap, headerMapValueDelimiter);
  } else if(  !Fs::isBinaryFile( uncModelFilePath ) ) {
    isModelHeader = _SnpModelConverter_TsvGetModelHeaders(isModelData, uncModelFilePath,  headerMap, headerMapValueDelimiter, derivedMap);
  }

  if(isModelData) {
    SnpModelConverter dummy;
    try {
      dummy.getArrayTypeAndAlgorithm(uncModelFilePath, arrayResults);
      if(! derivedMap.count(CHIP_TYPE)) {
        if(arrayResults[0].empty()) {
          derivedMap[CHIP_TYPE] = UNKNOWN_CHIP_TYPE;
        } else {
          derivedMap[CHIP_TYPE] = arrayResults[0];
        }
      }
      if(! derivedMap.count(ALGORITHM_NAME)) {
        derivedMap[ALGORITHM_NAME] = arrayResults[1];
      }
    }
    catch (...) {
    }
  }

  // VALIDATE HEADERS
  std::vector< SnpModelConverterValidation > validations;
  isModelHeader = isModelHeader && _SnpModelConverter_ValidateGetModelHeaders(uncModelFilePath, headerMap, derivedMap,  headerMapValueDelimiter, validations);

  // DERIVED HEADERS
  if(isModelHeader && (derivedType != SnpModelConverter::NONE)) {
    (*SMC_ADD_DERIVED_HEADERS_CALLBACK[derivedType].m_callback)(headerMap, validations, derivedMap, headerMapValueDelimiter);
  }

  Verbose::out(4, "SnpModelConverter::GetModelHeaders END");

  return isModelHeader && isModelData;
}
// end SnpModelConverter::GetModelHeaders
/*****************************************************************************/
///////////////////////////////////////////////////////////////////////////////
// SnpModelConverter::GetModelHeaders CODE END
///////////////////////////////////////////////////////////////////////////////
