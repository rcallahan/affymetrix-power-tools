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
 * @file   SnpModelConverter.h
 * @author Chuck Sugnet
 * @date   Mon Mar 15 10:31:04 2010
 * 
 * @brief  Class collecting utility methods for reading various formats of snp models like:
 * - Brlmm-p (aka labelz) used for snp5, snp6 and Axiom
 * - Birdseed priors used for snp 6
 * - Birdseed posteriors from snp 6
 * 
 * Additional details about file formats etc. are available on the internal wiki:
 * http://infowiki.ev.affymetrix.com/index.php/GenotypingModelFiles
 * 
 */


#ifndef SNP_MODEL_CONVERTER_H
#define SNP_MODEL_CONVERTER_H

#include "affy-pcre.h"
#include "label/snp.label.h"
#include "util/SQLite.h"
#include "file/TsvFile/TsvFile.h"
#include "chipstream/QuantLabelZ.h"
#include <string>

/**
 * Going forward we want to put:
 * 'official-affy-model' in files
 * 'chip-analysis-type' in file
 * 'total-model-count' in file
 * 'model-type' in file.
 */


extern const std::string FILE_INPUT;
extern const std::string TOTAL_MODEL_COUNT;
extern const std::string USED_FOR_PRIOR;
extern const std::string OFFICIAL_AFFY;
extern const std::string CHIP_ANALYSIS_TYPE;
extern const std::string CHIP_TYPE;
extern const std::string ARRAY_TYPE;
extern const std::string MODEL_TYPE;
extern const std::string ALGORITHM_NAME;

/**
 * @brief  Class collecting utility methods for reading various formats of snp models like:
 * - Brlmm-p (aka labelz) used for snp5, snp6 and Axiom
 * - Birdseed priors used for snp 6
 * - Birdseed posteriors from snp 6
 * 
 * Additional details about file formats etc. are available on the internal wiki:
 * http://infowiki.ev.affymetrix.com/index.php/GenotypingModelFiles
 * 
 */

extern const char *SnpModelFileType_str[];

// Possible type strings SnpModelConverter::GetModelHeaders.
// Correspends SnpDerivedHeaderType enum.
const std::string SNP_MODEL_DERIVED_HEADER_TYPE_STRINGS[] = {
  std::string("none"),
  std::string("GTC4.1"),
};

class SnpModelConverter {
    
public:

  enum SnpDerivedHeaderType {
    NONE,
    GTC4_1,
    LAST_DERIVED_TYPE,
  };
  
    enum SnpChipType {
        UNKNOWN_CHIP,
        SNP5,
        SNP6,
        AXIOM,
        LAST_CHP_TYPE,
    };

    enum SnpModelFileType {
        UNKNOWN_FILE,
        BIRDSEED_PRIOR,
        BIRDSEED_POSTERIOR,
        BRLMMP_MODEL,
        BRLMMP_MODEL_HDF5,
        SQLITE_MODEL_FILE,
        LAST_MODEL_TYPE
    };

    SnpModelConverter(enum SnpChipType chip=UNKNOWN_CHIP, enum SnpModelFileType model=UNKNOWN_FILE) : 
        m_ChipType(chip), m_FileType(model) {}

   static bool GetModelHeaders( const std::string & modelFilePath,
                                const std::string & headerMapType,
                                std::map< std::string, std::string *> & headerMap,
                                const char * headerMapValueDelimiter = NULL);
   

    bool vecContains(const std::string &needle, const std::vector<std::string> &haystack);

	enum SnpModelConverter::SnpChipType guessChipType(enum SnpModelConverter::SnpModelFileType modelType,
                                                      const std::vector<std::string> &chipTypes, 
                                                      const std::string &firstPsName);

    int guessNumProbesets(enum SnpModelConverter::SnpChipType type);

    enum SnpModelFileType determineFileType(const std::string &fileName);

    void convertToDbModel(const std::string &fileIn, const std::string &fileOut, const std::string &file5Path="");

    void createTables(SQLiteDatabase &db);

    void convertBrlmmPToDbModel(const std::string &fileIn, const std::string &fileOut);

    sqlite3_stmt *getInsertStatement(SQLiteDatabase &db);

    void writeSnpDistribution(const snp_labeled_distribution &p, sqlite3_stmt &stmt, SQLiteDatabase &db);

    void headersToMeta(affx::TsvFile &tsv, SQLiteDatabase &db);

    void convertBirdseedPriorToDbModel(const std::string &fileIn, const std::string &fileOut);

    void birdseedModelFromStrings(const std::string &name,
                                  int copynuber,
                                  const std::string &bb, 
                                  const std::string &ab, 
                                  const std::string &aa,
                                  char sep,
                                  snp_labeled_distribution &p);

    void birdseedStringToCluster(const std::string &s, char sep, cluster_data &d);

    void convertBirdseedPosteriorToDbModel(const std::string &fileIn, const std::string &fileOut);

    void convertBrlmmPFile5ToDbModel(const std::string &fileIn, const std::string &file5Path, const std::string &fileOut);

    enum SnpModelFileType guessTsvFileType(const std::string &file);

    enum SnpModelFileType guessFileType(const std::string &file);

    void addBirdseedPosteriorFromStrings(std::string &id, const std::string &aa, 
                                         const std::string &ab, const std::string &bb,
                                         snp_labeled_distribution &p);

	void getArrayTypeAndAlgorithm(const std::string &fileIn, std::vector<std::string> &result);

private:
    
    enum SnpChipType m_ChipType;
    enum SnpModelFileType m_FileType;

	std::string getHd5ArrayTypeAndAlgorithm(const std::string &fileIn, enum SnpModelFileType ft, std::string &alg);
	std::string getTsvArrayTypeAndAlgorithm(const std::string &fileIn, enum SnpModelFileType ft, std::string &alg);

};

#endif /* SNP_MODEL_CONVERTER_H */
