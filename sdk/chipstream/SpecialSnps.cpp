////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

/// @file SpecialSnps.cpp

//
#include "chipstream/SpecialSnps.h"
//
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "stats/stats.h"
#include "util/Fs.h"
#include "util/Util.h"

using namespace std;
using namespace affx;

void readCopyNumberFile(CopyNumberMap &copyNumberMap,
                        const std::string& copyNumberFilePath){

    // The format of the file that is read by this function is assumed to have the form: 
    /*
	  sample_id	<marker1>	<marker2>	etc.
	  <sample1>	<copyNumber>	<copyNumber>	etc.
	  <sample2> 	<copyNumber>	<copyNumber>	etc.	
	  .	    .		    .
	  .	    .  		    . 
      etc.	   etc.	           etc.
    */

//     affx::TsvFile tsv;
	affx::File5_File cn_file5;
    string sampleName;
	vector<string> markerNames;
    int markerCount;
    string label;

    /* sanity checks */
    if (cn_file5.open(copyNumberFilePath, affx::FILE5_OPEN) != FILE5_OK) {
        Err::errAbort("Couldn't open file: " + string(copyNumberFilePath) + " to read.");
    }
    string file_name = Fs::basename(copyNumberFilePath);
	Util::changeEnd(file_name, ".a5", "");
    affx::File5_Tsv *cn_file5_tsv = cn_file5.openTsv(file_name);
    cn_file5_tsv->getColumnName(0, 0, &label);
    if (label != "sample_id") {
		Err::errAbort("sample_id must be the first column header in the copyNumber file, instead of:'" + label + "'");
	}

    markerCount = cn_file5_tsv->getColumnCount(0); 

    for (int i=1; i<markerCount; ++i) {
		label.clear();
		cn_file5_tsv->getColumnName(0, i, &label);
		markerNames.push_back(label);
		vector<int> temp;
		copyNumberMap[label] = temp;
    }

    ///@todo cannot assume that the order is the same!
    while (cn_file5_tsv->nextLevel(0) == FILE5_OK) {
        for (int i=1; i<markerCount; i++) {
            int val;				 
            if (cn_file5_tsv->get(0, i, &val)==FILE5_OK) {
                copyNumberMap[markerNames[i-1]].push_back(val);
            }
            else {
                Err::errAbort("Unable to read column " + ToStr(i) + " from copyNumber file");
            }
        }
    }
    //
    cn_file5_tsv->close();
    delete cn_file5_tsv;
}


void readSpecialSnpsFile(SpecialSnpMap *specialSnps,
                         std::string *chipType,
                         std::string *specialSnpVersion,
                         const std::string& specialSnpsPath)
{
    affx::TsvFile tsv;
    string snpName;
    string chrom;
    int malecopy;
    int femalecopy;
    pair<int,int> duty;
    string specialType, chip;
    chipType->clear();
    specialSnpVersion->clear();
    tsv.bind(0, "probeset_id", &snpName, TSV_BIND_REQUIRED);
    tsv.bind(0, "chr", &chrom, TSV_BIND_REQUIRED);
    tsv.bind(0,"copy_male",&malecopy,TSV_BIND_REQUIRED);
    tsv.bind(0,"copy_female",&femalecopy,TSV_BIND_REQUIRED);
    /* sanity checks */
    if (tsv.open(specialSnpsPath) != TSV_OK) 
        Err::errAbort("Couldn't open file: " + string(specialSnpsPath) + " to read.");
    // set the special snps type.
    if (tsv.headersFindNext("special-snp-version", specialType) == TSV_OK) {
        *specialSnpVersion = specialType;
    }
    // read chipType from header
    while (tsv.headersFindNext("chip-type", chip) == TSV_OK) {
        *chipType = chip;
    }
    while (tsv.nextLevel(0) == TSV_OK) {
        duty.first=malecopy;
        duty.second=femalecopy;	
        (*specialSnps)[snpName] = duty;
    }
    tsv.close();
}

void readSpecialSampleSnpsFile(SpecialSampleSnpMap *specialSampleSnps,
                               SpecialSampleNameVector *specialSampleName,
                               const std::string& specialSampleSnpsPath)
{
    affx::TsvFile tsv;
    string snpName;
    vector<int> sample;
    int sampleCount;
    string specialType, chip;
    string label;

    /* sanity checks */
    if (tsv.open(specialSampleSnpsPath) != TSV_OK) 
        Err::errAbort("Couldn't open file: " + string(specialSampleSnpsPath) + " to read.");
    tsv.cidx2cname(0,0,label);
    if (label != "probeset_id")
        Err::errAbort("probeset_id must be the first column in the specialSampleSnps file");
    sampleCount = tsv.getColumnCount(0); 
    specialSampleName->resize(sampleCount);
    for (int i=1; i<sampleCount; ++i) {
        label.clear();
        tsv.cidx2cname(0, i, label);
        specialSampleName->push_back(label);
    }
    while (tsv.nextLevel(0) == TSV_OK) {
        snpName.clear();
        if (tsv.get(0,0,snpName)!=TSV_OK)
            Err::errAbort("Unable to read column " + ToStr(0) + " from specialSampleSnps file");
        sample.clear();
        for (int i=1; i<sampleCount; i++) {
            int val;
            if (tsv.get(0,i,val)==TSV_OK)
                sample.push_back(val);
            else
                Err::errAbort("Unable to read column " + ToStr(i) + " from specialSampleSnps file");
        }
        (*specialSampleSnps)[snpName] = sample;
    }
    tsv.close();
}

void makeSpecialSnpsFromChrXFile(SpecialSnpMap *specialSnps,const std::string& chrXPath)
{
    vector<string> chrXSNPNames;
    affx::TsvFile tsv;
    string snpName;

    tsv.bind(0, "all_chrx_no_par", &snpName, TSV_BIND_REQUIRED);
    /* sanity checks */
    if (tsv.open(chrXPath) != TSV_OK) 
        Err::errAbort("Couldn't open file: " + string(chrXPath) + " to read.");
    while (tsv.nextLevel(0) == TSV_OK) {
        chrXSNPNames.push_back(snpName);
    }
    pair<int,int> thePair(1, 2);
    for (vector<string>::const_iterator it = chrXSNPNames.begin();
         it != chrXSNPNames.end(); ++it) {
        (*specialSnps)[*it] = thePair;
    }
    tsv.close();
}

void haploidFromSpecial(SpecialSnpMap &SpecialSnps, map<string,bool> &haploidSnps ){
    SpecialSnpMapIter mapIter;
    /*for (mapIter = SpecialSnps.begin(); mapIter != SpecialSnps.end(); mapIter++) {
      printf("%s %d %d\n",mapIter->first.c_str(),mapIter->second.first,mapIter->second.second);
      }*/
    // grab off just the haploid snps
    // take a random sample of 10K (or actual size)
    vector<string> TmpName;
    for (mapIter = SpecialSnps.begin(); mapIter != SpecialSnps.end(); mapIter++) {
        if (mapIter->second.first==1 && mapIter->second.second==2) {
            // load up the hypothetical chrX snps
            TmpName.push_back(mapIter->first);
        }
    }
    int i;
    vector<string> SampleName;
    int Xsamplesize=10000;

    if (TmpName.size()>Xsamplesize) {
        SampleName.resize(Xsamplesize);
        reproducible_random_sample(TmpName.begin(), TmpName.end(),
                                   SampleName.begin(), SampleName.size(), 42);
    }
    else {
        SampleName.resize(TmpName.size());
        for (i=0; i<TmpName.size(); i++)
            SampleName[i]=TmpName[i];
    }
    haploidSnps.clear();
    for (i=0; i<SampleName.size(); i++) {
        haploidSnps[SampleName[i]]=true;
    }
}

void specialFromHaploid(SpecialSnpMap &SpecialSnps, map<string,bool> &haploidSnps){
    SpecialSnps.clear();
    map<string, bool >::iterator mapIter;
    pair<int,int> duty;
    for (mapIter = haploidSnps.begin(); mapIter != haploidSnps.end(); mapIter++) {

        duty.first=1;
        duty.second=2;
        SpecialSnps[mapIter->first]=duty;
    }
}

/**
 * Open a text file and read in snps that are possibly haploid. Currently
 * only chrX snps are supported.
 *
 * @param file        - chrX file to read
 * @param haploidSnps - Map of probeset names that are marked as chrX
 * @param layout      - 
 */
void fillInHaploidSnpsVersion1(string file, map<string,bool> &haploidSnps, ChipLayout *layout, bool permissive) {
    int cnt = 0;
    affx::TsvFile tsv;
    string snpName;
    tsv.bind(0, "all_chrx_no_par", &snpName, TSV_BIND_REQUIRED);
    /* sanity checks */
    if (tsv.open(file) != TSV_OK)
        Err::errAbort("Couldn't open file: " + file + " to read.");
    while (tsv.nextLevel(0) == TSV_OK) {
        if (permissive || (layout != NULL && layout->containsProbeSet(snpName.c_str())))
            haploidSnps[snpName] = true;
        else
            Err::errAbort("Can't find chrX snp: " + snpName + " in cdf file. Wrong chrX file?");
        cnt++;
    }
    tsv.close();
    Verbose::out(1,"Read " + ToStr(cnt) + " ChrX SNPs from chrX SNPs file " + file + " (v1)");
}

/**
 * Open a text file and read in snps that are possibly haploid. Currently
 * only chrX snps are supported.
 *
 * @param file        - chrX file to read
 * @param haploidSnps - Map of probeset names that are marked as chrX
 * @param layout      -
 * @param permissive  -
 */
void fillInHaploidSnpsVersion2(string file,
                               std::map<std::string,bool> &haploidSnps,
                               ChipLayout *layout,
                               bool permissive) {
    int cnt = 0;
    affx::TsvFile tsv;
    string snpName;
    string chrom;
    tsv.bind(0, "Probeset", &snpName, TSV_BIND_REQUIRED);
    tsv.bind(0, "chr", &chrom, TSV_BIND_REQUIRED);
    /* sanity checks */
    if (tsv.open(file) != TSV_OK)
        Err::errAbort("Couldn't open file: " + file + " to read.");
    while (tsv.nextLevel(0) == TSV_OK) {
        if (permissive || (layout != NULL && layout->containsProbeSet(snpName.c_str()))) {
            if (chrom == "X") {
                haploidSnps[snpName] = true;
                cnt++;
            }
        }
        else
            Err::errAbort("Can't find snp: " + snpName + " on chrom: " + chrom + " in cdf file. Wrong chrx file?");
    }
    tsv.close();
    Verbose::out(1,"Read " + ToStr(cnt) + " ChrX SNPs from chrX SNPs file " + file + " (v2)");
}

/**
 * Determine which version of a chrx snp have and open appropriately
 *
 * @param file        - chrX file to read
 * @param haploidSnps - Map of probeset names that are marked as chrX
 * @param layout      -
 * @param permissive  -
 */
void fillInHaploidSnps(string file,
                       std::map<std::string,bool> &haploidSnps,
                       ChipLayout *layout,
                       bool permissive) {
    affx::TsvFile tsv;
    string version;
    string key = "version";
    if (tsv.open(file) != TSV_OK)
        Err::errAbort("Couldn't open file: " + file + " to read.");
    // Is there an explicit setting for version?
    if (tsv.getHeader(key, version) == TSV_OK) {
        tsv.close();
        if (version == "2") {
            fillInHaploidSnpsVersion2(file, haploidSnps, layout, permissive);
        }
        else {
            Err::errAbort("In file: " + file + " Don't know how to handle version: " + version);
        }
    }
    // if not assume version one and move forward.
    else {
        tsv.close();
        fillInHaploidSnpsVersion1(file, haploidSnps, layout, permissive);
    }
}

bool isChrX(const SpecialSnpMapIter& specialSnp) {
    bool result = false;
    std::pair<int,int> maleFemaleCnPair = specialSnp->second;
    if (maleFemaleCnPair.first == 1 && maleFemaleCnPair.second == 2) {
        result = true;
    }
    return result;
}

bool isChrY(const SpecialSnpMapIter& specialSnp) {
    bool result = false;
    std::pair<int,int> maleFemaleCnPair = specialSnp->second;
    if (maleFemaleCnPair.first == 1 && maleFemaleCnPair.second == 0) {
        result = true;
    }
    return result;
}

bool isHaploid(const SpecialSnpMapIter& specialSnp) {
    bool result = false;
    if (isChrX(specialSnp) || isChrY(specialSnp)) {
        result = true;
    }
    return result;
}

bool isHaploid(const SpecialSnpMap& snpMap, const std::string& probeset_name) {
    bool result = false;
    SpecialSnpMapIter snpMapIter = snpMap.find(probeset_name);
    if (snpMapIter != snpMap.end() && isHaploid(snpMapIter)) {
        result = true;
    }
    return result;
}

/******************************************************************/
/**************************[END OF SpecialSnps.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
