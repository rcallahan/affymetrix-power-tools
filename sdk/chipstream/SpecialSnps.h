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

/// @file SpecialSnps.h

#ifndef _SPECIALSNPS_H
#define _SPECIALSNPS_H

#include "chipstream/ChipLayout.h"
//
#include <cstring>
#include <map>
#include <string>
//


typedef std::map<std::string, std::vector<int> > CopyNumberMap; 

// This function will read in a file which has copy number data for a set of samples. The first row
// will consist of the string "marker_id" followed by strings representing the names of the samples
// being analyzed.  The first column will be assumed to have a string value representing the marker 
// name. Subsequent columns will be integer values representing the copy number values for the 
// sample/marker.  
void readCopyNumberFile(CopyNumberMap &copyNumberMap,
						const std::string& copyNumberFilePath);


typedef std::map<std::string,std::pair<int, int> > SpecialSnpMap;
typedef std::map<std::string,std::pair<int, int> >::const_iterator SpecialSnpMapIter;
typedef std::map<std::string,std::vector<int> > SpecialSampleSnpMap;
typedef std::vector<std::string> SpecialSampleNameVector;

// Return a map where the key is the SNP name, and the value is a pair of ints:
// {# of copies for males, # of copies for females}
void readSpecialSnpsFile(SpecialSnpMap *specialSnps,
                         std::string *chipType,
                         std::string *specialSnpVersion,
                         const std::string& specialSnpsPath);

void readSpecialSampleSnpsFile(	SpecialSampleSnpMap *specialSampleSnps,
			 	SpecialSampleNameVector *specialSampleName,
			 	const std::string& specialSampleSnpsPath);

// For backward-compatibility with old chrX files
void makeSpecialSnpsFromChrXFile(SpecialSnpMap *specialSnps,
                                 const std::string& chrXPath);

// convert between in-memory formats
void specialFromHaploid(SpecialSnpMap &SpecialSnps,
                        std::map<std::string,bool> &haploidSnps);
void haploidFromSpecial(SpecialSnpMap &SpecialSnps,
                        std::map<std::string,bool> &haploidSnps );

// read in chrXsnps file format
void fillInHaploidSnps(std::string file,
                       std::map<std::string,bool> &haploidSnps,
                       ChipLayout *layout,
                       bool permissive);
void fillInHaploidSnpsVersion2(std::string file,
                               std::map<std::string,bool> &haploidSnps,
                               ChipLayout *layout,
                               bool permissive);
void fillInHaploidSnpsVersion1(std::string file, 
                               std::map<std::string,bool> &haploidSnps, 
                               ChipLayout *layout, 
                               bool permissive);

/**
 * @brief Determines if special snp is from chrX based on male/female cn
 */
bool isChrX(const SpecialSnpMapIter& specialSnp);

/**
 * @brief Determines if special snp is from chrY based on male/female cn
 */
bool isChrY(const SpecialSnpMapIter& specialSnp);

/**
 * @brief Determines if special snp is from chrX or chrY
 */
bool isHaploid(const SpecialSnpMapIter& specialSnp);

/**
 * @brief Determines if given probeset_name is a special snp and if it is from chrX or chrY
 */
bool isHaploid(const SpecialSnpMap& snpMap, const std::string& probeset_name);

#endif /* _SPECIALSNPS_H */

/******************************************************************/
/**************************[END OF SpecialSnps.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
