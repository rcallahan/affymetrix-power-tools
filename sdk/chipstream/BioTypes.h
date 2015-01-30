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
 * @file   BioTypes.h
 * @author Chuck Sugnet
 * @date   Mon Jul 10 10:05:42 2006
 * 
 * @brief Central location enumerations and other types that aren't
 * particular to a single class or application.
 * 
 */

#ifndef _BIOTYPES_H_
#define _BIOTYPES_H_

// The length of the longest string returned by getGenderString()
#define MAX_GENDER_STRING_LENGTH 10

namespace affx {

  typedef unsigned char GType; // 1B
  
  /**  Types of Genotype calls made. */
  enum GType_enum {
    AA, ///< Homozygous A version of allele. (0)
    AB, ///< Heterozygous, one copy of both A and B. (1)
    BB, ///< Homozygous B version of allele. (2)
    NN,  ///< No call, unknown genotype.  The internal is value is '3'; external is '-1', 
    // Many places use '4' for the total calls.
    // Should we define TOTAL as 4?
    INVALID
  };

  int   GType_to_int(GType gtype);

  GType GType_from_int(int gtype_int);

  int   GType_valid(int gtype);

  int   GType_called(int gtype);

  /** Convert the enumerant to a value used by the CHP files
   * @param call The gtype call value.
   * @return The CHP file call value.
   */
  int GTypeCallToChpValue(GType call);

  //////////
  /** Those pesky XY people need a different prior for all the SNPs on
      X that aren't on the psuedoautosomal region. */
  
  enum Gender {
    Female,  ///< Has XX chromosomes
    Male,    ///< Has XY chromosomes
    UnknownGender  ///< Don't know genotype, use XX model.
  };
  
  //
  const char* getGenderString(affx::Gender gender);
  
  affx::Gender flipMaleFemaleGender(affx::Gender gender);
  

};

#endif /* _BIOTYPES_H_ */
