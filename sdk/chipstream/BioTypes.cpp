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


/// @file   BioTypes.cpp
/// @brief  Functions to work on ChipStream BioTypes.

//
#include "chipstream/BioTypes.h"
//
#include "file/CHPFileData.h"
#include "util/Convert.h"
#include "util/Err.h"

//
using namespace std;
using namespace affx;

//////////

/*
 * Convert the enumerant to a value used by the CHP files
 */
int affx::GTypeCallToChpValue(affx::GType call)
{
    switch (call)	
    {
        case NN: return ALLELE_NO_CALL; break;
        case AA: return ALLELE_A_CALL; break;
        case AB: return ALLELE_AB_CALL; break;
        case BB: return ALLELE_B_CALL; break;
        default: Err::errAbort("Don't recognize call of type: " + ToStr(call));
    }
    // not reached, but return something for gcc.
    Err::errAbort("GTypeCallToChpValue: Internal error.");
    return ALLELE_NO_CALL;
}

/// @brief     Convert a GType enum to an external integer value
/// @param     gtype     The GType to convert
/// @return    the value to use for output
int
affx::GType_to_int(GType gtype)
{
  switch (gtype) 
    {
    case NN:
      return -1;
      break;
    case AA:
      return 0;
      break;
    case AB:
      return 1;
      break;
    case BB:
      return 2;
      break;
    default:
      Err::errAbort("Gtype_to_int: bad conversion");
      return 123; // prevent warning
      break;
    }
}

/// @brief     Convert an integer from a file to an internal GType
/// @param     gtype_int int to convert
/// @return    GType enum
GType
affx::GType_from_int(int gtype_int)
{
  switch (gtype_int) 
    {
    case -1:
      return NN;
      break;
    case 0:
      return AA;
      break;
    case 1:
      return AB;
      break;
    case 2:
      return BB;
      break;
    default:
      Err::errAbort("GType_from_int: bad conversion");
      return NN; // prevent warning
      break;
    }
}

/// @brief     Check that the gtype is AA,AB,BB,NN
/// @param     gtype     
/// @return    true if valid
int
 affx::GType_valid(int gtype)
{
  if ((gtype==AA)||(gtype==AB)||(gtype==BB)||(gtype==NN)) {
    return 1;
  }
  return 0;
}

/// @brief     Check that the gtype is AA,AB,BB  (not NN)
/// @param     gtype     
/// @return    true if called
int
affx::GType_called(int gtype)
{
  if ((gtype==AA)||(gtype==AB)||(gtype==BB)) {
    return 1;
  }
  return 0;
}

//////////

/// @brief     convert the BioType Gender to a string
/// @param     gender    an affx::Gender value
/// @return    pointer to a const string
const char*
affx::getGenderString(affx::Gender gender)
{
  // NOTE! Keep in sync with "enum Gender".
  const char* gender_string[]={"female", "male", "unknown"};
  //
  switch (gender) {
  case affx::Female:
    return gender_string[0];
  case affx::Male:
    return gender_string[1];
  case affx::UnknownGender:
    return gender_string[2];
  default:
    Err::errAbort("getGenderString: unhandled gender");
    // bogus return for gcc
    return NULL;
  }
}

/// @brief Change gender call Male <-> Female.  Unknown stays the same.  This
/// is used for chrZW gender calling.
/// @param     gender    an affx::Gender value
/// @return    gender    an affx::Gender value
affx::Gender affx::flipMaleFemaleGender (affx::Gender gender) {
    if (gender == affx::Male) {
        return affx::Female;
    }
    else if (gender == affx::Female) {
        return affx::Male;
    }
    else {
        return gender;
    }

//     return static_cast<affx::Gender>(((-1 * gender) + 4) % 3);

}
