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
 * @file   ProbeSet.cpp
 * @author Chuck Sugnet
 * @date   Mon May 23 13:28:33 2005
 * 
 * @brief  Classes and functions for dealing with probe sets.
 */

//
#include "chipstream/ProbeSet.h"
//
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <fstream>
#include <iostream>

using namespace std;

Probe::Probe()
{
  clear();
}

void Probe::clear()
{
  // @todo "probeid_t" should be "int"
  id=-1;
  type=BLANK;
  gcCount=NULLPROBEGC;
  // m_parentAtom=NULL;
  m_apid=-1;
}

bool Probe::isNull() {
  return (id==-1);
}

/** 
 * Constructor for creating a probe from a vector of words. 
 * @param vec - vector of strings containing probe information.
 */
Probe::Probe(const std::vector<std::string> &vec, unsigned int offSet) {
  if(vec.size() - offSet < numFields()) {
    Err::errAbort("Wrong number of fields to initialize Probe. Expecting: " +
                    ToStr( numFields() ) + " Got: " +
                    ToStr( vec.size() ));
  }
  id = Convert::toInt(vec[offSet++]);
  type = typeForString(vec[offSet++].c_str());
  gcCount = (unsigned char)Convert::toInt(vec[offSet++]);
  offSet+=2; // ignoring length and interrogation
  //  seq = string(vec[offSet++]);
  m_apid = -1;
}

/** 
 * Constructor for creating a probe from a vector of words. 
 * @param vec - vector of strings containing probe information.
 */
Probe::Probe(const std::vector<const char *> &vec, unsigned int offSet) {
  if(vec.size() - offSet < numFields()) {
    Err::errAbort("Wrong number of fields to initialize Probe. Expecting: " +
                    ToStr(numFields()) + " Got: " + ToStr(vec.size()) );
  }
  id = Convert::toInt(vec[offSet++]);
  type = typeForString(vec[offSet++]);
  gcCount = (unsigned char)Convert::toInt(vec[offSet++]);
  offSet+=2; // ignoring length and interrogation
  //  seq = string(vec[offSet++]);
  m_apid = -1;
}

const char *Probe::stringForType(Probe::ProbeType t) {
    switch(t)
      {
      case PMST :
        return "pm:st";
      case MMST :
        return "mm:st";
      case PMAT :
        return "pm:at";
      case MMAT :
        return "mm:at";
      case GENERICAT :
        return "generic:at";
      case GENERICST :
        return "generic:st";
      case JUMBOAT :
        return "jumbo-checkerboard:at";
      case JUMBOST :
        return "jumbo-checkerboard:st";
      case THERMOAT :
        return "thermo:at";
      case THERMOST :
        return "thermo:st";
      case TRIGRIDAT :
        return "trigrid:at";
      case TRIGRIDST :
        return "trigrid:st";
      case BLANK :
        return "blank";
      default:
        Err::errAbort("Don't recognize Probe::pType: " + ToStr(t));
      }
    return NULL;
  }

Probe::ProbeType Probe::typeForString(const std::string& cs) {
  // @todo order these by expected frequency.
  if (cs=="pm:st") {
    return PMST;
  }
  if (cs=="mm:st") {
    return MMST;
  }
  if ((cs=="pm:at") || (cs=="pm:target->at")) {
    return PMAT;
  }
  if ((cs=="mm:at") || (cs=="mm:target->at")) {
    return MMAT;
  }
  if ((cs=="generic:at") || (cs=="generic:target->at")) {
    return GENERICAT;
  }
  if ((cs=="generic:st") || (cs=="generic:target->st")) {
    return GENERICST;
  }
  if ((cs=="jumbo-checkerboard:at") || (cs=="edge:target->at")) {
    return JUMBOAT;
  }
  if ((cs=="jumbo-checkerboard:st") || (cs=="edge:target->st")) {
    return JUMBOST;
  }
  if ((cs=="thermo:at") || (cs=="thermo:target->at")) {
    return THERMOAT;
  }
  if ((cs=="thermo:st") || (cs=="thermo:target->st")) {
    return THERMOST;
  }
  if ((cs=="trigrid:at") || (cs=="ngrid:target->at")) {
    return TRIGRIDAT;
  }
  if ((cs=="trigrid:st") || (cs=="ngrid:target->st")) {
    return TRIGRIDST;
  }
  if ((cs=="blank") || (cs=="bar:target->st") || (cs=="bar:target->at")) {
    return BLANK;
  }
  //
  APT_ERR_ABORT("Cant convert ' "+cs+"' to Probe::ProbeType.");
  // return something for the compiler.
  return BLANK;
}
