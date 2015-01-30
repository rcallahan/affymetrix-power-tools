////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
// ~/Affy/apt2/trunk/BioSpecies.cpp ---
//
// $Id: BioSpecies.cpp,v 1.1 2009-10-28 00:06:03 harley Exp $
//

#include "BioSpecies.h"

BioSpecies::BioSpecies() {
}
BioSpecies::~BioSpecies() {
}

BioSpecies* BioSpecies::speciesFromName(const std::string& name)
{
  if (name=="human") {
    return new BioSpeciesHuman();
  }
  if (name=="chimp") {
    return new BioSpeciesChimp();
  }
  if (name=="dicty") {
    return new BioSpeciesDicty();
  }
  return NULL;
}

bool BioSpeciesHuman::isChrAutosome(int chr)
{
  return (chr<22);
}
    
bool BioSpeciesChimp::isChrAutosome(int chr)
{
  return (chr<22);
}
    
bool BioSpeciesDicty::isChrAutosome(int chr)
{
  return false;
}


    
    
