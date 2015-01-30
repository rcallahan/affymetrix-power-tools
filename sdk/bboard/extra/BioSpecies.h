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
// ~/Affy/apt2/trunk/BioSpecies.h ---
//
// $Id: BioSpecies.h,v 1.1 2009-10-28 00:06:03 harley Exp $
//

// Author:  harley <harley@mahalito.net>
// Keywords: 

// Commentary:
// * This is an example of an object which might appear on
//   an apt2 blackboard.

#ifndef _BIOSPECIES_H_
#define _BIOSPECIES_H_

//
#include <string>

class BioSpecies {
public:
  BioSpecies();
  virtual ~BioSpecies();

  static BioSpecies* speciesFromName(const std::string& name);

  virtual bool isChrAutosome(int chr)=0;
};

class BioSpeciesHuman: public BioSpecies {
  bool isChrAutosome(int chr);
};

class BioSpeciesChimp: public BioSpecies {
  bool isChrAutosome(int chr);
};

class BioSpeciesDicty: public BioSpecies {
  bool isChrAutosome(int chr);
};

#endif // _BIOSPECIES_H_
