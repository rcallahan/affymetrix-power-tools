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

#ifndef _KITAOCELLISTENER_H_
#define _KITAOCELLISTENER_H_

//
#include "chipstream/CelListener.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/GenderCelListener.h"
#include "chipstream/KitAODb.h"
//
#include <string>
#include <utility>
#include <vector>

class KitAOCelListener : public ChipSummary, public CelListener { 
public:
  enum {
    KIT_UNKNOWN=0,
    KIT_A=1,
    KIT_O=2,
  };

  KitAOCelListener();
  virtual ~KitAOCelListener();
  //
  void init();
  void clear();

  //
  void newChip(affymetrix_fusion_io::FusionCELData* cel);

  // 
  void computeKitAO();

  //
  void setDb(KitAODb* db);
  void setLayout(ChipLayout* layout);

  // # The cutoff value (VAL) which will be used to classify Cels.
  // # Less than this negative value is KitA.  (id=1  [-inf..-KitA])
  // # In between is unknown.                  (id=0  (-KitA..+KITO));
  // # More than this postitive value is KitO. (id=2  [+VAL..+inf])
  // The database file must have a line like this:
  // #%classifier-cutoff-a-value=-3.31
  // #%classifier-cutoff-o-value=-1.34
  double getClassiferAValue() const;
  double getClassiferOValue() const;

  // what is passed into computeKitAO as member vars
  // we dont own this memory.
  KitAODb*    m_db;
  ChipLayout* m_layout;
  //
  std::string m_celFileName;
  std::string m_exprAnalysisString;
  std::string m_sketchFileName;
  // what is returned from computeKitAO as member vars
  int m_return_kit_id;
  double m_return_kit_value;
};

#endif // _KITAOCELLISTENER_H_
