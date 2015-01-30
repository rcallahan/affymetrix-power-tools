////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
 * @file   DMET2Model.h
 * @author Mybrid Spalding
 * @date   Wed Jun  4 10:18:22 PDT 2008
 * @brief  Allele translation table regression class for all DMET2 input files.
 */

#ifndef TRANSLATION_REGRESSION_DMET2_MODEL_H
#define TRANSLATION_REGRESSION_DMET2_MODEL_H

//
#include "translation/TranslationInputTsvTableModel.h"
#include "translation/regression/ATDRegression.h"
//
#include <cstring>
#include <set>
#include <string>
//

using namespace std;

class DMET2MarkerModel : public TranslationInputTsvTableModel
{
public:

  //Experiment Gene ProbeSet A1 A2 Ref Var Call
  enum {
    EXPERIMENT,
    GENE,
    PROBE_SET,
    A1,
    A2,
    REF,
    VAR,
    CALL,
  };

  set< string > m_experiments;
  map< string, int>  m_egaIndex;

  DMET2MarkerModel() {}
  DMET2MarkerModel(RunTimeEnvironment *rte,
                   const string markerFile,
                   bool isDMET2);

};

class DMET2HaplotypeModel : public TranslationInputTsvTableModel
{
public:

  // Experiment Gene Call Call_Count Known_Count
  enum {
    EXPERIMENT,
    GENE,
    CALL,
    CALL_COUNT,
    KNOWN_COUNT,
  };

  set< string > m_experiments;
  map< string, int>  m_egcIndex;

  DMET2HaplotypeModel() {}
  DMET2HaplotypeModel(RunTimeEnvironment * rte,
                      const string haplotypeFile);

};


#endif /* TRANSLATION_REGRESSION_DMET2_MODEL_H */
