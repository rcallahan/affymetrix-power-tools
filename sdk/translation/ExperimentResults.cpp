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
 * @file   ExperimentResults.cpp
 * @author Mybrid Spalding
 * @date   Tue Sep 16 15:47:07 PDT 2008
 * @brief  Wrapper class for ExperimentGeneResults, all gene experiment results. There should ever only ever be one of these.
 */


#include "translation/ExperimentResults.h"

using namespace std;


/*****************************************************************************/
/**
 * ExperimentResults::Clear:
 * Synopsis:
 *
 *   Frees all the ExperimentGeneResults pointer, which in turn
 *   frees a subsequent list of CallResults pointers.
 *
 *
 */
/*****************************************************************************/
void ExperimentResults::Clear()
{
  std::map< std::string, ExperimentGeneResults *>::iterator it;


  for (it =  m_geneResults.begin(); it != m_geneResults.end(); it++) {
    delete it->second;
  }

  return;

}
// end ExperimentResults::Clear
/*****************************************************************************/
