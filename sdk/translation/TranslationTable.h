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
 * @file   TranslationTable.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  The primary business object class that manipulates the various input data models,
 *         creating other business objects appropriate for the translation algorithm and invokes translation. 
 */



#ifndef TRANSLATION_TRANSLATIONTABLE_H
#define TRANSLATION_TRANSLATIONTABLE_H

#include "translation/GeneCall.h"
//
#include <map>

class CallSet;
class CopyNumberTableModel;
class ExperimentGeneResults;
class RunTimeEnvironment;
class TranslationTableModel;

class TranslationTable
{
public:

  std::map<std::string, std::vector<GeneCall> > m_geneGroup;
  std::map<std::string, CallSet>                m_geneGroupCompleteCallSet;
  std::map<std::string, std::string>            m_geneExperimentCopyNumberCall;
  std::map<std::string, CallSet>                m_geneCopyNumberZeroCallSet;

  TranslationTable(const RunTimeEnvironment & rte, TranslationTableModel & ttm,
                   const CopyNumberTableModel & cntm);


  ~TranslationTable() {};

  bool                   hasGeneCall(const std::string & geneName);
  std::vector<GeneCall>  getGeneCallGroup(const std::string & geneName);

  // translate_experiment takes a set of experiment records filtered for just
  // the one experiment. For marker types this will mean just one record
  // or row.
  void translateExperimentGene(RunTimeEnvironment & rte,
                               GenotypeTableModel & gtm,
                               const TranslationTableModel & ttm,
                               const CopyNumberTableModel & cntm,
                               ExperimentGeneResults *egr);


private:

  GeneCall * _addGeneCallToGroup(const RunTimeEnvironment & rte,
                                 const std::string & gene,
                                 GeneCall * addGeneCall );

  void       _completeAlleleSets(GeneCall * gc);
  // Conveinance error messages

  bool       _duplicateCallSetCheck(const RunTimeEnvironment & rte);

  GeneCall * _findMarkerGeneCall(std::string gene, std::string probeSet);

  bool       _ignoreTTMRow(const RunTimeEnvironment & rte,
                           TranslationTableModel & ttm,
                           const int row);

  std::vector<CallSet>
  _invalidateExperimentGeneProbeSetBaseValues(
    const RunTimeEnvironment & rte, GenotypeTableModel &gtm,
    const TranslationTableModel & ttm);


  bool       _reconcileCopyNumberCallSet(const RunTimeEnvironment & rte,
                                         TranslationTableModel & ttm);

  bool       _validateTTMRowIsHeader(const RunTimeEnvironment & rte,
                                     TranslationTableModel & ttm,
                                     const int row);
};

#endif /* TRANSLATION_TRANSLATIONTABLE_H */

