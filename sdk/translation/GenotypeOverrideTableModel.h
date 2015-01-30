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
 * @file   GenotypeOverrideTableModel.h
 * @author Mybrid Spalding
 * @date   Tue Jul 22 08:33:25 PDT 2008
 * @brief  Class to supplant Genotype data (CHP) with user supplied data.
 */

#ifndef TRANSLATION_GENOTYPE_OVERRIDE_TABLE_MODEL_H
#define TRANSLATION_GENOTYPE_OVERRIDE_TABLE_MODEL_H

#include "translation/TranslationInputTsvTableModel.h"

extern const TittmColumnDefinition GOTM_COLUMN_DEFINTIONS[];

class GenotypeOverrideTableModel : public TranslationInputTsvTableModel
{

public:

  // DATA
  enum  {
    CHP_FILE_INDEX,
    GENE_INDEX,
    MARKER_NAME_INDEX,
    PROBE_SET_INDEX,
    BASECALL_INDEX,
    OVERRIDE_COMMENT_INDEX,
    REFERENCE_INDEX,
    VARIANT_INDEX,
    END_GOTM_FIELD_INDEX,
  };

  std::map< std::string, int >         m_chpFileProbeSetToRowIndex;
  std::map< std::string, std::string > m_chpFileProbeSetOriginalBasecall;
  std::map< std::string , std::vector< std::string > >
  m_chpFileProbeSetInvalidBasecallComment;

  GenotypeOverrideTableModel(const class RunTimeEnvironment & rte,
                             const std::string & genoOverrideFileName,
                             class TranslationTableModel & ttm);

  GenotypeOverrideTableModel() {};
  ~GenotypeOverrideTableModel() {};

  void                       describeVerbose(ADT_VERBOSE_ENUM level = ADT_VERBOSE_INPUT_FILES);

  std::vector< std::string > getOverrideAlleles(const std::string & chpFileName,
      const std::string & probeSetId, int row = -1);

  std::vector< std::string > getOverrideInfo(const std::string & chpFileName,
                                             const std::string & probeSetId,
                                             bool includeInvalidOriginalBaseCall = false);
  std::string                getZeroCopyNumberOverrideInfo( const std::string & chpFileName, const std::string & probeSetId );

private:

};

#endif /* TRANSLATION_GENOTYPE_OVERRIDE_TABLE_MODEL_H */
