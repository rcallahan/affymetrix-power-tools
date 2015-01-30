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
 * @file   TranslationCommonControl.h
 * @author Mybrid Spalding
 * @date   Tue Jul  8 08:55:25 PDT 2008
 * @brief  Controller code common to both command line and console controllers. 
 */


#ifndef _TRANSLATION_COMMON_CONTROL_H_
#define _TRANSLATION_COMMON_CONTROL_H_

#include <cstring>
#include <string>
//

extern const std::string HAPLOTYPE_REPORT_OPTION_FIRST;

extern const std::string HAPLOTYPE_REPORT_OPTION_DUPLICATE;


class TranslationCommonControl
{

public:

  static void defineCommonOptions(class TranslationEngine *engine);

  static bool initializeInputFileModels(class RunTimeEnvironment &rte,
                                        class MarkerListModel       **returnMlm,
                                        class CopyNumberTableModel  **returnCnrm,
                                        class GenotypeTableModel    **returnGtm,
                                        class TranslationTableModel **returnTtm,
                                        class SampleInfoTableModel  **returnSitm,
                                        class GenotypeOverrideTableModel
                                        **returnGotm);

  static void setCommonADTOptions(class ADTOptions & adtOpts, class TranslationEngine *engine);

  static void validateCommonOptions(class RunTimeEnvironment & rte , class TranslationEngine *engine);


};


#endif /* _TRANSLATION_COMMON_CONTROL_H_ */
