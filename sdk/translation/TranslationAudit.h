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
 * @file   TranslationAudit.h
 * @author Mybrid Spalding
 * @date   Thu Sep 25 13:56:24 PDT 2008
 * @brief  Business object that contains the annotation file and translation file audit requirements.
 */

#ifndef TRANSLATION_AUDIT_H
#define TRANSLATION_AUDIT_H
//
#include <cstring>
#include <sstream>
#include <string>
#include <vector>
//

class TranslationAudit
{

public:

  TranslationAudit() {};

  TranslationAudit(const class RunTimeEnvironment &rte);

  ~TranslationAudit() {};

  int audit(class RunTimeEnvironment & rte);

  int auditAnnotation(class RunTimeEnvironment & rte);

private:

  std::stringstream          m_invalidGeneProbeSetsSStr;
  std::stringstream          m_invalidSwitchDesignStrandSStr;
  std::stringstream          m_msgSStr;
  std::vector< std::string > m_missingProbeSets;
  int                        m_numInvalidGeneProbeSets;
  int                        m_numInvalidSwitchDesignStrands;

  int  _auditAnnotationReport(class RunTimeEnvironment & rte);
  int  _auditAnnotationReferenceVariant(class RunTimeEnvironment & rte);

};

class TranslationAuditType
{
public:

  std::string m_commandLine;
  int (TranslationAudit::*m_auditCallback)(class RunTimeEnvironment & rte);

};

#endif /* TRANSLATION_AUDIT_H */

