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
 * @file CNFamilialAnalysisMethodFactory.cpp
 *
 * @brief This file contains the CNFamilialAnalysisMethodFactory class members.
 */

#include "copynumber/CNFamilialAnalysisMethodFactory.h"
//
#include "copynumber/CNFamilialAnalysisMethodCNLossLOHConcordance.h"
#include "copynumber/CNFamilialAnalysisMethodCNNeutralLOHConcordance.h"
#include "copynumber/CNFamilialAnalysisMethodDenovoCopyNumber.h"
#include "copynumber/CNFamilialAnalysisMethodGenotypeConcordance.h"
#include "copynumber/CNFamilialAnalysisMethodGenotypeDiscordance.h"
#include "copynumber/CNFamilialAnalysisMethodHemizygousParentOfOrigin.h"
#include "copynumber/CNFamilialAnalysisMethodHeteroUPD.h"
#include "copynumber/CNFamilialAnalysisMethodIsoUPD.h"
#include "copynumber/CNFamilialAnalysisMethodLOD.h"
#include "copynumber/CNFamilialAnalysisMethodSegmentOverlap.h"
#include "copynumber/CNFamilialReporterFamilial.h"
//

/**
 * @brief Constructor. Registers the objects we know how to create and how
 * they describe themselves.
 */
CNFamilialAnalysisMethodFactory::CNFamilialAnalysisMethodFactory()
{
    registerMethod(CNFamilialAnalysisMethodLOD::explainSelf(), &CNFamilialAnalysisMethodLOD::newObject);
    registerMethod(CNFamilialAnalysisMethodSegmentOverlap::explainSelf(), &CNFamilialAnalysisMethodSegmentOverlap::newObject);
    registerMethod(CNFamilialAnalysisMethodGenotypeConcordance::explainSelf(), &CNFamilialAnalysisMethodGenotypeConcordance::newObject);
    registerMethod(CNFamilialAnalysisMethodGenotypeDiscordance::explainSelf(), &CNFamilialAnalysisMethodGenotypeDiscordance::newObject);
    registerMethod(CNFamilialAnalysisMethodCNNeutralLOHConcordance::explainSelf(), &CNFamilialAnalysisMethodCNNeutralLOHConcordance::newObject);
    registerMethod(CNFamilialAnalysisMethodCNLossLOHConcordance::explainSelf(), &CNFamilialAnalysisMethodCNLossLOHConcordance::newObject);
    registerMethod(CNFamilialAnalysisMethodHeteroUPD::explainSelf(), &CNFamilialAnalysisMethodHeteroUPD::newObject);
    registerMethod(CNFamilialAnalysisMethodIsoUPD::explainSelf(), &CNFamilialAnalysisMethodIsoUPD::newObject);
    registerMethod(CNFamilialAnalysisMethodDenovoCopyNumber::explainSelf(), &CNFamilialAnalysisMethodDenovoCopyNumber::newObject);
    registerMethod(CNFamilialAnalysisMethodHemizygousParentOfOrigin::explainSelf(), &CNFamilialAnalysisMethodHemizygousParentOfOrigin::newObject);

    registerMethod(CNFamilialReporterFamilial::explainSelf(), &CNFamilialReporterFamilial::newObject);
}

CNFamilialAnalysisMethodFactory::~CNFamilialAnalysisMethodFactory()
{
}

/**
 * @brief Create a pointer to a new CNFamilialAnalysisMethod object as described
 * in the string specification.
 * @param spec - Specification string
 * @return Pointer to new CNFamilialAnalysisMethod objects, must be deleted when
 * finished.
 */
CNFamilialAnalysisMethod *CNFamilialAnalysisMethodFactory::CNFamilialAnalysisMethodForString(const std::string &spec)
{
  CNFamilialAnalysisMethod *method = NULL;
  SelfCreate *create = NULL;
  create = SelfCreate::selfCreateFromString(spec, m_Docs, m_Creators, "CNFamilialAnalysisMethod", false);
  if (create == NULL) {throw(Except("SelfCreate::selfCreateFromString(...) failed using " + spec));}
  /* Check class type. */
  if(InstanceOf(create, CNFamilialAnalysisMethod))
  {
    method = static_cast<CNFamilialAnalysisMethod *>(create);
  }
  else {
    throw(Except("Class doesn't appear to be of type CNFamilialAnalysisMethod."));
  }
  return method;
}

/**
 * @brief Create a pointer to a new CNFamilialReporter object as described
 * in the string specification.
 * @param spec - Specification string
 * @return Pointer to new CNFamilialReporter objects, must be deleted when
 * finished.
 */
CNFamilialReporter* CNFamilialAnalysisMethodFactory::CNFamilialReporterForString(const std::string &spec)
{
  CNFamilialReporter* reporter = NULL;
  SelfCreate* create = NULL;
  create = SelfCreate::selfCreateFromString(spec, m_Docs, m_Creators, "CNFamilialReporter", false);
  if (create == NULL) {throw(Except("SelfCreate::selfCreateFromString(...) failed using " + spec));}
  /* Check class type. */
  if(InstanceOf(create, CNFamilialReporter))
  {
    reporter = static_cast<CNFamilialReporter*>(create);
  }
  else {
    throw(Except("Class doesn't appear to be of type CNFamilialReporter."));
  }
  return reporter;
}
