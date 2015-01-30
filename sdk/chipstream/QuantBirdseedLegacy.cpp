////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

/*
 * FILE QuantBirdseedLegacy.cpp
 */

#include "chipstream/QuantBirdseedLegacy.h"

QuantBirdseedLegacy::~QuantBirdseedLegacy() {};

/** Constructor, currently creates prior estimates from pre-calculated data from R. */
QuantBirdseedLegacy::QuantBirdseedLegacy(double confidenceThreshold, double correctionFactor) :
    QuantBirdseedv1(confidenceThreshold, correctionFactor)
{
  setupSelfDoc(*this);
}

/** Fill in our self documentation from current state. */
void QuantBirdseedLegacy::setupSelfDoc(SelfDoc &doc)
{
  doc.setDocName(QUANTBIRDSEEDLEGACY);
  doc.setDocDescription("Do genotyping calls using the Birdseed v1 algorithm. Legacy alias to birdseed-v1 method. You should use birdseed-v1 rather than birdseed.");
  doc.setDocOptions(getDefaultDocOptions());
}

/**
 * @brief What is the name of the quantification method?
 * @return name of adjuster.
 */
std::string QuantBirdseedLegacy::getType()
{
  return std::string(QUANTBIRDSEEDLEGACY);
}

SelfCreate * QuantBirdseedLegacy::newObject(std::map<std::string, std::string> &param)
{
  SelfDoc doc = explainSelf();
  double confThreshold;
  fillInValue(confThreshold, "conf-threshold", param, doc);
  double correctionFactor;
  fillInValue(correctionFactor, "correction-factor", param, doc);
  birdseed::v1::setSelfDocOptions(&doc, param);
  QuantBirdseedLegacy *obj = new QuantBirdseedLegacy(confThreshold, correctionFactor);

  /// @todo hack to get the self doc correctly set with user supplied values
  /// this should really go in the constructor
  birdseed::v1::setSelfDocOptions(obj, param);
  return obj;
}

SelfDoc QuantBirdseedLegacy::explainSelf()
{
  SelfDoc doc;
  setupSelfDoc(doc);
  return doc;
}


/******************************************************************/
/************************[END OF QuantBirdseed.cpp]****************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */

