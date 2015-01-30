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
 * @file   AnalysisStreamGType.cpp
 * @author Chuck Sugnet
 * @date   Fri Feb 24 15:49:35 2006
 * 
 * @brief AnalysisStream class specialized for quantification methods that are
 * used for determining genotype on SNP chips.
 */

#include "chipstream/AnalysisStreamGType.h"

inline QuantGTypeMethod *AnalysisStreamGType::getQuantGTypeMethod() {
  return m_QGTypeMethod; 
}

/** 
 * Set the object for summarizing an individual probe set.
 * @param qMethod - Quantification object.
 */
void AnalysisStreamGType::setQuantMethod(QuantMethod *qMethod) {
  QuantGTypeMethod *gMethod = dynamic_cast<QuantGTypeMethod *>(qMethod);
  if(gMethod == NULL)
    Err::errAbort("AnalysisStreamGType::setQuantMethod() - Can only set QuantGTypeMethods in AnalysisStreamGType.");
  m_QMethod = gMethod;
  m_QGTypeMethod = gMethod;
}

/** 
 * Fill in the information for Self documentation.
 * @param doc - Self documenter to be filled in.
 */
void AnalysisStreamGType::setupSelfDoc(SelfDoc &doc) {
  doc.setDocName("geno");
  doc.setDocDescription("Does genotype calls on probesets.");
}

/** 
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc AnalysisStreamGType::explainSelf() { 
  SelfDoc doc;
  setupSelfDoc(doc);
  return doc;
}
  
/** 
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 * 
 * @param param - Map of key/value pairs to initialize the object.
 * 
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate *AnalysisStreamGType::newObject(std::map<std::string,std::string> &param) {
  SelfDoc doc = explainSelf();
  AnalysisStreamGType *stream = new AnalysisStreamGType();
  return stream;
}


