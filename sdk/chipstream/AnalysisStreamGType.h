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
 * @file   AnalysisStreamGType.h
 * @author Chuck Sugnet
 * @date   Fri Feb 24 15:49:35 2006
 * 
 * @brief AnalysisStream class specialized for quantification methods that are
 * used for determining genotype on SNP chips.
 */

#ifndef _ANALYSISSTREAMGTYPE_H_
#define _ANALYSISSTREAMGTYPE_H_

//
#include "chipstream/AnalysisStream.h"
#include "chipstream/QuantGTypeMethod.h"

class AnalysisStreamGType : public AnalysisStream {

public:
  QuantGTypeMethod *getQuantGTypeMethod();

  /** 
   * Set the object for summarizing an individual probe set.
   * @param qMethod - Quantification object.
   */
  void setQuantMethod(QuantMethod *qMethod);

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc);

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf();
  
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
  static SelfCreate *newObject(std::map<std::string,std::string> &param);

private:
  QuantGTypeMethod *m_QGTypeMethod;
  
};

#endif /* _ANALYSISSTREAMGTYPE_H_ */
