////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
 * @file   QuantSea.h
 * @author Pete Klosterman
 * @date   Fri Jun 23 12:48:52 2006
 *
 * @brief Class for doing probe set quantification using the SEA, or
 * Simplified Expression Analysis, method.  This is a single chip
 * method which is a subset of the PLIER (Probe Logarithmic Error
 * Intensity Estimate) method.
 */

#ifndef _QUANTSEA_H_
#define _QUANTSEA_H_

//
#include "chipstream/QuantPlier.h"

/// String description of the sea method.
#define QUANT_SEA_STR "sea"

/**
 * Class for doing probe set quantification using the SEA, or
 * Simplified Expression Analysis, method.
 */
class QuantSea : public QuantPlier {

public:

  /** 
   * What version of the algorithm are we implementing?
   * @return version string.
   */
  std::string getVersion() {
    return "2.0";
  }

  /**
   * Constructor.
   */
  QuantSea(QuantPlierParams &plierParams);

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt>::iterator i;
    std::vector<SelfDoc::Opt> options = QuantPlierBase::getDefaultDocOptions();
    // Disable options provided only in the full PLIER method.
    for (i = options.begin(); i < options.end(); i++) {
      if (i->name == "optmethod") 
      {
        i->value = "1";
        i->descript = "Optimization method to use for plier 1 for SEA (Simplified Expression Analysis), 0 for full Plier optimization. [option forced to 1 for SEA]";
      }
      else if (i->name == "PlierFitFeatureResponse") 
      {
        i->value = "false";
        i->descript = "Fit Feature Response dynamically or don't update from initial values.  [option forced to false for SEA]";
      }
    }
    return options;
  }
  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    std::vector<SelfDoc::Opt> options;
    doc.setDocName(QUANT_SEA_STR);
    doc.setDocDescription("The SEA (Simplified Expression Analysis) method provides a simple signal estimate, using the initialization algorithm from the PLIER (Probe Logarithmic Error Intensity Estimate) method and omitting the PLIER parameter fitting. SEA is useful for single chip signal estimation. The version of PLIER used by SEA differs from the previous version by the addition of a SafteyZero, NumericalTolerance, and FixPrecomputed. These options are intended to improve the stability of PLIER results when using precomputed feature reponse values. To get the older PLIER behavior set SafetyZero to 0.0, NumericalTolerance to 0.0, and FixPrecomputed to false.");
    // Set common parameters and documentation also used by QuantPlier and QuantIterPlier.
    options = getDefaultDocOptions();
    doc.setDocOptions(options);
  }

  /**
   * @brief Provide explanatory documentation.
   *
   * Parameters are for the most part the same as those for
   * QuantPlier and QuantIterPlier, except that here
   * PlierFitFeatureResponse is hard coded to 0 and
   * PlierOptOptimizationMethod is hard coded to 1.
   *
   * @return SelfDoc
   */
  static SelfDoc explainSelf() {
    SelfDoc doc;
    setupSelfDoc(doc);
    return doc;
  }

  /**
   * @brief Set up parameters, construct QuantSea object.
   *
   * PlierFitFeatureResponse is hard coded to 0 and
   * PlierOptOptimizationMethod to 1.
   *
   * @param param - Map of key/value pairs to initialize the object.
   * @return Pointer to SelfCreate object.
   */
  static SelfCreate *newObject(std::map<std::string,std::string> &param) {
    QuantPlierParams plierParams;
    SelfDoc doc = explainSelf();
    fillPlierParams(param, plierParams, doc);

    plierParams.setPlierFitFeatureResponse(false);
    plierParams.setPlierOptOptimizationMethod(1);
    QuantSea *plier = new QuantSea(plierParams);
    plier->setDocOptions(doc.getDocOptions());
    return plier;
  }

};

#endif /* _QUANTSEA_H_ */
