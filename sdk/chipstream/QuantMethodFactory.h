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
 * @file   QuantMethodFactory.h
 * @author Chuck Sugnet
 * @date   Tue Oct 25 11:56:37 2005
 * 
 * @brief Factory class for making chip streams based on a string
 * representation.
 */

#ifndef _QUANTMETHODFACTORY_H_
#define _QUANTMETHODFACTORY_H_

//
#include "chipstream/QuantDabg.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantIterPlier.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantPlier.h"
#include "chipstream/QuantRma.h"
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
//
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

/**
 * Factory class for making chip streams based on a string representation. 
 */
class QuantMethodFactory {

public:
  
  enum QuantType { Expression, GenoType };

  /** 
   * @brief Constructor. Registers the objects we know how to create and how
   * they describe themselves.
   */    
  QuantMethodFactory();
    
  /** 
   * @brief Constructor. Registers the objects we know how to create and how
   * they describe themselves.
   * @param Type of Quantification methods to make available.
   */    
  QuantMethodFactory(enum QuantType type);

  virtual ~QuantMethodFactory();
  
  void registerMethod(SelfDoc &doc, SelfCreate::selfCreator maker) {
    m_Docs.push_back(doc);
    m_Creators.push_back(maker);
  }

  /** 
   * @brief Create a pointer to a new QuantMethod object as described
   * in the string specification.
   * @param spec - Specification string i.e. quant-norm.sketch=1000000
   * 
   * @return Pointer to new QuantMethod objects, must be deleted when
   * finished.
   */
  QuantMethod *qMethodForString(const std::string &spec);

  /** 
   * @brief Create a pointer to a new QuantMethod object as described
   * in the string specification.
   * @param spec - Specification string i.e. quant-norm.sketch=1000000
   * @param board - Blackboard with various state
   * @return Pointer to new QuantMethod objects, must be deleted when
   * finished.
   */
  QuantMethod *quantMethodForString(const string &spec, PsBoard &board, QuantType type);

  /** 
   * @brief Create a pointer to a new QuantMethod object as described
   * in the string specification.
   * @param spec - Specification string i.e. quant-norm.sketch=1000000
   * @param board - Blackboard with various state
   * @return Pointer to new QuantMethod objects, must be deleted when
   * finished.
   */
  QuantMethod *quantMethodForString(const string &spec, PsBoard &board);


  /** 
   * @brief Factory for creating QuantMethod from string representation.
   * 
   * @param spec - Name of QuantMethod.
   * @param layout - Probe and probe set info.
   * 
   * @return Quantification method requested.
   */
  QuantMethod *quantMethodForString(string &spec, ChipLayout &layout, QuantType type);

  /** 
   * @brief Factory for creating QuantMethod from string representation.
   * 
   * @param spec - Name of QuantMethod.
   * @param layout - Probe and probe set info.
   * 
   * @return Quantification method requested.
   */
  QuantExprMethod *quantExprMethodForString(string &spec, ChipLayout &layout, QuantType type);

  /** 
   * @brief Factory for creating genotyping QuantMethod from string
   * representation.
   * 
   * @param spec - Name of QuantMethod.
   * @param layout - Probe and probe set info.
   * 
   * @return Quantification method requested.
   */
  QuantGTypeMethod *quantGTypeMethodForString(string &spec, ChipLayout &layout, QuantType type);

  /** 
   * Set the probes to be used as gc background for quantification methods.
   * Note that memory for the probes is owned elsewhere.
   *
   * @param gcControlProbes - Probes that should be used for
   * background distribution.
   */
  void setGcControlProbes(const std::vector<Probe *> &gcControlProbes) {
    m_GcControlProbes = gcControlProbes;
  }

  /** 
   * Open a file and read the precomputed probe effects to use in plier
   * from it.
   * 
   * @param fileName - Path to file containing probe effects.
   * @param effects - Map of probe ids to probe effects to be filled in.
   */

  static int openPrecompFeatureEffects(       const std::string& fileName, 
                                              double* effects, 
                                              const ChipLayout &layout);

  int openPrecompFeatureEffects_A5(affx::File5_Group* group5,
                                           const std::string& name, 
                                           double* effects,
                                           const ChipLayout &layout);

  /** 
   * Read the feature responses from a file. File should have columns
   * for 'probeset_id', 'probe_id' and 'feature_response'.
   * 
   * @param fileName - text file with feature responses.
   */
  const void readPrecompFeatureEffectsFromFile(	int iProbeCount, 
                                                const std::string& fileName, 
                                                const ChipLayout &layout);

  void readPrecompFeatureEffectsFrom_A5(      int iProbeCount, 
                                              affx::File5_Group* group5,
                                              const std::string& name, 
                                              const ChipLayout &layout); 



  /**
   * Get the vector of documentation objects that this factory is aware of.
   * @return - vector of documentation objects.
   */
  std::vector<SelfDoc> getDocs() {
    return m_Docs;
  }

  void static setupQuantMethodFactory(
        QuantMethodFactory &qFactory,
        affx::File5_File * inputFile,
        affx::File5_Group * inputGroup,
        affx::File5_File * outputFile,
        affx::File5_Group * outputGroup,
        uint32_t probeCount,
        const std::string &outDir,
        const std::string &useFeatEff,
        bool a5FeatureEffectsInputGlobal,
        const std::string &a5FeatureEffectsInputFile,
        const std::string &a5FeatureEffectsInputName,
        const std::string &a5FeatureEffectsInputGroup,
        const std::string &a5InputGroup,
        const std::string &a5Group,
        const std::string &setAnalysisName,
        const std::string &quantMethodSpec,
        const ChipLayout &layout 
     );


  /** 
   * Can this factory be used to build from this specification?
   * @param spec - Specification in form "chipstream.key=value.key=value"
   * @return true if this factory knows how to build from this specification, false otherwise
   */
  bool canBuild(const std::string &spec) {
    return SelfCreate::canMake(spec, m_Docs);
  }
  

private:

  void demangleName(  const std::string inputString,
                      std::string &demangledProbeSetName,
                      int &allele_id,
                      int &context_id);

protected:

  // Feature effects of all probes. All data in feature effects file.
  double* m_PrecompFeatureEffects; 
  int m_iFeatureEffectSize;
  /// Background probes to use. Note that the memory is owned elsewhere.
  std::vector<Probe *> m_GcControlProbes;    
  /// Self documentation
  std::vector<SelfDoc> m_Docs;
  /// Self creation
  std::vector<SelfCreate::selfCreator> m_Creators;
};

#endif /* _QUANTMETHODFACTORY_H_ */
  
