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
 * @file   AdapterTypeNormTran.h
 *
 * @brief  Class for doing adapter type normalization
 */
#ifndef _ADAPTERTYPENORMTRAN_H_
#define _ADAPTERTYPENORMTRAN_H_

#define ADAPTER_TYPE_COUNT 10
//
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStream.h"
//
#include "util/AffxByteArray.h"
#include "util/AffxConv.h"

#include "util/AffxMultiDimensionalArray.h"
#include "util/Convert.h"
#include "util/Util.h"
//
#include <vector>
//

/// String describing median normalization algorithm
#define ADAPTERTYPENORMSTR "adapter-type-norm"

class CNProbeSetArray;
class CNProbeArray;
/**
 *  AdapterTypeNormTran -  Class for doing adapter type normalization
 */
class AdapterTypeNormTran : public ChipStream
{
public:

  /**
  * Constructor
  *
  */
  AdapterTypeNormTran(double dStyOnlySnpRatio, double dNspOnlySnpRatio, double dNspOnlyCnRatio,  double dNspBothSnpRatio, double dNspBothCnRatio);

  ~AdapterTypeNormTran();

  static std::string getDefaultOptions();

  /**
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions();

  void setup(const std::string& strAnnotationFileName, ChipLayout* pobjChipLayout, unsigned int uiCelFileCount) ;

  /**
   * @brief transform the intensity point supplied coming from a particular
   * probe in a particular microarray.
   *
   * @param iProbeId - Probe index from the cel file.
   * @param chipIx - chip index from sample.
   * @param intensity - Original intensity.
   */
  float transform(int iProbeId, int chipIx, float intensity);

  /**
   * @brief Process the annotation files to extract adapter types per probeset.
   */
  void loadAnnotation();

  /**
   * @brief Method for being passed a new cel file worth of data.
   * @param data - Vector of vectors of cel file data from same sample.
   */  
	void newDataSet(IntensityMart* iMart);

  /**
   * @brief Method for being passed a new cel file worth of data. Calculates
   * and stores the median (or average) of data supplied.
   * @param data - Vector of cel file data from same sample.
   */
  void initializeData(const std::vector<float> &data);

  void newChip(std::vector<float> &data);


  /**
   * @brief Signal that no more data is coming (i.e. newChip() will
   * not be called anymore. Currently the class builds up all the
   * summaries and figures out the normalization factores when this
   * function is called.  This makes it difficult to pass through data
   * to downstream.
   */
  void endDataSet();

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
  static SelfCreate *newObject(std::map<std::string, std::string> &param);

protected:
// void parseAnnotation(AffxString& strProbeSetName, AffxString& str, char& cNspAdapterType, char& cStyAdapterType);
  char parseAnnotation(AffxString& str, const AffxString& strEnzyme, AffxArray<AffxString>& vAdapters);
  double getNspOnlyMedian(const std::vector<float>& data, char cProbeSetType, char cAdapterType, AffxMultiDimensionalArray<float>& vIntensities);
  double getStyOnlyMedian(const std::vector<float>& data, char cProbeSetType, char cAdapterType, AffxMultiDimensionalArray<float>& vIntensities);
  double getNspBothMedian(const std::vector<float>& data, char cProbeSetType, char cAdapterType, AffxMultiDimensionalArray<float>& vIntensities);
  double getStyBothMedian(const std::vector<float>& data, char cProbeSetType, char cAdapterType, AffxMultiDimensionalArray<float>& vIntensities);

  bool isPm(int iIndex);
  bool isSnp(int iIndex);
  bool isCn(int iIndex);
  char getNspAdapterType(int iIndex);
  char getStyAdapterType(int iIndex);

protected:
  /// Annotation DB file name
  std::string m_strAnnotationFileName;
  /// SNP Annotation file name
  std::string m_strSnpAnnotationFileName;
  /// CN Annotation file name
  std::string m_strCnAnnotationFileName;
  /// Pointer to the ChipLayout object
  ChipLayout* m_pobjChipLayout;
  /// Grand Median for each chip seen.
  AffxMultiDimensionalArray<double> m_vGrandMedians;

  AffxMultiDimensionalArray<char> m_mxAdapterProbeTypes;
  AffxMultiDimensionalArray<char> m_mxAdapterPMTypes;
  AffxMultiDimensionalArray<char> m_mxAdapterNspTypes;
  AffxMultiDimensionalArray<char> m_mxAdapterStyTypes;

  AffxMultiDimensionalArray<double> m_mxNspOnlySnpMedians;
  AffxMultiDimensionalArray<double> m_mxNspOnlyCnMedians;
  AffxMultiDimensionalArray<double> m_mxStyOnlySnpMedians;
  AffxMultiDimensionalArray<double> m_mxStyOnlyCnMedians;
  AffxMultiDimensionalArray<double> m_mxNspBothSnpMedians;
  AffxMultiDimensionalArray<double> m_mxNspBothCnMedians;
  AffxMultiDimensionalArray<double> m_mxStyBothSnpMedians;
  AffxMultiDimensionalArray<double> m_mxStyBothCnMedians;

  double m_dStyOnlySnpRatio;
  double m_dNspOnlySnpRatio;
  double m_dNspOnlyCnRatio;
  double m_dNspBothSnpRatio;
  double m_dNspBothCnRatio;

  bool m_bAnnotationFilesProcessed;
  unsigned int m_uiCelFileCount;
  unsigned int m_uiCelFileIndex;

public:
  void setCelFileIndex(unsigned int ui);
};

#endif /* _ADAPTERTYPENORMTRAN_H_ */
