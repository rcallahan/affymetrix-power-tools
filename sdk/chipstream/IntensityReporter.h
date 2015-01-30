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
 * @file   IntensityReporter.h
 * 
 * @brief  Class for doing adapter type normalization
 */
#ifndef _INTENSITYREPORTER_H_
#define _INTENSITYREPORTER_H_

//
#include "chipstream/ChipStream.h"
//
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "util/AffxString.h"
#include "util/Guid.h"
#include "util/Util.h"
//
#include <vector>
//

/// String describing median normalization algorithm
#define INTENSITYREPORTER "intensity-reporter"

/**
 *  IntensityReporter -  Class for dumping intensities to a file.
 */
class IntensityReporter : public ChipStream 
{
public:
  static void setupSelfDoc(SelfDoc &doc) ;

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf() ;

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
  static SelfCreate *newObject(std::map<std::string,std::string> &param) ;
  
   /** 
	* Constructor
	* 
	*/
  IntensityReporter(const std::string& strFileName);

	/**
	 * @brief destructor
	 */
  virtual ~IntensityReporter();

  static std::string getDefaultOptions();

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions(); 

  /** 
   * @brief Method for being passed a new cel file worth of data. 
   * It simply counts the number of ceel files passed to the module.
   * @param data - Vector of vectors of cel file data from same sample.
   */  
/* 	void newChip(std::vector<float> &data)  */
/* 	{ */
/* 		m_iCelFileCount++;  */
/* 		chipStreamPassNewChip(data); */
/* 	} */

  /** 
   * @brief Signal that no more data is coming (i.e. newChip() will
   * not be called anymore. 
   * Define the columns for the HDF5 file.
   */
  void endChips() ;

   /** 
    * @brief transform the intensity point supplied coming from a particular
    * probe in a particular microarray.
    * Writes out the transformed intensities to the HDF5 file.
    * 
    * @param iProbeIndex - Probe index from the cel file.
    * @param chipIx - Set of chip indexes from same sample.
    * @param intensity - Original intensities. 
    * @param return - transformed intensities.
    */
  float transform(int iProbeIndex, int chipIx, float intensity);
  void transform(int chipIx, std::vector<float>& intensity);

  /** 
   * @brief Method for being passed a new cel file worth of data.
   * @param data - vector of vectors of cel file data from same sample.
   */
  void newDataSet(IntensityMart* iMart);


protected:
	AffxString m_strFileName;
	affx::File5_File m_file5;
	affx::File5_Group* m_group5;
	affx::File5_Tsv* m_tsv5;
	int m_iCelFileCount;
};

#endif /* _INTENSITYREPORTER_H_ */
