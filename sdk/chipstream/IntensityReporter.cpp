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

#include "chipstream/IntensityReporter.h"
//
#include "util/Fs.h"
//

void IntensityReporter::setupSelfDoc(SelfDoc &doc)
{
  doc.setDocName(INTENSITYREPORTER);
  doc.setDocDescription("Class for dumping intensity values to a file.");
  doc.setDocOptions(getDefaultDocOptions());
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc IntensityReporter::explainSelf()
{
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
SelfCreate * IntensityReporter::newObject(std::map<std::string, std::string> &param)
{
  SelfDoc doc = explainSelf();
  std::string strFileName;

  fillInValue(strFileName, "file", param, doc);
  IntensityReporter* p = new IntensityReporter(strFileName);
  return p;
}

/**
* Constructor
*
*/
IntensityReporter::IntensityReporter(const std::string& strFileName)
{
  m_iCelFileCount = 0;
  if (strFileName == "") {
    throw(Except("No file specified for intensity reporter."));
  }
  m_strFileName = strFileName;
  m_Type = INTENSITYREPORTER;
  setupSelfDoc(*this);

  bool bFileExists =  Fs::fileExists(strFileName);

  // @todo I think this could be "affx::FILE5_OPEN||affx::FILE5_CREATE" -jhg
  m_file5.open(strFileName, (bFileExists ? affx::FILE5_REPLACE : affx::FILE5_CREATE));
  m_group5 = m_file5.openGroup("IntensityReporter", affx::FILE5_REPLACE);
  m_tsv5 = m_group5->openTsv("Parameters", affx::FILE5_REPLACE);
  m_tsv5->defineColumn(0, 0, "Parameter", affx::FILE5_DTYPE_STRING, 1024);
  m_tsv5->set_string(0, 0, "#%guid=" + affxutil::Guid::GenerateNewGuid()); m_tsv5->writeLevel(0);
  m_tsv5->close();
  delete m_tsv5;
  m_tsv5 = m_group5->openTsv("Intensities", affx::FILE5_REPLACE);
}

/**
 * @brief destructor
 */
IntensityReporter::~IntensityReporter()
{
  m_tsv5->close();
  delete m_tsv5;
  m_group5->close();
  delete m_group5;
  m_file5.close();
}

std::string IntensityReporter::getDefaultOptions()
{
  return "";
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> IntensityReporter::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;
  SelfDoc::Opt strFileName = {"file", SelfDoc::Opt::String, "", "", "NA", "NA",
                              "The name of the file to write the intensity values to"
                             };
  opts.push_back(strFileName);
  return opts;
}

/**
 * @brief Method for being passed a new cel file worth of data.
 * It simply counts the number of ceel files passed to the module.
 * @param data - Vector of vectors of cel file data from same sample.
 */
/*  void newChip(std::vector<float> &data)  */
/*  { */
/*   m_iCelFileCount++;  */
/*   chipStreamPassNewChip(data); */
/*  } */

/**
 * @brief Signal that no more data is coming (i.e. newChip() will
 * not be called anymore.
 * Define the columns for the HDF5 file.
 */
void IntensityReporter::newDataSet(IntensityMart* iMart) {
  m_TransformedIMart = iMart;
  int dataset_count = m_TransformedIMart->getCelDataSetCount();
  std::vector<float> data;
  for (int d = 0; d < dataset_count; d++) {
    data = m_TransformedIMart->getCelData(d);
    assert(data.size() > 0);
    transform(d, data);
  }
  chipStreamPassNewChip(m_TransformedIMart);
}

void IntensityReporter::endChips()
{
  m_tsv5->defineColumn(0, 0, "probeid", affx::FILE5_DTYPE_INT);
  for (int iChipIndex = 0; (iChipIndex < m_iCelFileCount); iChipIndex++) {
    m_tsv5->defineColumn(0, (iChipIndex + 1), "Chip_" + ::getInt(iChipIndex + 1), affx::FILE5_DTYPE_FLOAT);
  }
  chipStreamEndChips();
}

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
float IntensityReporter::transform(int iProbeIndex, int chipIx, float intensity)
{
  /*    Verbose::out(1, "sample = " + ::getInt(chipIx[0]) + ", probe = " + ::getInt(iProbeId + 1)); */

  if (chipIx == 0) {
    m_tsv5->set_i(0, 0, (iProbeIndex + 1));
  }
  m_tsv5->set_f(0, (chipIx + 1), intensity);
  if ((chipIx + 1) == m_iCelFileCount) {
    m_tsv5->writeLevel(0);
  }

  return intensity;
}

void IntensityReporter::transform(int chipIx, std::vector<float>& intensity) {
  for (size_t probeIx = 0; probeIx < intensity.size(); probeIx++) {
    intensity[probeIx] = transform(probeIx, chipIx, intensity[probeIx]);
  }
}
