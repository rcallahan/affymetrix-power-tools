////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
//
#include "chipstream/DataStore.h"
//
#include "chipstream/PsBoard.h"
#include "chipstream/AnalysisStreamFactory.h"
#include "util/Fs.h"

using namespace affx;

#define DATA_ACCESSOR(TYPE,SUFFIX,DESCRIPTION,TYPE_FILE5,PROBEPREFIX)     \
bool DataStore::isSet##SUFFIX()                                           \
{                                                                         \
    int colIx = m_ColIndexes##PROBEPREFIX[SUFFIX];                        \
    return colIx >= 0;                                                    \
}                                                                         \
void DataStore::get##SUFFIX(std::vector<TYPE> &vec)                       \
{                                                                         \
    checkSetup(DESCRIPTION);                                              \
    int colIx = m_ColIndexes##PROBEPREFIX[SUFFIX];                        \
    APT_ERR_ASSERT(colIx >= 0, ToStr(DESCRIPTION) + " not set.");         \
    getTsv##PROBEPREFIX();                                                \
    File5_TsvColumn* column = m_Tsv##PROBEPREFIX->getColumnPtr(0,colIx);  \
    vector<TYPE> toFill(m_Order##PROBEPREFIX.size());                     \
    fill(toFill.begin(), toFill.end(), -1);                               \
    int read = column->read_array(0, toFill.size(), &toFill[0]);          \
    vec.resize(toFill.size());                                            \
    for(int i = 0; i < toFill.size(); i++) {                              \
      vec[i] = toFill[m_Map##PROBEPREFIX[i]];                             \
    }                                                                     \
    APT_ERR_ASSERT(read == m_Num##PROBEPREFIX, "bad read.");              \
}                                                                         \
void DataStore::set##SUFFIX(const std::vector<TYPE> &vec) {               \
    checkSetup(DESCRIPTION);                                              \
    int colIx = m_ColIndexes##PROBEPREFIX[SUFFIX];                        \
    if(colIx == -1) {                                                     \
      colIx = m_ColCount##PROBEPREFIX++;                                  \
      m_ColIndexes##PROBEPREFIX[SUFFIX] = colIx;                          \
    }                                                                     \
    getTsv##PROBEPREFIX();                                                \
    m_Tsv##PROBEPREFIX->defineColumn(0,colIx,DESCRIPTION,TYPE_FILE5);     \
    File5_TsvColumn* column = m_Tsv##PROBEPREFIX->getColumnPtr(0,colIx);  \
    vector<TYPE> data(vec.size());                                        \
    fill(data.begin(), data.end(), -1);                                   \
    column->resize(data.size());                                          \
    for(int i = 0; i < vec.size(); i++) {                                 \
      data[i] = vec[m_Order##PROBEPREFIX[i]];                             \
    }                                                                     \
    int written = column->write_array(0, data.size(), &data[0]);          \
    APT_ERR_ASSERT(written == data.size(), "Invalid write");              \
  }                                                                       \
TYPE DataStore::get##SUFFIX(int idx)                                      \
{                                                                         \
    getTsv##PROBEPREFIX();                                                \
    m_Tsv##PROBEPREFIX->gotoLine(m_Map##PROBEPREFIX[idx]);                \
    TYPE value = -1;                                                      \
    int colIx = m_ColIndexes##PROBEPREFIX[SUFFIX];                        \
    APT_ERR_ASSERT(colIx >= 0, ToStr(DESCRIPTION) +" not set.");          \
    if(!FILE5_OK == m_Tsv##PROBEPREFIX->get(0, colIx, &value)) {          \
      Err::errAbort("Error reading " + ToStr(idx)                         \
                    + " " + ToStr(DESCRIPTION));                          \
    }                                                                     \
   return value;                                                          \
}

DATA_ACCESSOR(char, ProbeGc,"ProbeGc", FILE5_DTYPE_CHAR, Probe);
DATA_ACCESSOR(char, ProbeGcBgrd,"ProbeGcBgrd", FILE5_DTYPE_CHAR, Probe);
DATA_ACCESSOR(char, ProbePm,"ProbePm", FILE5_DTYPE_CHAR, Probe);
DATA_ACCESSOR(char, ProbeMm,"ProbeMm", FILE5_DTYPE_CHAR, Probe);
DATA_ACCESSOR(int, ProbePmMm,"ProbePmMm", FILE5_DTYPE_INT, Probe);
DATA_ACCESSOR(char, ProbeInProbeset,"ProbeInProbeset", FILE5_DTYPE_CHAR, Probe);
DATA_ACCESSOR(int, ProbePmAlleleMatch,"ProbePmAlleleMatch", FILE5_DTYPE_INT, Probe);

void DataStore::closeResources() {

  if(m_OrderProbeStorage != NULL) {
    m_OrderProbeStorage->close();
    Freez(m_OrderProbeStorage);
  }

  if(m_Chips != NULL) {
    m_Chips->close();
    Freez(m_Chips);
  }

  if(m_TsvProbe != NULL) {
    m_TsvProbe->close();
    Freez(m_TsvProbe);
  }

  if(m_GroupProbe != NULL) {
    m_GroupProbe->close();
    Freez(m_GroupProbe);
  }
 
  if (m_ProbeSetOrderStorage != NULL) {
      m_ProbeSetOrderStorage->close();
      Freez(m_ProbeSetOrderStorage);
  }
  
  if (m_TsvProbeSet != NULL) {
      m_TsvProbeSet->close();
      Freez(m_TsvProbeSet);
  }
  
  if (m_GroupProbeSet != NULL) {
      m_GroupProbeSet->close();
      Freez(m_GroupProbeSet);
  }
   
  if(m_IntensityGroup != NULL) {
    m_IntensityGroup->close();
    Freez(m_IntensityGroup);
  }
  
  m_File5.close();
}

DataStore::~DataStore() {
  try {
    closeResources();
    Fs::rm(m_FileStorage);
  }
  catch (...) {
    Verbose::out(0, "****Error**** - Exception thrown in DataStore Destructor!");
  }
}

void DataStore::setValues(ChipLayout &layout, PsBoard &board) { 
  // Watchout - these vectors get reused for space efficency
  //  std::vector<int> probeReuse(m_NumProbe, 0);
  std::vector<bool> probeBool(m_NumProbe, false);    

  // ProbeGc
  
  vector<char> gc = layout.getGcProbeVec();
  if(gc.size() > 0) {
    setProbeGc(gc);
  }

  // ProbeGcBgrd
  vector<char> probeReuse(m_OrderProbe.size());
  if (board.getOptions()->isOptDefined("bgp-file") && board.getOptions()->getOpt("bgp-file") != "") {
    vector<Probe *> controlProbes;
    AnalysisStreamFactory::probeListFromBgpFile(board.getOptions()->getOpt("bgp-file"), controlProbes);
    fill(probeReuse.begin(), probeReuse.end(), 0);
    for(int i = 0; i < controlProbes.size(); i++) {
      probeReuse[controlProbes[i]->id] = 1;
    }
    setProbeGcBgrd(probeReuse);
  }

  // ProbePm
  probeBool = layout.getPmProbeMask();

  fill(probeReuse.begin(), probeReuse.end(), 0);
  for (int i = 0; i < probeBool.size(); i++) {
    if (probeBool[i]) {
      probeReuse[i] = 1;
    }
  }
  setProbePm(probeReuse);

  // ProbeMM
  probeBool = layout.getMmProbeMask();
  fill(probeReuse.begin(), probeReuse.end(), 0);
  for(int i = 0; i < probeBool.size(); i++) {
    if (probeBool[i]) {
      probeReuse[i] = 1;
    }
  }
  setProbeMm(probeReuse);

  // ProbePmMm
  vector<int> pmMm = layout.getPmMmVec();
  if (!pmMm.empty() ) {
    setProbePmMm(pmMm);
    pmMm.resize(0);
    vector<int>(pmMm).swap(pmMm);
  }

  // ProbeInProbeResueet
  probeBool = layout.getProbesetProbes();
  fill(probeReuse.begin(), probeReuse.end(), 0);
  for(int i = 0; i < probeBool.size(); i++) {
    if(probeBool[i]) {
      probeReuse[i] = 1;
    }
  }
  setProbeInProbeset(probeReuse);

  // Reusing the pmMm vector. bad style but keeps memory profile low
  pmMm = layout.getPmAlleleMatchVec();
  if (!pmMm.empty()) {
    setProbePmAlleleMatch(pmMm);
  }
}

void DataStore::setProbeOrder(const std::vector<int> &order) {
  if (!m_OrderProbe.empty()) {
    APT_ERR_ABORT("Can't reset order once it has been set.");
  }
  if(!m_OrderProbe.empty() && order.size() != m_OrderProbe.size()) {
    Err::errAbort("New order has: " + ToStr(order.size()) + " entries but orignial order has: " + ToStr(m_OrderProbe.size()));
  } 

  m_OrderProbe = order;
  affx::File5_Group *group = getIntensityGroup();
  m_OrderProbeStorage = group->openVector("/intensities/order", FILE5_DTYPE_INT, FILE5_OPEN_CREATE);
  m_OrderProbe = order;
  m_NumProbe = order.size();
  m_MapProbe.resize(m_OrderProbe.size());
  fill(m_MapProbe.begin(), m_MapProbe.end(), -1);
  int indexCount = 0;
  vector<bool> seen(m_OrderProbe.size(), false);
  for(int i = 0; i < m_OrderProbe.size(); i++) {
    m_MapProbe[m_OrderProbe[i]] = indexCount++;
    if(seen[m_OrderProbe[i]]) {
      Err::errAbort("Can't have dupes in the order!");
    }
    else {
      seen[m_OrderProbe[i]] = true;
    }
  }
  m_OrderProbeStorage->resize(m_OrderProbe.size());
  int written = m_OrderProbeStorage->write_array(0,m_OrderProbe.size(), &m_OrderProbe[0]);
  assert(written == m_OrderProbe.size());
  m_OrderProbeStorage->close();
  Freez(m_OrderProbeStorage);;
}

void DataStore::writeColumn(int chipIx, const std::string &colName, const std::vector<float> &data) {
  assert(data.size() == m_NumProbe);
  // Create a new column for our chip data in the channel matrix
  m_Chips->defineColumn(0,chipIx,colName,FILE5_DTYPE_FLOAT);
  affx::File5_TsvColumn* column = m_Chips->getColumnPtr(0,chipIx);
  //    column->setBufferSize(256000);
  //    m_CelNames.push_back(colName);
  // Write out the column data to HDF5
  column->resize(data.size());
  vector<float> ordered(data.size());
  fill(ordered.begin(), ordered.end(), -1.0f);
  for(int i = 0; i < data.size(); i++) {
    ordered[i] = data[m_OrderProbe[i]];
  }
  int written = column->write_array(0, ordered.size(), &ordered[0]);
  assert(written == data.size() );
}

void DataStore::newChip(affymetrix_fusion_io::FusionCELData *cel) {
    
  vector<float> data;
  string celName = cel->GetFileName();
  if(m_CelNames.empty()) {
    initFromCelFile(*cel);
  }
  std::vector<std::wstring> dataChannels = cel->GetChannels();
  if(dataChannels.size() > 0 && dataChannels.size() != m_NumChannels) {
    Err::errAbort("Expecting to get: " + ToStr(m_NumChannels) + " in cel file: " + 
                  m_CelNames[m_CelNames.size() - 1] + " but got: " + ToStr(dataChannels.size()) + " instead");
  }
  m_CelNames.push_back(celName);
  for(int channelIx = 0; channelIx < dataChannels.size() || channelIx == 0; channelIx++) {
    cel->SetActiveDataGroup(dataChannels[channelIx]);
      
    assert(m_NumProbe == cel->GetRows() * cel->GetCols());
      
    // Load the data up into our vector
    data.resize(m_NumProbe);
    fill(data.begin(), data.end(), -1.0f);
    cel->GetIntensities(0,data);
    writeColumn(m_NumChips++, celName, data);

  }
}

void DataStore::readCelFile(const std::string &celName) {
  FusionCELData cel;
  try {
    std::string tmp_unc_name=Fs::convertToUncPath(celName);
    cel.SetFileName(tmp_unc_name.c_str());
    if(!cel.Read()) {
      Err::errAbort("\nCan't read cel file: " + cel.GetFileName() + 
                    "\n>>> Error reported: " + StringUtils::ConvertWCSToMBS(cel.GetError()));
    }
    newChip(&cel);
    cel.Close();
  }
  catch(const Except &e) {
    Err::errAbort(ToStr("\n") + e.what());
  }
  catch(const std::bad_alloc &e) {
    Err::errAbort("\nRan out of memory when reading cel file: " + cel.GetFileName() + "\n");
  }
  catch(const std::exception &e) {
    Err::errAbort("\nException caught. Message is: " + ToStr(e.what()));
  }
  catch (...) {
    Err::errAbort("\nUnknown problem when reading cel file: " + cel.GetFileName() + "\n");
  }
}

/** 
 * @brief Return all of the intensities for given chipIx
 * @param chipIx - Chip Index number.
 * @return double - vector of all intensities in the DataStore for given chipIx.
 */
void DataStore::fillInCelData(chipid_t chipIx, std::vector<float> &data) const { 
  affx::File5_TsvColumn* column = m_Chips->getColumnPtr(0,chipIx);
  data.resize(column->reserved());
  vector<float> toFill(m_NumProbe, -1.0f);
  int read = column->read_array(0, toFill.size(), &toFill[0]);
  for(int i = 0; i < toFill.size(); i++) {
    data[i] = toFill[m_MapProbe[i]];
  }
  APT_ERR_ASSERT(read == m_NumProbe, "bad read.");
}

void DataStore::setGcControlProbes(const std::vector<Probe *> &vec) {
  checkSetup("setGcControlProbes");
  vector<char> isBgProbe(m_NumProbe);
  fill(isBgProbe.begin(), isBgProbe.end(), 0);
  for(int i = 0; i < vec.size(); i++) {
    isBgProbe[m_MapProbe[vec[i]->id]] = 1;
  }
  setProbeGcBgrd(isBgProbe);
}

void DataStore::getGcControlProbes(std::vector<int> &vec) {
  vector<char> isBgProbe(m_NumProbe);
  getProbeGcBgrd(isBgProbe);
  for(int i = 0; i < isBgProbe.size(); i++) {
    if(isBgProbe[i] == 1) {
      vec.push_back(i);
    }
  }
}

DataStore::DataStore() {
  initVariables();
}

DataStore::DataStore(const std::string &fileLocation) {
  m_FileStorage = fileLocation;
  int ret = m_File5.open(fileLocation, affx::FILE5_OPEN_CREATE | affx::FILE5_REPLACE);
  if(ret != FILE5_OK) {
    Err::errAbort("Couldn't open file: " + fileLocation);
  }
  initVariables();
}

void DataStore::initVariables() {
  // @todo - should there be a clear first?
  m_NumProbe = 0;
  m_NumChannels = 0;
  m_NumChips = 0;
  m_ProbeSetOrderStorage = NULL;
  m_IntensityGroup = NULL;
  m_OrderProbeStorage = NULL;
  m_ProbeSetOrderStorage = NULL;
  m_Chips = NULL;
  m_GroupProbe = NULL;
  m_GroupProbeSet = NULL;
  m_TsvProbe = NULL;
  m_TsvProbeSet = NULL;
  m_ColIndexesProbe.resize(ProbeLastEnum);
  fill(m_ColIndexesProbe.begin(), m_ColIndexesProbe.end(), -1);
  m_ColCountProbe = 0;
  m_ColCountProbeSet = 0;
}

void DataStore::initFromFile(const string &fileLocation) {
  m_FileStorage = fileLocation;
  int ret = m_File5.open(fileLocation, affx::FILE5_OPEN);
  if(ret != FILE5_OK) {
    Err::errAbort("Couldn't open file: " + fileLocation);
  }
  affx::File5_Group *group = getIntensityGroup();
  m_OrderProbeStorage = group->openVector("/intensities/order", FILE5_DTYPE_INT, FILE5_OPEN);
  int size = m_OrderProbeStorage->reserved();
  m_OrderProbe.resize(size);
  m_OrderProbeStorage->read_array(0,size,&m_OrderProbe[0]);
  m_NumProbe = size;
  // 
  m_Chips = group->openTsv(CEL_LABEL, FILE5_OPEN);    
  m_NumChannels = 1;
  m_NumChips = m_Chips->getColumnCount(0);
}

void DataStore::initFromDataStore(const DataStore &ds) {
  closeResources();
  int ret = m_File5.open(m_FileStorage, affx::FILE5_OPEN_CREATE | affx::FILE5_REPLACE);
  if(ret != FILE5_OK) {
    Err::errAbort("Couldn't open file: " + m_FileStorage);
  }
  initVariables();
  setProbeOrder(ds.getProbeOrder());
  initIntensity(ds.getProbeCount(), ds.getChannelCount());
  m_CelNames = ds.m_CelNames;
  m_NumChannels = ds.m_NumChannels;
  m_NumChips = ds.m_NumChips;
}

const std::vector<int> &DataStore::getProbeOrder() const {
  return m_OrderProbe;
}

void DataStore::slurpDataFromCelFiles(const std::vector<std::string> &celFiles) {
  for(int fileIx = 0; fileIx < celFiles.size(); fileIx++) {
    readCelFile(celFiles[fileIx]);
  }
}

int DataStore::getDataSetIx(int chipIx, int channelIx) const {
  return m_NumChannels*chipIx+channelIx; // So chip 0 channel 0 is 0 chip 1 channel 0 is 2, etc
}

std::string DataStore::getCelName(int dataSetIx) const {
  return m_CelNames[dataSetIx / m_NumChannels];
}

void DataStore::setBufferSize(int bytes) const {
    int perVector = bytes / m_Chips->getColumnCount(0);
    for (int i = 0; i < m_Chips->getColumnCount(0); i++) {
        affx::File5_TsvColumn *col = m_Chips->getColumnPtr(0, i);
        col->setBufferSize(perVector);
    }
}

float DataStore::getProbeIntensity(probeid_t probeIx, chipid_t chipIx, unsigned int channelIx) const {
  float value = -1;
    
  int dataSetIx = getDataSetIx(chipIx, channelIx);
  if(!FILE5_OK == m_Chips->getLine(0, dataSetIx, m_MapProbe[probeIx], &value)) {
    Err::errAbort("Error reading probe " + ToStr(probeIx) + " chip: " + ToStr(chipIx) + " channel: " + ToStr(channelIx));
  }
  return value;
}

void DataStore::setProbeIntensity(float value, probeid_t probeIx, chipid_t chipIx) {
  m_Chips->gotoLine(m_MapProbe[probeIx]);
  m_Chips->set_f(0, chipIx, value);
}

std::vector<float> DataStore::getCelData(int dataSetIx) { 
  std::vector<float> data; 
  fillInCelData(dataSetIx, data);
  return data; 
}

std::vector<float> DataStore::getCelData(chipid_t celIx, 
                                         unsigned int channelIx) {
  int dataSetIx = getDataSetIx(celIx, channelIx);
  std::vector<float> data = getCelData(dataSetIx);
  return data;      
}

void DataStore::setCelFileNames(const std::vector<std::string> &names) { 
  m_CelNames = names;
};

void DataStore::checkSetup(const std::string &msg) {
  if(m_OrderProbe.empty()) {
    APT_ERR_ABORT( "Must set probe order before: " + msg);
  }
}


void DataStore::getProbePm(std::vector<bool> &pmBool) {
  std::vector<char> pm;
  getProbePm(pm);
  pmBool.resize(pm.size());
  fill(pmBool.begin(), pmBool.end(), false);
  for(int i = 0; i < pm.size(); i++) {
    if(pm[i] > 0) {
      pmBool[i] = true;
    }
  }
}

void DataStore::getProbeGcVec(std::vector<int> &vec) {
  checkSetup("getProbeGcVec");
  int colIx = m_ColIndexesProbe[ProbeGc];
  APT_ERR_ASSERT(colIx >= 0, "ProbesGc not set.");
  affx::File5_TsvColumn* column = m_TsvProbe->getColumnPtr(0,colIx);
  vec.resize(m_NumProbe);
  fill(vec.begin(), vec.end(), -1);
  int read = column->read_vector(vec.size(), &vec);
  APT_ERR_ASSERT(read == m_NumProbe, "Didn't get a big enough read.");
}

DataStore* DataStore::copyMetaDataToEmptyMart() const {
  DataStore* newMart = new DataStore();
  return newMart;
}

affx::File5_Group *DataStore::getIntensityGroup() {
  if(m_IntensityGroup == NULL) {
    m_IntensityGroup = m_File5.openGroup("intensities", affx::FILE5_OPEN_CREATE);
  }
  return m_IntensityGroup;
}

affx::File5_Group *DataStore::getGroupProbe() {
  if(m_GroupProbe == NULL) {
    m_GroupProbe = m_File5.openGroup("probes", affx::FILE5_OPEN_CREATE);
  }
  return m_GroupProbe;
}

affx::File5_Group *DataStore::getGroupProbeSet() {
  if(m_GroupProbeSet == NULL) {
    m_GroupProbeSet = m_File5.openGroup("probe-sets", affx::FILE5_OPEN_CREATE);
  }
  return m_GroupProbeSet;
}

affx::File5_Tsv *DataStore::getTsvProbe() {
  if(m_TsvProbe == NULL) {
    affx::File5_Group *probeGroup = getGroupProbe();
    m_TsvProbe  = probeGroup->openTsv(PROBE_LABEL, affx::FILE5_OPEN_CREATE);
  }
  return m_TsvProbe;
}

affx::File5_Tsv *DataStore::getTsvProbeSet() {
  if(m_TsvProbe == NULL) {
    affx::File5_Group *probeSetGroup = getGroupProbeSet();
    m_TsvProbe  = probeSetGroup->openTsv(PROBESET_LABEL, affx::FILE5_OPEN_CREATE);
  }
  return m_TsvProbe;
}

void DataStore::initChannels(int numChannels) {
  // Allocate a tsv for each channel
  affx::File5_Group *group = getIntensityGroup();
  m_NumChannels = numChannels;
  m_Chips = group->openTsv(CEL_LABEL, FILE5_OPEN_CREATE);
}

void DataStore::initFromCelFile(const std::string &celName) {

  // How many channels are there?
  FusionCELData cel;
  std::string tmp_unc_name=Fs::convertToUncPath(celName);
  cel.SetFileName(tmp_unc_name.c_str());
  if(!cel.ReadHeader()) {
    Err::errAbort("\nCan't read cel file header: " + cel.GetFileName() + 
                  "\n>>> Error reported: " + StringUtils::ConvertWCSToMBS(cel.GetError()));
  }
  initFromCelFile(cel);
  cel.Close();

  initChannels(m_NumChannels);
}

void DataStore::initFromCelFile(affymetrix_fusion_io::FusionCELData &cel) {
  
  // How many probes are there
  m_NumProbe = cel.GetRows() * cel.GetCols();
  
  // How many channels are there?
  std::vector<std::wstring> dataChannels = cel.GetChannels();
  m_NumChannels = dataChannels.size() > 0 ? dataChannels.size() : 1;
  
  initChannels(m_NumChannels);
}

void DataStore::initIntensity(int numProbes, int numChannels) {
  m_NumProbe = numProbes;
  initChannels(numChannels);
}
