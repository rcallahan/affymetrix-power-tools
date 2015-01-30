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

#ifndef _PROBESETSUMMARIZESTAGEENGINE_H
#define _PROBESETSUMMARIZESTAGEENGINE_H

//
#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/SpfReader.h"
#include "chipstream/AptTypes.h"
#include "chipstream/CelStatListener.h"
#include "chipstream/MetaProbeset.h"
#include "chipstream/ProbeListFactory.h"
#include "AnalysisStage.h"
#include "chipstream/QuantDataTransform.h"
#include "util/BaseEngine.h"
//
#include <cstring>
#include <ctime>
#include <ostream>
#include <string>
#include <vector>

class ProbesetSummarizeStageEngine : public BaseEngine {

public:

  virtual std::string getEngineName() { return ProbesetSummarizeStageEngine::EngineName(); }
  static const std::string EngineName() { return "ProbesetSummarizeStageEngine"; }

  /**
   * Constructor
   */
  ProbesetSummarizeStageEngine();

  /**
   * Destructor
   */
  ~ProbesetSummarizeStageEngine();

  /*! A class to register the summary engine. */
  class Reg : public EngineReg
  {
  public:
    /*! Constructor - register the summary engine. */
    Reg() : EngineReg(ProbesetSummarizeStageEngine::EngineName())
    {
    }

    /*! Creates an object.
     * @return The object.
     */
    BaseEngine *MakeObject() { return new ProbesetSummarizeStageEngine; }
  };

  /*! The one and only registration object. */
  static Reg reg;

  /*! Converts the type to the summary engine type.
   * @param chip The pointer to the base engine object.
   * @return The summary engine type or NULL if not compatible.
   */
  static ProbesetSummarizeStageEngine * FromBase(BaseEngine *engine);

private:

  void printStandardMethods(std::ostream &out);
  void extraHelp();
  void explain();
  void defineOptions();
  void defineStates();
  void defineStdMethods();
  void checkOptionsImp();
  void fillInLayout(ChipLayout **cLayout, string &spfFile, 
                    string &spfFileOrig, std::vector<MetaProbeset *> &metaSets,
                    std::vector<const char *> &probesetNames);
  void runImp();
  void writeHeaders(IntensityMart &iMart,
                    QuantMethod *qMethod,
                    std::vector<QuantMethodReport *> &reporters,
                    AnalysisInfo &info);
  void mergeProbeSubset(const std::vector<ProbeListPacked> &plVec, 
                        std::vector<bool> &probeSubset);

  AnalysisInfo makeAnalysisInfo(ChipLayout &layout, 
                                const std::string &name,
                                vector<MetaProbeset *> &metaSets,
                                const std::string &guid,
                                const std::string &qSpec,
                                QuantMethodFactory::QuantType type,
                                PsBoard &board);
                                

  void setupReporters(std::vector<QuantMethodReport *> &reporters, 
                      AnalysisInfo &info,
                      QuantMethod *eMethod, 
                      vector<MetaProbeset *> &metaSets,
                      CelStatListener *stats);

/*   void setupReporters( */
/*                       AnalysisStreamExpression *as,  */
/*                       QuantExprMethod *eMethod,  */
/*                       std::vector<MetaProbeset *> &metaSets,  */
/*                       CelStatListener *stats); */
  void setupCelStats(
                     CelStatListener &celStats, 
                     const std::vector<Probe *> &controlProbes, 
                     ChipLayout &layout);
  void loadCdfLayout(
                     ChipLayout &layout, 
                     std::vector<MetaProbeset *> &metaToLoad,
                     const std::vector<Probe *> &probesToLoad, 
                     std::vector<const char *> *probesetNames,
                     std::vector<bool> &probeSubset, 
                     bool justStats, 
                     bool doAll,
                     probeidmap_t &killList);
  void loadPgfLayout(
                     ChipLayout &layout, 
                     std::vector<MetaProbeset *> &metaToLoad,
                     const std::vector<Probe *> &probesToLoad, 
                     std::vector<const char *> *probesetNames,
                     std::vector<bool> &probeSubset, 
                     probeidmap_t &killList,
                     bool justStats, 
                     bool doAll);
  void loadChipLayout(
                      ChipLayout &layout, 
                      //const std::vector<ProbeListPacked> &plVec,
                      std::vector<MetaProbeset *> &metaToLoad,  
                      probeidmap_t &killList,
                      std::vector<const char *> *probesetNames,
                      const std::vector<Probe *> &probesToLoad, 
                      std::vector<bool> &probeSubset, 
                      bool justStats, 
                      bool doAll);
  ProbeSetGroup *makeProbesetGroupFromMeta(
                                           MetaProbeset &meta, 
                                           ChipLayout &layout);
  void doSummaries(SpfReader &reader, IntensityMart &iMart,
                   AnalysisStreamExpression *analysis,
                   int numPsSets);

  void makeSpfFile(std::string inSpf, std::string outSpf, vector<MetaProbeset *> metaSets);
  void addProbeSet(ProbeSet *basicPs, ProbeSet *currentPs);
  void determineDesiredOrder(
                             ChipLayout &layout, 
                             vector<int> &desiredOrder, 
                             bool haveSubset,
                             vector<MetaProbeset *> &metaSets);

  void fillInAnalysisInfo(AnalysisInfo &info, 
                          QuantMethod *qMethod,
                          std::string prefix = "apt-");
  
  void runChipStreamStage(AnalysisStage &stage, PsBoard &board, 
                          DataStore &in, DataStore &out);

  void runQuantificationStage(AnalysisStage &stage, PsBoard &board, 
                              DataStore &in,
                              std::vector<QuantMethodReport *> &reporters,
                              AnalysisInfo &info,
                              QuantMethodFactory::QuantType type,
                              vector<MetaProbeset *> &metaSets,
                              CelStatListener *celStats,
                              const std::string &spfFile);
  void initializeAsFactory(AnalysisStreamFactory &asFactory, ChipLayout *layout);
  void setupMask(ChipLayout *layout, map<string, vector<bool> > masks, DataStore &probeInfo);
  void readMetaProbesets(std::vector<MetaProbeset *> &metaSets);
  void closeGlobalA5();
  void addReporters_Expression(affx::TsvReport::TsvReportFmt_t report_format,
                               std::vector<QuantMethodReport *> &reporters, 
                               AnalysisInfo& info,
                               QuantMethod *qMethod,
                               int precision);

private:
  /// List of aliases for standardized methods
  std::map<std::string,std::string> m_stdMethods; 
  /// Global A5 Input File
  affx::File5_File* m_a5_global_input_file;
  affx::File5_Group* m_a5_global_input_group;
  /// Global A5 Output File
  affx::File5_File* m_a5_global_output_file;
  affx::File5_Group* m_a5_global_output_group;
};

#endif
