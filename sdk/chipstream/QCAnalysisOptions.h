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

#ifndef _QCANALYSISOPTIONS_H_
#define _QCANALYSISOPTIONS_H_

//
#include "chipstream/AdjHomHiLoCelListener.h"
#include "chipstream/CelListener.h"
#include "chipstream/DmListener.h"
#include "chipstream/EmGenderCelListener.h"
#include "chipstream/HomHiLoCelListener.h"
#include "chipstream/MultiChannelHomHiLoCelListener.h"
#include "chipstream/MultiChannelNonOverlapCelListener.h"
#include "chipstream/QCProbesetOptions.h"
#include "chipstream/SignalBackgroundCelListener.h"
#include "chipstream/VariabilityScoreCelListener.h"
//
#include <cstring>
#include <map>
#include <string>
#include <utility>
#include <vector>
//

/// A simple class which contains the basic specification for a QC
/// analysis method. Typically populated from a row in the QCA file.
class QCMethod {

    public:
        /// Enumeration of possible analysis types
        enum AnalysisType {
            DM,             ///< DM Call Rate
            GENDER,         ///< EM ChrX Het Rate Gender Call
            HOMHILO,        ///< Min Diff between Hom Peak and Valley
            ADJHOMHILO,     ///< Adjusted Hom Min Hi Low used for SNP6
            //
            MULTICHANNEL_DISHQC,
            MULTICHANNEL_HOMHILO,
            MULTICHANNEL_SIGNALBACKGROUND,
            MULTICHANNEL_VARIABILITY,
        };

        /// Virtual Destructor
        virtual ~QCMethod() { }

        std::string analysisName; ///< The name of the analysis
        std::string groupName; ///< The group of probesets to use
        AnalysisType analysis; ///< The analysis type
        std::map<std::string,std::string> options; ///< Analysis options

        ///@todo it would be nice to handle more sophisticated
        ///      options (ie not just doubles). Something
        ///      like ChipStream and QuantMethod analysis
        ///      string configurations

        /**
         * Create a CelListener for the encapsulated method
         * @param layout - ChipLayout info w/ QC probesets loaded
         * @param celFiles - vector of cel file names to process
         * @param psOpts - QCProbesetOptions instance with info on QC probesets
         * @return pointer to CelListener. caller needs to free.
         */
        CelListener *createCelListener(ChipLayout *layout,
                                       vector<string> &celFiles,
                                       QCProbesetOptions &psOpts);

        /**
         * Create a DmListener from an DM analysis method
         * @param layout - ChipLayout info w/ QC probesets loaded
         * @param probesetNames - vector of QC probesets to use
         * @param celFiles - vector of cel file names to process
         * @return pointer to DmListener. caller needs to free.
         */
        DmListener *createDmListener(ChipLayout *layout,
                                     const std::vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                     const std::vector<std::string>& celFiles);

        /**
         * Create a EmGenderCelListener from a GENDER analysis method
         * @param layout - ChipLayout info w/ QC probesets loaded
         * @param probesetNames - vector of QC probesets to use
         * @param celFiles - vector of cel file names to process
         * @return pointer to EmGenderCelListener. caller needs to free.
         */
        EmGenderCelListener *createEmGenderCelListener(
                                     ChipLayout *layout,
                                     const std::vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                     const std::vector<std::string>& celFiles);

        /**
         * Create a HomHiLoCelListener from a homhilo analysis method
         * @param layout - ChipLayout info w/ QC probesets loaded
         * @param probesetNames - vector of QC probesets to use
         * @param celFiles - vector of cel file names to process
         * @return pointer to HomHiLoCelListener. caller needs to free.
         */
        HomHiLoCelListener *createHomHiLoCelListener(
                                     ChipLayout *layout,
                                     const std::vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                     const std::vector<std::string>& celFiles);

  /**
   * Create a AdjHomHiLoCelListener from a adj-homhilo analysis method
   * @param layout - ChipLayout info w/ QC probesets loaded
   * @param psOpts - QCProbesetOptions
   * @param celFiles - vector of cel file names to process
   * @return pointer to AdjHomHiLoCelListener. caller needs to free.
   */
  AdjHomHiLoCelListener*
  createAdjHomHiLoCelListener(ChipLayout *layout,
                              QCProbesetOptions &psOpts,
                              std::vector<std::string> &celFiles);

  
  // MULTICHANNEL_DISHQC
  MultiChannelNonOverlapCelListener*
  createMultiChannelNonOverlapCelListener(ChipLayout *layout,
                                          const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                          QCProbesetOptions &psOpts);
  // MULTICHANNEL_HOMHILO
  MultiChannelHomHiLoCelListener*
  createMultiChannelHomHiLoCelListener(ChipLayout* layout,
                                       const std::vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                       QCProbesetOptions& psOpts);
  // MULTICHANNEL_SIGNALBACKGROUND
  SignalBackgroundCelListener*
  createSignalBackgroundCelListener(ChipLayout *layout,
                                    const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                    QCProbesetOptions &psOpts);
  // MULTICHANNEL_VARIABILITY
  VariabilityScoreCelListener*
  createVariabilityScoreCelListener(ChipLayout *layout,
                                    const vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                                    QCProbesetOptions &psOpts);

  private:
  /**
   * Create a vector of probeLists
   */
  void fillProbeListVector(std::vector<ProbeListPacked>& probeSets,
                           const std::vector<QCProbesetOptions::nameLigationBase_t>& probesetNames,
                           ChipLayout* layout);

  void fillProbeListVector(std::vector<ProbeListPacked>& probeSets,
                           const std::vector<std::string>& probesetNames,
                           ChipLayout* layout);

};

/// A simple class to contain the multiple analysis specified
/// in a QCA file.
class QCAnalysisOptions {
    public:
        std::string defaultAnalysisName; ///< The default analysis name
        std::vector<QCMethod> methods; ///< Vector of analysis methods

        /**
         * Read in a QCA File and populate a QCAnalysisOptions instance
         * @param qcaFileName - the name of the qca file to load
         * @param opts - the QCAnalysisOptions instance to populate
         */
        static void readQCAFile(const std::string &qcaFileName,
                                QCAnalysisOptions &opts);

        /// Virtual Destructor
        virtual ~QCAnalysisOptions() { }

        static float ToFloat(std::string valString);

    private:
        /**
         * Private utility method to split up the options string
         * @param str - the string to parse
         * @param val - the parsed out key=val pair to be filled in
         * @param sep1 - separator between key and val (default '=')
         * @param sep2 - separator between options (default ',')
         * @param return - bool indicating if there are more options
         */
        static bool SplitOptions(std::string& str,
                                 std::pair<std::string,std::string>& val,
                                 const char sep1='=',
                                 const char sep2=',');
};

#endif // _QCANALYSISOPTIONS_H_
