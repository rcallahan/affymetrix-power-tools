////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
 * @file   ADTOptions.h
 * @author Mybrid Spalding
 * @date   Thu Oct 30 15:21:21 2008
 *
 * @brief  Console and command line options.
 *
 *
 */


#ifndef TRANSLATION_ADT_OPTIONS_H
#define TRANSLATION_ADT_OPTIONS_H

//
#include <cstring>
#include <string>
#include <vector>
//

//! CONTROLLER: type implemented for the engine. 
const unsigned char C_CONSOLE = 0x01;
//! CONTROLLER: type implemented for the engine. 
const unsigned char C_CMDLINE = 0x02;

//! CHP or TSV file streamed by experiment to conserve memory. 
enum ADT_EXPERIMENT_STREAM_TYPE_ENUM {
  ADT_EXPERIMENT_STREAM_TYPE_NULL,
  ADT_EXPERIMENT_STREAM_TYPE_TSV,
  ADT_EXPERIMENT_STREAM_TYPE_CHP,
};

//! DMET2 and DMET3 are differentiated by column names
enum ADT_TRANSLATION_TABLE_TYPE_ENUM {
  ADT_TRANSLATION_TABLE_TYPE_INVALID,
  ADT_TRANSLATION_TABLE_TYPE_DMET2,
  ADT_TRANSLATION_TABLE_TYPE_DMET3
};


class ADTOptions
{

public:




  // PROGRAM OPTIONS \cond INGORE_OBVIOUS
  std::string m_audit;
  bool        m_dmet2Calling;
  bool        m_enforceCompleteHaplotypeGroup;
  bool        m_ignoreUnknownAlleles;
  bool        m_ignoreReportAllele;
  std::string m_progName;
  bool        m_useFirstDupAlleleDef;
  std::string m_useFirstDupAlleleDefHeaderText;

  // INPUT DATA OPTIONS
  bool m_prototypeCHPFiles;
  std::string m_inputAnnotationFile;
  std::string m_inputDir;
  std::string m_inputTTableFile;
  ADT_TRANSLATION_TABLE_TYPE_ENUM
  m_inputTTableType;
  std::string m_inputGenoFile;
  std::string m_inputGenotypeOverrideFile;

  std::string m_inputCopyFile;
  std::string m_inputMarkerListFile;
  std::vector< std::string > m_inputProbeSetVector;
  std::string m_inputSampleFile;
  int         m_inputSampleTable;
  std::vector< std::vector < std::string > >
  m_sampleTable;


  /// INPUT DATA, EXPERIMENT LIST OPTIONS
  ADT_EXPERIMENT_STREAM_TYPE_ENUM
  m_streamType;
  std::string m_inputExperimentList;
  std::string m_inputExperimentListFile;
  std::vector< std::string >
  m_inputExperimentFiles;
  std::vector< std::string >  m_inputExperimentListVector;


  // OUTPUT DATA OPTIONS

  std::string m_outputDir;
  std::string m_logFile;
  std::string m_outputReportPrefix;

  // REPORTING OPTIONS
  bool m_markerReport;
  bool m_profile;
  int  m_regression;
  bool m_summaryReportSort;
  int  m_verbosity;
  bool m_uncalledReportAllMarkers;
  
  ADTOptions() {
    /* Some defaults */
    m_dmet2Calling                   = false;
    m_enforceCompleteHaplotypeGroup  = false;
    m_ignoreUnknownAlleles           = false;
    m_ignoreReportAllele             = false;
    m_markerReport                   = false;
    m_profile                        = false;
    m_prototypeCHPFiles              = false;
    m_streamType                     = ADT_EXPERIMENT_STREAM_TYPE_NULL;
    m_summaryReportSort              = true;
    m_useFirstDupAlleleDef           = false;
    m_verbosity                      = 1;
    m_uncalledReportAllMarkers       = false;
    
  };

  ~ADTOptions() {
  };

  //! \endcond
  
};


#endif /* TRANSLATION_ADT_OPTIONS_H */
