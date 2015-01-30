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
 * @file   ExperimentReport.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:03:01 PDT 2008
 * @brief  Factory class for all translation reports to sub-class.
 */

#ifndef TRANSLATION_EXPERIMENT_REPORT_H
#define TRANSLATION_EXPERIMENT_REPORT_H

#include <cstring>
#include <map>
#include <string>
//

//

typedef enum _ExperimentReportTypeEnum {
  NULLREPORT,
  MARKER,
  HAPLOTYPE,
  REGRESSION,
  TRANSLATION,
} ExperimentReportTypeEnum;


class ExperimentReport
{
public:

  ExperimentReport() {};

  virtual ~ExperimentReport() {
    close(true);
  };

  /**
   * ExperimentReport::generate
   * Synopsis:
   * The report generation method to override.
   *
   * This gets called once per experiment so as to stream
   * output one experiment at a time.
   *
   * @param rte  - the single run time environment instance
   * @param ttm  - the single translation table model
   * @param er   - container for all CallResults for one experiment
   *
   * @return
   */
  virtual bool                      generate(class RunTimeEnvironment & rte,
      class TranslationTableModel & ttm,
      std::map< std::string, class ExperimentResults*> & er) {
    return false;
  }

  /**
   * ExperimentReport::name
   * Synopsis:
   *
   * Simply returns a constant std::string of the sub-classed report type.
   *
   *
   * @return name - the report name or type
   */
  virtual std::string               name() const {
    return std::string("ExperimentReport base class");
  }

  /**
   * ExperimentReport::report_type
   * Synopsis:
   *
   * Registering a new report requires adding an ENUM to the
   * report type enum in this file. This method returns that ENUM
   * specified.
   *
   *
   * @return ExperimentReportTypeEnum - NULLREPORT by default.
   */
  virtual ExperimentReportTypeEnum  report_type() {
    return NULLREPORT;
  }

  /**
   * ExperimentReport::close
   * Synopsis:
   *
   * Clean up memory, typically called in the destructor, but
   * can be called for reuse of report object.
   *
   * !!!! WARNING !!!! DO NOT USE Verbose:: if the 'abort'
   * parameter is passed. The abort indicates that
   * the Windows triggered an exception and the Verbose::
   * object gets comprised under these conditions.
   * All Verbose:: actions must be conditional upon abort.
   *
   *
   * @param abort - see warning above. Consoole sent an abort.
   *
   * @return true - close happened ok.
   
   */
  virtual bool                      close(bool abort) {
    return false;
  };


};


#endif /* TRANSLATION_EXPERIMENT_REPORT_H */
