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
 * @file   QuantMethodRunReport.cpp
 * @author Chuck Sugnet
 * @date   Fri Jul 28 12:12:05 2006
 * 
 * @brief Summarize and output basic statistics on the performance of
 * a particular summarization or detection method across a number of
 * cel files.
 */

//
#include "chipstream/QuantMethodRunReport.h"
//
#include "file/TsvFile/TsvFile.h"
#include "stats/stats.h"
#include "util/Fs.h"

using namespace std;
using namespace affx;

/** Constructor. */
QuantMethodRunReport::QuantMethodRunReport(const std::vector<std::string> &chipNames) {
  m_ChipNames = chipNames;
}

bool QuantMethodRunReport::prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
    return true;
}

/** 
 * After every probeset computation this function is called an is an opportunity
 * to query the quantification method for results, residuals, etc.
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * @return true if success, false otherwise.
 */  
bool QuantMethodRunReport::report(ProbeSetGroup &psGroup, QuantMethod &qMethod,
                                  const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                                  PmAdjuster &pmAdjust) {
  // Nothing to do as we cull stats from the ChipSummary interface
  return true;
}

/** Print out our stats as a tab delmited text file. */
bool QuantMethodRunReport::finish() {
  // There are two passes here: 1) determine what columns we are going to 
  // print and outputting their headers. 2) Loop through all of the 
  // cel files and output their statistics.
  
  // Here are the column headers.

  // Column for cel file names.
  int cidx = 0;
  defineStringColumn(0, cidx++, "cel_files",TSVREPORT_CELFILE_STRLEN);
  
  // Column for each metric
  for(int source=0; source<m_ChipSummaries.size(); source++) {
    addHeaderComments(m_ChipSummaries[source]->getHeaderComments());
    ChipSummary::metricDefVec_t metricdefs=m_ChipSummaries[source]->getMetricDefs();
    for(int m=0; m<metricdefs.size(); m++) {
      if(metricdefs[m].m_type == ChipSummary::Metric::Integer) {
        defineColumn(0, cidx, metricdefs[m].m_name,affx::FILE5_DTYPE_INT);
      }
      else if (metricdefs[m].m_type == ChipSummary::Metric::Double) {
        defineColumn(0, cidx, metricdefs[m].m_name,affx::FILE5_DTYPE_DOUBLE);
        int p=metricdefs[m].m_precision;
        if (p>0) {
          setPrecision(0, cidx, p);
        }
      } 
      else if (metricdefs[m].m_type == ChipSummary::Metric::String) {
        defineStringColumn(0, cidx, metricdefs[m].m_name,TSVREPORT_METRIC_STRLEN);
      }
      else {
        Err::errAbort("QuantMethodRunReport: Unable to handle unknown type: " + ToStr(metricdefs[m].m_type) );
      }
      // next col
      cidx++;
    }
  }

  // Meta Info and Comments
//  std::vector<std::string>::iterator keysIx;
//  std::vector<std::string>::iterator valsIx;
//  for(keysIx = m_HeaderKeys.begin(), valsIx = m_HeaderVals.begin();
//      keysIx != m_HeaderKeys.end() && valsIx != m_HeaderVals.end();
//      ++keysIx, ++valsIx) {
//    addHeader(*keysIx, *valsIx);
//  }
//  for(valsIx = m_OtherLines.begin(); valsIx != m_OtherLines.end(); ++valsIx) {
//          tsv.addHeaderComment(*valsIx);
//  }

  // Initialize the file
  writeTsv_v1();
  
  // Print the statistics for each individual chip.
  for(uint32_t chipIx = 0; chipIx < (uint32_t)m_ChipNames.size(); ++chipIx) {
    int count = 0;
    std::vector<CumulativeStats<double> >::iterator celProbeStatsIx;
    
    // cel file name
    set_string(0, count++, Fs::basename(m_ChipNames[chipIx]));

    // Chip Summary metrics
    for(int source=0; source<m_ChipSummaries.size(); source++) {
      std::vector<ChipSummary::Metric> metrics = m_ChipSummaries[source]->getMetrics(chipIx);
      for(int m=0; m<metrics.size(); m++) {
        if (metrics[m].m_Type == ChipSummary::Metric::Integer) {
          set_i(0, count, metrics[m].m_Integer);
        } 
        else if (metrics[m].m_Type == ChipSummary::Metric::Double) {
          set_d(0, count, metrics[m].m_Double);
        } 
        else if(metrics[m].m_Type == ChipSummary::Metric::String) {
          set_string(0, count, metrics[m].m_String);
        } else {
          Err::errAbort("QuantMethodRunReport: Unable to handle unknown type: " + ToStr(metrics[m].m_Type) );
        }
        // next col
        count++;
      }
    }
    writeLevel(0);
  }
  close();

  return true;
}

void QuantMethodRunReport::registerChipSummary(ChipSummary *summary) {
      m_ChipSummaries.push_back(summary);
}

/* Catch meta header info from engine */
/*
void QuantMethodRunReport::out(const std::string &s, bool comment) {
      std::vector<std::string> words;
      Util::chopString(s, '=', words);
      if(words.size() == 2) {
          m_HeaderKeys.push_back(words[0]);
          m_HeaderVals.push_back(words[1]);
      }
      else if(s.find("=") == s.length() - 1) {
          std::string temp = s;
          Util::trimString(temp, "=");
          m_HeaderKeys.push_back(temp);
          m_HeaderVals.push_back("");
      }
      else {
          m_OtherLines.push_back(s);
      }
}
*/
