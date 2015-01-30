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

#ifndef _CUMULATIVESTATS_H_
#define _CUMULATIVESTATS_H_

#include "portability/affy-base-types.h"
#include "util/Err.h"
//
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
//

//#include "stats/stats.h"

/* Utility class for keeping track of summary statistics. */
template <class T, class LargeNumber = double>
class CumulativeStats {

  public:
  
  /** Constructor. */
  CumulativeStats() : m_Name(""), m_Counts(0),
  m_Max(0),     m_Min(0),
  m_Mean(0),    m_S(0), m_Full(false),
  m_Median(0), m_MedianTwo(0), m_IQR(0)
  {std::vector< LargeNumber > m_List;}
  
  /** Constructor. */
  CumulativeStats(const std::string &name) : m_Name(name), m_Counts(0),
  m_Max(0),     m_Min(0),
  m_Mean(0),    m_S(0), m_Full(false),
  m_Median(0), m_IQR(0)
  {std::vector< LargeNumber > m_List;}
  
  
  /** 
   * Add a vector of data, keeping a running tab of the mean, and
   * variance statistics. Algorithm from: Donald E. Knuth
   * (1998). The Art of Computer Programming, volume 2:
   * Seminumerical Algorithms, 3rd edn., p. 232. Boston:
   * Addison-Wesley.
   *
   * @param dat - New vector of data with one data point per chip.
   */
  void addData(T dat) {
    if(m_Counts == 0) 
    m_Min = m_Max = dat;
    else {
      m_Min = min(m_Min, dat);
      m_Max = max(m_Max, dat);
    }
    m_Counts++;
    T delta = dat - m_Mean;
    m_Mean += delta / m_Counts;
    m_S += delta * (dat - m_Mean);
    if (m_Full == true){
      m_List.push_back(dat);
    }
  }

  /* Complains if full method is not specified
   * Complains if less than 20 data points are provided
   * Computes median and IQR
   */
  void setUpFullMethod() {
    // wasnt set up to produce median and IQR.
    if (m_Full != true) {
      Err::errAbort("CumulativeStats::setUpFullMethod(): "
                    "You need to set your CumulativeStats object "
                    "to full method for getting median");
    }
    // to few points to be meaningful.
    if (m_List.size() < 20) {
      Err::errAbort("CumulativeStats::setUpFullMethod(): "
                    "You need to have at least 20 numbers "
                    "("+ToStr(m_List.size())+" given) "
                    "to compute fullMethod statistics (median, IQR, etc.) ");
    }
    // ...to many to handle.
    if (m_List.size() > 50000) {
      Err::errAbort("CumulativeStats::setUpFullMethod: "
                    "You have more than 50000 points,"
                    "is this too much? "+ToStr(m_List.size()));
    }

    
    std::sort(m_List.begin(), m_List.end());
    int quartile_ix = (m_List.size() / 4);
    int third_ix = quartile_ix + quartile_ix;
    m_IQR = m_List[third_ix] - m_List[quartile_ix];
    int odd = m_List.size() % 2;
    if (odd == 0) { //even
      m_Median = (m_List[(m_List.size() /2)] +m_List[(m_List.size()/2)-1])/2;
    }
    else {
      m_Median = m_List[(m_List.size() /2)];
    }
  }
  
  /* Accessors. */
  inline LargeNumber getStdev() { return static_cast<LargeNumber>(sqrt(m_S / (m_Counts - 1))); } 
  inline LargeNumber getMean() { return static_cast<LargeNumber>(m_Mean); } 
  inline T getMax() { return m_Max; }
  inline T getMin() { return m_Min; }
  inline void setFull(){ m_Full = true;}
  inline uint64_t getCount() { return m_Counts; }
  inline std::string getName() { return m_Name; }
  inline LargeNumber getMedian() { return static_cast<LargeNumber>(m_Median);}
  inline LargeNumber getIQR() { return static_cast<LargeNumber>(m_IQR);}

  private:
  
  std::string m_Name; ///< Name for this cumulative statistic.
  uint64_t m_Counts;  ///< How many examples have been seen.
  T m_Max;            ///< Maximum value seen.
  T m_Min;            ///< Minimum value seen.
  LargeNumber m_Mean; ///< Running mean calculation.
  LargeNumber m_S;    ///< Running variance calculation, divide by (n-1) for actual variance.
  bool m_Full;          ///< boolean for using full method or sparse
  std::vector<T> m_List;
  LargeNumber m_Median;
  LargeNumber m_MedianTwo;
  LargeNumber m_IQR;
};

#endif /* _CUMULATIVESTATS_H_ */
