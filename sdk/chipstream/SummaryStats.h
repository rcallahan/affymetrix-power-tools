////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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
 * @file   SummaryStats.h
 * @author Chuck Sugnet
 *
 * @brief  Utility class for keeping track of summary statistics.
 *
 */

#ifndef SUMMARYSTATS_H
#define SUMMARYSTATS_H

//
#include <cfloat>
#include <cmath>
#include <iostream>
#include <vector>
//
//

class SummaryStats
{

public:

  /**
   * Add a vector of data, keeping a running tab of the mean, and
   * variance statistics. Algorithm from: Donald E. Knuth
   * (1998). The Art of Computer Programming, volume 2:
   * Seminumerical Algorithms, 3rd edn., p. 232. Boston:
   * Addison-Wesley.
   *
   * @param dat - New vector of data with one data point per chip.
   */
  void addData(const std::vector<double> &dat);

  void resize(int size);

  double getStdev(int index);

  double getMean(int index);

  void reportIndex(std::ofstream &out, int index);

  void report(std::ofstream &out, const std::string &name);

private:
  std::vector<int> m_Counts; ///< How many examples have been seen.
  std::vector<double> m_Means; ///< Running mean calculation.
  std::vector<double> m_S; ///< Running variance calculation, divide by (n-1) for actual variance.
};

#endif /* SUMMARYSTATS_H */
