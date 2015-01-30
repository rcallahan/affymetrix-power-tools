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
 * @file   SummaryStats.cpp
 * @author Chuck Sugnet
 *
 * @brief  Utility class for keeping track of summary statistics.
 *
 */

//
#include "chipstream/SummaryStats.h"
//
#include "util/Verbose.h"
//
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
//
//
//

void SummaryStats::addData(const std::vector<double> &dat)
{
  assert(m_Counts.size() == m_Means.size());
  assert(m_Counts.size() == m_S.size());
  assert(m_Counts.size() == dat.size());
  for (unsigned int i = 0; i < m_Counts.size(); i++) {
    m_Counts[i]++;
    double delta = dat[i] - m_Means[i];
    m_Means[i] += delta / m_Counts[i];
    m_S[i] += delta * (dat[i] - m_Means[i]);
  }
}

void SummaryStats::resize(int size)
{
  m_Counts.clear();
  m_Counts.resize(size);
  m_Means.clear();
  m_Means.resize(size);
  m_S.clear();
  m_S.resize(size);
}

double SummaryStats::getStdev(int index)
{
  assert(index < m_Counts.size());
  if (m_Counts[index] <= 1) {
    Verbose::out(1, "Warning! SummaryStats::getStdev() - Must have counts > 1 for calculating stdev.");
    return 0;
  }
  return sqrt(m_S[index] / (m_Counts[index] - 1));
}

double SummaryStats::getMean(int index)
{
  return m_Means[index];
}

void SummaryStats::reportIndex(std::ofstream &out, int index)
{
  out << "\t" << m_Means[index];
  out << "\t" << getStdev(index);
}

void SummaryStats::report(std::ofstream &out, const std::string &name)
{
  out << name << "_mean";
  for (unsigned int i = 0; i < m_Counts.size(); i++) {
    out << "\t" << m_Means[i];
  }
  out << std::endl;
  out << name << "_stdev";
  for (unsigned int i = 0; i < m_Counts.size(); i++) {
    out << "\t" << getStdev(i);
  }
  out << std::endl;
}

