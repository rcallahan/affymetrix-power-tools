////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
 * @file   QuantAvgDiff.cpp
 * @author Pete Klosterman
 * @date   Mon Jun 26 11:31:23 2006
 *
 * @brief  Class for doing MAS 4 average difference signal estimation.
 */

//
#include "chipstream/QuantAvgDiff.h"

using namespace std;

/**
 * @brief Constructor.
 * @param param Map of key/value pairs to set user-defined parameters.
 */
QuantAvgDiff::QuantAvgDiff (std::map<std::string,std::string>& param)
{
  setupSelfDoc(*this);
  // Set the string descriptor of this method.
  m_Type = getDocName();
  m_ChipCount = 0;
  m_ProbeCount = 0;
}

/**
 * @brief Set the number of probes and number of chips to be used for estimation.
 * @param numProbes Probe count.
 * @param numChips Microarray count.
 */
void QuantAvgDiff::setBounds (unsigned int numProbes, unsigned int numChips)
{
  m_ProbeCount = numProbes;
  m_ChipCount = numChips;
  // Set things up the first time through.
  if (m_PM.empty())
  {
    m_SetSignal.resize (numChips);
    while (m_PM.size() < numChips)
    {
      m_PM.push_back (vector<float>(numProbes));
      m_MM.push_back (vector<float>(numProbes));
    }
  }
  // Resize the number of columns (probes).
  if (numProbes != m_PM[0].size())
  {
    const unsigned int PMSize = m_PM.size();
    for (unsigned int i = 0; i < PMSize; ++i)
    {
      m_PM[i].resize (numProbes);
      m_MM[i].resize (numProbes);
    }
  }
  // Increase the number of rows (chips).
  if (numChips > m_PM.size())
  {
    const unsigned int PM0Size = m_PM[0].size();
    m_SetSignal.resize (numChips);
    while (m_PM.size() < numChips)
    {
      m_PM.push_back (vector<float>(PM0Size));
      m_MM.push_back (vector<float>(PM0Size));
    }
  }
}

/**
 * @brief Set up the quantification method given all the data about the probe
 * set, chip and data.
 *
 * @param psGroup Probes to be used.
 * @param iMart Raw intensities.
 * @param iTrans Transformations to be applied to data before use.
 * @param pmAdjust How to estimate background or MM probe.
 * @return bool True if setup sucessful, false otherwise.
 */
bool QuantAvgDiff::setUp (ProbeSetGroup& psGroup, const IntensityMart& iMart, 
                          std::vector<ChipStream *>& iTrans, PmAdjuster& pmAdjust)
{
  const unsigned int chipCount = iMart.getCelFileCount();
  const unsigned int probeSetsSize = psGroup.probeSets.size();
  assert (probeSetsSize > 0);

  unsigned int atomPmCount = psGroup.countPmProbes();
  if (atomPmCount == 0)
    return false;
  // Prepare the quantification method for this many atoms and chips.
  setBounds (atomPmCount, chipCount);

  atomPmCount = 0;
  m_Probes.clear();
  // Loop through and fill in the data as coming from the intensity mart.
  for (unsigned int psIx = 0; psIx < probeSetsSize; ++psIx)
  {
    const ProbeSet *ps = psGroup.probeSets[psIx];
    if (ps == NULL)
    {
      Verbose::out (1, "No probeset for index: " + ToStr (psIx));
      continue;
    }
    const unsigned int atomsSize = ps->atoms.size();
    for (unsigned int atomIx = 0; atomIx < atomsSize; ++atomIx)
    {
      Atom& atom = *(ps->atoms[atomIx]);
      const unsigned int probesSize = atom.probes.size();
      unsigned int channelIx = atom.getChannelCode();
      for (unsigned int probeIx = 0; probeIx < probesSize; ++probeIx)
      {
        Probe *p = atom.probes[probeIx];
        if (p->type == Probe::PMST || p->type == Probe::PMAT)
	{
          const unsigned int probeIndex = p->id;
          // Add probe to our probes and fill in the intensity data for each chip.
          m_Probes.push_back (p);
          for (unsigned int chipIx = 0; chipIx < chipCount; ++chipIx)
	  {
            // PM transformation.
            float intensity = transformPrimaryData (probeIndex, chipIx, iMart, iTrans, channelIx);
            float mmIntensity = 0;
            // MM transformation.
	    ///@todo use channel to grab correct adjustment 
            pmAdjust.pmAdjustment (probeIndex, chipIx, iMart, iTrans, intensity, mmIntensity);
            setPMDataAt (atomPmCount, chipIx, intensity);
            setMMDataAt (atomPmCount, chipIx, mmIntensity);
          }
          ++atomPmCount;
        }
      }
    } // atoms.
  } // probe sets.

  return true;
}

/**
 * @brief Create a new AvgDiff object.
 *
 * @param param Map of key/value pairs to initialize the object.
 * @return Pointer to SelfCreate object.
 */
SelfCreate* QuantAvgDiff::newObject (std::map<std::string,std::string>& param)
{
  if(param.size() != 0) { Err::errAbort(ToStr("No parameters for ") + ToStr(QUANT_AVGDIFF_STR)); }
  return new QuantAvgDiff (param);
}

/**
 * @brief Perform the estimation.
 *
 * Calculate the average pm - mm difference.
 */
void QuantAvgDiff::computeEstimate()
{
  for (unsigned int chipIx = 0; chipIx < m_ChipCount; ++chipIx)
  {
    const unsigned int pmSize = m_PM[chipIx].size();
    assert (pmSize);
    assert (pmSize == m_MM[chipIx].size());
    double diff = 0;
    for (unsigned int i = 0; i < pmSize; ++i)
      diff += m_PM[chipIx][i] - m_MM[chipIx][i];
    m_SetSignal[chipIx] = diff / pmSize;
  }
}
