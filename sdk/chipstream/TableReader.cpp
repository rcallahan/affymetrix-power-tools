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

//
#include "chipstream/TableReader.h"
//
#include "util/Convert.h"
#include "util/Err.h"
#include "util/TableFile.h"
#include "util/Verbose.h"
//
#include <vector>
//

using namespace std;
/** 
 * @brief Read the file that contains the intensity data.
 * @return true if opened sucessfully, false otherwise.
 */
bool TableReader::readFiles() {
  unsigned int clfIx = 0, colIx = 0, index = 0, rowIx = 0;
  /* Use this to pass data to streams and intensity marts. */
  vector<float> chip(m_ClfProbeCount);
  vector<int> clfMap(m_ClfProbeCount);
  
    /* Sanity check. */
  if(m_FileNames.size() != 1) {
    Err::errAbort("Can't have more than one text file.");
  }

  /* Read into a tablefile. This is inefficient, but easy for now
     as already coded. */
  TableFile tf;
  if(!m_RowColNames) {
    tf.setUseColNames(false);
    tf.setUseRowNames(false);
  }
  tf.open(m_FileNames[0].c_str());
  /* Make the map once (rather than looking up for each chip. */
  for(clfIx = 0; clfIx < clfMap.size(); clfIx++) {
    clfMap[clfIx] = tf.rowIndex(ToStr(clfIx).c_str());
  }
  Verbose::progressBegin(3, "Reading " + ToStr(m_FileNames.size()) + " Files", m_FileNames.size(), 0, m_FileNames.size());
  //Verbose::setDot(tf.numRows());
  for(rowIx = 0; rowIx < tf.numRows(); rowIx++) {
    Verbose::progressStep(2);
    for(colIx = 0; colIx < tf.numCols(); colIx++) {
      float intensity = Convert::toFloat(tf.getData(rowIx, colIx).c_str());
      chip[colIx] = intensity;
    }
    /* Load data into our intensity marts. */
    for(index = 0; index < m_IntenMarts.size(); index++) 
      m_IntenMarts[index]->setProbeIntensity(rowIx, chip);
    /* Push data through our streams, have to make a new copy for
       each stream as it may get modified along the way. */
    for(index = 0; index < m_Streams.size(); index++) {
      std::vector<float> copy(chip);
      m_Streams[index]->newChip(copy);
    }
  }
  Verbose::progressEnd(2, "Done.");
  for(index = 0; index < m_Streams.size(); index++) {
    m_Streams[index]->endDataSet();
  } 
  m_Read = true;
  m_Passed = true;
  return true;
}     
