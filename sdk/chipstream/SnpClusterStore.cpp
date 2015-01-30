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

#include "chipstream/QuantLabelZIO.h"
#include "chipstream/SnpClusterStore.h"
#include "chipstream/SnpModelConverter.h"
#include "util/Err.h"
#include "util/Fs.h"
//
#include <algorithm>
#include <map>
//
using namespace std;
using namespace affx;

SnpClusterStore::SnpClusterStore(const std::string &fileName, const std::string &tmpFileName) {
    if (SnpModelDb::isSqlModelFile(fileName)) {
        m_SnpDb.open(fileName);
    }
    else {
        if (tmpFileName.empty()) {
            TmpFileFactory tmpFac;
            m_TmpFile = tmpFac.genFilename_basic("SnpClusterStore", ".tmp");
        }
        else {
            m_TmpFile = tmpFileName;
        }
        Verbose::out(1, "Converting db to file: " + m_TmpFile);
        SnpModelConverter conv;
        conv.convertToDbModel(fileName, m_TmpFile);
        m_SnpDb.open(m_TmpFile);
    }
}

SnpClusterStore::~SnpClusterStore() {
    m_SnpDb.close();
    if (!m_TmpFile.empty()) {
        Fs::rm(m_TmpFile, false);
    }
}

bool SnpClusterStore::snpClusterExists(const std::string &snpName, int copyNumber) const {
    bool found = m_SnpDb.getSnpDistribution(snpName, copyNumber, m_DistCache);
    return found;
}

snp_distribution SnpClusterStore::getSnpCluster(const std::string &snpName, int copyNumber) const {
  bool found = false;
  string name = copyNumber == 1 ? snpName + ":1" : snpName;
  if (m_DistCache.probeset_id == name) {
      found = true;
      return m_DistCache.Dist;
  }
  else {
      found = m_SnpDb.getSnpDistribution(snpName, copyNumber, m_DistCache);

  }
  APT_ERR_ASSERT(found, "No model for: " + snpName + " copynumber: " + ToStr(copyNumber));
  return m_DistCache.Dist;
}

void SnpClusterStore::getHeaderValue(const std::string &key, std::vector<std::string> &values) const {
    m_SnpDb.getMetaValue(key, values);
}
