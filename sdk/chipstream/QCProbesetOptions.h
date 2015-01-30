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

#ifndef _QCPROBESETOPTIONS_H_
#define _QCPROBESETOPTIONS_H_

//
#include <cstring>
#include <map>
#include <string>
#include <utility>
#include <vector>
//

/// A simple class which contains various QC probeset information.
/// Typically populated from the QCC file.
class QCProbesetOptions {
public:

  QCProbesetOptions();
  /// Indicates if the primary key is probeset name or ID
  bool useNames;

  /// A map keyed by the groupName(label) of probeset name (or IDs) vector
  // the extra char is for the ligation base.
  typedef std::pair<std::string,std::string> nameLigationBase_t;
  std::map<std::string, std::vector<nameLigationBase_t> > probesetGroups;

  bool m_hasLigationBase;

  /// A vector of probesets flagged for inclusion in CHP header
  std::vector<std::string> headerProbesets;

	/* need to carry the following forward */
	//bool multichannel;
	//std::vector<string> metaCelFileNames;
	//vector<vector<string> > metaCelFiles; //channel then file index
  std::string m_sketchInFile;
  std::string m_channelFile;
	std::string m_exprAnalysisString;

  // where the DiskIntensityMart should put its temp files.
  std::string m_tempdir;

  /**
   * Read in a QCC File and populate a QCProbesetOptions instance
   * @param qccFileName - the name of the qcc file to load
   * @param opts - the QCProbesetOptions instance to populate
   */
  static void readQCCFile(const std::string& qccFileName,
                          QCProbesetOptions& opts);

  /// Virtual Destructor
  virtual ~QCProbesetOptions() {
  };

private:
  /**
   * Private utility method to split up the group column
   * @param str - the string to parse
   * @param val - the parsed out group to be filled in
   * @param sep - separator between groups (default ',')
   * @param return - bool indicating if there are more groups
   */
  static bool SplitGroups(std::string& str, std::string& val,const char sep=',');
};

#endif // _QCPROBESETOPTIONS_H_
