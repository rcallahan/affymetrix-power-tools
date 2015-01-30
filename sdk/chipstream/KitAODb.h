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
 * @file   AdapterTypeNormTran.cpp
 *
 * @brief  Class for doing adapter type normalization
 */

#ifndef _KITAODB_H_
#define _KITAODB_H_

//
#include <map>
#include <set>
#include <string>


class KitAODbEntry {
public:
  int m_linenum;
  double m_pc1;
  double m_means;
  //
  KitAODbEntry();
  KitAODbEntry(int linenum,float pc1,float m_means);
};

class KitAODb {
public:
  //
  KitAODb();
  ~KitAODb();

  //
  void clear();

  //
  void addEntry(const std::string& name,KitAODbEntry* entry);
  KitAODbEntry* findEntry(const std::string& name);

  //
  int nameCount() const;
  int entryCount() const;
  int addedCount() const;

  //
  void readTsv(const std::string& filename);
  void writeTsv(const std::string& filename);

  //
  double getClassiferAValue() const;
  double getClassiferOValue() const;

  /// the filename which was read.
  std::string m_filename;
  /// the number of added entries.
  int m_added_count;
  /// The cutoff value read from the file headers.
  double m_classiferAValue;
  double m_classiferOValue;

  /// The number of probesets which are supposed to be in the db.
  int m_file_probeset_count;

  typedef std::map<std::string,KitAODbEntry*>::iterator name2entry_iter_t;
  typedef std::pair<std::string,KitAODbEntry*> name2entry_val_t;
  std::map<std::string,KitAODbEntry*> m_name2entry;

  typedef std::set<KitAODbEntry*>::iterator entry_iter_t;
  std::set<KitAODbEntry*> m_entries;

};

#endif // _KITAODB_H_
