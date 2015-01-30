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
 * @file   IdxGroup.h
 * @author Ray Wheeler
 * @date   Mon Feb  2 16:24:08 PST 2009
 * 
 * @brief  Named partitions of numbers (presumably indexes)
 */



#include "chipstream/IdxGroup.h"
//
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//


/**
 * @brief Insert group of indexes into grouping of given name
 * @param groupname - name of grouping
 * @param group - integers that belong to a group
 */
unsigned int IdxGroup::addGroup(std::string grouping_name, std::vector<unsigned int> group) {
  unsigned int grouping_idx;
  if (groupingNameExists(grouping_name)) {
    grouping_idx = m_grouping_name_to_grouping_idx[grouping_name];
    m_groups[grouping_idx].push_back(group);

  }
  else {
    m_grouping_name_to_grouping_idx.insert(std::pair<std::string,unsigned int>(grouping_name, m_grouping_count++));
    std::vector<std::vector<unsigned int> > new_grouping_vec(1, group);
    m_groups.push_back(new_grouping_vec);
    std::vector<unsigned int> new_element_idx_vec;
    m_element_idx_to_group_idx.push_back(new_element_idx_vec);
    grouping_idx = m_grouping_count - 1;
  }
  
  unsigned int new_group_idx = m_groups[grouping_idx].size() - 1;
  for (int i = 0; i < group.size(); i++) {
    if (group[i] == m_element_idx_to_group_idx[grouping_idx].size()) {
      m_element_idx_to_group_idx[grouping_idx].push_back(new_group_idx);
    }
    else {
      if (group[i] > m_element_idx_to_group_idx[grouping_idx].size()) {
	m_element_idx_to_group_idx[grouping_idx].resize(group[i] + 1, 0);
      }
      else { // if (group[i] < m_element_idx_to_group_idx[grouping_idx].size())
	if (m_element_idx_to_group_idx[grouping_idx][group[i]] != 0) {
	  Verbose::warn(1, "IdxGroup::addGroup -- Adding element index:" + ToStr(group[i]) + " more than once to same grouping name:" + grouping_name);
	}
      }
      m_element_idx_to_group_idx[grouping_idx][group[i]] = new_group_idx;
    }
  }
  return new_group_idx;
}


/**
 * @brief Get count of groups with given grouping name
 * @param groupname - name of grouping
 */
int IdxGroup::getGroupCount(std::string grouping_name) {
  int rtn_val = 0;
  if (groupingNameExists(grouping_name)) {
    rtn_val =  m_groups[m_grouping_name_to_grouping_idx[grouping_name]].size();
  }
  return rtn_val;
}

/**
 * @brief Get vector of vectors constituting the groups of indexes for given grouping name
 * @param groupname - name of grouping
 */
std::vector<std::vector<unsigned int> > IdxGroup::getGroupingVec(std::string grouping_name) {
  if (groupingNameRequire(grouping_name)) {
    return m_groups[m_grouping_name_to_grouping_idx[grouping_name]];
  }
  else {
    std::vector<std::vector<unsigned int> > empty;
    return empty;
  }
}

/**
 * @brief Get group at given index with given grouping name
 * @param groupname - name of grouping
 * @param group_idx - index of group in grouping vector
 */
std::vector<unsigned int> IdxGroup::getGroup(std::string grouping_name, unsigned int group_idx) {
  if (groupIdxRequire(grouping_name, group_idx)) {
    return m_groups[m_grouping_name_to_grouping_idx[grouping_name]][group_idx];
  }
  else {
    std::vector<unsigned int> empty;
    return empty;
  }
}

/**
 * @brief Get group at given index with given grouping name
 * @param groupname - name of grouping
 * @param group_idx - index of group in grouping vector
 */
// std::vector<unsigned int> IdxGroup::getGroup(std::string grouping_name, unsigned int group_idx) {

// }

/**
 * @brief Get group containing given index with given grouping name
 * @param groupname - name of grouping
 * @param idx - element inside group of interest
 */
std::vector<unsigned int> IdxGroup::getGroupContainingIdx(std::string grouping_name, unsigned int idx) {
  if (elementIdxRequire(grouping_name, idx)) { 
    unsigned int grouping_idx = m_grouping_name_to_grouping_idx[grouping_name];
    return m_groups[grouping_idx][m_element_idx_to_group_idx[grouping_idx][idx]];
  }
  else {
    std::vector<unsigned int> empty;
    return empty;
  }
}

/**
 * @brief Get other elements of a group containing given index with given grouping name
 * @param groupname - name of grouping
 * @param idx - element inside group of interest
 */
std::vector<unsigned int> IdxGroup::getOthersGroupContainingIdx(std::string grouping_name, unsigned int idx) {
  std::vector<unsigned int> rtn_vec;
  if (elementIdxRequire(grouping_name, idx)) {
    unsigned int grouping_idx = m_grouping_name_to_grouping_idx[grouping_name];
    rtn_vec = m_groups[grouping_idx][m_element_idx_to_group_idx[grouping_idx][idx]];
    std::vector<unsigned int>::iterator it = rtn_vec.begin();
    while (it != rtn_vec.end() && idx != *it) {
      it++;
    }
    if (it == rtn_vec.end()) {
      Err::errAbort("IdxGroup::getOthersGroupContainingIdx -- Data integrity error.  Cannot find expected index element: " + ToStr(idx) + " in grouping name: " + grouping_name + " at group index: " + ToStr(m_element_idx_to_group_idx[grouping_idx][idx]));
    }
    else {
      rtn_vec.erase(it);
    }
  }
  return rtn_vec;
}


/**
 * @brief Verify that the grouping name exists
 * @param groupname - name of grouping
 */
bool IdxGroup::groupingNameRequire(std::string grouping_name) {
  if (m_grouping_name_to_grouping_idx.count(grouping_name)) {
    return true;
  }
  else {
    Verbose::warn(0, "IdxGroup::groupingNameRequire -- Given group_name: " + grouping_name + " is not a valid groupname");
    return false;
  }
}

/**
 * @brief Verify that the given group index is valid for the given grouping name
 * @param groupname - name of grouping
 * @param group_idx - index of group in grouping vector
 */
bool IdxGroup::groupIdxRequire(std::string grouping_name, unsigned int group_idx) {
  if (groupingNameRequire(grouping_name)) {
    if (group_idx < m_groups[m_grouping_name_to_grouping_idx[grouping_name]].size()) {
      return true;
    }
    else {
      Verbose::warn(0, "IdxGroup::groupIdxRequire -- Given group_idx: " + ToStr(group_idx) + " is not a valid index for given grouping_name: " + grouping_name);
    }
  }
  return false;
}

/**
 * @brief Verify that the given group index is valid for the given grouping name
 * @param groupname - name of grouping
 * @param group_idx - index of group in grouping vector
 */
bool IdxGroup::elementIdxRequire(std::string grouping_name, unsigned int idx) {
  if (groupingNameRequire(grouping_name)) {
    unsigned int grouping_idx = m_grouping_name_to_grouping_idx[grouping_name];
    if (idx < m_element_idx_to_group_idx[grouping_idx].size()) {
      return true;
    }
    else {
      Verbose::warn(0, "IdxGroup::elementIdxRequire -- Given idx element: " + ToStr(idx) + " has not been entered for given grouping_name: " + grouping_name);
    }
  }
  return false;
}
