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

#ifndef _IDXGROUP_H_
#define _IDXGROUP_H_

#include <cstring>
#include <map>
#include <string>
#include <vector>
//

class IdxGroup {

  //  typedef unsigned int element_idx;
  
public:
  
  /** 
   * @brief Constructor
   */
  IdxGroup() : m_grouping_count(0) {}
  
  /**
   * @brief Insert group of indexes into grouping of given name
   * @param groupname - name/category of grouping
   * @param group - integers that belong to a group
   */
  unsigned int addGroup(std::string grouping_name, std::vector<unsigned int> group);
  
  /**
   * @brief Get count of groups with given grouping name
   * @param groupname - name/category of grouping
   */
  int getGroupCount(std::string grouping_name);
  
  /**
   * @brief Get vector of vectors constituting the groups of indexes for given grouping name
   * @param groupname - name/category of grouping
   */
  std::vector<std::vector<unsigned int> > getGroupingVec(std::string grouping_name);// {return return m_groups[m_grouping_name_to_grouping_idx[grouping_name]];}
  
  /**
   * @brief Get group at given index with given grouping name
   * @param groupname - name/category of grouping
   * @param group_idx - index of group in grouping vector
   */
  std::vector<unsigned int> getGroup(std::string grouping_name, unsigned int group_idx);// {return m_groups[m_grouping_name_to_grouping_idx[grouping_name]][group_idx];}
  
  /**
   * @brief Get group containing given index with given grouping name
   * @param groupname - name/category of grouping
   * @param idx - element inside group of interest
   */
  std::vector<unsigned int> getGroupContainingIdx(std::string grouping_name, unsigned int idx);  
  
  /**
   * @brief Get other elements of a group containing given index with given grouping name
   * @param groupname - name/category of grouping
   * @param idx - element inside group of interest
   */
  std::vector<unsigned int> getOthersGroupContainingIdx(std::string grouping_name, unsigned int idx);
  
  /**
   * @brief Verify that the grouping name exists
   */
  bool groupingNameExists(std::string grouping_name) {return m_grouping_name_to_grouping_idx.count(grouping_name)!=0;} 
  
  /**
   * @brief Verify that the grouping name exists
   * @param groupname - name of grouping
   */
  bool groupingNameRequire(std::string grouping_name);  
  
  
  /**
   * @brief Verify that the given group index is valid for the given grouping name
   * @param groupname - name of grouping
   * @param group_idx - index of group in grouping vector
   */
  bool groupIdxRequire(std::string grouping_name, unsigned int group_idx);
  
  /**
   * @brief Verify that the given index element has been entered for the given grouping name
   * @param groupname - name of grouping
   * @param idx - index that is an element of some group
   */
  bool elementIdxRequire(std::string grouping_name, unsigned int idx);
  
private:
  /// map of groupname to group index
  std::map<std::string, unsigned int> m_grouping_name_to_grouping_idx;
  /// vector of indexes into group
  // group_idx = m_element_idx[grouping_idx][element_idx]
  std::vector<std::vector<unsigned int> > m_element_idx_to_group_idx;
  /// vector of groups
  std::vector<std::vector<std::vector<unsigned int> > > m_groups;
  /// grouping count
  unsigned int m_grouping_count;
};

#endif /* _IDXGROUP_H_ */
