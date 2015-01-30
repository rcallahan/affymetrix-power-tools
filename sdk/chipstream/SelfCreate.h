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
 * @file   SelfCreate.h
 * @author Chuck Sugnet
 * @date   Tue Oct 25 17:10:41 2005
 * 
 * @brief Small interface for functions that know how to make an instance of
 * themselves given a map of key, value parameters. Idea is that you can call
 * this class's static methods generically and cast appropriately.
 * 
 */
#ifndef _SELFCREATE_H_
#define _SELFCREATE_H_

//
#include "chipstream/SelfDoc.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

/**
 * Small interface for functions that know how to make an instance of themselves
 * given a map of key, value parameters. Idea is that you can call this class's
 * static methods generically and cast appropriately.
 */
class SelfCreate { 

public:

  /**
   * @brief virtual destructors for virtual class.
   */
  virtual ~SelfCreate() {};

  /** 
   * @brief Take a string of from "name.param1=value1.param2=value2" and populate the
   * name and paramMap objects with it.
   * @param spec String description of name and parameters.
   * @param name Will contain name after call.
   * @param paramMap Parameter key/value pairs will be filled in.
   */
  static void fillInNameParam(const std::string &spec, std::string &name, 
                              std::map<std::string, std::string> &paramMap);

  /** 
   * @brief Check the proposed parameters in query against the allowed
   * in gold. Specifically checks to make sure all of the keys in
   * allowed are keys in gold.
   * 
   * @param query Proposed key/value pairs for parameters.
   * @param gold Allowed key/description pairs for parameters.
   * @param badParam Vector filled with keys for param found to be bad.
   * 
   * @return true if all keys in query are keys in gold, false otherwise
   */
  static bool validParam(const std::map<std::string,std::string> &query,
                         SelfDoc &doc,
                         std::vector<std::string> &badParam);

  /// Makes for some nice syntactic sugar. 
  typedef SelfCreate *(*selfCreator)(std::map<std::string,std::string> &param);

  /** 
   * @brief This static function should be overridden by child classes
   * to return an object of the correct type initialized correctly
   * with the parameters in the string, string map. All objects
   * created this way should be deleted when finished using.
   * 
   * @param param - Map of key/value pairs to initialize the object.
   * 
   * @return Pointer toCreate object, this should be sub casted as necessary.
   */
  static SelfCreate *newObject(std::map<std::string,std::string> &param) {
#ifdef _MSC_VER
    (param); /* unused. */
#endif
    return NULL;
  }

  /** 
   * Is there a match in this SelfDoc vector for a the specification provided? In other
   * words should this set of SelfDoc objects be trying to construct this specification.
   * 
   * @param spec - String specification of object requested.
   * @param docs - Vector of selfdoc objects to look in for name of object to be made.
   * 
   * @return true if there is a valid match for this specification in these docs, false otherwise
   */
  static bool canMake(const std::string &spec, std::vector<SelfDoc> &docs);

  /** 
   * @brief Create an object requested by the specification string by matching
   * it against the self documenting classes and then calling the matching self
   * creator function. Resulting object must be delete'd when done.
   * @param spec - String specification of object requested.
   * @param docs - Vector of selfdoc objects to look in for name of object to be made.
   * @param creators - Vector of creation functions to use for creating object.
   * @param typeDescription - What type of object are we attempting to make?
   * @param error - If true call errAbort when don't recognize value, otherwise return NULL.
   * @return Pointer to new SelfCreate object, delete when finished.
   */
  static SelfCreate *selfCreateFromString(const std::string &spec, std::vector<SelfDoc> &docs,
                                          std::vector<SelfCreate::selfCreator> &creators, 
                                          const std::string &typeDescription, bool error=true);
  
  /** 
   * Chop up a string into a vector of words separated by '.' while trying to
   * respect floating point digits.
   * 
   * @param s - string of interest.
   * @param words - vector to put words into, will be cleared then filled.
   */
  static void chopStringWithFloats(const std::string &s, std::vector<std::string> &words);

  static void setValue(std::string name, std::map<std::string,std::string> &param, SelfDoc &doc) {
    std::map<std::string,std::string>::iterator iter;
    SelfDoc::Opt &opt = doc.getDocOption(name);
    iter = param.find(name);
    if(iter != param.end()) {
      opt.value = iter->second;
    }
	  if (!opt.type == SelfDoc::Opt::String)
	  {
		  // Strip surrounding single quotes.
		  if ((opt.value.length() > 2) && (opt.value[0] == '\'') && (opt.value[opt.value.length() - 1] == '\''))
		  {
			  opt.value.assign(opt.value.substr(1, opt.value.length() - 2));
		  }
	  }
    if(!opt.checkOption()) {
      Err::errAbort("SelfCreate::setValue() - '" + opt.value + "' is not a valid value for parameter: '" + name + "'. The specified range is " + opt.minVal + " to " + opt.maxVal);
    }
  }

  static void fillInValue(int &value, std::string name, std::map<std::string,std::string> &param, SelfDoc &doc) {
    setValue(name, param, doc);
    SelfDoc::Opt &opt = doc.getDocOption(name);
    value = opt.asInt();
  }

  static void fillInValue(double &value, std::string name, std::map<std::string,std::string> &param, SelfDoc &doc) {
    setValue(name, param, doc);
    SelfDoc::Opt &opt = doc.getDocOption(name);
    value = opt.asDouble();
  }

  static void fillInValue(float &value, std::string name, std::map<std::string,std::string> &param, SelfDoc &doc) {
    setValue(name, param, doc);
    SelfDoc::Opt &opt = doc.getDocOption(name);
    value = opt.asFloat();
  }

  static void fillInValue(bool &value, std::string name, std::map<std::string,std::string> &param, SelfDoc &doc) {
    setValue(name, param, doc);
    SelfDoc::Opt &opt = doc.getDocOption(name);
    value = opt.asBool();
  }

  static void fillInValue(std::string &value, std::string name, std::map<std::string,std::string> &param, SelfDoc &doc) {
    setValue(name, param, doc);
    SelfDoc::Opt &opt = doc.getDocOption(name);
    value = opt.asString();
  }
    
};

#endif /* _SELFCREATE_H_ */
