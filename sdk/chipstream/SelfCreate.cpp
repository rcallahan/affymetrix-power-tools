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
 * @file   SelfCreate.cpp
 * @author Chuck Sugnet
 * @date   Tue Nov  8 08:56:26 2005
 * 
 * @brief Small interface for functions that know how to make an instance of
 * themselves given a map of key, value parameters. Idea is that you can call
 * this class's static methods generically and dynamic_cast appropriately.
 */

//
#include "chipstream/SelfCreate.h"
//
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
//

using namespace std;
/** 
 * Chop up a string into a vector of words separated by '.' while trying to
 * respect floating point digits.
 * 
 * @param s - string of interest.
 * @param words - vector to put words into, will be cleared then filled.
 */
void SelfCreate::chopStringWithFloats(const std::string &s, std::vector<std::string> &words) {
  int len = 0, start = 0, next = 0, peek = 0, eqPos = 0;
  char delim = '.', pairDelim = '=';
  words.clear();
  len = s.length();
  bool bQuoted = false;
  while(start < len) {
    bool doFloat = true;
    int pos = 0;
    /* Look for the next '.' character and peek at the
       one after that to see if we are doing a float. */
	bQuoted = false;
    next = (int)s.find(delim, start);
    peek = (int)s.size();
    if(next == -1) { 
      next = (int)s.size(); 
    }
    else {
		peek = (int)s.find('\'', start+1);
		if (peek != -1) 
		{
			bQuoted = true;
			next = (int)s.find('\'', peek+1);
			if (next == -1) {next = (int)s.size();}
		}
		else {peek = (int)s.find(delim, next+1);}
      if(peek == -1)
	  {
        peek = (int)s.size();
	  }
    }
	if (bQuoted)
	{
		std::string str = s.substr(start, next - start);
		int iPosition = str.find('\'');
		str = str.substr(0, iPosition) + str.substr(iPosition+1);
		words.push_back(str);
		next++;
	}
	else
	{		
     
    /* The doFloat checks are to handle cases like 
       .myparam=2.0. */
    /* Check to see if we are all digits past the '=' separator... */
    eqPos = (int)s.find(pairDelim, start);
    if(eqPos == -1 || eqPos >= peek)
      doFloat = false;
    for(pos = eqPos+1; pos < next && doFloat; pos++) {
      if( !( isdigit(s[pos]) || s[pos] == '-' || s[pos] == '+' ) ) {
        doFloat = false;
        break;
      }
    }
    /* check to see that second half of float is also digits. */
    for(pos = next+1; pos < peek && doFloat; pos++) {
      if( !isdigit(s[pos]) ) {
        doFloat = false;
        break;
      }
    }
    if(doFloat)
      next = peek;
    words.push_back(s.substr(start, next - start));
    }
    start = next+1;
  }
}

/** 
 * @brief Take a string of from "name.param1=value1.param2=value2" and populate the
 * name and paramMap objects with it.
 * @param spec String description of name and parameters.
 * @param name Will contain name after call.
 * @param paramMap Parameter key/value pairs will be filled in.
 */
void SelfCreate::fillInNameParam(const std::string &spec, std::string &name, 
                                 std::map<std::string, std::string> &paramMap) {
  vector<string> words;
  string param;
  paramMap.clear();
  string::size_type position = spec.find(".");
  if(position == string::npos) {
    name = spec;
  }
  else {
    unsigned int i = 0;
    vector<string> params;
    name = spec.substr(0, position);
    chopStringWithFloats(spec.substr(position+1), params);
    for(i = 0; i < params.size(); i++) {
		
      string::size_type ePos = params[i].find("=");
      if(ePos == string::npos) {
		Verbose::out(1, "*");
        for(int j = 0; j < params.size(); j++) {
          Verbose::out(1, "Parsed parameter " + ToStr(j) + " is " + params[j]);
        }		
		Verbose::out(1, "*");
		Verbose::out(1, "Try surrounding the value with single quotes.");
		Verbose::out(1, "*");
        Err::errAbort("Must specify a key=value pair in parameter: " + spec);
      }
      paramMap[params[i].substr(0,ePos)] = params[i].substr(ePos+1);
    }
  }
}

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
bool SelfCreate::validParam(const std::map<std::string,std::string> &query,
                            SelfDoc &doc,
                            std::vector<std::string> &badParam) {
  std::map<std::string,std::string>::const_iterator iter;
  bool allOk = true;
  badParam.clear();
  for(iter = query.begin(); iter != query.end(); iter++) {
    if(!doc.validParam(iter->first)) {
      badParam.push_back(iter->first);
      allOk = false;
    }
  }
  return allOk;
}

bool SelfCreate::canMake(const std::string &spec, std::vector<SelfDoc> &docs) {
  string name;
  map<string, string> param;
  unsigned int i = 0;
  SelfCreate::fillInNameParam(spec, name, param);
  for (i = 0; i < docs.size(); i++) {
    if (docs[i].getDocName() == name) { 
      return true;
    }
  }
  return false;
}

/** 
 * @brief Create an object requested by the specification string by matching
 * it against the self documenting classes and then calling the matching self
 * creator function. Resulting object must be delete'd when done.
 * 
 * @param spec - String specification of object requested.
 * @param docs - Vector of selfdoc objects to look in for name of object to be made.
 * @param creators - Vector of creation functions to use for creating object.
 * @param typeDescription - What type of object are we attempting to make?
 * @param error - If true call errAbort when don't recognize value, otherwise return NULL.
 * @return Pointer to new SelfCreate object, delete when finished.
 */
SelfCreate *SelfCreate::selfCreateFromString(const std::string &spec, std::vector<SelfDoc> &docs,
                                             std::vector<SelfCreate::selfCreator> &creators, 
                                             const std::string &typeDescription, bool error) {
  string name;
  map<string, string> param;
  vector<string> badParam;
  unsigned int i = 0;
  SelfCreate *create = NULL;
  SelfCreate::fillInNameParam(spec, name, param);
  /* Loop through known types. */
  for(i = 0; i < docs.size(); i++) {
    if(docs[i].getDocName() == name) {
      /* Check the parameters. */
      if(SelfCreate::validParam(param, docs[i], badParam)) {
        create = (*creators[i])(param);
      }
      else { /* Bad param, tell user. */
        string bad;
        unsigned int badIx = 0;
        for(badIx = 0; badIx < badParam.size(); badIx++) {
          bad += badParam[badIx];
          bad += ",";
        }
        Err::errAbort("Invalid parameters for " + name + ": " + bad);
      }
      return create;
    }
  }
  if(error)
    Err::errAbort("Don't know how to make object for name: " + name + " of type: " + typeDescription);
  return NULL;
}
