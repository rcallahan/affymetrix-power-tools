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
 * @file   SelfDoc.h
 * @author Chuck Sugnet
 * @date   Tue Oct 25 11:57:11 2005
 * 
 * @brief Small interface to for algorithmic classes that can explain
 * themselves.
 */
#ifndef _SELFDOC_H_
#define _SELFDOC_H_

//
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <cstring>
#include <map>
#include <string>
//
/** 
 * Small interface to for algorithmic classes that can explain
 * themselves.
 */
class SelfDoc {

public:

  /**
   * Description of one possible option/parameter for a self describing class.
   */ 
  class Opt {
  public: 
    /** Types of options supported. */
    enum OptType {
      String,
      Double,
      Float,
      Integer,
      Boolean
    };

    std::string name;           ///< Name of option.
    enum SelfDoc::Opt::OptType type; ///< What data type is our option.
    std::string value;          ///< Current value for option.
    std::string defaultVal;     ///< Default value for option.
    std::string minVal;         ///< Minimum value acceptable for numeric types, or NA if not used.
    std::string maxVal;         ///< Maximum value acceptable for numeric types, or NA if not used.
    std::string descript;       ///< Short description of our option and what it is used for.
    
    /** 
     * Check to see that an option's value is valid.
     * @return true if option is valid - false otherwise.
     */

    bool checkOption() const;
    /* Conversion utilities. */
    double asDouble() const;
    float asFloat() const;
    int asInt() const;
    bool asBool() const;
    std::string asString() const;

  };

  /// Typedef to make syntactic sugar nicer.
  typedef SelfDoc *explainer ();

  /**
   * @brief virtual destructors for virtual class.
   */
  virtual ~SelfDoc() {};

  /** 
   * @brief Setter method for name appearing in doc.
   * @param s - short name or identifier.
   */
  void setDocName(const std::string &s) { m_DocName =  s; }

  /** 
   * @brief Getter method for documentation name.
   * @return - Short name or identifier
   */
  std::string getDocName() { return m_DocName;}

  /** 
   * @brief Setter method for description appearing in doc.
   * @param s - description of why/what/how we do.
   */
  void setDocDescription(const std::string &s) { m_DocDescription = s; }

  /** 
   * @brief Getter method for documentation description.
   * @return - description of why/what/how we do.
   */
  std::string getDocDescription() { return m_DocDescription; }

  /** 
   * Get a string representation of the state of options in this object.
   * @return - string representation of the kind used by SelfCreate to
   * construct objects. Looks like "name.opt1=val1.opt2=val2"
   */
  virtual std::string getState();

  /** 
   * @brief Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  std::vector<SelfDoc::Opt> getDocOptions() { return m_Options; }

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt> dummy;
    return dummy;
  }
  
  /** 
   * Determine if the name provided is the identifier for a valid option.
   * @param name - possible identifier.
   * @return true if valid option identifier, false otherwise.
   */
  bool validParam(const std::string &name) {
    return m_OptionsIndex.find(name) != m_OptionsIndex.end();
  }

  /** 
   * @brief Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  SelfDoc::Opt & getDocOption(const std::string &name);

  /** 
   * Look up option with name and set it to value.
   * @param name - identifier for option.
   * @param value - string to set value to.
   */
  void setOptValue(const std::string &name, const std::string &value);

  /** 
   * Specialized version for booleans as ToStr() appears to print them
   * as '1' and '0' instead of 'true' and 'false'.
   * @param name - identifier for option.
   * @param value - boolean to set value to.
   */
  void setOptValue(const std::string &name, const bool &value);

  /** 
   * @brief Setter method for parameters and their documentation.
   * @param param - Map of parameter names to descriptions.
   */
  void setDocOptions(const std::vector<SelfDoc::Opt> &param);

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf();

  /** 
   * Write out a self documenter name, description, and all of its parameters.
   * 
   * @param doc - Self documenting class.
   * @param out - Stream to write the description to.
   */
  static void printExplanation(SelfDoc &doc, std::ostream &out);

protected:
  /// Name of class or algorithm.
  std::string m_DocName;
  /// Description of why/what/how we do.
  std::string m_DocDescription;
  /// Possible parameters.
  std::vector<SelfDoc::Opt> m_Options;
  /// Possible parameters map.
  std::map<std::string,int> m_OptionsIndex;

};


#endif /* _SELFDOC_H_ */
