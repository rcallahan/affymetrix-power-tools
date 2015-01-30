////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
* @file   SnpModelConverterEngine.h
* @author Mybrid Spalding
*
*/
#ifndef _SNPMODELCONVERTER_ENGINE_H_
#define _SNPMODELCONVERTER_ENGINE_H_

//
#include "util/BaseEngine.h"
//
#include <string>
#include <vector>
//

class SnpModelConverterEngine : public BaseEngine
{

public:

  virtual std::string getEngineName() {
    return SnpModelConverterEngine::EngineName();
  }
  static const std::string EngineName() {
    return "SnpModelConverterEngine";
  }

  /**
   * Constructor
   */
  SnpModelConverterEngine();

  /**
   * Destructor
   */
  virtual ~SnpModelConverterEngine();

  /*! A class to register the summary engine. */
  class Reg : public EngineReg
  {
  public:
    /*! Constructor - register the summary engine. */
    Reg() : EngineReg(SnpModelConverterEngine::EngineName()) {
    }

    /*! Creates an object.
     * @return The object.
     */
    BaseEngine *MakeObject() {
      return new SnpModelConverterEngine;
    }
  };

  /*! The one and only registration object. */
  static Reg reg;

  /*! Converts the type to the engine type.
   * @param chip The pointer to the base engine object.
   * @return The summary engine type or NULL if not compatible.
   */
  static SnpModelConverterEngine * FromBase(BaseEngine *engine);

private:
  std::vector< std::string > m_modelFiles;
  void defineOptions(SnpModelConverterEngine *engine);
  void checkOptionsImp();
  void runImp();

  void runImpDumpHeaders(const std::map< std::string, std::string *> & headers);
  void runImpConvertToDb(const std::string & modelFilePath );
  
};

#endif /* _SNPMODELCONVERTER_ENGINE_H_ */

