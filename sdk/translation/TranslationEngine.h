////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
/*
 * @author Alan Williams
 * @date   Mon Jun 23 14:57:34 PDT 2008
 *
 * @brief Factory design pattern, the TranslationEngine for DMET 3.0 that implements BaseEngine. 
 */
#ifndef _TRANSLATINE_ENGINE_H_
#define _TRANSLATION_ENGINE_H_

//
#include "translation/RunTimeEnvironment.h"
//
#include "util/BaseEngine.h"
//
#include <cstring>
#include <string>
#include <vector>
//


/**
 * @brief Class to run the DMET analysis
 */
class TranslationEngine  : public BaseEngine
{

public:

	virtual std::string getEngineName() { return TranslationEngine::EngineName(); }
	static const std::string EngineName() { return "TranslationEngine"; }

	/*! A class to register the summary engine. */
	class Reg : public EngineReg
	{
	public:
		/*! Constructor - register the summary engine. */
		Reg() : EngineReg(TranslationEngine::EngineName())
		{
		}

		/*! Creates an object.
		 * @return The object.
		 */
		BaseEngine *MakeObject() { return new TranslationEngine; }
	};

	/*! The one and only registration object. */
	static Reg reg;

	/*! Converts the type to the summary engine type.
	 * @param chip The pointer to the base engine object.
	 * @return The summary engine type or NULL if not compatible.
	 */
	static TranslationEngine * FromBase(BaseEngine *engine);



  unsigned char      m_controllerMask;
  RunTimeEnvironment m_rte;

  /**
   * Constructor
   */
  TranslationEngine(char * argv[] = NULL,
                    unsigned char controllerMask = C_CONSOLE,
                    void (*messageHandlerCallback)(ADTOptions &adtOpts)
                    = NULL);

  /**
   * Destructor
   */
  ~TranslationEngine() {};


  /**
   * Do the analysis
   * Will call Verbose::out() and Err::errAbort().
   */
  void runImp();

  void checkOptionsImp();


private:

  void (*m_messageHandlerCallback)(class ADTOptions & adtOpts);

  void _setADTOptions(class ADTOptions & adtOpts,
                      void (*messageHandlerCallback)(ADTOptions *adtOpts)
                      = NULL);
  
};

#endif /* _TRANSLATION_ENGINE_H_ */

