////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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

#ifndef _ANNOTATIONCONVERTERENGINE_H_
#define _ANNOTATIONCONVERTERENGINE_H_

//
#include "util/BaseEngine.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

#include "sqlite3.h"

#include "Configuration.h"
#include "AnnotationTable.h"
#include "LocalizationTable.h"
#include "ChromosomeTable.h"
#include "CdfInformationTable.h"

#include "CsvHeadersToInformationTable.h"

/**
 */
class AnnotationConverterEngine : public BaseEngine {
    public:
			
		virtual std::string getEngineName() { return AnnotationConverterEngine::EngineName(); }
		static const std::string EngineName() { return "AnnotationConverterEngine"; }

		/**
		* Constructor
		*/
		AnnotationConverterEngine();
		

		/**
		* Destructor
		*/
		~AnnotationConverterEngine();
		

		/*! A class to register the engine. */
		class Reg : public EngineReg
		{
		public:
			/*! Constructor - register the summary engine. */
			Reg() : EngineReg(AnnotationConverterEngine::EngineName())
			{
			}

			/*! Creates an object.
			 * @return The object.
			 */
			BaseEngine *MakeObject() { return new AnnotationConverterEngine; }
		};

		/*! The one and only registration object. */
		static Reg reg;

		/*
		 */
		static AnnotationConverterEngine * FromBase(BaseEngine *engine);



        /** 
         * @brief Compute the stats given the input options
         */
       //$$ void ComputeAnnotationConverter();

       
        string programPath;
    private:
        void checkOptionsImp();
        void runImp();
        void defineOptions();
        void defineStates();
		void CheckChipType();
        void clear();
        void checkDiskSpaceImp();
        
        // My Methods
        void CheckSingleRequiredSwitch(const char* SwitchName, std::string& SwitchValue);
        void LoadConfiguration(string FileName); 

        void AnalyzeTsvFiles(); 
            void AnalyzeHeaders();                                 
        
        sqlite3* OpenDb(std::string dbFileName);
        void CloseDbAndRemoveFile();
        void BeginTransaction();
        void CommitTransaction();
        void RollbackTransaction();
        void ExecuteNonQuery(std::string SqlStatement);
                
        void GetTsvFiles();
        // std::string GetErrorMsg();
        
        // Members
        sqlite3* m_annotationDb;        
        Configuration m_configuration; 
        CsvHeadersToInformationTable* m_csvHeadersToInformationTablePtr;
        AnnotationTable* m_annotationTable;
        LocalizationTable* m_localizationTable;
        ChromosomeTable* m_chromosomeTable;                
        CdfInformationTable* m_cdfInformationTable;
        
        static void ValueMismatchReport(string VarName, map<string, string> FileValueMap);
        static void MissingValueReport(string VarName, vector<string> Files);
        static vector<string> m_valueMismatchReportContent;
        static vector<string> m_missingValueReportContent;
        
        //Options
        std::string m_dbFileName;
        std::string m_dbTemplate;
        std::string m_tsvFileName;
        std::string m_tsvListFileName;
        std::string m_arraySpecificConfig; 
		std::string m_arraySetName;
        ///////////////
        
        std::vector<std::string> m_listOfTsvFiles;
           
};

#endif // _SNPSUMMARY_H_
