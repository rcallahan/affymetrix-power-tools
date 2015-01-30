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

#ifndef _CNWorkflowEngine_H_
#define _CNWorkflowEngine_H_
/**
 * @file CNWorkflowEngine.h
 *
 * @brief This header contains the Log2RatioEngine class definition.
 */

#include "chipstream/apt-geno-qc/GenoQC.h"
#include "chipstream/apt-probeset-genotype/ProbesetGenotypeEngine.h"
#include "copynumber/CNLog2RatioData.h"
#include "util/BaseEngine.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
#include "util/Verbose.h"
//

using namespace std;

/**
 * @brief  The engine for controlling the copynumber workflows.
 */
class CNWorkflowEngine : public BaseEngine
{
private:
    /// List of aliases for standardized methods
    std::map<std::string,std::string>* m_pstdMethods;
    ProbesetGenotypeEngine *m_apgEngine;
    void checkOptionsImp();
    void checkDiskSpaceImp();
    void runImp();
    void defineOptions();
    void defineStates();
    void printStandardMethods(std::ostream &out);
    void extraHelp();
    void explain();

    /**
     * Print out a warning message
     */
    void warning(const std::string& strMessage) {
        Verbose::out(1, "WARNING: " + strMessage);
    }

    /**
     * Error out by calling clear() and Err::errAbort()
     */
    void error(const std::string& strMessage) {
        clear();
        Err::errAbort(strMessage);
    }

public:
        void clear();
        virtual std::string getEngineName() { return CNWorkflowEngine::EngineName(); }
        static const std::string EngineName() { return "CNWorkflowEngine"; }
        /**
         * Constructor
         */
        CNWorkflowEngine();

        /**
         * Destructor
         */
        ~CNWorkflowEngine();

        /*! A class to register the engine. */
        class Reg : public EngineReg
        {
        public:
            /*! Constructor - register the engine. */
            Reg() : EngineReg(CNWorkflowEngine::EngineName())
            {
            }

            /*! Creates an object.
             * @return The object.
             */
            BaseEngine *MakeObject() { return new CNWorkflowEngine; }
        };

        /*! The one and only registration object. */
        static Reg reg;

        /*! Converts the type to the engine type.
         * @param chip The pointer to the base engine object.
         * @return The engine type or NULL if not compatible.
         */
        static CNWorkflowEngine * FromBase(BaseEngine *engine);



    void addCel(const std::string& strCelFileName, const std::string& strSampleFileName = "", const std::string& strResultFileName = "");

    bool getAnalysisName(const AffxString& strFileName, AffxString& strAnalysisName);
    bool getAnnotationFileName(const AffxString& strFileName,AffxString& strAnnotationFileName);
    // This method is called by GTC without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    bool getNormalizationTypes(const AffxString& strFileName, int& iNormalizationType, bool& bAdapterNormlization, int& iSampleCount);
    bool getReferenceGenderCounts(const AffxString& strFileName, int& iReferenceMaleCount, int& iReferenceFemaleCount, int& iReferenceUnknowCount);

    // This method is called by GTC without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationParameter(const AffxString& strFileName, const AffxString& strParamterName);
    // This method is called by GTC without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationDBFileName(const AffxString& strFileName) {return Fs::basename(getAnnotationParameter(strFileName, "affymetrix-algorithm-param-apt-opt-annotation-file"));}
    // This method is called by GTC without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationSNPFileName(const AffxString& strFileName) {return Fs::basename(getAnnotationParameter(strFileName, "affymetrix-algorithm-param-apt-opt-netaffx-snp-annotation-file"));}
    // This method is called by GTC without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationCNFileName(const AffxString& strFileName) {return Fs::basename(getAnnotationParameter(strFileName, "affymetrix-algorithm-param-apt-opt-netaffx-cn-annotation-file"));}

    // This method is called by GTC without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    bool isCN5Reference(const AffxString& strReferenceFileName);
    bool isCyto2Reference(const AffxString& strReferenceFileName);

protected:
    void defineStdMethods();
    void reportBasics();

    void setupParameters();
    void setupReferenceRunParameters();
    void setupCopyNumberRunParameters();
    void callGenoQCEngine();
    void callProbesetGenotypeEngine();
    void callCNReferenceEngine();
    void callCNLog2RatioEngine();
    void buildCNAnalysisStrings();
    void calculateGenders();
    void checkChipType();
};


#endif


