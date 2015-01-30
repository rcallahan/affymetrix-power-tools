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

#ifndef _CNCytoEngine_H_
#define _CNCytoEngine_H_
/**
 * @file CNCytoEngine.h
 *
 * @brief This header contains the Log2RatioEngine class definition.
 */

#include "copynumber/CNLog2RatioData.h"
//
#include "chipstream/apt-geno-qc/GenoQC.h"
#include "util/BaseEngine.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include<string>



using namespace std;

/**
 * @brief  The engine for controlling the copynumber Cytos.
 */
class CNCytoEngine : public BaseEngine
{
private:

    /// List of aliases for standardized methods
    std::map<std::string,std::string>* m_pstdMethods;
    std::vector<string> m_vCNReferenceHeader;
    void checkOptionsImp();
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
        void checkParameters();
        virtual std::string getEngineName() { return CNCytoEngine::EngineName(); }
        static const std::string EngineName() { return "CNCytoEngine"; }
        /**
         * Constructor
         */
        CNCytoEngine();

        /**
         * Destructor
         */
        ~CNCytoEngine();

        /*! A class to register the engine. */
        class Reg : public EngineReg
        {
        public:
            /*! Constructor - register the engine. */
            Reg() : EngineReg(CNCytoEngine::EngineName())
            {
            }

            /*! Creates an object.
             * @return The object.
             */
            BaseEngine *MakeObject() { return new CNCytoEngine; }
        };

        /*! The one and only registration object. */
        static Reg reg;

        /*! Converts the type to the engine type.
         * @param chip The pointer to the base engine object.
         * @return The engine type or NULL if not compatible.
         */
        static CNCytoEngine * FromBase(BaseEngine *engine);

    void addCel(const std::string& strCelFileName, const std::string& strSampleFileName = "", const std::string& strResultFileName = "");


    void loadReferenceHeader(std::string strReferenceFileName);

    AffxString getAnnotationFileName(const AffxString& strFileName);
    // This method is called by ChAS without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    bool getNormalizationTypes(const AffxString& strFileName, int& iNormalizationType, bool& bAdapterNormlization, int& iSampleCount);
    bool getReferenceGenderCounts(const AffxString& strFileName, int& iReferenceMaleCount, int& iReferenceFemaleCount, int& iReferenceUnknowCount);

    bool isDualNormalizeReference(const std::string& strReferenceFileName);

    static bool isCN5Reference(const AffxString& strReferenceFileName);
    // This method is called by ChAS without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static bool isCyto2Reference(const AffxString& strReferenceFileName);

    static bool isValidReferenceFile(const AffxString& strReferenceFileName);
    // This method is called by ChAS without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationParameter(const AffxString& strFileName, const AffxString& strParamterName);
    // This method is called by ChAS without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationArrayType(const AffxString& strFileName) {return getAnnotationParameter(strFileName, "affymetrix-array-type");}
    // This method is called by ChAS without having called checkOptions.
    // In short, best not to refer to Engine options/state in this method
    static AffxString getAnnotationVersion(const AffxString& strFileName) {return getAnnotationParameter(strFileName, "affymetrix-algorithm-param-netaffx-build");}

protected:
    void createReference();
    void useReference(bool bAnalyis);
    void defineStdMethods();

    void callGenoQCEngine();
    void calculateGenders(CNLog2RatioData& data);

    AffxString runReferencePart1(CNLog2RatioData& data);
    void runReferencePart2(CNLog2RatioData& data, const AffxString& strTempFileName);

    bool writeOutputText(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets);
    bool writeOutputHdf5(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets);

    void checkChipType();
    bool isCyto2SnpReference(const AffxString& strReferenceFileName);


    bool predictedIntensitiesExist(const std::string strInputReferenceFileName);
    bool checkSnpReferenceInputFile(const std::string strSnpReferenceFileName);
    bool checkGenotypeCallOverrideFile(    const std::string &strGenotypeCallOverrideFile,
                                        const std::vector<std::string> &celFiles,
                                        const std::string  &strSnpReferenceInputFileName );
    bool checkGenderOverrideFile(       const std::string &strGenderOverrideFileName,
                                        const std::vector<std::string> &celFiles);

    void checkConsistency(const vector<string>& analysis);
	
	void loadProbes(CNLog2RatioData& data, const string& chrXProbeFile, const string& chrYProbeFile,
                  std::vector< std::vector<probeid_t> >& chrXProbes, 
                  std::vector< std::vector<probeid_t> >& chrYProbes);

    std::string saveProbes(const std::string& basename, std::vector< std::vector<probeid_t> >& chrProbes);

    void removeTempProbeFiles();
};


class StringLessThan
{
    public:
       bool operator()(const std::string &string1, const std::string &string2)
                      {return string1 < string2;}
};


#endif


