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
 * @file CNGenderEngine.cpp
 *
 * @brief This file contains the CNGenderEngine class members.
 */
#include "chipstream/BioTypes.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/CnProbeGenderCelListener.h"
#include "copynumber/CNGenderEngine.h"
//
#include "util/Fs.h"
#include "util/Verbose.h"

using namespace std;

 /*
 * @brief Constructor
 */
CNGenderEngine::CNGenderEngine()
{
    defineOptions();
}

/**
 * @brief Destructor
 */
CNGenderEngine::~CNGenderEngine()
{

}

void CNGenderEngine::defineOptions() 
{
    defineOptionSection("Input Options");
    defineOption("", "chrX-probes", PgOpt::STRING_OPT, "File containing probes on X chromosome for gender computation.", "");
    defineOption("", "chrY-probes", PgOpt::STRING_OPT, "File containing probes on Y chromosome for gender computation.", "");
    defineOption("", "cel-files", PgOpt::STRING_OPT, "Text file specifying cel files to process, one per line with the " "first line being 'cel_files'.", "");
    defineOption("", "male-gender-ratio-cutoff", PgOpt::DOUBLE_OPT, "Male gender ratio cutoff", "1.3"); 
    defineOption("", "female-gender-ratio-cutoff", PgOpt::DOUBLE_OPT, "Female gender ratio cutoff", "1.0");


  defineOption("", "check-for-invalid-cel-files", PgOpt::BOOL_OPT,
                     "Does an upfront check on the intensities in the input cel files and warns if there are any peculiar values .",
                     "false");


    defineOptionSection("Engine Options (Not used on command line)");
    defOptMult(	"", "cels", PgOpt::STRING_OPT, "CEL files to process.", "");

}


/**
 * @brief Make sure that our options are sane. Call Err::errAbort if not.
 */
void CNGenderEngine::checkOptionsImp() 
{
    std::vector<std::string> celFiles;
    EngineUtil::getCelFiles(celFiles, this);
    if(celFiles.size() == 0)
        Err::errAbort("No cel files specified.");
    setOpt("cels",celFiles);

}

bool CNGenderEngine::checkForInvalidCelFiles(const std::vector<std::string>& celFiles)
{
    bool returnValue=true;
    affymetrix_fusion_io::FusionCELData cel;
    for (unsigned int uiCelFileIndex = 0; (uiCelFileIndex < celFiles.size()); uiCelFileIndex++)
    {
cout <<  celFiles[uiCelFileIndex].c_str() << endl; 
        cel.SetFileName(celFiles[uiCelFileIndex].c_str());
        if (!cel.Read()) {Err::errAbort("In CNGenderEngine cannot read cel file for gender calling: " + celFiles[uiCelFileIndex]);}

        unsigned int uiRowCount = cel.GetRows();
        unsigned int uiColCount = cel.GetCols();
        unsigned int uiCount = uiRowCount * uiColCount;
        int maxWarnings = 3;
        for (int iProbeIndex = 0; iProbeIndex < uiCount; iProbeIndex++)
        {
            float fIntensity = cel.GetIntensity(iProbeIndex);
            if (fIntensity == 0.0)
            {
                 if(maxWarnings > 0)
                 {
                     Verbose::out(1,"CEL file contains zero intensity values:\t" + celFiles[uiCelFileIndex]);
                     maxWarnings--;
                     returnValue = false;
                 }
                 continue;
            }
            else if (fIntensity < 0)
            {
                Verbose::out(1,"CEL file contains negative intensity values:\t" + celFiles[uiCelFileIndex]);
                returnValue = false;
            }
            else if ((double)fIntensity != (double)fIntensity)
            {
                Verbose::out(1, "CEL file contains NaN intensity values:\t" + celFiles[uiCelFileIndex]);
                returnValue = false;
            }
            else if (!Util::isFinite(fIntensity))
            {
                Verbose::out(1,"CEL file contains infinite intensity values:\t" + celFiles[uiCelFileIndex]);
                returnValue = false;
            }
        }
    } 
    return returnValue;
}

/**
 * @brief Run the specified analysis
 */
void CNGenderEngine::runImp() {

        std::vector< std::vector<probeid_t> > chrXProbes;
        std::vector< std::vector<probeid_t> > chrYProbes;

        EngineUtil::readProbeFile(chrXProbes,getOpt("chrX-probes")); 
        EngineUtil::readProbeFile(chrYProbes,getOpt("chrY-probes")); 

        std::vector<std::string> celFiles = getOptVector("cels");

        if(getOptBool("check-for-invalid-cel-files"))
        {
            checkForInvalidCelFiles(celFiles); 
        }


        CnProbeGenderCelListener objGenderCaller(	chrXProbes, 
                                                        chrYProbes, 
                                                        getOptDouble("female-gender-ratio-cutoff"), 
                                                        getOptDouble("male-gender-ratio-cutoff"));

        affymetrix_fusion_io::FusionCELData cel;

        affx::TsvFile *tsv = new affx::TsvFile;
        tsv->writeTsv(Fs::join(getOpt("out-dir"),"genderCalls.txt"));
        tsv->defineColumn(    0,0,"Sample Name", affx::TSV_TYPE_INT);
        tsv->defineColumn(    0,1,"GenderCall", affx::TSV_TYPE_DOUBLE);
        tsv->defineColumn(    0,2,"Intensity Ratio", affx::TSV_TYPE_STRING);

        for (unsigned int uiCelFileIndex = 0; uiCelFileIndex < celFiles.size(); uiCelFileIndex++)
        {
            cel.SetFileName(celFiles[uiCelFileIndex].c_str());
            if (!cel.Read()) {error("Cannot read cel file for gender calling: " + celFiles[uiCelFileIndex]);}

            objGenderCaller.newChip(&cel); 
            int iMaleCount =0;
            int iFemaleCount =0;
            int iUnknownCount =0;

            float fRatio = objGenderCaller.getLastRatio();

            tsv->set(0,0,Fs::basename(celFiles[uiCelFileIndex]));
            tsv->set(0,2,fRatio);
            switch(objGenderCaller.getLastGender())
            {
                case affx::Male: 
                    iMaleCount++; 
                    Verbose::out(1, "Sample " + ::getUnsignedInt(uiCelFileIndex + 1) 
                                              + ":\t" 
                                              + Fs::basename(celFiles[uiCelFileIndex]) 
                                              + "\tMale\t" 
                                              + ::getDouble(fRatio, 6)); 
                    tsv->set(0,1,"Male");
                    break;


                case affx::Female: 
                    iFemaleCount++; 
                    Verbose::out(1, "Sample " + ::getUnsignedInt(uiCelFileIndex + 1) 
                                              + ":\t" 
                                              + Fs::basename(celFiles[uiCelFileIndex]) 
                                              + "\tFemale\t" 
                                              + ::getDouble(fRatio, 6)); 
                    tsv->set(0,1,"Female");
                    break;


                default: 
                    iUnknownCount++; 
                    Verbose::out(1, "Sample " + ::getUnsignedInt(uiCelFileIndex + 1) 
                                              + ":\t" 
                                              + Fs::basename(celFiles[uiCelFileIndex]) 
                                              + "\tUnknown\t" 
                                              + ::getDouble(fRatio, 6));
                    tsv->set(0,1,"UnknownGender");
            }
            tsv->writeLevel(0);
        }
        delete tsv;
}


