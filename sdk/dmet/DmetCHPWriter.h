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

#include "calvin_files/data/src/CHPMultiDataData.h"
#include "calvin_files/utils/src/GenoCallCoder.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
#include "chipstream/AnalysisInfo.h"
#include "chipstream/ChipLayout.h"
#include "file5/File5_File.h"
#include "file5/File5_Tsv.h"
#include "util/BaseEngine.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>

#define DMET_INVALID_INTEGER (-2)
#define DMET_INVALID_DOUBLE  (-2.0)
#define DMET_INVALID_UINT    (255)

class DmetCHPWriter : public BaseEngine {

  public:

	virtual std::string getEngineName() { return DmetCHPWriter::EngineName(); }
	static const std::string EngineName() { return "DmetCHPWriter"; }

	/**
	* Constructor
	*/
	DmetCHPWriter();

	/**
	* Destructor
	*/
	~DmetCHPWriter(){
      delete m_CallCoder;
    };

	/*! A class to register the summary engine. */
	class Reg : public EngineReg
	{
	public:
		/*! Constructor - register the summary engine. */
		Reg() : EngineReg(DmetCHPWriter::EngineName())
		{
		}

		/*! Creates an object.
		 * @return The object.
		 */
		BaseEngine *MakeObject() { return new DmetCHPWriter; }
	};

	/*! The one and only registration object. */
	static Reg reg;

	/*! Converts the type to the summary engine type.
	 * @param chip The pointer to the base engine object.
	 * @return The summary engine type or NULL if not compatible.
	 */
	static DmetCHPWriter * FromBase(BaseEngine *engine);

    void setAnalysisInfo(const AnalysisInfo &info) {
        m_Info = info;
    }

  private:
    void defineOptions(); 
    void defineStates();
    void runImp();
    void checkOptionsImp();

    void loadLayout(ChipLayout &layout, std::vector<std::string> &probesetNames, 
            std::vector<int> &probesetAlleles);

    void parseRegionsFile(std::vector<std::string> &copyNumberRegionNames, int &regionCount);

    void setHeader(affymetrix_calvin_io::CHPMultiDataData  *data, ParameterNameValueTypeList &params, ParameterNameValueTypeList &metrics);

    void setCopyNumberData(affymetrix_calvin_io::CHPMultiDataFileWriter &writer, std::vector<std::string> copyNumberRegionNames, std::string celFileName);

    void setGenotypeData(affymetrix_calvin_io::CHPMultiDataFileWriter &writer,
                                std::vector<std::string> const &probesetNames,
                                std::vector<int> const &probesetAlleles,
                                std::map<std::string, int> &alleleMap,
                                std::string celFileName,
                                ChipLayout &layout,
                                std::map<std::string, int> &calls,
                                std::map<std::string, int> &forcedCalls,
                                std::map<std::string, double> &confidences,
                                std::map<std::string, double> &summaries,
                                std::map<std::string, std::string> &contexts);

    void getCelNamesFromSummaryFile(std::vector<std::string> &cellFileNames); 

    template <typename T>
    void fillVectorFromA5(const std::string &fileName, const std::string &tsvName, const std::string &sampleName, std::map<std::string,T> &fillMap) {
        affx::File5_File* file5=new affx::File5_File();
        if(file5 == NULL)
            Err::errAbort("Unable to create File5_File pointer.");
        file5->open(fileName,affx::FILE5_OPEN);
        affx::File5_Tsv* tsv5=file5->openTsv(tsvName);
        if(tsv5 == NULL)
            Err::errAbort("Unable to get tsv '"+tsvName+"' in a5 file '"+fileName+"'.");
        
        std::string probeSetName;
        T data;
        while(tsv5->nextLevel(0)==affx::FILE5_OK) {
            tsv5->get(0,0,&probeSetName);
            tsv5->get(0,sampleName,&data);
            fillMap[probeSetName] = data;
        } 

        tsv5->close();
        delete tsv5;
        file5->close();
        delete file5;
    }
    void fillInCallMetrics(     std::map<std::string,int> &calls, 
                                ParameterNameValueTypeList &metrics);

    void fillInReportMetrics(   const std::string &reportFile, 
                                const std::string &celFile, 
                                ParameterNameValueTypeList &metrics);
    AnalysisInfo      m_Info;
    GenoCallCoder*    m_CallCoder;  

}; 




