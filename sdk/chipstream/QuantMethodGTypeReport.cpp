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

//
#include "chipstream/QuantMethodGTypeReport.h"
//
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantLabelZMulti.h"
#include "file5/File5.h"
#include "util/Fs.h"

QuantMethodGTypeReport::QuantMethodGTypeReport(bool outputForcedCalls, bool outputContext, bool outputProbabilities, int probFileSampleCount) {
    affx::TsvReport::init();
    m_is_header_buffer=true;
    m_outputForcedCalls = outputForcedCalls;
    m_outputContext = outputContext;
    m_outputProbabilities = outputProbabilities;
    m_ProbeSetsToReport = NULL;
    m_File5 = NULL;
    m_F5Group = NULL;
    m_MaxNameLength = TSVREPORT_PROBESET_STRLEN;
    m_DoCompact = false;
    m_probFileSampleCount = probFileSampleCount;
} 

/** Destructor. */
QuantMethodGTypeReport::~QuantMethodGTypeReport() {
  if (m_F5Group != NULL) {
    m_F5Group->close();
    Freez(m_F5Group);
  }
  if (m_File5 != NULL) {
    m_File5->close();
    Freez(m_File5);
  }
};

void QuantMethodGTypeReport::setCompactFile5Format(int maxNameLength) {
    std::string fileName = Fs::join(m_dir_path,m_file_prefix + ".genotypes.a5");
    m_File5 = TsvReport::openA5File(fileName, true);
    m_group_name ="/";
    m_MaxNameLength = maxNameLength;
    m_F5Group = m_File5->openGroup(m_group_name,affx::FILE5_CREATE|affx::FILE5_OPEN);
    m_a5_shared_group = m_F5Group;
    m_CallsOutTsv.setA5SharedGroup(m_F5Group);
    m_ConfsOutTsv.setA5SharedGroup(m_F5Group);
    m_ContextOutTsv.setA5SharedGroup(m_F5Group);
    m_ForcedCallsOutTsv.setA5SharedGroup(m_F5Group);
    for (std::vector<affx::TsvReport>::iterator it = m_ProbabilitiesOutTsvVec.begin(); it != m_ProbabilitiesOutTsvVec.end(); ++it) {
        it->setA5SharedGroup(m_F5Group);
    }
    m_file_prefix="/";
    m_format = TsvReport::FMT_A5;
    m_DoCompact = true;
  }

/** 
 * Get set up for a run of reporting probesets. Often used to open file
 * streams and print headers to files etc.
 * 
 * @param qMethod - Quantification method to be used.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeReport::prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
    // printf("### QuantMethodGTypeReport::prepare!\n");

    // change col names "/a/b/c.cel" => "c.cel"
    std::vector<std::string> col_names=Fs::basename(iMart.getCelFileNames());
    Verbose::out(2,"Maximum name length is: " + ToStr(m_MaxNameLength));
    string delim = ".";
    if (m_DoCompact) {
      delim = "";
    }

    // Calls
    m_CallsOutTsv.setDirPath(m_dir_path);
    m_CallsOutTsv.setFilename(m_file_prefix + delim + "calls");
    m_CallsOutTsv.setA5SharedGroup(m_a5_shared_group);
    m_CallsOutTsv.setFormat(m_format);
    m_CallsOutTsv.addHeaderComment("Calls: -1=NN, 0=AA, 1=AB, 2=BB");
    m_CallsOutTsv.defineStringColumn(0,0,"probeset_id",m_MaxNameLength);
    m_CallsOutTsv.setPrecision(0);
    if (m_DoCompact) {
      m_CallsOutTsv.defineColumns(col_names,affx::FILE5_DTYPE_CHAR,0);
    }
    else {
      m_CallsOutTsv.defineColumns(col_names,affx::FILE5_DTYPE_INT,0);
    }
    
    m_CallsOutTsv.addHeadersFrom(m_header_buffer);
    m_CallsOutTsv.writeTsv_v1();
  
    // Confidences
    m_ConfsOutTsv.setDirPath(m_dir_path);
    m_ConfsOutTsv.setFilename(m_file_prefix + delim + "confidences");
    m_ConfsOutTsv.setA5SharedGroup(m_a5_shared_group);
    m_ConfsOutTsv.setFormat(m_format);
    m_ConfsOutTsv.defineStringColumn(0,0,"probeset_id",m_MaxNameLength);
    if (m_DoCompact) {
      m_ConfsOutTsv.defineColumns(col_names,affx::FILE5_DTYPE_FLOAT,0);
    }
    else {
      m_ConfsOutTsv.defineColumns(col_names,affx::FILE5_DTYPE_DOUBLE,0);
    }
    m_ConfsOutTsv.addHeadersFrom(m_header_buffer);
    m_ConfsOutTsv.writeTsv_v1();
    
    // Forced Calls
    if (m_outputForcedCalls) {
        m_ForcedCallsOutTsv.setDirPath(m_dir_path);
        m_ForcedCallsOutTsv.setFilename(m_file_prefix + delim + "forced-calls");
        m_ForcedCallsOutTsv.setA5SharedGroup(m_a5_shared_group);
        m_ForcedCallsOutTsv.setFormat(m_format);
        m_ForcedCallsOutTsv.defineStringColumn(0,0,"probeset_id",m_MaxNameLength);
        m_ForcedCallsOutTsv.defineColumns(col_names,affx::FILE5_DTYPE_CHAR,0);
        m_ForcedCallsOutTsv.addHeadersFrom(m_header_buffer);
        m_ForcedCallsOutTsv.writeTsv_v1();
    }
    
    // Allele/Wobble Context
    QuantLabelZMulti *qlzm = dynamic_cast<QuantLabelZMulti *>(&qMethod);
    if (qlzm == NULL)
        m_outputContext = false;
    if (m_outputContext) {
        m_ContextOutTsv.setDirPath(m_dir_path);
        m_ContextOutTsv.setFilename(m_file_prefix + delim + "context");
        m_ContextOutTsv.setA5SharedGroup(m_a5_shared_group);
        m_ContextOutTsv.setFormat(m_format);
        m_ContextOutTsv.defineStringColumn(0,0,"probeset_id",m_MaxNameLength);
        for (int i=0; i<col_names.size(); i++) 
            m_ContextOutTsv.defineStringColumn(0,1+i,col_names[i],20); ///@todo should not hard code size
        m_ContextOutTsv.addHeadersFrom(m_header_buffer);
        m_ContextOutTsv.writeTsv_v1();
    }
    // Allele probabilities
    QuantLabelZ* qlz = dynamic_cast<QuantLabelZ*>(&qMethod);
    if (qlz == NULL) {
      m_outputProbabilities = false;
    }
    if (m_outputProbabilities) {
      m_ProbabilitiesOutTsvVec.resize(std::ceil((double)col_names.size()/m_probFileSampleCount));
      for (std::vector<affx::TsvReport>::iterator it = m_ProbabilitiesOutTsvVec.begin(); it != m_ProbabilitiesOutTsvVec.end(); ++it) {
        
        it->setDirPath(m_dir_path);
        int fileIndex = it -  m_ProbabilitiesOutTsvVec.begin();
        std::string suffix = "";
        if (m_ProbabilitiesOutTsvVec.size() > 1) {
          suffix = ToStr(fileIndex);
        }
        it->setFilename(m_file_prefix + delim + "probabilities" + suffix);
        it->setA5SharedGroup(m_a5_shared_group);
        it->setFormat(m_format);
        it->defineStringColumn(0,0,"probeset_id",m_MaxNameLength);
        int start = m_probFileSampleCount*(fileIndex);
        int count = Min(m_probFileSampleCount, (int)col_names.size()-start);
        for (int i=0; i<count; ++i) {
          it->defineStringColumn(0,1+i,col_names[start+i],48); // @todo should not hard code size
        }
        it->addHeadersFrom(m_header_buffer);
        it->writeTsv_v1();
      }
    }
    return true;
}

/** 
 * After every probeset computation this function is called an is an opportunity
 * to query the quantification method for results, residuals, etc.
 * 
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeReport::report(ProbeSetGroup& psGroup,
                                    QuantMethod& qMethod, 
                                    const IntensityMart& iMart, 
                                    std::vector<ChipStream *>& iTrans, 
                                    PmAdjuster& pmAdjust) {
    QuantGTypeMethod *gMethod = dynamic_cast<QuantGTypeMethod *>(&qMethod);
    if (gMethod == NULL) {
      Err::errAbort("Can only use a QuantMethodGTypeReport with QuantGTypeMethods.");
    }
    std::string name = gMethod->getProbeSetName();
    
    if (m_ProbeSetsToReport != NULL)
        if (m_ProbeSetsToReport->find(name.c_str())==m_ProbeSetsToReport->end())
            return false;

    // "0" is the probeset_id
    m_CallsOutTsv.set_string(0,0,name);
    m_ConfsOutTsv.set_string(0,0,name);
    if (m_outputForcedCalls)
        m_ForcedCallsOutTsv.set_string(0,0,name);
    if (m_outputContext)
        m_ContextOutTsv.set_string(0,0,name);
    if (m_outputProbabilities) {
      for (std::vector<affx::TsvReport>::iterator it = m_ProbabilitiesOutTsvVec.begin(); it != m_ProbabilitiesOutTsvVec.end(); ++it) {
        it->set_string(0,0,name);
      }
    }
    // printf("%s : ",name.c_str());     // JHG -debug

    QuantLabelZMulti *qlzm = dynamic_cast<QuantLabelZMulti *>(&qMethod);
    if ((qlzm == NULL) && (m_outputContext == true))
        Err::errAbort("Context output set to true but quant method is not a LabelZMulti!");
    QuantLabelZ *qlz = dynamic_cast<QuantLabelZ *>(&qMethod);
    if ((qlz == NULL) && (m_outputProbabilities == true)) {
        Err::errAbort("allele-probabilities output set to true but quant method is not a QuantLabelZ or QuantLabelZMulti!");
    }
    for (unsigned int i = 0; i < gMethod->getNumCalls(); i++) {
      //
      if (m_DoCompact) {
        m_ConfsOutTsv.set_f(0,i+1,(float)gMethod->getConfidence(i));
        m_CallsOutTsv.set_c(0,i+1,(char)(gMethod->getCall(i)));
      }
      else {
        m_ConfsOutTsv.set_d(0,i+1,gMethod->getConfidence(i));
        m_CallsOutTsv.set_i(0,i+1,(gMethod->getCall(i)));
      }
      //
      if (m_outputForcedCalls)
        m_ForcedCallsOutTsv.set_i(0,i+1,(gMethod->getForcedCall(i)));
      if (m_outputContext)
        m_ContextOutTsv.set_string(0,i+1,qlzm->getContext(i));
      if (m_outputProbabilities) {
        int file_idx = (std::floor)((float)i/m_probFileSampleCount);
        int col_idx = i-(file_idx*m_probFileSampleCount)+1;
        m_ProbabilitiesOutTsvVec[file_idx].set_string(0,col_idx,qlz->getCallProbabilities(i));
      }
    }
    // printf("\n"); // JHG -debug

    m_CallsOutTsv.writeLevel(0);
    m_ConfsOutTsv.writeLevel(0);
    if (m_outputForcedCalls)
        m_ForcedCallsOutTsv.writeLevel(0);
    if (m_outputContext)
        m_ContextOutTsv.writeLevel(0);
    if (m_outputProbabilities) {
        for (std::vector<affx::TsvReport>::iterator it = m_ProbabilitiesOutTsvVec.begin(); it != m_ProbabilitiesOutTsvVec.end(); ++it) {      
            it->writeLevel(0);
        }
    }

    return true;
}
    
/** 
 * No more probesets will be processed, this is a chance to finish outputting
 * results and clean up.
 * 
 * @param qMethod - Quantification method that was used.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeReport::finish(QuantMethod &qMethod) {
    m_CallsOutTsv.close();
    m_ConfsOutTsv.close();
    if (m_outputForcedCalls)
        m_ForcedCallsOutTsv.close();
    if (m_outputContext)
        m_ContextOutTsv.close();
    if (m_outputProbabilities) {
        for (std::vector<affx::TsvReport>::iterator it = m_ProbabilitiesOutTsvVec.begin(); it != m_ProbabilitiesOutTsvVec.end(); ++it) {
            it->close();
        }
    }
    return true;
}
