////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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
 * @file CNFamilialReporterFamilial.cpp
 *
 * @brief This file contains the CNFamilialReporterFamilial class members.
 */

//
#include "copynumber/CNFamilialReporterFamilial.h"
//
#include "util/Fs.h"
#include "util/AffxSTL.h"

#include <iterator>


/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNFamilialReporterFamilial::explainSelf()
{
    SelfDoc doc;
    doc.setDocName(CNFamilialReporterFamilial::getType());
    doc.setDocDescription(CNFamilialReporterFamilial::getDescription());
    doc.setDocOptions(CNFamilialReporterFamilial::getDefaultDocOptions());
    return doc;
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> CNFamilialReporterFamilial::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)

  return opts;
}

/**
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 *
 * @param param - Map of key/value pairs to initialize the object.
 *
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate* CNFamilialReporterFamilial::newObject(std::map<std::string,std::string>& params)
{
    //SelfDoc doc = explainSelf();
    //std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNFamilialReporterFamilial* pMethod = new CNFamilialReporterFamilial();
    //std::string strPrefix = getPrefix();

    return pMethod;
}

/**
 * @brief Constructor
 */
CNFamilialReporterFamilial::CNFamilialReporterFamilial()
{
}

/**
 * @brief Destructor
 */
CNFamilialReporterFamilial::~CNFamilialReporterFamilial()
{
}

/**
 * @brief Run the analysis.
 */
void CNFamilialReporterFamilial::run()
{
    if (m_pEngine == NULL) {return;}
    AffxString famFileName = m_pEngine->getOpt("familial-file");
    AffxString strOutDirName = m_pEngine->getOpt("out-dir");
    AffxString strFileName = Fs::join(strOutDirName,famFileName);
    Verbose::out(1, "Writing familial-file: " + strFileName);

    affymetrix_calvin_io::GenericData data;
    setupHeader(strFileName, data);

    // write the file.
    affymetrix_calvin_io::GenericFileWriter* writer = new affymetrix_calvin_io::GenericFileWriter(&data.Header());
    writer->WriteHeader();
    writeDataGroupAndDataSetHeaders(writer);
    writeDataSetData(writer);

    delete writer;

    // Convert from calvin to text format?
    if (m_pEngine->getOptBool("text-output"))
    {
        CalvinToText converter;
        converter.run(strFileName, strFileName + ".txt", true, true, true);
    }
}

/**
 * @brief Setup the file header in memory before writing if out
 * @param const AffxString& - The file name
 * @param affymetrix_calvin_io::GenericData& - The object to load the data into
 */
void CNFamilialReporterFamilial::setupHeader(const AffxString& strFileName, affymetrix_calvin_io::GenericData& data)
{
    addParentHeaders(data);
    data.Header().SetFilename(strFileName);
    data.Header().AddDataGroupHdr(affymetrix_calvin_io::DataGroupHeader(L"Familial"));
    //data.Header().AddDataGroupHdr(affymetrix_calvin_io::DataGroupHeader(L"Segments"));

    // Setup the header.
    affymetrix_calvin_io::GenericDataHeader* pHeader = data.Header().GetGenericDataHdr();
    pHeader->SetFileCreationTime(DateTime::GetCurrentDateTime().ToString());
    addStandardHeader(pHeader, m_pEngine);
    ///@todo add in headers
    //addPgOptions(pHeader, m_pEngine->getOptions(), "affymetrix-algorithm-param-option-");
    //addPgOptions(pHeader, m_pEngine->getStates(), "affymetrix-algorithm-param-state-");
    for (int iIndex = 0; (iIndex < (int)CNFamilialAnalysisMethod::getParams().size()); iIndex++)
    {
        pHeader->AddNameValParam(CNFamilialAnalysisMethod::getParams().at(iIndex));
    }

    // set up the data set headers.
    int iGroupIndex = 0;
    //JPB: setupDataSetHeaderSamples() and setupDataSetHeaderPaternityTesting() may end up inside the for loop
    setupDataSetHeaderSamples(data, iGroupIndex);
    setupDataSetHeaderPaternityTesting(data, iGroupIndex);
    setupDataSetHeaderMIE(data, iGroupIndex);
    for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++)
    {
        CNFamilialAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
        if (pMethod->isFamilialAnalysis())  {continue;}
        else if (pMethod->isSegmentOverlapAnalysis()) {setupDataSetHeaderSegmentOverlaps(data, pMethod, iGroupIndex);}
        else if (pMethod->isSegmentTypeAnalysis()) {continue;}
        else {throw(Except("CNFamilialReporterFamilial::run(...) Unknown CNFamilialAnalysisMethod: " + pMethod->getName()));}
    }
    iGroupIndex++;
    for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++)
    {
        CNFamilialAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
        if (pMethod->isFamilialAnalysis())  {continue;}
        else if (pMethod->isSegmentOverlapAnalysis()) {continue;}
        else if (pMethod->isSegmentTypeAnalysis()) {setupDataSetHeaderSegments(data, pMethod, iGroupIndex);}
        else {throw(Except("CNFamilialReporterFamilial::run(...) Unknown CNFamilialAnalysisMethod: " + pMethod->getName()));}
    }
}

/**
 * @brief Add parent file header to the header for this file
 * @param affymetrix_calvin_io::GenericData& - The object to load the data into
 */
void CNFamilialReporterFamilial::addParentHeaders(affymetrix_calvin_io::GenericData& data)
{
    CNCychp& cychpIndex = getCychpIndex();
    CNCychp& cychpMother = getCychpMother();
    CNCychp& cychpFather = getCychpFather();

    // Copy over headers from the cychp files.
    try
    {
        affymetrix_calvin_io::GenericData genericData;
        if (cychpIndex.getFileName() != "")
        {
            affymetrix_calvin_io::GenericFileReader reader;
            reader.SetFilename(cychpIndex.getFileName().c_str());
            reader.Open(genericData);
            data.Header().GetGenericDataHdr()->AddParent(*genericData.Header().GetGenericDataHdr());
            genericData.Clear();
        }
        if (cychpMother.getFileName() != "")
        {
            affymetrix_calvin_io::GenericFileReader reader;
            reader.SetFilename(cychpMother.getFileName().c_str());
            reader.Open(genericData);
            data.Header().GetGenericDataHdr()->AddParent(*genericData.Header().GetGenericDataHdr());
            genericData.Clear();
        }
        if (cychpFather.getFileName() != "")
        {
            affymetrix_calvin_io::GenericFileReader reader;
            reader.SetFilename(cychpFather.getFileName().c_str());
            reader.Open(genericData);
            data.Header().GetGenericDataHdr()->AddParent(*genericData.Header().GetGenericDataHdr());
            genericData.Clear();
        }
    } catch(...) {Verbose::out(1, "WARNING: Load parent header failed.");}

}

/**
 * @brief Add the standard file header information
 * @param affymetrix_calvin_io::GenericDataHeader* - The header to load
 * @param BaseEngine* - The engine to get parameter information from
 */
void CNFamilialReporterFamilial::addStandardHeader(affymetrix_calvin_io::GenericDataHeader* pHeader, BaseEngine* pEngine)
{
    affymetrix_calvin_parameter::ParameterNameValueType param;

    pHeader->SetFileTypeId("affymetrix-multi-data-type-analysis");

    param.SetName(L"affymetrix-algorithm-param-file-name");
    param.SetValueText(StringUtils::ConvertMBSToWCS(Fs::basename(m_pEngine->getXMLParameterFileName())));
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-param-file-guid");
    param.SetValueText(StringUtils::ConvertMBSToWCS(m_pEngine->getXMLParameterFileGuid()));
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-name");
    param.SetValueText(L"CYTO2");
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-version");
    param.SetValueText(L"1.0.0");
    pHeader->AddNameValParam(param);

    param.SetName(L"program-name");
    param.SetValueText(StringUtils::ConvertMBSToWCS(pEngine->getOpt("program-name")));
    pHeader->AddNameValParam(param);

    param.SetName(L"program-version");
    param.SetValueText(StringUtils::ConvertMBSToWCS(pEngine->getOpt("program-version")));
    pHeader->AddNameValParam(param);

    param.SetName(L"program-company");
    param.SetValueText(StringUtils::ConvertMBSToWCS("Affymetrix"));
    pHeader->AddNameValParam(param);

    param.SetName(L"create-date");
    param.SetValueText(StringUtils::ConvertMBSToWCS(Util::getTimeStamp()));
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-param-file-name");
    param.SetValueText(StringUtils::ConvertMBSToWCS(Fs::basename(pEngine->getXMLParameterFileName())));
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-param-file-guid");
    param.SetValueText(StringUtils::ConvertMBSToWCS(pEngine->getXMLParameterFileGuid()));
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-param-familial-file-version");
    param.SetValueInt32(2);
    pHeader->AddNameValParam(param);
}

/**
 * @brief Dump the specified PgOptions into the header
 * @param affymetrix_calvin_io::GenericDataHeader* - The header to load
 * @param PgOptions& - The options to dump
 * @param const AffxString& - The prefix to use when naming the options
 */
void CNFamilialReporterFamilial::addPgOptions(affymetrix_calvin_io::GenericDataHeader* pHeader, PgOptions& opt, const AffxString& strPrefix)
{
    affymetrix_calvin_parameter::ParameterNameValueType param;
    for(int i=0; i< opt.m_option_vec.size(); i++) {
        if ((int)opt.m_option_vec[i]->m_values.size() > 1)
        {
            for(int j=0; j< opt.m_option_vec[i]->m_values.size(); j++)
            {
                loadParam(strPrefix + opt.m_option_vec[i]->m_longName + "-" + ::getInt(j+1), opt.m_option_vec[i]->m_type, opt.m_option_vec[i]->m_values[j], param);
                pHeader->AddNameValParam(param);
            }
        }
        else
        {
            loadParam(strPrefix + opt.m_option_vec[i]->m_longName, opt.m_option_vec[i]->m_type, opt.m_option_vec[i]->getValue(), param);
            pHeader->AddNameValParam(param);
        }
    }
}

/**
 * @brief Load a parameter from it's components
 * @param const AffxString& - The name of the parameter
 * @param PgOpt::PgOptType - The type of the parameter
 * @param const AffxString& - The value of the parameter
 * @param affymetrix_calvin_parameter::ParameterNameValueType& - The parameter to load
 */
void CNFamilialReporterFamilial::loadParam(const AffxString& strName, PgOpt::PgOptType type, const AffxString& strValue, affymetrix_calvin_parameter::ParameterNameValueType& param)
{
    param.SetName(StringUtils::ConvertMBSToWCS(strName));
    switch(type)
    {
    case PgOpt::BOOL_OPT: param.SetValueInt8(((strValue == "true") ? (char)1 : (char)0)); break;
    case PgOpt::DOUBLE_OPT: param.SetValueFloat((float)::getDouble(strValue)); break;
    case PgOpt::INT_OPT: param.SetValueInt32(::getInt(strValue)); break;
    case PgOpt::STRING_OPT:
        if ((strName.indexOf("command-line") != -1) || (strName.indexOf("analysis") != -1) || (strName.indexOf("program-cvs-id") != -1) || (strName.indexOf("version-to-report") != -1) || (strName.endsWith("-dir")))
        {
            param.SetValueText(StringUtils::ConvertMBSToWCS(strValue));
        }
        else
        {
            param.SetValueText(StringUtils::ConvertMBSToWCS(Fs::basename(strValue)));
        }
        break;
    default: throw(Except("Cannot find PgOpt type for: " + strName));
    }
}

/**
 * @brief Setup the Samples data set
 * @param affymetrix_calvin_io::GenericData& - The object to load the data into
 * @param int - The group index
 */
void CNFamilialReporterFamilial::setupDataSetHeaderSamples(affymetrix_calvin_io::GenericData& data, int iGroupIndex)
{
    CNCychp& cychpMother = getCychpMother();
    CNCychp& cychpFather = getCychpFather();

    affymetrix_calvin_io::DataSetHeader dsHeader;
    if ((cychpMother.getFileName() != "") && (cychpFather.getFileName() != ""))
    {
        dsHeader.SetRowCnt(3);
    }
    else
    {
        dsHeader.SetRowCnt(2);
    }
    dsHeader.SetName(L"Samples");
    dsHeader.AddUIntColumn(L"SampleKey");
    dsHeader.AddUnicodeColumn(L"CHPFilename", 256);
    dsHeader.AddAsciiColumn(L"CHPID", 55);
    dsHeader.AddAsciiColumn(L"Role", 21);
    data.Header().GetDataGroup(iGroupIndex).AddDataSetHdr(dsHeader);
}

/**
 * @brief Setup the PaternityTesting data set
 * @param affymetrix_calvin_io::GenericData& - The object to load the data into
 * @param int - The group index
 */
void CNFamilialReporterFamilial::setupDataSetHeaderPaternityTesting(affymetrix_calvin_io::GenericData& data, int iGroupIndex)
{
    CNCychp& cychpMother = getCychpMother();
    CNCychp& cychpFather = getCychpFather();

    affymetrix_calvin_io::DataSetHeader dsHeader;
    if ((cychpMother.getFileName() != "") && (cychpFather.getFileName() != ""))
    {
        dsHeader.SetRowCnt(3);
    }
    else
    {
        dsHeader.SetRowCnt(1);
    }
    dsHeader.SetName(L"PaternityTesting");
    dsHeader.AddUByteColumn(L"AnalysisType");
    dsHeader.AddUIntColumn(L"ReferenceSampleKey");
    dsHeader.AddUIntColumn(L"FamilialSampleKey");
    dsHeader.AddUByteColumn(L"RoleValidity");
    dsHeader.AddFloatColumn(L"RoleIndexScore");
    data.Header().GetDataGroup(iGroupIndex).AddDataSetHdr(dsHeader);
}

void CNFamilialReporterFamilial::setupDataSetHeaderMIE(affymetrix_calvin_io::GenericData& data, int iGroupIndex)
{
    affymetrix_calvin_io::DataSetHeader dsHeader;
    dsHeader.SetRowCnt(23);     // # autosomes + X chromosome

    dsHeader.SetName(L"MIE");
    dsHeader.AddUByteColumn(L"Chromosome");
    dsHeader.AddAsciiColumn(L"Display", 4);
    dsHeader.AddUIntColumn(L"MarkerCount");
    dsHeader.AddIntColumn(L"MIE-Trio");
    dsHeader.AddIntColumn(L"MIE-Mat");
    dsHeader.AddIntColumn(L"MIE-Pat");
    data.Header().GetDataGroup(iGroupIndex).AddDataSetHdr(dsHeader);
}

/**
 * @brief Setup the Segment Overlaps data set
 * @param affymetrix_calvin_io::GenericData& - The object to load the data into
 * @param CNFamilialAnalysisMethod* - The analysis method to load data from
 * @param int - The group index
 */
void CNFamilialReporterFamilial::setupDataSetHeaderSegmentOverlaps(affymetrix_calvin_io::GenericData& data, CNFamilialAnalysisMethod* pMethod, int iGroupIndex)
{
    if (pMethod == NULL) {return;}
    int iSegmentCount = ((pMethod == NULL) ? 0 : pMethod->getSegmentOverlapCount());
    affymetrix_calvin_io::DataSetHeader dsHeader;
    dsHeader.SetRowCnt(iSegmentCount);
    dsHeader.SetName(L"SegmentOverlaps");
    dsHeader.AddAsciiColumn(L"SegmentType", 21);
    dsHeader.AddUIntColumn(L"ReferenceSampleKey");
    dsHeader.AddAsciiColumn(L"ReferenceSegmentID", 55);
    dsHeader.AddUIntColumn(L"FamilialSampleKey");
    dsHeader.AddAsciiColumn(L"FamilialSegmentID", 55);
    data.Header().GetDataGroup(iGroupIndex).AddDataSetHdr(dsHeader);
}

/**
 * @brief Setup the Segments data set
 * @param affymetrix_calvin_io::GenericData& - The object to load the data into
 * @param CNFamilialAnalysisMethod* - The analysis method to load data from
 * @param int - The group index
 */
void CNFamilialReporterFamilial::setupDataSetHeaderSegments(affymetrix_calvin_io::GenericData& data, CNFamilialAnalysisMethod* pMethod, int iGroupIndex)
{
    if (pMethod == NULL) {return;}
    int iSegmentCount = ((pMethod == NULL) ? 0 : pMethod->getSegmentCount());
    affymetrix_calvin_io::DataSetHeader dsHeader;
    dsHeader.SetRowCnt(iSegmentCount);
    dsHeader.SetName(StringUtils::ConvertMBSToWCS(pMethod->getSegmentTypeString()));
    dsHeader.AddAsciiColumn(L"SegmentID", 55);
    dsHeader.AddUIntColumn(L"ReferenceSampleKey");
    dsHeader.AddUIntColumn(L"FamilialSampleKey");
    dsHeader.AddUByteColumn(L"Chromosome");
    dsHeader.AddUIntColumn(L"StartPosition");
    dsHeader.AddUIntColumn(L"StopPosition");
    dsHeader.AddUByteColumn(L"Call");
    dsHeader.AddFloatColumn(L"Confidence");
    dsHeader.AddIntColumn(L"MarkerCount");
    dsHeader.AddFloatColumn(L"Homozygosity");
    dsHeader.AddFloatColumn(L"Heterozygosity");
    data.Header().GetDataGroup(iGroupIndex).AddDataSetHdr(dsHeader);
}

/**
 * @brief Write the data group and data set headers
 * @param affymetrix_calvin_io::GenericFileWriter* - The generic calvin file writer
 * @return int - The offset to the samples data set
 */
void CNFamilialReporterFamilial::writeDataGroupAndDataSetHeaders(affymetrix_calvin_io::GenericFileWriter* writer)
{
    int iGroupIndex = 0;
    int iSetIndex = 0;

    // write out the familial analysis header
    writer->GetDataGroupWriter(iGroupIndex).WriteHeader();

    // write out the segement overlap analysis header
    for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++)
    {
        CNFamilialAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
        if (pMethod->isFamilialAnalysis())
        {
            pMethod->setDataSetOffsetVec(setupDataSet(writer, iGroupIndex, iSetIndex));
            iSetIndex++;
            pMethod->setDataSetOffsetVec(setupDataSet(writer, iGroupIndex, iSetIndex));
            iSetIndex++;
            pMethod->setDataSetOffsetVec(setupDataSet(writer, iGroupIndex, iSetIndex));
            iSetIndex++;
        }
        else if (pMethod->isSegmentOverlapAnalysis())
        {
            pMethod->setDataSetOffsetVec(setupDataSet(writer, iGroupIndex, iSetIndex));
            iSetIndex++;
        }
        else if (pMethod->isSegmentTypeAnalysis())  {continue;}
        else {throw(Except("CNFamilialReporterFamilial::run(...) Unknown CNFamilialAnalysisMethod: " + pMethod->getName()));}
    }

    // write out the segment type analysis header
    //iGroupIndex++;
    //iSetIndex = 0;
    //writer->GetDataGroupWriter(iGroupIndex).WriteHeader();
    //for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++)
    //{
    //    CNFamilialAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
    //    if (pMethod->isFamilialAnalysis())  {continue;}
    //    else if (pMethod->isSegmentOverlapAnalysis()) {continue;}
    //    else if (pMethod->isSegmentTypeAnalysis())
    //    {
    //        pMethod->setDataSetOffset(setupDataSet(writer, iGroupIndex, iSetIndex));
    //        iSetIndex++;
    //    }
    //    else {throw(Except("CNFamilialReporterFamilial::run(...) Unknown CNFamilialAnalysisMethod: " + pMethod->getName()));}
    //}
    //return iFamilialOffset;
}

/**
 * @brief Setup the data set
 * @param affymetrix_calvin_io::GenericFileWriter* - The generic calvin file writer
 * @param int - The group index
 * @param int - The set index
 * @return int - The offset to the data set
 */
int CNFamilialReporterFamilial::setupDataSet(affymetrix_calvin_io::GenericFileWriter* writer, int iGroupIndex, int iSetIndex)
{
    affymetrix_calvin_io::DataSetWriter& setWriter = writer->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    setWriter.WriteHeader();
    int iOffset = writer->GetFilePos();
    writer->SeekFromCurrentPos(setWriter.GetDataSetSize() + 1);
    setWriter.UpdateNextDataSetOffset();
    writer->GetDataGroupWriter(iGroupIndex).UpdateNextDataGroupPos();
    return iOffset;
}

/**
 * @brief Write the data set data
 * @param affymetrix_calvin_io::GenericFileWriter* - The generic calvin file writer
 * @param int - The offset to the Samples data set
 */
void CNFamilialReporterFamilial::writeDataSetData(affymetrix_calvin_io::GenericFileWriter* writer)
{
    int iGroupIndex = 0;
    int iSetIndex = 0;

    // write out the familial and segment overlap analysis data
    for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++)
    {
        CNFamilialAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
        if (pMethod->isFamilialAnalysis())
        {
            writeDataSetSamples(writer, pMethod, iGroupIndex, iSetIndex);
            iSetIndex++;
            writeDataSetPaternityTesting(writer, pMethod, iGroupIndex, iSetIndex);
            iSetIndex++;
            writeDataSetMIE(writer, pMethod, iGroupIndex, iSetIndex);
            iSetIndex++;
        }
        else if (pMethod->isSegmentOverlapAnalysis())
        {
            writeDataSetSegmentOverlaps(writer, pMethod, iGroupIndex, iSetIndex);
            iSetIndex++;
        }
        else if (pMethod->isSegmentTypeAnalysis())  {continue;}
        else {throw(Except("CNFamilialReporterFamilial::run(...) Unknown CNFamilialAnalysisMethod: " + pMethod->getName()));}
    }

    // write out the segment type analysis data
    iGroupIndex++;
    iSetIndex = 0;
    for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++)
    {
        CNFamilialAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
        if (pMethod->isFamilialAnalysis())  {continue;}
        else if (pMethod->isSegmentOverlapAnalysis()) {continue;}
        else if (pMethod->isSegmentTypeAnalysis())
        {
            writeDataSetSegments(writer, pMethod, iGroupIndex, iSetIndex);
            iSetIndex++;
        }
        else {throw(Except("CNFamilialReporterFamilial::run(...) Unknown CNFamilialAnalysisMethod: " + pMethod->getName()));}
    }
}

/**
 * @brief Write the Samples data set data
 * @param affymetrix_calvin_io::GenericFileWriter* - The generic calvin file writer
 * @param CNFamilialAnalysisMethod* - The analysis method to load data from
 * @param int - The group index
 * @param int - The set index
 */
void CNFamilialReporterFamilial::writeDataSetSamples(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex)
{
    CNCychp& cychpIndex = getCychpIndex();
    CNCychp& cychpMother = getCychpMother();
    CNCychp& cychpFather = getCychpFather();

    writer->SeekFromBeginPos(pMethod->getDataSetOffsetVec().at(iSetIndex));
    affymetrix_calvin_io::DataSetWriter& setWriter = writer->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    setWriter.Write(0);
    setWriter.Write(StringUtils::ConvertMBSToWCS(Fs::basename(cychpIndex.getFileName())), 256);
    setWriter.Write(cychpIndex.getCychpHeader().getChpID(), 55);
    setWriter.Write(cychpIndex.getFamilialType(), 21);

    if (cychpMother.getFileName() != "")
    {
        setWriter.Write(1);
        setWriter.Write(StringUtils::ConvertMBSToWCS(Fs::basename(cychpMother.getFileName())), 256);
        setWriter.Write(cychpMother.getCychpHeader().getChpID(), 55);
        setWriter.Write(cychpMother.getFamilialType(), 21);
    }

    if (cychpFather.getFileName() != "")
    {
        setWriter.Write(2);
        setWriter.Write(StringUtils::ConvertMBSToWCS(Fs::basename(cychpFather.getFileName())), 256);
        setWriter.Write(cychpFather.getCychpHeader().getChpID(), 55);
        setWriter.Write(cychpFather.getFamilialType(), 21);
    }
}

/**
 * @brief Write the PaternityTesting data set data
 * @param affymetrix_calvin_io::GenericFileWriter* - The generic calvin file writer
 * @param CNFamilialAnalysisMethod* - The analysis method to load data from
 * @param int - The group index
 * @param int - The set index
 */
void CNFamilialReporterFamilial::writeDataSetPaternityTesting(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex)
{
    CNCychp& cychpMother = getCychpMother();
    CNCychp& cychpFather = getCychpFather();

    writer->SeekFromBeginPos(pMethod->getDataSetOffsetVec().at(iSetIndex));
    affymetrix_calvin_io::DataSetWriter& setWriter = writer->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);

    if (cychpMother.getFileName() != "" && cychpFather.getFileName() != "")
    {
        setWriter.Write((unsigned char)0);    // paternity with maternity assumed
        setWriter.Write(0);                   // child sample ID
        setWriter.Write(2);                   // father sample ID
        setWriter.Write((unsigned char)cychpFather.getFamilialCall());
        setWriter.Write((float)cychpFather.getFamilialConfidence());
    }

    if (cychpMother.getFileName() != "")
    {
        setWriter.Write((unsigned char)1);    // maternity
        setWriter.Write(0);                   // child sample ID
        setWriter.Write(1);                   // mother sample ID
        setWriter.Write((unsigned char)cychpMother.getFamilialCall());
        setWriter.Write((float)cychpMother.getFamilialConfidence());
    }

    if (cychpFather.getFileName() != "")
    {
        setWriter.Write((unsigned char)2);    // paternity
        setWriter.Write(0);                   // child sample ID
        setWriter.Write(2);                   // father sample ID
        setWriter.Write((unsigned char)cychpFather.getFamilialDuoCall());
        setWriter.Write((float)cychpFather.getFamilialDuoConfidence());
    }
}

void CNFamilialReporterFamilial::writeDataSetMIE(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex)
{
    const int iXChromosome = m_pEngine->getOptInt("xChromosome");
    const int iYChromosome = m_pEngine->getOptInt("yChromosome");
    const int lastAutosome = 22;

    CNCychp& cychpMother = getCychpMother();
    CNCychp& cychpFather = getCychpFather();
    const bool motherPresent = cychpMother.getFileName() != "";
    const bool fatherPresent = cychpFather.getFileName() != "";

    writer->SeekFromBeginPos(pMethod->getDataSetOffsetVec().at(iSetIndex));
    affymetrix_calvin_io::DataSetWriter& setWriter = writer->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);

    std::vector<unsigned char> vChromosomes;

    // copy all keys of CNCychp::m_mMarkerCount to vChromosomes
    std::transform(
                CNCychp::m_mMarkerCount.begin(),
                CNCychp::m_mMarkerCount.end(), 
                std::back_inserter(vChromosomes),
                affxstl::select1st<std::map<unsigned char, unsigned int>::value_type>()
                );

    std::sort(vChromosomes.begin(), vChromosomes.end());

    for (unsigned char chr = 1; chr < iYChromosome; chr++)
    {
        if (chr > lastAutosome && chr < iXChromosome) {
            continue;
        }
        AffxString displayStr;
        if (chr == iXChromosome) {
            displayStr = "  X";
        } else {
            displayStr = StringUtils::ConvertWCSToMBS(StringUtils::ToString(chr, 4, L' '));
        }
        if (std::binary_search(vChromosomes.begin(), vChromosomes.end(), chr))
        {
            // data available
            setWriter.Write(chr);                               // chromosome number
            setWriter.Write(displayStr, 4);                     // chromosome number as string
            setWriter.Write(CNCychp::m_mMarkerCount[chr]);      // number of genotypeable markers
            if (motherPresent && fatherPresent)
            {
                setWriter.Write(CNCychp::m_mMIE_Trio[chr]);         // MIE-Trio
                setWriter.Write(CNCychp::m_mMIE_Mat[chr]);          // MIE-Mat
                setWriter.Write(CNCychp::m_mMIE_Pat[chr]);          // MIE-Pat
            }
            else if (motherPresent)
            {
                setWriter.Write(-1);                                // MIE-Trio
                setWriter.Write(CNCychp::m_mMIE_Mat[chr]);          // MIE-Mat
                setWriter.Write(-1);                                // MIE-Pat
            }
            else if (fatherPresent)
            {
                setWriter.Write(-1);                                // MIE-Trio
                setWriter.Write(-1);                                // MIE-Mat
                setWriter.Write(CNCychp::m_mMIE_Pat[chr]);          // MIE-Pat
            }
        }
        else {
            // data not available
            setWriter.Write(chr);                               // chromosome number
            setWriter.Write(displayStr, 4);                     // chromosome number as string
            setWriter.Write(0);                                 // number of genotypeable markers
            setWriter.Write(-1);                                // MIE-Trio
            setWriter.Write(-1);                                // MIE-Mat
            setWriter.Write(-1);                                // MIE-Pat
        }
    }
}

/**
 * @brief Write the SegmentOverlaps data set data
 * @param affymetrix_calvin_io::GenericFileWriter* - The generic calvin file writer
 * @param CNFamilialAnalysisMethod* - The method to get the data from
 * @param int - The group index
 * @param int - The set index
 */
void CNFamilialReporterFamilial::writeDataSetSegmentOverlaps(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex)
{
    affymetrix_calvin_io::DataSetWriter& setWriter = writer->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    if (pMethod != NULL)
    {
        writer->SeekFromBeginPos(pMethod->getDataSetOffsetVec().at(iSetIndex));
        for (int i = 0; (i < pMethod->getSegmentOverlaps().getCount()); i++)
        {
            CNSegmentOverlap* p = pMethod->getSegmentOverlaps().getAt(i);
            setWriter.Write(p->SegmentType, 21);
            setWriter.Write(p->ReferenceSampleKey);
            setWriter.Write(p->ReferenceSegmentID, 55);
            setWriter.Write(p->FamilialSampleKey);
            setWriter.Write(p->FamilialSegmentID, 55);
        }
    }
}

/**
 * @brief Write the Segments data set data
 * @param affymetrix_calvin_io::GenericFileWriter* - The generic calvin file writer
 * @param CNFamilialAnalysisMethod* - The method to get the data from
 * @param int - The group index
 * @param int - The set index
 */
void CNFamilialReporterFamilial::writeDataSetSegments(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex)
{
    if (pMethod != NULL)
    {
        writer->SeekFromBeginPos(pMethod->getDataSetOffsetVec().at(iSetIndex));
        for (int i = 0; (i < pMethod->getIndexSegments().getCount()); i++)
        {
            CNSegment* p = pMethod->getIndexSegments().getAt(i);
            writeStandardSegment(writer, p, iGroupIndex, iSetIndex);
        }
        for (int i = 0; (i < pMethod->getMotherSegments().getCount()); i++)
        {
            CNSegment* p = pMethod->getMotherSegments().getAt(i);
            writeStandardSegment(writer, p, iGroupIndex, iSetIndex);
        }
        for (int i = 0; (i < pMethod->getFatherSegments().getCount()); i++)
        {
            CNSegment* p = pMethod->getFatherSegments().getAt(i);
            writeStandardSegment(writer, p, iGroupIndex, iSetIndex);
        }
    }
}

/**
 * @brief Write the Segments data set data
 * @param affymetrix_calvin_io::GenericFileWriter* - The generic calvin file writer
 * @param CNSegment* - The segment to write out
 * @param int - The group index
 * @param int - The set index
 */
void CNFamilialReporterFamilial::writeStandardSegment(affymetrix_calvin_io::GenericFileWriter* writer, CNSegment* p, int iGroupIndex, int iSetIndex)
{
    affymetrix_calvin_io::DataSetWriter& setWriter = writer->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    setWriter.Write(p->getSegmentName(), 55);
    setWriter.Write(0);
    setWriter.Write(p->getFamilialSampleKey());
    setWriter.Write((unsigned char)p->getChromosome());
    setWriter.Write((unsigned int)p->getStartPosition());
    setWriter.Write((unsigned int)p->getEndPosition());
    setWriter.Write((unsigned char)p->getCall());
    setWriter.Write((float)p->getConfidence());
    setWriter.Write(p->getMarkerCount());
    setWriter.Write((float)p->getHomozygosity());
    setWriter.Write((float)p->getHeterozygosity());
}
