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
#include "ConfigFileHandler.h"
#include "Configuration.h"
#include "ResourceStrings.h"
#include "StringUtilities.h"

#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/XMLString.hpp>

#include <sstream>


static const char* NODE_NAME_CONVERSION_TABLE = "ConversionTable";
static const char* NODE_NAME_CONVERSION = "Conversion";
static const char* NODE_NAME_COLUMNS_CLONING = "ColumnsCloning";
static const char* NODE_NAME_DBCLONE = "DbClone";
static const char* NODE_NAME_CSVCOLUMN = "CsvColumn";
static const char* NODE_NAME_COLUMN_DEFINITIONS = "ColumnDefinitions";
static const char* NODE_NAME_COLUMN = "Column";
static const char* NODE_NAME_ALTERNATIVE_COLUMN_MAPPINGS = "AlternativeColumnMappings";
static const char* NODE_NAME_ALTERNATIVE_COLUMN_MAPPING = "AlternativeColumnMapping";
static const char* NODE_NAME_CDF_INFORMATION_TABLE_ROW = "CdfInformationTableRow";
static const char* NODE_NAME_CDF_INFORMATION_TABLE = "CdfInformationTable";
static const char* NODE_NAME_PREFIX =  "Prefix";

static const char* HEADER_VARIABLE_TO_MATCH = "HeaderVariableToMatch";
static const char* HEADER_VARIABLES_TO_MATCH = "HeaderVariablesToMatch";
static const char* RESERVED_HEADER_VARIABLE = "ReservedHeaderVariable";
static const char* RESERVED_HEADER_VARIABLES = "ReservedHeaderVariables";


// Note: This relies on w_char being implemented as 16bit ints - NOT PORTABLE
static const char* ATTRIB_NAME_COLUMN_NAME = "columnName";
static const char* ATTRIB_NAME_CSVVALUE = "csvValue";
static const char* ATTRIB_NAME_DBVALUE =  "dbValue";

static const char* ATTRIB_NAME_REQUIREMENT_TYPE = "requirementType";
static const char* ATTRIB_NAME_SEVERITY =  "severity";
static const char* ATTRIB_NAME_MESSAGE =  "message";
static const char* ATTRIB_NAME_PREFIX =  "prefix";

static const char* ATTRIB_NAME_NAME = "name";
static const char* ATTRIB_NAME_VALUE = "value";
static const char* ATTRIB_NAME_CSV_COLUMN_NAME = "csvColumnName";
static const char* ATTRIB_NAME_DB_COLUMN_NAME = "dbColumnName";
static const char* ATTRIB_NAME_DELETE = "delete";
static const char* ATTRIB_NAME_COLUMN_TYPE = "columnType";

ConfigFileHandler::ConfigFileHandler(Configuration* Config)
{    
    m_inCloning = false;    
    m_configuration = Config;
    
    m_cdfInformationTableRows.push_back("filename");
    m_cdfInformationTableRows.push_back("guid");
    m_cdfInformationTableRows.push_back("checksum");
    m_cdfInformationTableRows.push_back("version");

}

std::string ConfigFileHandler::GetAttrib(const Attributes& Attrs, const char* AttribName)
{
    return GetAttrib(Attrs, AttribName, true);
}
std::string ConfigFileHandler::GetAttrib(const Attributes& Attrs, const char* AttribName, bool ReturnEmptyStringIfNull)
{
    // Convert char* to XMLCh*
    XMLCh* xAttribName = XMLString::transcode(AttribName);
    // get XMLCh* result
    const XMLCh* valuePtr = Attrs.getValue(xAttribName);
    // release the attr name
    XMLString::release(&xAttribName);
    // check for missing attr
     if (NULL == valuePtr)
     {
        if (ReturnEmptyStringIfNull)
        {
            return "";
        }
        else
        {
            return NULL;            
        }
     }          
     return  ToStdString(valuePtr);     
}

bool ConfigFileHandler::IsDeleteAtribSet(const Attributes& Attrs)
{     
     if (GetAttrib(Attrs, ATTRIB_NAME_DELETE) == "yes")
     {
        return true;
     }
     return false;
}
void ConfigFileHandler::AddOrRemoveRequiredColumn(std::map<string, RequiredColumn> &ColumnDefinitions, const Attributes& Attrs)
{
    bool deleteValue = IsDeleteAtribSet(Attrs);
    m_lastColumnName = GetAttrib(Attrs, ATTRIB_NAME_COLUMN_NAME);
           
    std::map<string, RequiredColumn>::iterator found = ColumnDefinitions.find(m_lastColumnName);
        
    if (deleteValue)
    {
        if (found != ColumnDefinitions.end())
        {
            ColumnDefinitions.erase(found);
        }
    }
    else
    {// add or update           
        SeverityLevel severity = Error;
        
        string severityAttribValue = GetAttrib(Attrs, ATTRIB_NAME_SEVERITY);
        if (severityAttribValue == "warning")
        {
            severity = Warning;
        }
        else if (severityAttribValue == "information")
        {
            severity = Information;
        }                
        
        if (found == ColumnDefinitions.end())
        {// Add  
            RequiredColumn requiredColumn;
            requiredColumn.Message = GetAttrib(Attrs, ATTRIB_NAME_MESSAGE);
            requiredColumn.Severity = severity;
            ColumnDefinitions.insert(std::make_pair(m_lastColumnName,requiredColumn));
        }
        else
        {// update
            found->second.Message = GetAttrib(Attrs, ATTRIB_NAME_MESSAGE);            
            found->second.Severity = severity;
        }
    }     
}
void ConfigFileHandler::AddOrRemoveValue(std::vector<string> &Vertor, const Attributes& Attrs, const char* AttribName, bool FixHeaderVar)
{     
     
     bool deleteValue = IsDeleteAtribSet(Attrs);
     
     // Get AttribName Value     
     XMLCh* xAttribName = XMLString::transcode(AttribName);
     const XMLCh* valuePtr = Attrs.getValue(xAttribName);
     // release the attr name
     XMLString::release(&xAttribName);

     string value;
          
     if (valuePtr != NULL)
     {
        value = ToStdString(valuePtr);
        std::vector<string>::iterator found =  find(Vertor.begin(),Vertor.end(), value);
            
        if (deleteValue)
        {
            if (found != Vertor.end())
            {
                Vertor.erase(found);
            }
        }
        else
        {   
            if (found == Vertor.end())
            {
                if (FixHeaderVar)
                {
                    StringUtilities::FixHeaderVarName(value);
                }
                Vertor.push_back(value);
            }
        }
     }
     //          
}
void ConfigFileHandler::startElement(const XMLCh* const uri, const XMLCh* const localname,
                                     const XMLCh* const qname, const Attributes& attrs)
{
    std::string nodeName = ToStdString(qname);
    
    if (nodeName == NODE_NAME_CONVERSION)
    {
        std::string columnName = GetAttrib(attrs, ATTRIB_NAME_COLUMN_NAME);
        std::string csvValue = GetAttrib(attrs, ATTRIB_NAME_CSVVALUE);
        std::string dbValue = GetAttrib(attrs, ATTRIB_NAME_DBVALUE);
        
        CsvToDbConversionMap::iterator pos = m_configuration->m_csvToDbConversionMap.find(columnName);
        if ( pos == m_configuration->m_csvToDbConversionMap.end())
        {            
            // Add new map for this column
            map<string,string> conversionMap;
            conversionMap.insert(make_pair(csvValue,dbValue));
            m_configuration->m_csvToDbConversionMap.insert(make_pair(columnName, conversionMap));            
        }
        else
        {
            if (IsDeleteAtribSet(attrs))
            {
                m_configuration->m_csvToDbConversionMap.erase(pos);
            }
            else
            {
                m_configuration->m_csvToDbConversionMap[columnName].insert(make_pair(csvValue,dbValue));
            }
        }
    }
    else if (nodeName == NODE_NAME_CONVERSION_TABLE)          
    {
        if (IsDeleteAtribSet(attrs))
        {
           m_configuration->m_csvToDbConversionMap.clear();
        }        
    }
    else if (nodeName == HEADER_VARIABLES_TO_MATCH)
    {
        if (IsDeleteAtribSet(attrs))
        {
            m_configuration->m_varsToMatch.clear();
        }
    }
    else if (nodeName == HEADER_VARIABLE_TO_MATCH)
    {
        AddOrRemoveValue(m_configuration->m_varsToMatch, attrs, ATTRIB_NAME_NAME, true);        
    }    
    else if (nodeName == RESERVED_HEADER_VARIABLES)
    {
        if (IsDeleteAtribSet(attrs))
        {
            m_configuration->m_reservedVars.clear();
        }
    }
    else if (nodeName == RESERVED_HEADER_VARIABLE)
    {
        AddOrRemoveValue(m_configuration->m_reservedVars, attrs, ATTRIB_NAME_NAME, false);        
    }  
    else if (nodeName == NODE_NAME_COLUMNS_CLONING)
    {
        m_inCloning = true;
        if (IsDeleteAtribSet(attrs))
        {
            m_configuration->m_columnsCloning.clear();
        }
    }    
    else if (nodeName == NODE_NAME_CSVCOLUMN && m_inCloning)
    {        	    
	    m_columnNameToClone = GetAttrib(attrs, ATTRIB_NAME_COLUMN_NAME);
	    if (IsDeleteAtribSet(attrs))
	    {
	        ColumnsCloningMap::iterator columnsIterator = m_configuration->m_columnsCloning.find(m_columnNameToClone);
	        if (columnsIterator != m_configuration->m_columnsCloning.end())
	        {
	            m_configuration->m_columnsCloning.erase(columnsIterator); 	        
	        }   	       
	    }
    }
    else if (nodeName == NODE_NAME_DBCLONE && m_inCloning)
    {               
        string cloneName = GetAttrib(attrs, ATTRIB_NAME_COLUMN_NAME);
        string cloneType = GetAttrib(attrs, ATTRIB_NAME_COLUMN_TYPE);
        if (cloneType.empty())
        {
            cloneType = DB_TYPE_TEXT;
        }
            
        std::vector<string> columns; 	    
	    if (false == m_configuration->IsCloneableColumn(m_columnNameToClone)) // Check if it wasn't added
	    {	        
	        m_configuration->m_columnsCloning.insert(std::make_pair(m_columnNameToClone , columns));
	    }
	    
	    ColumnsCloningMap::iterator columnsIterator = m_configuration->m_columnsCloning.find(m_columnNameToClone);
	    if (columnsIterator != m_configuration->m_columnsCloning.end())
	    {	                    
            columnsIterator->second.push_back(cloneName);                        
	    }	    
	    
	    map<string, string>::iterator columnTypes = m_configuration->m_cloneColumnTypes.find(cloneName);        
        if (columnTypes == m_configuration->m_cloneColumnTypes.end())
        {            
            m_configuration->m_cloneColumnTypes.insert(std::make_pair(cloneName, cloneType));
        }
            
    }    
    else if (nodeName == NODE_NAME_COLUMN_DEFINITIONS)          
    {
        if (IsDeleteAtribSet(attrs))
        {
            m_configuration->m_requiredColumnsPresence.clear();
            m_configuration->m_requiredColumnsSomeValue.clear();
            m_configuration->m_requiredColumnsNoEmptyValues.clear();
            m_configuration->m_requiredColumnsPrefixCheck.clear();
        }        
    }
    else if (nodeName == NODE_NAME_COLUMN)          
    {
        
        string columnType = GetAttrib(attrs, ATTRIB_NAME_COLUMN_TYPE);
        string columnName = GetAttrib(attrs, ATTRIB_NAME_COLUMN_NAME);
        if (columnName.empty())
        {
            throw Except("One of the 'Column' nodes has missing or empty 'columnName' attribute!");
        }
        if (false == columnType.empty())
        {
            if (false == (DB_TYPE_TEXT == columnType || DB_TYPE_REAL == columnType || DB_TYPE_INTEGER == columnType))
            {
                std::stringstream message;
                message << " Column Type '" << columnType << "' is incorrect! ";
                message << "Must be '" << DB_TYPE_REAL << "', '" << DB_TYPE_INTEGER << "' or '" << DB_TYPE_TEXT << "'.";
                throw Except(message.str().c_str()); 
            }
            
            string columnName = GetAttrib(attrs, ATTRIB_NAME_COLUMN_NAME);
            if (false == columnName.empty())
            {
                m_configuration->m_columnTypes.insert(std::make_pair(columnName,columnType));
            }   
        }
        string requirementTypeAttribValue = GetAttrib(attrs, ATTRIB_NAME_REQUIREMENT_TYPE);
        
        if (requirementTypeAttribValue == "Presence")
        {
            AddOrRemoveRequiredColumn(m_configuration->m_requiredColumnsPresence, attrs);
        }
        else if (requirementTypeAttribValue == "SomeValue")
        {
            AddOrRemoveRequiredColumn(m_configuration->m_requiredColumnsSomeValue, attrs);
        }
        else if (requirementTypeAttribValue == "NoEmptyValues")
        {
            AddOrRemoveRequiredColumn(m_configuration->m_requiredColumnsNoEmptyValues, attrs);            
        }
        else if (requirementTypeAttribValue == "PrefixCheck")
        {
            AddOrRemoveRequiredColumn(m_configuration->m_requiredColumnsPrefixCheck, attrs);            
        }
    }    
    else if (nodeName == NODE_NAME_PREFIX)
    {
        if (m_configuration->m_requiredColumnsPrefixCheck.size() > 0)
        {
            m_configuration->m_requiredColumnsPrefixCheck[m_lastColumnName].Prefixes.push_back(GetAttrib(attrs, ATTRIB_NAME_VALUE));
        }
    }
    else if (nodeName == NODE_NAME_ALTERNATIVE_COLUMN_MAPPINGS)
    {
        if (IsDeleteAtribSet(attrs))
        {
            m_configuration->m_alternativeColumnMappings.clear();
        }
    }    
    else if (nodeName == NODE_NAME_ALTERNATIVE_COLUMN_MAPPING)          
    {                
        string csvColumnName = GetAttrib(attrs, ATTRIB_NAME_CSV_COLUMN_NAME);
        string dbColumnName = GetAttrib(attrs, ATTRIB_NAME_DB_COLUMN_NAME);
        
        map<string, string>::iterator pos = m_configuration->m_alternativeColumnMappings.find(csvColumnName);
        if (pos != m_configuration->m_alternativeColumnMappings.end())
        { // found
            if (IsDeleteAtribSet(attrs))
            {
                m_configuration->m_alternativeColumnMappings.erase(pos);
            }
            else
            { // Update
                pos->second = dbColumnName;
            }
        }
        else
        {
            if (false == IsDeleteAtribSet(attrs))
            {
                m_configuration->m_alternativeColumnMappings.insert(std::make_pair(csvColumnName, dbColumnName));
            }
        }
    }    
    else if (nodeName == NODE_NAME_CDF_INFORMATION_TABLE)
    {
         if (IsDeleteAtribSet(attrs))
         {
            m_configuration->m_cdfInformationTable.clear();
         }
    }
    else if (nodeName == NODE_NAME_CDF_INFORMATION_TABLE_ROW)
    {   
        map<string, string> cdfInformationTableRow;
        for(unsigned int i = 0; i < m_cdfInformationTableRows.size(); i++)
        {
            std::string a = m_cdfInformationTableRows[i];
            cdfInformationTableRow.insert(std::make_pair(a, GetAttrib(attrs, m_cdfInformationTableRows[i])));
        }  
                 
        m_configuration->m_cdfInformationTable.push_back(cdfInformationTableRow);
    }    
    
};
 
 void ConfigFileHandler::endElement (const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname)
 {
    std::string nodeName = ToStdString(qname);
    if (nodeName == NODE_NAME_COLUMNS_CLONING)
    {
        m_inCloning = false;
    }
 };
	


std::string ConfigFileHandler::ToStdString(const XMLCh* const in) 
{
    char* p = XMLString::transcode(in);
    std::string str = p;
    XMLString::release(&p);
    return str;
}
  
