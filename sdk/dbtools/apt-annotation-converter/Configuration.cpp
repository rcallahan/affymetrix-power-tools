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
#include "Configuration.h"
#include "ConfigFileHandler.h"
#include "ResourceStrings.h"

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/util/PlatformUtils.hpp>

XERCES_CPP_NAMESPACE_USE;

#include "util/Util.h"
#include "util/Fs.h"
#include "util/Verbose.h"


static const char* XERCES_XML_INITIALIZATION_ERROR = "XERCES XML initialization error: ";
static const char* UNABLE_TO_PARSE_CONFIGURATION_FILE = "Unable to parse configuration file: ";

void Configuration::LoadConfigurationFile(string FileName)
{
    m_par1CsvColumnName = "ChrX pseudo-autosomal region 1";
    m_par2CsvColumnName = "ChrX pseudo-autosomal region 2";
    m_parDbColumnName = "ChrX_PAR";
    m_parDbColumnType = DB_TYPE_INTEGER;
    std::string err_message = 	UNABLE_TO_PARSE_CONFIGURATION_FILE; 	
    err_message += FileName.c_str();
    err_message += ". ";
    
    try
	{
		XMLPlatformUtils::Initialize();
	}

	catch (const XERCES_CPP_NAMESPACE::XMLException& xmlException)
	{	
	    std::string message = 	XERCES_XML_INITIALIZATION_ERROR; 	
	    message += ConfigFileHandler::ToStdString(xmlException.getMessage());
	    
        throw Except(message.c_str());					
	}
	
	//
	//  Create a SAX parser object to use and create our SAX event handlers
	//  and plug them in.
	//
	SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
	ConfigFileHandler* const configFileHandler = new ConfigFileHandler(this);
	parser->setContentHandler(configFileHandler);
	parser->setErrorHandler(configFileHandler);
	parser->setFeature(XMLUni::fgSAX2CoreValidation, false);
    parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, false);
		    
    try
	{
		parser->parse(FileName.c_str());
		int errorCount = parser->getErrorCount();
		if (errorCount != 0)
		{	
		    	throw Except(err_message.c_str());
		}
	}
	catch (SAXParseException& parseException)
	{
		std::stringstream message; 	
		message << err_message.c_str();
	    message << ConfigFileHandler::ToStdString(parseException.getMessage());
	    message << ". Line " << parseException.getLineNumber() << ".";
	    	    
        throw Except(message.str().c_str());			
	}
	
	catch (std::exception& exception)
	{
	    err_message += exception.what();
	    
        throw Except(err_message.c_str());			
	}
	catch (...)
	{	
            err_message += "Unknown exception.";
	    throw Except(err_message.c_str());	
	}
	delete configFileHandler;
	delete parser;
	XMLPlatformUtils::Terminate();
}

bool Configuration::IsCloneableColumn(const string& ColumnName)
{
    ColumnsCloningMap::iterator cloneableColumn  = m_columnsCloning.find(ColumnName);
    if (cloneableColumn != m_columnsCloning.end() && cloneableColumn->second.size() > 0)
    {        
      return true;                
    }
    
    return false;
}

std::vector<string> Configuration::GetColumnClones(const string& ColumnName)
{
    ColumnsCloningMap::iterator cloneableColumn  = m_columnsCloning.find(ColumnName);
    if (cloneableColumn != m_columnsCloning.end())
    {
        return cloneableColumn->second;
    }
    vector<string> foo;
    return foo;
}    

void Configuration::ConvertTsvToDbValue(const string& ColumnName, string& ValueToConvert)
{
    CsvToDbConversionMap::iterator columnsIterator;
    columnsIterator = m_csvToDbConversionMap.find(ColumnName);
    if (columnsIterator != m_csvToDbConversionMap.end())
    {
       map<string, string> csvValueDbValueMap = columnsIterator->second;// map<TsvValue, DbValue >
       map<string, string>::iterator tsvToDbIterator;
       tsvToDbIterator = csvValueDbValueMap.find(ValueToConvert);
       if (tsvToDbIterator != csvValueDbValueMap.end())
       {
            ValueToConvert = tsvToDbIterator->second;
       }
    }
}   

void Configuration::Clear()
{
    CsvToDbConversionMap::iterator columnConversions;
    for (columnConversions = m_csvToDbConversionMap.begin(); columnConversions != m_csvToDbConversionMap.end(); columnConversions++) 
    {
        columnConversions->second.clear();        
    }
    
    m_csvToDbConversionMap.clear();
}

Configuration::~Configuration()
{
    Clear();
}
