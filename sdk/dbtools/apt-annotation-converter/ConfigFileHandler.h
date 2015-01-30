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
#ifndef CONFIGFILEHANDLER_H
#define CONFIGFILEHANDLER_H

#define XML_LIBRARY
#define XERCES_STATIC_LIBRARY
#include <xercesc/sax2/DefaultHandler.hpp>
XERCES_CPP_NAMESPACE_USE;

#include <string>
#include <algorithm>

#include "Configuration.h"
#include "../../util/Except.h"

class ConfigFileHandler : public DefaultHandler
{
	public:
	ConfigFileHandler(Configuration* Config);
	
    virtual void startElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs);
    virtual void endElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname);
	
    
    // Utils
    static std::string ToStdString(const XMLCh* const in);
    
    private:   
    void AddOrRemoveRequiredColumn(std::map<string, RequiredColumn> &ColumnDefinitions, const Attributes& Attrs);
    void AddOrRemoveValue(std::vector<string> &Vertor, const Attributes& Attrs, const char* AttribName, bool FixHeaderVar);
    bool IsDeleteAtribSet(const Attributes& Attrs);
    std::string GetAttrib(const Attributes& Attrs, const char* AttribName);
    std::string GetAttrib(const Attributes& Attrs, const char* AttribName, bool ReturnEmptyStringIfNull);
    string m_lastColumnName;
    
    // Members
    std::vector<const char*> m_cdfInformationTableRows;    
    bool m_inCloning;
    std::string m_columnNameToClone;    
    Configuration* m_configuration;
    
    
};

#endif
