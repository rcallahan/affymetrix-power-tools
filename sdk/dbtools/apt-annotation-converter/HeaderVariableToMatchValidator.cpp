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
#include "HeaderVariableToMatchValidator.h"

HeaderVariableToMatchValidator::HeaderVariableToMatchValidator()
{
    m_status = Ok;
}
void HeaderVariableToMatchValidator::AddVar(string FileName, string VarValue)
{    
    m_fileValueMap.insert(std::make_pair(FileName, VarValue));
}

ValidationResult HeaderVariableToMatchValidator::Validate()
{
    map<std::string, std::string>::iterator valueIterator;
        
    for (valueIterator = m_fileValueMap.begin(); valueIterator != m_fileValueMap.end(); valueIterator++)
    {
        string fileName = valueIterator->first;
        string varValue = valueIterator->second;
        
        if (varValue.empty())
        {
          m_filesWithEmptyValues.push_back(fileName);
          m_status = MissingValue;
        }
        else if (m_value.empty())
        {
            m_value = varValue;
        }
        else if (m_value != varValue)
        {
            m_status = ValueMismatch;
            break;            
        }
    }    
    return m_status;
}
void HeaderVariableToMatchValidator::CopyValueToFilesWithEmptyValues()
{
    map<std::string, std::string>::iterator valueIterator;
    
    if (m_status != ValueMismatch && false == m_value.empty())
    {
        for (valueIterator = m_fileValueMap.begin(); valueIterator != m_fileValueMap.end(); valueIterator++)
        {
            valueIterator->second = m_value;
        }
    }
}
