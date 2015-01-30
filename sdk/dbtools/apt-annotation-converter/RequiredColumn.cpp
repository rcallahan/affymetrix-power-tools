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
#include "RequiredColumn.h"

RequiredColumn::RequiredColumn()
{
    Severity = Error;    
    m_verified = false;
    m_error = false;    
}
RequiredColumn::RequiredColumn(const RequiredColumn &requiredColumn)
{
    Severity = requiredColumn.Severity;
    Message = requiredColumn.Message;    
    m_verified = requiredColumn.m_verified;
    m_error = requiredColumn.m_error;
}

void RequiredColumn::AddRequiredInEachFile(std::map<string, RequiredColumn> RequiredColumns,  std::vector<string> &RequiredInEachFile)
{
    std::map<string, RequiredColumn>::iterator requiredColumnIterator;
    for (requiredColumnIterator = RequiredColumns.begin(); requiredColumnIterator != RequiredColumns.end(); requiredColumnIterator++)
    {            
        if (requiredColumnIterator->second.Severity == Error)
        {
            std::vector<string>::iterator columnToFind;
            columnToFind = find(RequiredInEachFile.begin(), RequiredInEachFile.end(), requiredColumnIterator->first);
            if (columnToFind == RequiredInEachFile.end())
            {
                  RequiredInEachFile.push_back(requiredColumnIterator->first/*Column Name*/);
            }
        }
    }
}

void RequiredColumn::AddRequiredRecommended(std::map<string, RequiredColumn> RequiredColumns,  std::vector<string> &MandatoryColumns, 
                                                                                          std::vector<string> &RecommendedColumns,
                                                                                          std::vector<string> &InformationColumns)
{
    std::map<string, RequiredColumn>::iterator requiredColumnIterator;
    for (requiredColumnIterator = RequiredColumns.begin(); requiredColumnIterator != RequiredColumns.end(); requiredColumnIterator++)
    {
       std::vector<string>::iterator columnToFind;
       switch(requiredColumnIterator->second.Severity)
       {
            case Error:
            columnToFind = find(MandatoryColumns.begin(),MandatoryColumns.end(), requiredColumnIterator->first);
            if (columnToFind == MandatoryColumns.end())
            {
                  MandatoryColumns.push_back(requiredColumnIterator->first/*Column Name*/);
            }                
            break;
            case Warning:
            columnToFind = find(RecommendedColumns.begin(),RecommendedColumns.end(), requiredColumnIterator->first);
            if (columnToFind == RecommendedColumns.end())
            {
                  RecommendedColumns.push_back(requiredColumnIterator->first/*Column Name*/);
            }                
            break;
            case Information:
            columnToFind = find(InformationColumns.begin(),InformationColumns.end(), requiredColumnIterator->first);
            if (columnToFind == InformationColumns.end())
            {
                  InformationColumns.push_back(requiredColumnIterator->first/*Column Name*/);
            }                 
            break;
       }       
    }
}

int RequiredColumn::ErrorCounter = 0;
int RequiredColumn::WarningCounter = 0;
