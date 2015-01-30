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
#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <map>
#include <string>
#include <vector>
using namespace std;

#include "RequiredColumn.h"

typedef map<string, map<string,string > > CsvToDbConversionMap; // map<ColumnName, map<TsvValue, DbValue > >
typedef map<string, vector<string> > ColumnsCloningMap; // map<CSV ColumnName, List of columns in DB>


class Configuration
{
 friend class ConfigFileHandler;
public:
    void LoadConfigurationFile(string FileName);    
    virtual ~Configuration();
    void ConvertTsvToDbValue(const string& ColumnName, string& ValueToConvert);
    
    /// CSVTOSQLITE-74
    bool IsCloneableColumn(const string& ColumnName);    
    std::vector<string> GetColumnClones(const string& ColumnName);
    ///
    
    /// Vars for the second, hardcoded part of the CSVTOSQLITE-74
    string m_par1CsvColumnName;
    string m_par2CsvColumnName;
    string m_parDbColumnName;                
    string m_parDbColumnType;
    //
        
          
    std::vector<string> m_varsToMatch;
    std::vector<string> m_reservedVars;
    
    map<string, string> m_alternativeColumnMappings; // <csvColumnName, dbColumnName>
    
    std::map<string, RequiredColumn> m_requiredColumnsPresence; // <csvColumnName, Column object>
    std::map<string, RequiredColumn> m_requiredColumnsSomeValue; // <csvColumnName, Column object>
    std::map<string, RequiredColumn> m_requiredColumnsNoEmptyValues; // <csvColumnName, Column object>
    std::map<string, RequiredColumn> m_requiredColumnsPrefixCheck; // <csvColumnName, Column object>
    std::map<string, string> m_columnTypes; // <csvColumnName, Column Type>
    std::map<string, string> m_cloneColumnTypes; // <dbColumnName, Column Type>
    
    std::vector<map<string, string> > m_cdfInformationTable; // vector<map<Column Name,Column Value> >
    private:
    CsvToDbConversionMap m_csvToDbConversionMap;
    ColumnsCloningMap m_columnsCloning;
    
    void Clear();    
};

#endif
