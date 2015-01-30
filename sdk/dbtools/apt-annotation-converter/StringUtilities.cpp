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
#include "StringUtilities.h"

#include "util/Util.h"

    //CSVTOSQLITE-240
    void StringUtilities::FixHeaderVarName(string& VarName)
    {
        if (VarName == "dbsnp-date")
        {
            VarName = "dbsnp_date";
        }
        else if (VarName == "dbsnp-version")
        {
            VarName = "dbsnp_version";
        }
    }
    
    void StringUtilities::CheckColumnName(const string& ColumnName, int ColumnIndex)
    {
         if (ColumnName.empty())
        {
            std::stringstream message;
            message << "Column name at position " << (ColumnIndex + 1) <<  " is empty!";            
            throw Except(message.str().c_str()); 
        }
    }
    string StringUtilities::MakeColumnName(string OriginalName)
    {
        Util::trimString(OriginalName);
        
        if (OriginalName.empty())
        {
           throw Except("Columns name is empty!"); 
        }
        string columnName;
        for (unsigned int i = 0; i < OriginalName.length(); i++)
        {
            char c = OriginalName[i];
            if (('a' <=  c && c <= 'z') || ('A' <=  c && c <= 'Z') || ('0' <=  c && c <= '9') )
            {
                columnName += c;
            }
            else if ((c == ' ') || (c == '.') || (c == '-')  || (c == '/') || (c == '\\'))
            {
                columnName += '_';
            }
        }
        
        if (columnName.length() == 0)
        {
            string message = "Column Name '";
            message += OriginalName;
            message += "' has no basic ASCII characters!";
            throw Except(message.c_str());
        }
        
        while (columnName.find("__") != std::string::npos)
        {
            Util::replaceString(columnName, "__","_");
        }
        
        return columnName;
    }
    
    bool StringUtilities::LooseCompare(string Str1, string Str2)
    {
        return true;
    }

