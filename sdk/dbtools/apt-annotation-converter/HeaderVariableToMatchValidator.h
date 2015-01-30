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
#ifndef HEADERVARIABLETOMATCHVALIDATOR_H
#define HEADERVARIABLETOMATCHVALIDATOR_H

#include <string>
#include <map>
#include <vector>

using namespace std;

typedef enum {Ok, MissingValue/*Waring*/, ValueMismatch/*Error*/} ValidationResult;

class HeaderVariableToMatchValidator
{
    public:
    HeaderVariableToMatchValidator();
    void AddVar(string FileName, string VarValue);
    ValidationResult GetStatus() { return m_status;};
    map<std::string, std::string> GetFileValueMap() { return m_fileValueMap;};
    std::vector<string> GetFilesWithEmptyValues() {return m_filesWithEmptyValues;};
    string GetValue(){return m_value;};
    
    ValidationResult Validate();
    void CopyValueToFilesWithEmptyValues();
    private:
    string m_value;
    map<std::string, std::string> m_fileValueMap; // <File Name, Var Value>  
    std::vector<string> m_filesWithEmptyValues;  
    ValidationResult m_status;    
};

#endif
