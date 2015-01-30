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

#include <string>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

typedef enum 
{
    Error,
    Warning,
    Information
} SeverityLevel;

class RequiredColumn
{
    public:
    RequiredColumn();
    RequiredColumn(const RequiredColumn &requiredColumn);
          
    SeverityLevel Severity;
    string Message;
    vector<string> Prefixes;
    
    bool m_verified;
    bool m_error;
   
    static const int MaxNumberOfMessages = 20;// Do not put too many errors in the log file
    static int ErrorCounter;
    static int WarningCounter;
    
    static void AddRequiredInEachFile(std::map<string, RequiredColumn> RequiredColumns,  std::vector<string> &RequiredInEachFile);
    static void AddRequiredRecommended(std::map<string, RequiredColumn> RequiredColumns,  std::vector<string> &MandatoryColumns, 
                                                                                          std::vector<string> &RecommendedColumns,
                                                                                          std::vector<string> &InformationColumns);
    
};
