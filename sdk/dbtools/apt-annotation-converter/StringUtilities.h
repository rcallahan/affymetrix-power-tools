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
using namespace std;

class StringUtilities
{
    public:    
    static string MakeColumnName(string OriginalName);// Makes name suitable for SQLite Column     
    static void CheckColumnName(const string& ColumnName, int ColumnIndex);
    static void FixHeaderVarName(string& VarName);// CSVTOSQLITE-240
    static bool LooseCompare(string Str1, string Str2);
};
