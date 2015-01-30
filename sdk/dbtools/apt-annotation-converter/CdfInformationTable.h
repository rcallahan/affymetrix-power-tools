////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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

#ifndef CDFINFORMATIONTABLE_H
#define CDFINFORMATIONTABLE_H

#include <map>
#include <string>
#include <algorithm>
using namespace std;
#include "sqlite3.h"

#include "util/BaseEngine.h"

#include "Configuration.h"

class CdfInformationTable
{
    public:
    CdfInformationTable(sqlite3* AnnotationDb, Configuration* ConfigPtr);
    void Update();
    
    private:
    sqlite3* m_annotationDb;     
    Configuration* m_configPtr;
};

#endif
