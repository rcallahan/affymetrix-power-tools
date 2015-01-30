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

#include "AnnotationTsv.h"

AnnotationTsv::AnnotationTsv()
{
    _myColumnCount = -1;
    m_optCheckFormatOnOpen=false;
}
/////CSVTOSQLITE-249 to fix the problem in the CSV file: It has the comma at the end of the column headers line:
                          // ,"Fragment Enzyme Type Length Start Stop",
int AnnotationTsv::GetColumnCount()
{
    if (_myColumnCount == -1)
    {
        _myColumnCount = getColumnCount(0);
        if (_myColumnCount > 0)
        {
            while (_myColumnCount > 0 && this->getColumnName(0, _myColumnCount -1).empty())
            {
                _myColumnCount--;
            }
        }        
    }
    return _myColumnCount;
}
