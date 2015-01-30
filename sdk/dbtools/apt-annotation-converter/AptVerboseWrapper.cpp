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
#include "AptVerboseWrapper.h"

#include "util/Verbose.h"

int AptVerboseWrapper::Delay = 0;

void AptVerboseWrapper::out(int level, const std::string &s, bool nl)
{
    Verbose::out(level, s, nl);
    if (Delay > 0)
    {
//        Sleep(Delay);
    }
}

void AptVerboseWrapper::progressStep(int verbosity)
{
    Verbose::progressStep(verbosity);
    if (Delay > 0)
    {
//        Sleep(Delay);
    }
}

void AptVerboseWrapper::progressEnd(int verbosity, const std::string &msg)
{
    if (Delay > 0)
    {
//        Sleep(Delay);
    }
    
    Verbose::progressEnd(verbosity, msg);
    
    if (Delay > 0)
    {
//        Sleep(Delay);
    }
}

void AptVerboseWrapper::progressBegin(int verbosity, const std::string &msg, int total, int dotMod, int maxCalls)
{     
    Verbose::progressBegin(verbosity, msg, total, dotMod, maxCalls);    
    if (Delay > 0)
    {
//        Sleep(Delay);
    }
}
