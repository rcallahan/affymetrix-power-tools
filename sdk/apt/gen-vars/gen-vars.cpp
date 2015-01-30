////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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


#pragma once
#define WIN32_LEAN_AND_MEAN
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <tchar.h>
//

using namespace std;

int main(int argc, char* argv[])
{
	string path = argv[0];
	path.resize(path.length()-17);
	cout << "APT Install Location: " << path << '\n';
	string file = path + "\\bin\\apt-vars.bat";
	cout << "Generating Var File: " << file << '\n';
	ofstream varFile(file.c_str());
	varFile << "@echo off\n";
	varFile << "@set APTROOT=" << path << "\n";
	varFile << "@set PATH=%APTROOT%\\bin;%PATH%\n";
	varFile << "echo Configuring Affymetrix Power Tools (APT) Command Prompt\n";
	varFile << "echo   APT Root: %APTROOT%\n";
	varFile << "echo   APT Version: NON-OFFICIAL-RELEASE\n";
	varFile << "@cd %userprofile%\\Desktop\n";
	varFile << "%comspec% /k\n";
	return 0;
}

