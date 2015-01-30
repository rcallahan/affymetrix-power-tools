////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////
#include "StdAfx.h"
#include "COMStringUtils.h"
#include <comdef.h>

using namespace std;

COMStringUtils::COMStringUtils()
{
}

BSTR COMStringUtils::ConvertString(const string &str)
{
	_bstr_t b(str.c_str());
	return b.copy();
}

BSTR COMStringUtils::ConvertWString(const wstring &str)
{
	_bstr_t b(str.c_str());
	return b.copy();
}

string COMStringUtils::ConvertString(BSTR str)
{
	_bstr_t bstr(str);
	return string(bstr);
}

wstring COMStringUtils::ConvertWString(BSTR str)
{
	_bstr_t bstr(str);
	return wstring(bstr);
}
