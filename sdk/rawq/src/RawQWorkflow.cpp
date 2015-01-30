////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

//
#include "rawq/src/RawQWorkflow.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "util/Fs.h"

using namespace affymetrix_fusion_io;

/////////////////////////////////////////////////////////////////////////////

CRawQWorkflow::CRawQWorkflow()
{
}

/////////////////////////////////////////////////////////////////////////////

CRawQWorkflow::~CRawQWorkflow()
{
	Clear();
}

/////////////////////////////////////////////////////////////////////////////

void CRawQWorkflow::Clear()
{
	m_Cel.Clear();
	m_Cdf.Close();
	m_ArrayType = "";
}

/////////////////////////////////////////////////////////////////////////////

bool CRawQWorkflow::ReadCELFile(const char *celFile)
{
	// Read the CEL file
	m_Cel.SetFileName(celFile);
	if (m_Cel.Exists() == false)
	{
		return false;
	}
	try
	{
		if (m_Cel.Read() == false)
		{	
			return false;
		}
	}
	catch(...)
	{
		return false;
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////

bool CRawQWorkflow::LoadCDFFile(const char *libPath)
{
	// Return if already loaded.
	std::string arrayType = StringUtils::ConvertWCSToMBS(m_Cel.GetChipType());
	if (arrayType == m_ArrayType)
		return true;

	// Read the CDF file.
	std::string cdfFile = Fs::join(libPath,arrayType+".CDF");
	m_Cdf.SetFileName(cdfFile.c_str());
	if (m_Cdf.Exists() == false)
	{
		return false;
	}
	if (m_Cdf.Read() == false)
	{
		return false;
	}

	// Store the array type for the next attempt to load the CDF file
	m_ArrayType = arrayType;

	return true;
}

/////////////////////////////////////////////////////////////////////////////

bool CRawQWorkflow::IsExpressionArray()
{
	// Check CDF file for all expression units.
	int iunit;
	int n = m_Cdf.GetHeader().GetNumProbeSets();
	for (iunit=0; iunit<n; iunit++)
	{
		if (m_Cdf.GetProbeSetType(iunit) != affxcdf::ExpressionProbeSetType)
		{
			return false;
		}
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////

CRawQWorkflow::RawQError CRawQWorkflow::ComputeRawQ(const char *celFile, const char *libPath)
{
	// Read the data files and check if expression array type
	if (ReadCELFile(celFile) == false)
		return UnableToReadCELFile;

	if (LoadCDFFile(libPath) == false)
		return UnableToReadCDFFile;

	if (IsExpressionArray() == false)
		return NotExpressionArray;

	// Compute the raw Q value.
	m_RawQ = m_RawQAlg.ComputeRawQ(m_Cel, m_Cdf);
	return NoError;
}

/////////////////////////////////////////////////////////////////////////////
