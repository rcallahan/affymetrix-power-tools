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

/*! \file RawQWorkflow.h This file provides the ability to compute a RawQ value from a CEL file. */

#ifndef _RAWQ_WORKFLOW_HEADER_
#define _RAWQ_WORKFLOW_HEADER_

/////////////////////////////////////////////////////////////////////////////

#include "rawq/src/RawQ.h"
//
#include "calvin_files/fusion/src/FusionCDFData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
//

/////////////////////////////////////////////////////////////////////////////

/*! This class provides interfaces to compute the RawQ value from a CEL file. */
class CRawQWorkflow
{
public:
	/*! The error codes for the RawQ function. */
	typedef enum {
		NoError,
		UnableToReadCELFile,
		UnableToReadCDFFile,
		NotExpressionArray
	} RawQError;

protected:
	/*! The probe array type. */
	std::string m_ArrayType;

	/*! The CDF file object. */
	affymetrix_fusion_io::FusionCDFData m_Cdf;

	/*! The CEL file object. */
	affymetrix_fusion_io::FusionCELData m_Cel;

	/*! The computed RawQ value. */
	float m_RawQ;

	/*! The RawQ algorithm interfaces. */
	CRawQ m_RawQAlg;

	/*! Reads the CEL file into memory.
	 * @param celFile The name of the CEL file.
	 * @return True if sucessful.
	 */
	bool ReadCELFile(const char *celFile);

	/*! Reads the CDF file into memory.
	 * @param libPath The path to the CDF files.
	 * @return True if sucessful.
	 */
	bool LoadCDFFile(const char *libPath);

	/*! Checks that the files are of the right tile type.
	 * @return True if the array type is an expression array type.
	 */
	bool IsExpressionArray();

public:
	/*! Constructor */
	CRawQWorkflow();

	/*! Destructor */
	~CRawQWorkflow();

	/*! Resets the default parameters. */
	void SetDefaults() { m_RawQAlg.SetDefaults(); }

	/*! Sets the number of vertical zones to use.
	 * @param v The number of vertical zones for the background calculation.
	 */
	void SetVerticalZones(int v) { m_RawQAlg.SetVerticalZones(v); }

	/*! Sets the number of horizontal zones to use.
	 * @param h The number of horizontal zones for the background calculation.
	 */
	void SetHorizontalZones(int h) { m_RawQAlg.SetHorizontalZones(h); }

	/*! Sets the percentage of background probes to use in each zone.
	 * @param p The percent of probes to use for the background calculation.
	 */
	void SetPercentBG(float p) { m_RawQAlg.SetPercentBG(p); }

	/*! Gets the number of vertical zones to use.
	 * @return The number of vertical zones for the background calculation.
	 */
	int GetVerticalZones() const { return  m_RawQAlg.GetVerticalZones(); }

	/*! Gets the number of horizontal zones to use.
	 * @return The number of horizontal zones for the background calculation.
	 */
	int GetHorizontalZones() const { return  m_RawQAlg.GetHorizontalZones(); }

	/*! Gets the percentage of background probes to use in each zone.
	 * @return The percent of probes to use for the background calculation.
	 */
	float GetPercentBG() const { return  m_RawQAlg.GetPercentBG(); }

	/*! Clear memory associated with the class. */
	void Clear();

	/*! Calculates the RawQ value.
	 * @param celFile The full path of the CEL file.
	 * @param libPath The path to where the CDF files live.
	 * @return The status of the calculation.
	 */
	RawQError ComputeRawQ(const char *celFile, const char *libPath);

	/*! Gets the RawQ value.
	 * @return The RawQ value.
	 */
	float GetRawQ() { return m_RawQ; }
};

/////////////////////////////////////////////////////////////////////////////

#endif
