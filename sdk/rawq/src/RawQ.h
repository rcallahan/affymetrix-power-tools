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

/*! \file RawQ.h The file contains classes to compute the RawQ value. */

#ifndef _RAWQ_HEADER_
#define _RAWQ_HEADER_

/////////////////////////////////////////////////////////////////////////////

#include "calvin_files/fusion/src/FusionCDFData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
//

/////////////////////////////////////////////////////////////////////////////

/*! A class to compute the RawQ value of an expression array. */
class CRawQ
{
public:
	/*! Constructor */
	CRawQ();

	/*! Destructor */
	~CRawQ();

private:
	/*! The number of vertical zones for the background calculation. */
	int m_VertZones;

	/*! The number of horizontal zones for the background calculation. */
	int m_HorzZones;

	/*! The percentage of probes in a zone to use for the background calculation. */
	float m_PercentBGCells;

public:
	/*! Sets the default values. */
	void SetDefaults();

	/*! Sets the number of vertical zones for the background calculation.
	 * @param v The number of vertical zones for the background calculation.
	 */
	void SetVerticalZones(int v) { m_VertZones = v; }

	/*! Gets the number of vertical zones for the background calculation.
	 * @return The number of vertical zones for the background calculation.
	 */
	int GetVerticalZones() const { return m_VertZones; }

	/*! Sets the number of horizontal zones for the background calculation.
	 * @param h The number of horizontal zones for the background calculation.
	 */
	void SetHorizontalZones(int h) { m_HorzZones = h; }

	/*! Gets the number of horizontal zones for the background calculation.
	 * @return The number of horizontal zones for the background calculation.
	 */
	int GetHorizontalZones() const { return m_HorzZones; }

	/*! Sets the percentage of probes in a zone to use for the background calculation.
	 * @param p The percentage of probes in a zone to use for the background calculation.
	 */
	void SetPercentBG(float p) { m_PercentBGCells = p; }

	/*! Gets the percentage of probes in a zone to use for the background calculation.
	 * @return The percentage of probes in a zone to use for the background calculation.
	 */
	float GetPercentBG() const { return m_PercentBGCells; }

	/*! Computes the RawQ value. The CDF file must be of an expression array.
	 * @param cell The CEL file intensities.
	 * @param cdf The CDF file data.
	 * @return The RawQ value.
	 */
	float ComputeRawQ(affymetrix_fusion_io::FusionCELData &cell, affymetrix_fusion_io::FusionCDFData &cdf);
};

/////////////////////////////////////////////////////////////////////////////

#endif
