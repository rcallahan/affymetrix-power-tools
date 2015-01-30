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

#ifndef _ProbeSetStats_HEADER_
#define _ProbeSetStats_HEADER_

/*! \file ProbeSetStats.h Defines a class to store statistics about the probe set results. */

namespace ExpressionReport
{

/*! Stores statistics about the probe set results. */
class ProbeSetStats
{
protected:
	/*! Present call stats. */
	DetectionStats present;

	/*! Absent call stats. */
	DetectionStats absent;

	/*! Marginal call stats. */
	DetectionStats marginal;

	/*! Total number of probe sets. */
	int numSets;

public:
	/*! Add one to the probe set count. */
	void AddSet() { ++numSets; }

	/*! The total number of probe sets. */
	int &NumSets() { return numSets; }

	/*! Present call stats. */
	DetectionStats &PresentCalls() { return present; }

	/*! Absent call stats. */
	DetectionStats &AbsentCalls() { return absent; }

	/*! Marginal call stats. */
	DetectionStats &MarginalCalls() { return marginal; }

	/*! Constructor - clears all members. */
	ProbeSetStats() { Clear(); }

	/*! Clears all data. */
	void Clear() { numSets=0; present.Clear(); absent.Clear(); marginal.Clear(); }
};

}

#endif
