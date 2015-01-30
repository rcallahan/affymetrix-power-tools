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

#ifndef _DetectionStats_HEADER_
#define _DetectionStats_HEADER_

/*! \file DetectionStats.h Defines a class to store information about a detection value. */

namespace ExpressionReport
{

/*! Stores information about a detection value. */
class DetectionStats
{
private:
	/*! The number of detection values. */
	int count;

	/*! The total signal for a detection value. */
	double signal;

public:
	/*! Add one to the detection count. */
	void IncrementCount() { ++count; }

	/*! Add to the signal value.
	 * @param s The signal to add.
	 */
	void AddSignal(double s) { signal += s; }

	/*! The total number of detection calls. */
	int &Count() { return count; }	

	/*! The total signal. */
	double &Signal() { return signal; }

	/*! Constructor */
	DetectionStats() { Clear(); }

	/*! Clears members. */
	void Clear() { count=0; signal=0.0f; }
};

}

#endif
