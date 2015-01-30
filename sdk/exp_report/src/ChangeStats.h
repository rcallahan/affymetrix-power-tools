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

#ifndef _ChangeStats_HEADER_
#define _ChangeStats_HEADER_

/*! \file ChangeStats.h Defines a class to store statistics about a change type results. */

namespace ExpressionReport
{

/*! Stores statistics about a change type results. */
class ChangeStats
{
public:
	/*! The number of bins in the fold change counter. */
	#define NUMBER_FOLD_CHANGE_BINS 5

private:
	/*! The number of calls of a type. */
	int changeCount;

	/*! The number of detection calls that changed from baseline to experiment. */
	int detectionChangeCount;

	/*! The number of detection calls that stayed the same at present. */
	int detectionPresentCount;

	/*! The number of detection calls that stayed the same at absent. */
	int detectionAbsentCount;

	/*! The number of moderate calls. */
	int moderateCount;

	/*! Bin counts of fold change range. */
	int foldChangeCount[NUMBER_FOLD_CHANGE_BINS];

public:
	/*! Initialize the members. */
	void Clear()
	{
		changeCount=0;
		detectionChangeCount=0;
		detectionPresentCount=0;
		detectionAbsentCount=0;
		moderateCount=0;
		for (int i=0; i<NUMBER_FOLD_CHANGE_BINS; i++)
			foldChangeCount[i] = 0;
	}

	/*! The number of calls of a type. */
	int ChangeCount() const { return changeCount; }

	/*! The number of detection calls that changed from baseline to experiment. */
	int DetectionChangeCount() const { return detectionChangeCount; }

	/*! The number of detection calls that stayed the same at present. */
	int DetectionPresentCount() const { return detectionPresentCount; }

	/*! The number of detection calls that stayed the same at absent. */
	int DetectionAbsentCount() const { return detectionAbsentCount; }

	/*! The number of moderate calls. */
	int ModerateCount() const { return moderateCount; }

	/*! Bin count of fold change range.
	 * @param bin The bin number.
	 */
	int FoldChangeCount(int bin) const { return foldChangeCount[bin]; }

	/*! Increment the number of calls of a type. */
	void IncrementChangeCount() { ++changeCount; }

	/*! Increment the number of detection calls that changed from baseline to experiment. */
	void IncrementDetectionChangeCount() { ++detectionChangeCount; }

	/*! Increment the number of detection calls that stayed Increment the same at present. */
	void IncrementDetectionPresentCount() { ++detectionPresentCount; }

	/*! Increment the number of detection calls that stayed Increment the same at absent. */
	void IncrementDetectionAbsentCount() { ++detectionAbsentCount; }

	/*! Increment the number of moderate calls. */
	void IncrementModerateCount() { ++moderateCount; }

	/*! Increment the bin count of fold change range.
	 * @param bin The bin number.
	 */
	void IncrementFoldChangeCount(int bin) { ++foldChangeCount[bin]; }

	/*! Constructor */
	ChangeStats() { Clear(); }

	/*! Destructor */
	~ChangeStats() { }

	/*! Get the bin ranges for the fold change counter. */
	static void GetBinRange(int bin, float &lower, float &upper)
	{
		float ranges[] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 3.402823466e+38F};
		lower = ranges[bin];
		upper = ranges[bin+1];
	}
};

}

#endif
