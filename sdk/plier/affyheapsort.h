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

/*
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The PLIER (Probe Logarithmic Error Intensity Estimate) method produces
an improved signal by accounting for experimentally observed patterns 
in probe behavior and handling error at the appropriately at low and 
high signal values.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/

#ifndef _AFFX_HEAP_SORT_
#define _AFFX_HEAP_SORT_

//////////////////////////////////////////////////////////////////////

inline long LeftChild(long a)
{
	return(2*a+1);
}

//////////////////////////////////////////////////////////////////////

inline long RightChild(long a)
{
	return(2*a+2);
}

//////////////////////////////////////////////////////////////////////

inline void Heap(double *Vector, long *TieBreaker, long *Rank, long ListSize, long start)
{
	long notdone;
	long left, right, largest, current;
	long tmprank;
	// object at start has been put on the appropriate subheap
	// propagate it down the heap as required
	current = start;
	notdone = 1;
	largest = current;
	while (notdone)
	{
		left = LeftChild(current);
		if (left<ListSize)
		{
			if (Vector[Rank[largest]]<Vector[Rank[left]])
				largest = left;
		}

		right = RightChild(current);

		if (right<ListSize)
		{
			if (Vector[Rank[largest]]<Vector[Rank[right]])
				largest = right;
		}

		if (largest==current)
			notdone = 0; // stop now with heap property
		else
		{
			tmprank = Rank[current];
			Rank[current] = Rank[largest];
			Rank[largest] = tmprank;
			current = largest; // just swapped into this position
		}
	}
}

//////////////////////////////////////////////////////////////////////

inline void HeapIndex(double *Vector, long *TieBreaker, long *Rank, long ListSize)
{
	long i;
	long tmprank;
	// provides a rank vector to go with an array by using heapsort
	// all vectors must be allocated beforehand
	for (i=0; i<ListSize; i++)
		Rank[i] = i; // put pointer to Vector
	// now sort
	for (i=ListSize/2; i>-1; i--)
		Heap(Vector, TieBreaker, Rank, ListSize, i); // do a heap-building operation
	for (i=ListSize-1; i>-1; i--)
	{
		tmprank = Rank[i]; // sort index
		Rank[i] = Rank[0]; // pop off heap to current end of list
		Rank[0] = tmprank; // put element at top of heap
		Heap(Vector, TieBreaker, Rank, i, 0); // build heap again
	}
}

//////////////////////////////////////////////////////////////////////

inline int CompareV(double *VA, double *VB, long VLength)
{
	long i;
	// Vlength entries
	for (i=0; i<VLength; i++)
	{
		if (VA[i]<VB[i])
			return(1);
		if (VA[i]>VB[i])
			return(0);
		// if equal, try the next one in line
	}
	return(0); // not less than, all mysteriously equal!
}

//////////////////////////////////////////////////////////////////////

inline void HeapMatrix(double **Vector, long *Rank, long ListSize, long VLength, long start)
{
static	long notdone;
static	long left, right, largest, current;
static	long tmprank;
	// object at start has been put on the appropriate subheap
	// propagate it down the heap as required
	current = start;
	notdone = 1;
	largest = current;
	while (notdone)
	{
		left = LeftChild(current);

		if (left<ListSize)
		{
			if (CompareV(Vector[Rank[largest]],Vector[Rank[left]], VLength))
				largest = left;
		}

		right = RightChild(current);

		if (right<ListSize)
		{
			if (CompareV(Vector[Rank[largest]],Vector[Rank[right]], VLength))
				largest = right;
		}

		if (largest==current)
			notdone = 0; // stop now with heap property
		else
		{
			tmprank = Rank[current];
			Rank[current] = Rank[largest];
			Rank[largest] = tmprank;
			current = largest; // just swapped into this position
		}
	}
}

//////////////////////////////////////////////////////////////////////

inline int HeapIndexMatrix(double **Matrix, long *Rank, long ListSize, long VLength)
{
static	long i;
static	long tmprank;
	// provides a rank vector to go with an array by using heapsort
	// all vectors must be allocated beforehand
	for (i=0; i<ListSize; i++)
		Rank[i] = i; // put pointer to entries in Matrix
	// now sort
	for (i=ListSize/2; i>-1; i--)
		HeapMatrix(Matrix, Rank, ListSize, VLength, i); // do a heap-building operation
	for (i=ListSize-1; i>-1; i--)
	{
		tmprank = Rank[i]; // sort index
		Rank[i] = Rank[0]; // pop off heap to current end of list
		Rank[0] = tmprank; // put element at top of heap
		HeapMatrix(Matrix, Rank, i, VLength, 0); // build heap again
	}
	return(1);
}

//////////////////////////////////////////////////////////////////////

#endif
