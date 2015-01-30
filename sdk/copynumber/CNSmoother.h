////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
/**
 * @file CNSmoother.h
 *
 * @brief This header contains the CNSmoother.h class definition.
 */

#ifndef COPY_NUMBER_SMOOTHER_H
#define COPY_NUMBER_SMOOTHER_H

// The shortest that a segment can be is 1 marker in a row
#ifndef CNSMOOTHER_SEGMENT_LENGTH
#define CNSMOOTHER_SEGMENT_LENGTH 1
#endif


#include <algorithm>
#include <queue>
#include <valarray>
//

// Basic information about a copy number segment
struct segmentInfo {
    segmentInfo() : position(0),size(0),value(0),prev(0),next(0) { }
    segmentInfo(int pos, unsigned int sz, unsigned int val) :
        position(pos), size(sz), value(val), prev(0), next(0) { }
    int position;            // number of markers from the first
    unsigned int size;       // number of markers in the segment
    unsigned int value;      // usually translates to copy number
    struct segmentInfo * prev;  // doubly linked list
    struct segmentInfo * next;
    };

// give it a name
typedef struct segmentInfo segment_info_t;


/*
 * @brief class that wraps the segment info struct
 */

class CNSmootherSegment {
public:

CNSmootherSegment() : m_segment(0) { }
CNSmootherSegment(segment_info_t * seg_ptr) : m_segment(seg_ptr) { }

// Minimum lenght of a segment
static int get_min_length() { return m_min_length; }
static void set_min_length(int length) { m_min_length = length; }

// Number of markers from the start, number of markers and their common value
int position() const { return m_segment->position; }
unsigned int size() const { return m_segment->size; }
unsigned int value() const { return m_segment->value; }

// Set or get the previous marker in the doubly linked list
void set_prev(segment_info_t * prev) { m_segment->prev = prev; }
segment_info_t * get_prev() const { return m_segment->prev; }

// Set or get the next marker in the doubly linked list
void set_next(segment_info_t * next) { m_segment->next = next; }
segment_info_t * get_next() const { return m_segment->next; }

// True if the segment has too few markers to be kept
bool too_short() const { return (m_segment->size < m_min_length); }

// True if the segment and its neighbours will be kept
bool fixed() const;

// Grow the segment one marker at either end that can steal from a
// neighbour that is too_short().
void expand();


private:
// The basic information about a segment (start,size,value)
segment_info_t * m_segment;

// Used to expand a segment by stealing markers fro its neighbour
void take_from_prev();
void take_from_next();

// minimum segment length is a global variable
static int m_min_length;
};


// For sorting segments by number of markers
class segment_sort_size {
public:
bool operator()(const CNSmootherSegment & S1,
            const CNSmootherSegment & S2) const
    {
    return (S1.size() < S2.size());
    }
};


class CNSmoother {
public:
CNSmoother(const int mlength) { CNSmootherSegment::set_min_length(mlength); }

// The one function that actually does something
void smooth(std::valarray<int> & input_values);
};

#endif
