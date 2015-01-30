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
#include "copynumber/CNSmoother.h"
//

// Any copy number segments shorter than m_min_length are replace with
// values conforming to their neighbors

using namespace std;

// need to do initialize a class variable in a cpp file
int CNSmootherSegment::m_min_length = CNSMOOTHER_SEGMENT_LENGTH;



// if a segment at least minimum length and so are its neighbors its
// membership of markers will be fixed
bool CNSmootherSegment::fixed() const
    {
    if (m_segment->size < m_min_length) return false;
    if (m_segment->prev) if (m_segment->prev->size < m_min_length) return false;
    if (m_segment->next) if (m_segment->next->size < m_min_length) return false;
    return true;
    }


// expansion of a segment consists of taking a marker from any neighbor
// with fewer than the minimum number of markers
void CNSmootherSegment::expand()
    {
    // these will do nothing in the cases that neighbor has enough markers
    take_from_prev();
    take_from_next();
    }


//  Take a marker from the previous member in the doubly linked list.
void CNSmootherSegment::take_from_prev()
    {
    // If this is the beginning of the list do nothing
    if (!m_segment->prev) return;

    // If the previous member is of minimum size do nothing
    if (m_segment->prev->size >= m_min_length) return;

    // Push the start back one and add to the marker count
    m_segment->position -= 1;
    m_segment->size += 1;

    // The previous segment has one marker fewer
    m_segment->prev->size -= 1;

    // Check to see the previous segment has at least one member
    if (m_segment->prev->size) return;

    // Otherwise the previous segment empty so skip it for a neighbor
    m_segment->prev = m_segment->prev->prev;

    // If the new neighbor is not the terminus make the skip reciprocal
    if (m_segment->prev) m_segment->prev->next = m_segment;
    }


// Exact same as take_from_prev but instead looking at the segment following
void CNSmootherSegment::take_from_next()
    {
    // If this is the end of the list do nothing
    if (!m_segment->next) return;

    // If the following member is of minimum size do nothing
    if (m_segment->next->size >= m_min_length) return;

    // Add to the marker count
    m_segment->size += 1;

    // Push the start of the following segment ahead one marker
    m_segment->next->position += 1;

    // The following segment has one marker fewer
    m_segment->next->size -= 1;

    // Check to see the previous segment has at least one member
    if (m_segment->next->size) return;

    // Otherwise the following segment empty so skip it for a neighbor
    m_segment->next = m_segment->next->next;

    // If the new neighbor is not the terminus make the skip reciprocal
    if (m_segment->next) m_segment->next->prev = m_segment;
    }


// This does all the work.  Values is a list of copy numbers that gets
// summarized by segments.  Short segments are then consumed by their
// larger neighbors one marker at a time from their ends.  This is
// accomplished by markers at the end of short segments to neighboring
// larger segments.
void CNSmoother::smooth(valarray<int> & values)
    {
    // If the list of marker values is empty do nothing.
    if (!values.size()) return;

    // values is summarised as segments here
    vector<segment_info_t> stored_segments;

    // Initialize the first segment
    int start=0;               // the start position of the segment
    int seg_value=values[0];   // the segment copy number value

    // Assemble the segments one marker at a time
    for (unsigned int end=1; end<values.size(); end++)
        {
        // roll along until copy number changes
        if (values[end] == seg_value) continue;

        // something changed so push back the segment
        segment_info_t segment(start, end - start, seg_value);
        stored_segments.push_back(segment);

        // one past the end is becomes the start of an empty segment
        start = end;
        seg_value = values[start];
        }

    // Close up the last segment and push it back on the list
    segment_info_t segment(start,values.size() - start,seg_value);
    stored_segments.push_back(segment);


    // Initial segments are all in a row.  Need them in a double linked list
    vector<CNSmootherSegment> cns_vec(stored_segments.size());

    // Initialize the values in the doubly linked list
    for (unsigned int i=0; i<stored_segments.size(); i++)
        {
        cns_vec[i] = CNSmootherSegment(&stored_segments[i]);
        }

    // Hook up the reciprocal links in the list
    for (unsigned int i=1; i<stored_segments.size(); i++)
        {
        cns_vec[i].set_prev(&stored_segments[i-1]);
        cns_vec[i-1].set_next(&stored_segments[i]);
        }

    // Order them largest to smallest.  This does not affect linkage
    sort(cns_vec.begin(),cns_vec.end(),segment_sort_size());

    // A queue of segments that are not fixed.
    queue<CNSmootherSegment>cns_queue;

    // Load the queue with segments that are larger than the minimum size
    // but also have at least one neighbor that is shorter than minimum size
    for (vector<CNSmootherSegment>::reverse_iterator iter=cns_vec.rbegin();
            iter != cns_vec.rend(); iter++)
        {
        // Because of the sorting the first segment found too short signals
        // that so too are all remaining segments
        if (iter->too_short()) break;

        // Anything segment with no neighbors to acquire markers from does
        // not belong in the queue
        if (iter->fixed()) continue;

        // Otherwise push it back
        cns_queue.push(*iter);
        }

    // Done when the queue is empty
    while (cns_queue.size())
        {
        // If its neighbors are all of minimum size toss it from the queue
        if (cns_queue.front().fixed())
            {
            cns_queue.pop();
            continue;
            }

        // Otherwise acquire a marker from any neighbor with fewer markers
        // than the minimum.
        cns_queue.front().expand();

        // Put it at the back of the queue.
        cns_queue.push(cns_queue.front());
        cns_queue.pop();
        }


    // The doubly linked list of segments has been maintained.  The starts
    // and lengths of segments has changed so that they all have the minimum
    // number of markers.  Any segments starting with fewer than the
    // minimum numbers of markers are now empty with a marker count of zero.
    // Overwrite the input copy numbers.
    for (unsigned int i=0; i<stored_segments.size(); i++)
        {
        if (!stored_segments[i].size) continue;
        int pos = stored_segments[i].position;
        unsigned int size = stored_segments[i].size;
        unsigned int value = stored_segments[i].value;
        for (unsigned int j=0; j<size; j++) values[pos + j] = value;
        }
    }
