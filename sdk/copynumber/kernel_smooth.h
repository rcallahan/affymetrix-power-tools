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
#ifndef KERNEL_SMOOTH_H
#define KERNEL_SMOOTH_H

/** A template function for kernel smoothing.  At the time this was written,
 *  the primary data types of interest were vector<T> and valarray<T>.
 *  Since valarray does not provide iterators, traversing the data from
 *  begin() to end() was not possible.  However both supply a size() method
 *  and operator[](unsigned int).
 *
 *  Neither the kernel nor kernel weights need to be symmetrical.
 *
 *  typename D data, are the array of data expected to be of type double
 *  or float.  These data are overwritten by their smoothed values.  The
 *  container of the data array must provide a method size() to determine
 *  the number of elements of the array and operator[](unsigned int) to
 *  dereference indices into the container.
 *
 *  typename T, is the kernel for smoothing.  The kernel is instantiated
 *  and configured before the algorithm is called.  It operates by returning
 *  smoothed value to replace the array value between lead and lag.  The
 *  The intended application in mind has the smoothed value as a weighted
 *  value of array elements in the span of the kernel.  However, other
 *  transformations are possible but better implemented using <algorithm>
 *  from the standard template library.  This implementation takes care
 *  of end points in a way that STL was not designed to do.
 *
 *  The kernel requires 3 methods in its interface.
 *  unsigned int lead() returns how many data points are considered in its lead.
 *  unsigned int lag() returns how many data points are in its lag.
 *  double step(double) accepts a new value.
 *
 *  To understand how the kernel and the algorithm operate together it is
 *  best to consider it in the middle of the data, then the ends can be
 *  treated as special cases.
 *
 *  In figure 1, a kernel with a symmetric span of 7 data points is show
 *  sitting over an array of data.  The center of the kernel is sitting
 *  over the ^.  The algorithm has just replaced the array value at ^ with
 *  the weighted average given by the kernel.  The next value to feed the
 *  kernel is the one over ~.  When the kernel receives the value over ~
 *  using its method step(double) the algorithm will be able to replace
 *  the value to the right of ^ with the kernel smoothed value.  This is
 *  how the kernel is expected to be able to step through the data.  The
 *  algorithm expects the implementation of the kernel to hold a private
 *  copy of data that it is sitting over in order to avoid averaging over
 *  already smoothed values.
 *
 *                 _______
 *  fig 1:    xxxxxxxxxxxxxxxxxxxx
 *                    ^   ~
 *
 *  Now the end points need to be considered as special cases.  This
 *  implementation of the algorithm assumes that augmenting each end of
 *  the data with average values near the ends will suffice.  This
 *  assumption may not satisfy all situations in which case a new method
 *  should be written.  Almost all the code here deals with the ends.
 *
 *  One can, with no loss in accuracy of calculations, imagine that the
 *  starts well to the left of the data and keeps moving to the right using
 *  step(double).  The picture in figure 2 shows this.  The imaginary data
 *  values way off to the left are indicated by ........ and once the
 *  kernel gets close enough to the beginning of the real data, imputed
 *  values need to be provided to smooth the transition between imaginary
 *  data that matters.  The ooo are data values that will to be imputed
 *  and they matter because the kernel will use them to compute the first
 *  smoothed value.  The left-most imputed o in figure 2 will affect the
 *  first smoothed data value.
 *
 *                _______
 *  fig 2:    ...........oooxxxxxxxxx
 *                   ^   ~
 *
 *
 *  In figure 3, the kernel has marched so that its center is just to the
 *  left of where the first observation is.  When the kernel is fed the
 *  observation over the ~ using the method step(double), it will return
 *  the smoothed value computed to replace the first observation.
 *                      _______
 *  fig 3:        .......oooxxxxxxxxx
 *                         ^   ~
 *
 *
 *  In figure 4, the kernel object has been moved by feeding it the
 *  observation over the ~ in figure 3.  It has replaced the value
 *  over the ^ and is ready to repeat the process one step at a time.
 *  It is the burden of the implementer of the kernel object to manage
 *  the data being fed to it and return a smooth value for over the ^.
 *  The data can be thought of as being fed to the kernel object from
 *  the leading end as it steps.  The middle is easy to code.
 *                   _______
 *  fig 4:    .......oooxxxxxxxxx
 *                      ^   ~
 *
 *
 *  In figure 5, the leading edge of the kernel object has run up to the
 *  end of the array of data.  In order for the kernel to step all the way
 *  to the end so that the ^ is under the last element of the array,
 *  imputed data must be fed to the kernel.  Once again the imputed data
 *  are indicated by ooo and imaginary data that does not matter by .....
 *                        _______
 *  fig 5:         xxxxxxxxxxxxxxooo...............
 *                           ^   ~
 *
 *
 *  In figure 6, the leading edge of the kernel object has reached the
 *  end of the imputed data and the center is over the last element
 *  of the data.  The algorithm is finished at this point.
 *                           _______
 *  fig 6:         xxxxxxxxxxxxxxooo
 *                              ^   ~
 *
 *
 *  Finally in figure 7, the Macarena.
 *
 *  fig 7:
 *
 *           o      o     o    o     o    <o     <o>    o>    o
 *          .|.    \|.   \|/   //    X     \      |    <|    <|>
 *           /\     >\   /<    >\   /<     >\    /<     >\   /<
 *
 *
 */


template<typename D,typename T>
void kernel_smooth(D & data, T & kernel)
    {
    // If there is no data return
    if (!data.size()) return;

    // This does not assume that the kernel span or weights are symmetrical.
    // See the long winded ascii art description above for details of what
    // the kernel and the algorithm are responsible for.

    // The lag of the kernel will hang off the left hand side of the array.
    // Assign to the lag an average of the span of the kernel at the left
    // hand side.  The user may have submitted an array shorter than the lead.
    unsigned int right_span = kernel.lead() + 1;
    if (data.size() < right_span) right_span = data.size();

    double left_avg = 0.0;
    for (unsigned int i=0; i<right_span; i++) left_avg += data[i];
    left_avg /= right_span;


    // Do the same as above but at the right hand side of the array.
    unsigned int left_span = kernel.lag() + 1;
    if (data.size() < left_span) left_span = data.size();

    double right_avg = 0.0;
    for (unsigned int i=0; i<left_span; i++)
        right_avg += data[data.size() - 1 - i];
    right_avg /= left_span;

    // Here goes the left hand side.  Stuff the left average into the lag.
    for (unsigned int i=0; i<kernel.lag(); i++) kernel.step(left_avg);

    // Step the kernel into the data, stopping the center one short.
    // Make sure that the kernel lead does not extend past the data.  If it
    // does extend past then pad it with the imputed average on the right.
    for (unsigned int i=0; i<kernel.lead(); i++)
        {
        if (i < data.size()) kernel.step(data[i]);
        else kernel.step(right_avg);
        }

    // Here goes the middle. Pad on the right when the lead goes pas the end.
    int lead =  kernel.lead();
    int data_size = data.size();
    for (unsigned int i=0; i<data_size; i++)
        {
        unsigned int lead_pos = i + lead;
        if (lead_pos < data_size) data[i] = kernel.step(data[lead_pos]);
        else data[i] = kernel.step(right_avg);
        }
    }

#endif
