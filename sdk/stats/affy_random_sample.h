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

////////////////////////////////////////////////////////////////
//
// Original Copyright and Permission Notice from SGI and HP
// for the random_sample code that this file is based on:
//
// Copyright © 1996-1999
// Silicon Graphics Computer Systems, Inc.
//
// Permission to use, copy, modify, distribute and sell this 
// software and its documentation for any purpose is hereby 
// granted without fee, provided that the above copyright 
// notice appears in all copies and that both that copyright 
// notice and this permission notice appear in supporting 
// documentation. Silicon Graphics makes no representations 
// about the suitability of this software for any purpose. 
// It is provided "as is" without express or implied 
// warranty.
//
// Copyright © 1994
// Hewlett-Packard Company
//
// Permission to use, copy, modify, distribute and sell this 
// software and its documentation for any purpose is hereby 
// granted without fee, provided that the above copyright 
// notice appears in all copies and that both that copyright 
// notice and this permission notice appear in supporting 
// documentation. Hewlett-Packard Company makes no 
// representations about the suitability of this software 
// for any purpose. It is provided "as is" without express 
// or implied warranty. 
//
////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
//
// Commentary:
//
// a replacement random_sample function 
// suitable for use with systems which dont have it.
//
// This code was taken from the SGI STL.
// and modified for it to work.
// The modification was to replace "__random_number(__t) with
// "rand()%__t" as the former is an internal SGI-only STL call.
//
// You will want something like the following in your file:
// #ifdef _WIN32
//   #include "affy_random_sample.h"
// #else
//   #include <ext/algorithm>
// #endif
//
/////////////////////////////////////////////////////////////////

#ifndef __AFFY_RANDOM_SAMPLE_H_
#define __AFFY_RANDOM_SAMPLE_H_

#include <cstdio>
#include <cstdlib>
//

template <class _InputIter, class _RandomAccessIter, class _Distance>
_RandomAccessIter __random_sample(_InputIter __first, _InputIter __last,
                                  _RandomAccessIter out,
                                  const _Distance __n)
{
  _Distance __m = 0;
  _Distance __t = __n;
  for ( ; __first != __last && __m < __n; ++__m, ++__first) 
    out[__m] = *__first;

  while (__first != __last) {
    ++__t;
    //_Distance __M = __random_number(__t);
    _Distance __M = (rand()%__t);
    if (__M < __n)
      out[__M] = *__first;
    ++__first;
  }

  return out + __m;
}

template <class _InputIter, class _RandomAccessIter>
inline _RandomAccessIter
affy_random_sample(_InputIter __first, _InputIter __last,
              _RandomAccessIter __out_first, _RandomAccessIter __out_last) 
{
  //__STL_REQUIRES(_InputIter, _InputIterator);
  //__STL_REQUIRES(_RandomAccessIter, _Mutable_RandomAccessIterator);
  return __random_sample(__first, __last,
                         __out_first, __out_last - __out_first);
}

// Now we make it the replacement
#define random_sample affy_random_sample

#endif // __AFFY_RANDOM_SAMPLE_H_
