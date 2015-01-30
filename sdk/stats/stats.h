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

/// @file   stats.h
/// @brief  Stats functions written for the transcriptome project.

#ifndef __STATS_H_
#define __STATS_H_

#include "stats/stats-util.h"
//
#include "portability/affy-base-types.h"
#include "util/Util.h"
//
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//
#ifdef WIN32
#define trunc(x) (floor(x))
#ifndef HAVE_RANDOM_SAMPLE // Even in windows some STLs (like stlport) have random sample
#include "affy_random_sample.h"
#endif
#else // g++ keeps random_sample in ext/algorithm.
#include <ext/algorithm>
#endif
#ifdef sparc
#define trunc(x) (floor(x))
#endif

using namespace std;

#undef _RANDOM_SAMPLE_ERROR_

/// @param     start     beginning iterator
/// @param     stop      end iterator
/// @param     percentage which precentile to return (0-100)
/// @param     sampleSize (optional) size of sample to take from the input data.
/// @brief     Return the value at a given percentile
/// @remark This does not modify the input data.  If the
/// sampleSize is larger than the inputsize, all the input
/// is used.
template <class RAIterator>
typename iterator_traits<RAIterator>::value_type  
percentile(RAIterator start, RAIterator stop, double percentage, size_t sampleSize = ULONG_MAX)
{
  vector<typename iterator_traits<RAIterator>::value_type> tmp;

#ifndef _RANDOM_SAMPLE_COMPILE_ERROR_
    // tmp.resize(min(sampleSize, stop - start));   // original
    // JHG SIZE_T DEBUG
    assert(stop>start);
    size_t tmp_size;
    tmp_size=stop - start;
    tmp_size=min(sampleSize,tmp_size);
    // Check for negative values or a truly large allocation 
    assert(0<tmp_size);

    // since we shouldnt modify the data, we make a copy.
    tmp.resize(tmp_size);
    if (sampleSize > tmp_size) { //stop - start) {
      copy(start, stop, tmp.begin());
    }
    else {
      random_sample(start, stop, tmp.begin(), tmp.end());
    }
#else
    tmp.resize(stop - start);
    copy(start, stop, tmp.begin());
#endif

    return percentile_in_place(tmp.begin(), tmp.end(), percentage);
}

/// @brief     Find the value for a given precentile, modifying the data is ok.
/// @param     start     beginning iterator
/// @param     stop      end iterator
/// @param     percentage the percentile to find (0-100)
/// @return    the value found
/// @remark Work with the data "in_place" -- we may shuffle
/// the values around as we work. No copy needed.
template <class RAIterator>
typename iterator_traits<RAIterator>::value_type  
percentile_in_place(RAIterator start, RAIterator stop, double percentage) {
  // be sure percentage is in range.
  assert((0.0<=percentage)&&(percentage<=100.0));
  //
  double bucket_index = (stop - start - 1)* percentage / 100.0;
  return (bucket_index == trunc(bucket_index)?
	  RandomizedSelect(start, stop, (int) bucket_index):
	  (RandomizedSelect(start, stop, (int) floor(bucket_index)) +
	   RandomizedSelect(start, stop, (int) ceil(bucket_index)))/2);
};

template <class RAIterator>
class Median: public binary_function<RAIterator, RAIterator, typename iterator_traits<RAIterator>::value_type> {
 public:
  long sampleSize;
  Median(long ss):sampleSize(ss) {};
    typename iterator_traits<RAIterator>::value_type
      operator()(RAIterator start, RAIterator stop) {
      return median(start, stop, sampleSize);
    }
};

template <class RAIterator>
typename iterator_traits<RAIterator>::value_type  
median(RAIterator start, RAIterator stop, size_t sampleSize=ULONG_MAX) 
{
  //typename iterator_traits<RAIterator>::difference_type sampleSize = LONG_MAX) {
  //printf("MEDIAN: ss=%lu\n",sampleSize); // JHG 
  assert(sampleSize>0); // JHG SIZE_T DEBUG
  return percentile(start, stop, 50, sampleSize);
}

template <class RAIterator>
typename iterator_traits<RAIterator>::value_type  
median_in_place(RAIterator start, RAIterator stop) {
  return percentile_in_place(start, stop, 50);
}

template <class RAIterator, class LargeNumber = double>
class Average: public binary_function<RAIterator, RAIterator, double> {
public:
  double operator()(RAIterator start, RAIterator stop) {
    int size = stop - start;
    adder<LargeNumber> sum = 
      for_each(start, stop, adder<LargeNumber>());
    return sum.result/size;
  }
};

template <class RAIterator>
double average(RAIterator start, RAIterator stop) {
  return Average<RAIterator, typename iterator_traits<RAIterator>::value_type>()(start,stop);
};

template <class RAIterator>
double pseudo_median(RAIterator start, RAIterator stop) {
  vector<double> averages;
  averages.reserve((start-stop)*(start-stop));
  for (RAIterator i = start; i!=stop; i++) {
    for (RAIterator j = start; j!=stop; j++) {
      averages.push_back((*i + *j)/2);
    }
  }
  return median(averages.begin(), averages.end());
};

template <class RAIterator, class LargeNumber = double>
class Variance: public binary_function<RAIterator, RAIterator, double> {
public:
  double operator()(RAIterator start, RAIterator stop) {
    int size = stop - start;
    square_adder<LargeNumber> square_sum = 
      for_each(start, stop, square_adder<LargeNumber>());
    double avg = Average<RAIterator, LargeNumber>()(start, stop);
    double result = square_sum.result/size - avg * avg;
    return (result <= 0 ? 0: result);
  }
};

template<class RAIterator>
double variance (RAIterator start, RAIterator stop) {
  return Variance<RAIterator, typename iterator_traits<RAIterator>::value_type>()(start, stop);
}

template <class RAIterator, class LargeNumber = double>
class UnbiasedVariance: public binary_function<RAIterator, RAIterator, double> {
public:
  double operator()(RAIterator start, RAIterator stop) {
    int size = stop - start;
    return Variance<RAIterator,LargeNumber>()(start, stop) * size / (size - 1);
  }
};

template <class RAIterator>
double unbiased_variance(RAIterator start, RAIterator stop) {
  //int size = stop - start;
  return UnbiasedVariance<RAIterator, typename iterator_traits<RAIterator>::value_type>()(start, stop);
}


template <class RAIterator>
double varcoeff(RAIterator start, RAIterator stop) {
  return (stddev(start, stop) / average(start, stop));
}

template <class RAIterator>
double unbiased_stddev(RAIterator start, RAIterator stop) {
  return sqrt(unbiased_variance (start, stop));
}

template <class RAIterator>
double stddev(RAIterator start, RAIterator stop) {
  return sqrt(variance (start, stop));
}

template <class RAIterator>
typename iterator_traits<RAIterator>::value_type  
interquartile_range(RAIterator start, RAIterator stop, long sampleSize = LONG_MAX) {
  return percentile(start, stop, 75, sampleSize) - percentile(start, stop, 25, sampleSize);
}

template <class RAIterator, class LargeNumber = double>
class Covariance: public binary_function<RAIterator, RAIterator, double> {
public:
  double operator()(RAIterator start1, RAIterator stop1, RAIterator start2) {
    int size = stop1 - start1;
    vector <LargeNumber> mult_tmp(size);
    transform(start1, stop1, start2, mult_tmp.begin(), multiplies<LargeNumber>());
    return Average<typename vector<LargeNumber>::iterator, LargeNumber>()(mult_tmp.begin(), mult_tmp.end()) -  
           Average<RAIterator,LargeNumber>()(start1, stop1) * 
           Average<RAIterator,LargeNumber>()(start2, start2 + size);
  }
};

template <class RAIterator, class LargeNumber = double>
class UnbiasedCovariance: public binary_function<RAIterator, RAIterator, double> {
public:
  double operator()(RAIterator start1, RAIterator stop1, RAIterator start2) {
    int size = stop1 - start1;
    return Covariance<RAIterator,LargeNumber>()(start1, stop1, start2) * size / (size -1);
  }
};

template <class RAIterator>
double unbiased_covariance(RAIterator start1, RAIterator stop1, RAIterator start2) {
  int size = stop1 - start1;  
  return (covariance(start1, stop1, start2) * size) /(size - 1);
}

template <class RAIterator>
double covariance(RAIterator start1, RAIterator stop1, RAIterator start2) {
  return Covariance<RAIterator, typename iterator_traits<RAIterator>::value_type>()(start1, stop1, start2);
};

template <class RAIterator>
double correlation_coeff(RAIterator start1, RAIterator stop1, RAIterator start2) {
   return (covariance(start1, stop1, start2)/sqrt(variance(start1, stop1) * variance(start2, start2 + (stop1 - start1))));
};

template <class RAIterator>
double variation_coeff(RAIterator start, RAIterator stop) {
  return stddev(start, stop)/
    (DBL_MIN + fabs(average(start,stop)));
}

template <class RAIterator>
double robust_variation_coeff(RAIterator start, RAIterator stop) {
  return interquartile_range(start, stop)/
    (DBL_MIN + fabs(median(start,stop)));
}

template <class RAIterator>
RAIterator reproducible_random_sample(RAIterator start, RAIterator  stop,
                                      RAIterator out, uint32_t size, int32_t seed=0)
{
  // Random number generator can't be seeded with zero, our clue to call time.
  if(seed == 0)
    seed = (uint32_t) time(NULL);
  
  uint32_t end = size;
  uint32_t index = 0;

  /* Copy initial values. */
  for( ; start != stop && index < size; ++index, ++start) 
    out[index] = *start;

  /* Give every value a chance to overwrite the previous versions. */
  while(start != stop) {
    ++end;
    uint32_t offset = Util::schrageRandom(&seed) % end;
    if (offset < size)
      out[offset] = *start;
    ++start;
  }
  
  return out + index;
}

#endif // __STATS_H_
