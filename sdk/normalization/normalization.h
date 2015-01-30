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

#ifndef __NORMALIZATION_H_
#define __NORMALIZATION_H_

#include "normalization/vectormap.h"
//
#include "portability/affy-base-types.h"
#include "stats/stats-util.h"
#include "stats/stats.h"
//
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//

// random_sample is and extension to STL and not in MS.
// use the one we copied if ext/algorithm isnt available.
#ifdef WIN32
#include "stats/affy_random_sample.h"
#else
#include <ext/algorithm>
#endif

using namespace std;

/**
@page normalization Architecture: Normalize Template Code

In general microarray data must be normalized before comparisons and
other summaries can be made. The normalization module consists of a
number of template that implement a variety of normalization
algorithms including median normalization (MedianNormalization) which
scales the data linearly to have the same median and quantile
normalization (QuantileNormalization) which modifies the data
non-linearly to ensure that each microarray has the same
distribution. There are also methods to do an approximate (sketch)
quantile normalization using just a subset of the data which can
greatly reduce the memory impact usually associated with a full
quantile normalization.
 
@section quantileNorm Quantile Normalization

Classes for quantile normalizing the data that a set of iterators
point to. Basic algorithm is that data is sorted and the average value
across all iterators is used as the normalized data value for all of
the original data. This ensures that the distribution is identical for
the data pointed to by each iterator. Ties can occur when a particular
iterator points to data where there are multiple values that are
identical. In these cases the ranking of a particular data value in
the sort can be variable. These ties are resolved by taking the
average value of the normalized distribution for all of the ties. These
packages also support another popular way of resolving ties as
implemented by the bioconductor (http://www.bioconductor.org/) "affy"
package where the value of the middle position of the ties is used for
all the ties.

<b>Warning:</b> To resolve the ties efficiently with the average value
a vector of partial sums is used. Be sure the LargeNumber type has
sufficient range (i.e. doesn't overflow) and has sufficient precision
(i.e. enough bits to represent mantissa). Doubles seem to do well.
 
A sample usage would be:
@code
vector<vector<float> > data;
vector<vector <float>::iterator> dataIterVec;
// Note that the second type is double to ensure that no overflow and no loss of 
// precision in partial sums averaging.
QuantileNormalization<vector<vector<float>::iterator >::iterator, double> qn(false);  
for(int i = 0; i < data.size(); i++) {
  dataIterVec.push_back(data[i].begin());
}
qn(dataIterVec.begin(), dataIterVec.end(), *(dataIterVec.begin()) + data[0].size());
@endcode

This would modify the matrix data consisting of (chips are rows, probes are columns):
@verbatim
1 2 3  3  3  4  5  6  7
2 4 6  8  8  8  8 12 14
3 6 9 12 12 12 15 18 21
@endverbatim
to a normalized data matrix with the following values:
@verbatim
2 4 7.1 7.1 7.1 8   9.3 12 14
2 4 6   8.2 8.2 8.2 8.2 12 14
2 4 6   7.8 7.8 7.8 9.3 12 14
@endverbatim
note that using the bioconductor manner of resolving ties would result in
slightly different values:
@verbatim
2 4 7.7 7.7 7.7 8   9.3 12 14
2 4 6   7.7 7.7 7.7 7.7 12 14
2 4 6   7.7 7.7 7.7 9.3 12 14
@endverbatim
As seen using the following commands:
@code
library(affy);
dat <- read.table('CPPTest/input/doc-example.txt')
# Note that normalize.quantiles expects chips to be columns 
# and probes to be rows so transpose with t().
ndat <- normalize.quantiles(t(dat))
t(round(ndat, digits=1))
@endcode

@section sketchNorm Sketch Quantile Normalization

It is possible to achieve a very good approximation to the quantile
normalization using much less memory by using a sketch, or subset, of
the data and extrapolating the normalization of all data values from the normalization of
the sketch. As memory can be a limiting factor for many analysis with
the higher density Affymetrix arrays and the difference between full
quantile normalization and a sketch quantile normalization (or just
"sketch" for short) is very small the sketch normalization has become
the preferred manner for normalizing microarray data.

There are three phases to a sketch normalization:
 -# Extract the sketch (or subsample) of the data using the ExtractSketch class.
 -# Find the average distribution of all the sketches using AverageSketch.
 -# Create an interpolator (InterpolateSketch object) per microarray that will use linear interpolation to determine the normalized value of any intensity for that microarray.

A sample usage is:
@code
// Declarations
vector<vector<float> > data;
// Load data...
int sketchSize = data[0].size()/10;
// Note again that LargeNumber type is double
vector< InterpolateSketch< vector<float>::iterator, double> > m_Interpolators;
vector< float > m_AverageSketch;
vector< vector<float>::iterator > m_Sketches;
ExtractSketch<vector<float>::iterator > extractor(sketchSize);
bool doBioC = false; // Not being compatible with bioconductor.
// For each microarrays worth of data;
for(int i = 0; i < data.size(); i++) {
  // Allocate sketch. Don't forget to delete this later.
  std::vector<float> *pVec = new vector<float>(sketchSize);
  // Extract sketch.
  extractor(data[i].begin(), data[i].end(), pVec->begin());
  // Save sketch for later.
  m_Sketches.push_back(pVec->begin());
}
// Average sketches together.
AverageSketch< vector< vector<float>::iterator >::iterator, vector<float>::iterator, float> as;
m_AverageSketch.resize(sketchSize);
as(m_Sketches.begin(), m_Sketches.end(), (*m_Sketches.begin()) + sketchSize, m_AverageSketch.begin());
// Make interpolator for each chip.
for(int i = 0; i < m_Sketches.size(); i++) {
 m_Interpolators.push_back(InterpolateSketch<vector<float>::iterator, double>(m_Sketches[i], m_Sketches[i] + sketchSize,
                                                                              m_AverageSketch.begin(), 
                                                                              m_AverageSketch.end(), 
                                                                              doBioc));
  }

// Use interpolators to get normalized values.
for(int i = 0; i < data.size(); i++) {
  for(int j = 0; j < data[i].size(); j++) {
     float normValue = m_Interpolators[i](data[i][j]);
     // Do something with our normalized value.
  }
}
@endcode
*/

/** first parameter is a nested pointer **p
 * is a number *p is a sample; the second parameter should be a
 * numeric type that can contain the sum of all the numbers in all the
 * samples 
 * Also see: @ref quantileNorm
 */
//! NumberStarStar is a nested pointer. (**Number) 
//! LargeNumber is a type which can contain the sum of all the samples.
template <class NumberStarStar, class LargeNumber>
class QuantileNormalization {
  typedef typename iterator_traits<NumberStarStar>::value_type NumberStar;
  typedef typename iterator_traits<NumberStar>::value_type Number;
  //! Be compatible/comparable with the bioconductor.
  bool affyBioconductor;
public:

  /// @brief Create a new quantile normalization function object
  /// @param ab be compatible with bioconductor? (resolve ties via R rank() function)
  QuantileNormalization(bool ab = false): affyBioconductor(ab) {}

  /// @brief     quantile normalization of vectors of vectors
  /// @param     start    pointer to start of pointers to data
  /// @param     stop     pointer to end   of pointers to data
  /// @param     istop    pointer to end   of data in the first row
  void operator()(NumberStarStar start, NumberStarStar stop, NumberStar istop) {
    vector<LargeNumber> normVal; // for normalized values

    // Change this line to use the old "multimap" code...
    // (we dont use this because of the space requirements)
    // typedef multimap<Number,long> Perm_t;
    // ...use the new "vectormap" code instead.
    typedef vectormap<Number,int> Perm_t;
    vector <Perm_t> perm; // for permutations

    size_t nSamples = stop - start;
    NumberStar istart = *start;
    size_t sampleSize = istop - istart;
    assert(sampleSize>0); // JHG SIZE_T DEBUG

    normVal.resize(sampleSize);
    for(NumberStarStar currentStart = start; currentStart != stop;
        currentStart++) { //for each sample
      istart = *currentStart;
      istop = istart + sampleSize;
      vector <Number> sortedVal(istart, istop);
      sort(sortedVal.begin(),sortedVal.end()); //sort it
      perm.push_back(Perm_t());
      for (typename vector<Number>::iterator i = sortedVal.begin();
           i != sortedVal.end(); i++) {
        perm[currentStart - start].insert(pair<Number, int>(*i,(int)(i - sortedVal.begin())));
      } //save sort permutation
      transform(normVal.begin(), normVal.end(), sortedVal.begin(),
                normVal.begin(), plus<LargeNumber>()); // add up to normalized values
    }
    istart = *start; /// @todo move division out of loop
    if (!affyBioconductor) {
      //compute partial sums (reusing vector) (1,2,3)->(1,3,6)
      transform(normVal.begin(), normVal.end(),normVal.begin(), adder<LargeNumber>());
    }
    //divide by number of samples
    transform(normVal.begin(), normVal.end(),normVal.begin(), UnaryDivides<LargeNumber>((int)nSamples));
    // for each sample
    for(NumberStarStar currentStart=start; currentStart!=stop; currentStart++) {
      istart = *currentStart;
      // for each data point
      for(istop = istart + sampleSize; istart != istop; istart++) {
        // find all ranks with this value
        pair<typename Perm_t::iterator, typename Perm_t::iterator>
          range = perm[currentStart - start].equal_range(*istart);
        //
        int pos1 = range.first->second;
        int pos2 = (--range.second)->second;
        if (!affyBioconductor) {
          pos1--;
          // to work with partial sums one needs to subtract partial sum up to excluding pos1
          *istart=(Number)((normVal[pos2] - (pos1 < 0 ? 0 :normVal[pos1]))/(pos2 - pos1));
          //compute average value for these ranks
          // NOTE: This cast to "(Number)" should be left here.
          //       We want to convert the returned value to Number and avoid the warning.
          assert(pos1<pos2); // otherwise (pos2-pos1) is negative
        }
        else {
          *istart=(Number) normVal[(pos2 + pos1)/2];
          // median of even size, slightly off for odd ones
        }
      }
    }
  }
};


template <class NumberStarStar>
class MedianNormalization {
  typedef typename iterator_traits<NumberStarStar>::value_type NumberStar;
  typedef typename iterator_traits<NumberStar>::value_type Number;
public:
  /// @brief     Median normalization of vectors of vectors.
  /// @param     start     pointer to start of pointers of data
  /// @param     stop      pointer to end   of pointers of data
  /// @param     istop     pointer to end   of data in the first row
  /// @param     target    target to normalize to
  /// @param     subSampleSize subsample the data when computing the median of the data row.
  /// @remark    Like quantile norm, but with a target value.
  void operator()(NumberStarStar start, NumberStarStar stop, NumberStar istop, Number target, size_t subSampleSize = ULONG_MAX) {
    size_t sampleSize = istop - *start; // length of first vec
    assert(sampleSize>0); // JHG SIZE_T DEBUG
    //printf("MEDIAN SAMPLESIZE=%u   SUBSAMPLESIZE=%u\n",(int32_t)sampleSize,(int32_t)subSampleSize); // DEBUG
    for(; start < stop; start ++) {
      Number sampleMedian = median(*start, *start + sampleSize, subSampleSize);
      UnaryMultiplies<double, Number> umult(double(target)/double(sampleMedian));
      transform(*start, *start + sampleSize, *start, umult);
    }
  }
};


// Sketch normalization and supporting functions from here
//
// Sketches are a means of working with the data "out of core".  The
// idea is that we can extract a subset of the data (the sketch) and
// use these sketches to compute the what the normalization targets
// should be.  Once this has been figured out, the resulting sketch
// can be applied to the input data.
//


/**
 * Functor class for extracting a subset of data from a larger set.
 * See also: @ref sketchNorm
 */
template <class NumberStar1, class NumberStar2=NumberStar1>
class ExtractSketch {
  typedef typename iterator_traits<NumberStar1>::value_type Number1;
public:
  unsigned int requestedSketchSize;
  /// @brief     Creates a functor which creates sketches of a given size
  /// @param     sks  Size of the sketch to create
  ExtractSketch (int sks = INT_MAX): requestedSketchSize(sks){}
  /// @brief     Extracts a "sketch" from
  /// @param     start     beginning of data to sample
  /// @param     stop      end of data to sample
  /// @param     sketchStart  output iterator
  void operator()(NumberStar1 start, NumberStar1 stop, NumberStar2 sketchStart) {
    unsigned int sampleSize = stop - start;
    /// @todo this changes the slot
    unsigned int sketchSize = min (sampleSize, requestedSketchSize);
    // JHG: make a copy of the vector -- this one will be modified in place
    // we save 'sampleSize'-1 copies by doing this once.
    vector<Number1> sorted_vec(start, stop);
    // sort the data once.
    sort(sorted_vec.begin(),sorted_vec.end());
    // create the sketch of the data across the percentiles. Note that
    // we must use 64bit int to ensure that we don't overflow in integer math.
    // Also note that we use integer math to prevent rounding errors confounding
    // ceil and floor.
    for (uint64_t i = 0; i < sketchSize; i++, sketchStart++) {
      // if percentile falls in between two elements take average
      uint64_t id1 = (i * (sampleSize - 1)) / (sketchSize -1);
      uint64_t id2 = ((i * (sampleSize - 1)) + (sketchSize - 2)) / (sketchSize -1);
      *sketchStart = (sorted_vec[(unsigned int)id1] + sorted_vec[(unsigned int)id2])/2;
    }
  }
};

/// @brief     Combine all the sketches into a summary by averaging them.
/// See also: @ref sketchNorm
template <class NumberStarStar, class NumberStarOut, class LargeNumber>
class AverageSketch {
  typedef typename iterator_traits<NumberStarStar>::value_type NumberStar;
  typedef typename iterator_traits<NumberStar>::value_type Number;
public:
  /// @brief     The function operator
  /// @param     start     Start of vector of sketches
  /// @param     stop      End of vector of sketches
  /// @param     istop     End of the first sketch in the vector of sketches
  /// @param     result    Where the place the generated sketch
  void operator()(NumberStarStar start, NumberStarStar stop, NumberStar istop, NumberStarOut result) {
    for (int i = 0; i < istop - *start;  i++) {
      LargeNumber tmpAcc = 0;
      for (NumberStarStar sampleIt = start; sampleIt != stop; sampleIt ++) {
        // cast to LargeNumber so we can do "+=" without needing a operator defined.
        tmpAcc += (LargeNumber)(*(*(sampleIt) + i));
      }
      *(result + i) = (Number)(tmpAcc / (stop - start));
    }
  }
};

/**
 * Function to interpolate a value in the container with sStart, sEnd
 * into either the average sketch container or avgerage partial sums
 * vector depending on if using the middle values during ties or doing
 * true quantile normalization via partial sums container. useMiddle
 * corresponds to the method used by affy package in bioconductor. The
 * avgStart container is the same size as the sStart (sketch start)
 * while the avgPartSum has an additional 0 element at the beginning.
 */
template <class Number, class NumberStar, class LargeNumberStar>
Number
interpolate_qnorm(Number x , NumberStar sStart, NumberStar sEnd, 
                  NumberStar avgStart, NumberStar avgEnd, LargeNumberStar avgPartSumStart,
                  LargeNumberStar avgPartSumEnd, bool useMiddle = false, Number minVal=0.0, 
                  bool hardMin=true) {
  pair<NumberStar, NumberStar> sketchRange = equal_range(sStart, sEnd, x);
  Number delta = 0;
  double theta = 0;
  /* Two ways to handle the minimum extrapolating from sketch or
     supplying an artificial minimum value for
     interpolation. Mathmatically the best thing to do would be to
     extrapolate from the first two points in the sketch to the value
     requested. Unfortunately, this can lead to negative values which
     downstream methods have problems with. Thus, we allow the user to
     specify a "minVal" which will be can be tacked onto the sketch
     invisibly to ensure that we never go below that minimum value
     (usually 0). Essentially this interpolates between minVal and the
     smallest value in the sketch rather than extrapolate from the
     second smallest two values. */
  if(sketchRange.second == sStart && !hardMin) {
    assert(x < *sStart);
    // extrapolate
    delta = *sStart - x;
    theta = (*(avgStart + 1) - *avgStart)/(double)(*(sStart + 1) - *sStart);
    return (Number)(*avgStart - delta * theta);
  }
  else if(sketchRange.second == sStart && hardMin) {
    assert(x < *sStart);
    // interpolate between minVal and the smallest value in the sketch.
    delta = x - minVal;
    assert(minVal < (*avgStart) && minVal < (*sStart) && 
           "normalization::interpolate_qnorm() - minVal must be smaller than items in sketch already.");
    theta = (*(avgStart) - minVal) / (double)(*sStart - minVal);
    return (Number)(minVal + delta * theta);
  }
  // extrapolate off the right (greater) end...
  else if (sketchRange.first == sEnd) {
    assert(x > *(sEnd -1));
    delta = x - *(sEnd - 1);
    theta = (*(avgEnd - 1) - *(avgEnd - 2))/(double)((*(sEnd - 1) - *(sEnd - 2)));
    return (Number)(*(avgEnd - 1) + delta * theta);
  }
  // a single value, interpolate between two points...
  else if (sketchRange.first == sketchRange.second) {
    delta = x - *(sketchRange.first - 1);
    int offset = sketchRange.first - sStart;
    theta = (*(avgStart + offset) - *(avgStart + offset -1))/
      (double)(*(sStart + offset) - *(sStart + offset -1));
    return (Number)(*(avgStart + offset - 1)  + delta * theta);
  }
  // a range of buckets...
  // Here we have to handle the ties in the data appropriately.
  else {
    int offsetStart = sketchRange.first - sStart;
    int offsetStop = sketchRange.second - sStart;
    // result depends on how ties are resolved
    return (useMiddle ?
            // the casts are inside the "?" so both results have the same type.
            (Number)( *(avgStart + (int)((offsetStop + offsetStart - 1)/2)) ) :
            (Number)((*(avgPartSumStart + offsetStop) - *(avgPartSumStart + offsetStart))/(offsetStop - offsetStart))
            );
  }
}


/**
 * Functor class to determine normalized value given raw data point and the
 * vector of sketch values for that chip and average values for all sketches.
 * See also: @ref sketchNorm
 */
template <class NumberStar, class LargeNumber>
class InterpolateSketch{
  typedef typename iterator_traits<NumberStar>::value_type Number;
  vector<LargeNumber> avgSketchPartialSums;
public:
  NumberStar sketchStart, sketchStop, avgSketchStart, avgSketchStop;
  bool affyBioconductor; // Bioconductor compatibility mode
  
  /// @brief
  /// @param     sksta
  /// @param     sksto
  /// @param     avgsta
  /// @param     avgsto
  /// @param     ab
  /// @return
  /// @remarks
  /// The
  InterpolateSketch (NumberStar sksta, NumberStar sksto, NumberStar avgsta, NumberStar avgsto, bool ab = false):
    sketchStart(sksta), sketchStop(sksto), avgSketchStart(avgsta), avgSketchStop(avgsto), affyBioconductor(ab) {
    if(!ab) {
      avgSketchPartialSums.resize(avgSketchStop - avgSketchStart + 1);
      avgSketchPartialSums[0] = 0;
      transform(avgSketchStart, avgSketchStop, avgSketchPartialSums.begin() + 1, adder<LargeNumber>());
    }
  }
  Number operator()(Number x) {
    return interpolate_qnorm(x, sketchStart, sketchStop, avgSketchStart, avgSketchStop, 
                                     avgSketchPartialSums.begin(), avgSketchPartialSums.end(), affyBioconductor);
  }
};

/// @brief     A convenience function for making InterpolateSketch objects
template <class NumberStar, class LargeNumber>
InterpolateSketch<NumberStar, LargeNumber> makeInterpolateSketch (NumberStar sksta, NumberStar sksto, NumberStar avgsta, NumberStar avgsto,bool affyBioconductor = false) {
  return InterpolateSketch<NumberStar, LargeNumber>(sksta, sksto, avgsta, avgsto, affyBioconductor);
}

template <class NumberStarSample, class NumberStarSketch, class NumberStarAverageSketch, class LargeNumber>
class ApplySketch{
  typedef typename iterator_traits<NumberStarSample>::value_type Number;
  bool affyBioconductor;
public:
  NumberStarAverageSketch avgSketchStart, avgSketchStop;
  ApplySketch(NumberStarSketch sta, NumberStarSketch sto, bool ab = false): affyBioconductor(ab), avgSketchStart(sta), avgSketchStop(sto) {}
  NumberStarSample operator()(NumberStarSample start, NumberStarSample stop, NumberStarSketch skstart, NumberStarSketch skstop) {
    return transform(start, stop, start, makeInterpolateSketch<NumberStarSketch, LargeNumber>(skstart, skstop, avgSketchStart, avgSketchStop, affyBioconductor));
  }
};

template <class NumberStarSample, class NumberStarSketch, class NumberStarAverageSketch, class LargeNumber>
ApplySketch<NumberStarSample, NumberStarSketch, NumberStarAverageSketch, LargeNumber> makeApplySketch(NumberStarSketch sta, NumberStarSketch sto) {
  return ApplySketch<NumberStarSample, NumberStarSketch, NumberStarAverageSketch, LargeNumber> (sta,sto);
}

template <class NumberStarStar, class LargeNumber>
class SketchNormalization {
  typedef typename iterator_traits<NumberStarStar>::value_type NumberStar;
  typedef typename iterator_traits<NumberStar>::value_type Number;
public:

  /// @brief     in-core normalization with sketches.
  /// @param     start     start of vector of vectors
  /// @param     stop      end of vector of vectors
  /// @param     istop     end of first vector
  /// @param     sketchSizeMax the size of sketches to use
  /// @return    the data normalized in place
  /// @remarks
  ///   This provides the same interface as the other norm functions.
  ///   Note that this function is for reference as it does not save any memory.
  ///   (But it can save time.)  To save memory, you need to discard the celfile
  ///   input as they are converted into sketches.  Then read them back in
  ///   apply the sketch and write it out.

  NumberStarStar operator()(NumberStarStar start, NumberStarStar stop, NumberStar istop, int sketchSizeMax = INT_MAX) {
    typedef vector <typename vector<Number>::iterator > VectOfVectIter;
    typedef typename VectOfVectIter::iterator VectOfVectIterIter;

    VectOfVectIter sketches;            // Vector of iterators pointing to sketches.
    int sampleSize = istop - *start;    // Size of input vector.
    int sketchSize = min (sampleSize, sketchSizeMax);  // Size of sketch

    // For each chip extract the sketch.
    for(NumberStarStar i = start; i != stop; i++) {
      // Create a vector for sketch
      vector<Number> * tmpSketch = new vector<Number>(sketchSize);
      // Create extract sketch functor.
      ExtractSketch<NumberStar, typename vector<Number>::iterator > es(sketchSize);
      // Execute extract sketch functor
      es(*i, *i + sampleSize, tmpSketch->begin());
      // Save this sketch.
      sketches.push_back(tmpSketch->begin());
    }
    vector<Number> averageSketch(sketchSize);
    AverageSketch<VectOfVectIterIter, typename vector<Number>::iterator, LargeNumber> as;
    as(sketches.begin(), sketches.end(), (*sketches.begin()) + sketchSize, averageSketch.begin());
    ApplySketch<NumberStar, typename vector<Number>::iterator,typename vector<Number>::iterator, LargeNumber> appSketch(averageSketch.begin(), averageSketch.end());

    for(int i = 0;start + i != stop; i ++) {
      appSketch(*(start+i), *(start + i) + sampleSize,
                *(sketches.begin()+ i),*(sketches.begin()+ i) + sketchSize);
    }
    return start;
  }
};

#endif // __NORMALIZATION_H_
