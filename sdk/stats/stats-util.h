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

#ifndef __STATS_UTIL_H_
#define __STATS_UTIL_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
//
// "values.h" is not needed by the code which follows.
// Also, it has been replaced by "features.h" and "limits.h"
//#ifndef CYGWIN
//#include <values.h>
//#endif

using namespace std;

template <class T>
class UnaryPlus: public plus<T> {
 public:
  UnaryPlus(const T& x) {
    plusConst = x;
  }
  T operator()(const T& x) {
    return x + plusConst;
  }
  T plusConst;
};

template <typename ArgType, typename ResultType = ArgType>
  class UnaryMultiplies: public unary_function<ResultType, ResultType> {
 public:
  UnaryMultiplies(const ArgType& x) {
    multConst = x;
  }
  ResultType operator()(const ResultType& x) {
    // the cast to double is to force "operator*(sat_int16,double) to be used.
    return (ResultType) (double(x) * double(multConst));
  }
  ArgType multConst;
};

template <class T>
class UnaryDivides: public divides<T> {
 public:
  UnaryDivides(const T& x) {
    divConst = x;
  }
  T operator()(const T& x) {
    return x / divConst;
  }
  T divConst;
};

template <class RAIterator>
typename iterator_traits<RAIterator>::value_type
max_value(RAIterator start, RAIterator stop) {
  return *max_element(start, stop);
}

template <class RAIterator>
void ratio(RAIterator start, RAIterator stop) {
  UnaryMultiplies<typename iterator_traits<RAIterator>::value_type> div(1 / *start);
  transform(start, stop, start, div);
}

template <class Container>
Container & 
difference(Container & first, Container & second) {
  transform(first.begin(), first.end(), second.begin(), first.begin(), minus<typename Container::value_type>());
  return first;
}

template <class RAIterator>
class Ratio: public binary_function<RAIterator, RAIterator, typename iterator_traits<RAIterator>::value_type>  {
 public:
  typename iterator_traits<RAIterator>::value_type denominator;
  Ratio(typename iterator_traits<RAIterator>::value_type x): denominator(x) {};
  void operator()(RAIterator start, RAIterator stop) {
    transform(start, stop, start, UnaryMultiplies<typename iterator_traits<RAIterator>::value_type>(1/denominator));
  };
};

template <class RAIterator>
class Add: public binary_function<RAIterator, RAIterator, typename iterator_traits<RAIterator>::value_type>  {
 public:
  typename iterator_traits<RAIterator>::value_type constant;
  Add(typename iterator_traits<RAIterator>::value_type x): constant(x) {};
  void operator()(RAIterator start, RAIterator stop) {
    transform(start, stop, start, UnaryPlus<typename iterator_traits<RAIterator>::value_type>(constant));
  };
};

template <class RAIterator1, class RAIterator2>
typename iterator_traits<RAIterator1>::value_type
scalar_product(RAIterator1 start1, RAIterator1 stop1, RAIterator2 start2 ) {
  typename iterator_traits<RAIterator1>::value_type sum = 0;
  for (;start1 < stop1; start1++, start2++) {
    sum += *start1 * *start2;
  }
  return sum;
};

template <class Container1, class Container2>
Container1 insert(Container1 into, Container2 from) {
  into.insert(into.end(), from.begin(), from.end());
  return into;
} 


// template <class Container>
// class SubSeqErase: public unary_function<Container, void> {
// public:
//   vector<int> selection;
//   SubSeqErase(typename vector<int>::iterator start, typename vector<int>::iterator stop): selection(start, stop) {};
//   void operator()(Container C){
//     for(i=selection.begin(); i!=selection.end(); i++) {
//       if(C.begin()+*i < C.end()){
// 	C.erase(C.begin() + *i);
//       }
//     }
//   }
// };


template <class RAIterator1, class RAIterator2>
class ScalarProduct: public binary_function<RAIterator1, RAIterator1, typename iterator_traits<RAIterator1>::value_type >  {
 public:
  vector<typename iterator_traits<RAIterator1>::value_type> v2;
  ScalarProduct(RAIterator2 start, RAIterator2 stop) :v2(start, stop){};
  typename iterator_traits<RAIterator1>::value_type 
    operator()(RAIterator1 start, RAIterator2 stop) {
    return scalar_product(start, stop, v2.begin());
  }
};

template <class T>
class adder : public unary_function<T, void> {
public:
  adder() : result(0) {}
  T result;
  T operator()(T x) { 
    T last = result;
    result += x; 
    assert(((x >= 0 && result >= last) || (x <= 0 && result <= last)) && "Possible overflow in adder<>");
    return result;
  }
};

template <class T>
T square(T x) {
  return x * x;
}

template <class T>
class square_adder : public unary_function<T, void> {
public:
  square_adder() : result(0) {}
  T result;
  void operator()(T x) { 
    T last = result;
    result += square(x); 
    assert( result >= last && "Possible overflow in square_adder<>");
  }
};

template <class T>
class resizer : public unary_function<T, void> {
public:
  long size;
  resizer(long s) {
    size = s;
  }
  void operator()(T x) { x.resize(size); }
};

template <class T>
class deleter : public unary_function<T*, void> {
 public:
  void operator()(T* ptr) {
    delete ptr;
  }
};

template <class RAIterator>
void exp(RAIterator start, RAIterator stop) {
  typedef typename iterator_traits<RAIterator>::value_type Real;
  transform(start, stop, start, ptr_fun((Real (*)(Real)) exp));
}

template <class RAIterator>
void log(RAIterator start, RAIterator stop) {  
typedef typename iterator_traits<RAIterator>::value_type Real;
  transform(start, stop, start, (Real (*)(Real)) log);
}

template <class Number>
Number pos_log(Number x) {
  return log(max(x, (Number)1.0));
}

template <class RAIterator>
void positive_log(RAIterator start, RAIterator stop) {  
typedef typename iterator_traits<RAIterator>::value_type Real;
  transform(start, stop, start, (Real (*)(Real)) pos_log);
}

/*class exponential : public unary_function<double, double> {
public:
  exponential() {}
  double operator()(double x) { 
    return exp(x); 
  }
};

class cexponential: 
public unary_function<vector<double>, vector<double> > {
 public:
  cexponential() {};
  vector<double> operator() (vector <double> x) {
    vector<double> result (x.size());
    transform(x.begin(), x.end(), result.begin(), exponential());
    return result;
  }
};

class logarithm : public unary_function<double, void> {
 public:
  logarithm() {}
  double operator()(double x) { 
    return log((x > MINDOUBLE ? x: MINDOUBLE));
  }
};
*/

template<class OutputIterator>
class NumericSequence {
  typedef typename iterator_traits<OutputIterator>::value_type Number;
public:
  OutputIterator operator()(Number offset, Number increment, int size, OutputIterator result) {
    *result = offset; result ++;
    for (; result < result + size; result ++) {
      *result = *(result - 1) + increment;
    }
    return result;
  }
};


/*
template<class T> 
vector<T> sequence(long size, T start, T inc = 1) {
  vector <T> result;
  for (long i = 0; i < size; i++) {
    result.insert(start, result.end);
    start += inc;
  }
  return result;
}
*/

template<class T> 
vector<T> intSequence(long size, T start = 0, T inc = 1, T mod = 0, int rep = 1) {
  vector <T> result;
  for (long i = 0; i < size/rep; i++) {
    for ( int j = 0; j < rep; j++) {
      result.push_back(mod == 0 ? start: start % mod);
    }
    start += inc;
  }
  return result;
}

template <class RandomAccessIterator>
string join (RandomAccessIterator start, RandomAccessIterator stop, string sep = "\t") {
  if (start < stop) {
    ostringstream oneoSstream;
    for(; start < stop-1; start++) {
      oneoSstream << *start << sep;
    }
    oneoSstream << *start;
    return oneoSstream.str();
  }
  else return "";
}
 
template<class T> 
class Print : public unary_function<T, void> {
 public:
  string sep;
  ostream & oneOstream;
  Print(ostream & ost = cout) : sep("\t"), oneOstream(ost) {}
  Print(string s, ostream & ost=cout): sep(s), oneOstream(ost) {
    //    sep = s;
    //oneOstream = ost;
  }
  void operator()(const T & x) { 
    oneOstream << x << sep;
  }
  void operator() (const T & x, ostream & ost) {
    ost << x << sep;
  }
  void operator()(const T & x, const string & filename) {
    ostream ost(filename.c_str());
    ost << x;
  }
};

template <class RandomAccessIterator>
void shuffle(RandomAccessIterator first, RandomAccessIterator last,
	     vector < typename iterator_traits<RandomAccessIterator>::difference_type> permutation) {
  vector <typename iterator_traits<RandomAccessIterator>::value_type> temp (first, last);
  typename vector<typename iterator_traits<RandomAccessIterator>::difference_type>::difference_type i = 0;
  for (RandomAccessIterator current = first; current != last; current ++, i++) {
    *current = temp[permutation[i]]; 
  }
}  

template <class T>
T RandomizedPartition(T start, T stop) {
typename iterator_traits<T>::value_type pivot = 
    * (start + (typename iterator_traits<T>::difference_type)
       (((stop - start)-1) * ((rand() + 0.0)/RAND_MAX)));
 //--start; // keep this! will be incremented once before use. 
 bool first = true; /* the "first" flag is a hack as the microsoft 
					   stl doesn't like decrementing start outside of 
					   range in debug mode (even if not referenced). */
 while (true) {
    do {
      --stop;
    } 
    while (*stop > pivot);
    do {
		if(first) 
			first = false;
		else 
	 	  ++start;
    }
    while (*start < pivot);
    if (start < stop) {
      swap (*start, *stop);
    }
    else {
      return stop;
    }
  }
}

template <class T>
typename iterator_traits<T>::value_type
RandomizedSelect(T start, T stop, 
		 typename iterator_traits<T>::difference_type i) {
  /*  vector<typename iterator_traits<T>::value_type> tmp(start,stop);
  sort(tmp.begin(), tmp.end());
  copy(start, stop, ostream_iterator<double>(cerr, " "));
  cerr << "\npos "<< i << " value " << tmp[i] << " size " << stop - start << endl;*/
  if (start == stop - 1) {
    return *start ;
  }
  T pivot = RandomizedPartition(start, stop);
  //   copy(start, stop, ostream_iterator<double>(cerr, " "));
  //   cerr << "\n";
  typename iterator_traits<T>::difference_type k = pivot - start;
  //  cerr << "i: " << i << " k: " << k << " pivot: " << *pivot << endl;
  if (i <= k) {
    return RandomizedSelect(start, pivot + 1, i);
  }
  else {
    return RandomizedSelect(pivot + 1, stop, i - k -1);
  }
}


template <class T>
T next(T x) {
  T tmp(x);
  tmp++;
  return  tmp;
}

template <class T>
T prev(T x) {
  T tmp(x);
  tmp--;
  return  tmp;
}

template <class RAit>
class skip_iterator: public RAit {
 public:
  skip_iterator(RAit it, int skip): RAit(it), n(skip) {} 
  skip_iterator<RAit> operator ++() {
    *this += n;
    return *this;
  }
  skip_iterator<RAit> operator ++(int) {
    skip_iterator<RAit> tmp = *this;
    *this += n;
    return tmp;
  }
  skip_iterator<RAit> operator --() {
    *this -= n;  
    return *this;
  }
  skip_iterator<RAit> operator --(int) {
    skip_iterator<RAit> tmp = *this;
    *this -= n;  
    return tmp;
  }
  skip_iterator<RAit> operator + (long jump) {
    return (RAit::operator+ (n * jump));
  } 
  skip_iterator<RAit> operator - (long jump) {
    return (*this - n * jump);
  } 
  long operator - (skip_iterator<RAit> sub) {
    return (RAit(*this) - RAit(sub))/n;
  }
  private :
    int n;
};

template <class RAit>
skip_iterator<RAit> make_skip_iterator(RAit it, int skip) {
  return (skip_iterator<RAit>(it,skip));
}


// for_each.  Apply a function to every element of a range.
// this is my binary version. Should be in the STL to start with
template <class _InputIter1, class _InputIter2, class _Function>
_Function for_each(_InputIter1 __first1, _InputIter1 __last1, _InputIter2 __first2, _Function __f) {
  for ( ; __first1 != __last1; ++__first1, ++__first2)
    __f(*__first1, *__first2);
  return __f;
}


/// @todo complete multiple_range_iterator class and use it for zero-copy implementaton of bcelnorm 
/* 
template class<Iterator>
class multiple_range_iterator: 
  public: boost::iterator_adaptor<
  multiple_range_iterator<Iterator::value_type>,
  Iterator> {
private:
  struct enabler {};
  typedef boost::iterator_adaptor<
    multiple_range_iterator<Iterator::value_type>,
    Iterator> super_t ;
public:
  multiple_range_iterator
*/

#endif // __STATS_UTIL_H_
