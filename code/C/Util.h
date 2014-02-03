/* 
 * File:   Util.h
 * Author: Robby McKilliam
 *
 * Created on 31 January 2014, 9:48 PM
 */

#ifndef UTIL_H
#define	UTIL_H

#include <complex>

static constexpr double pi = 3.141592653589793238463;    

typedef std::complex<double> complexd;

/** Class for sorting and storing permuations */
class IndexedReal
{

public:
  IndexedReal( double value, int index ) : v(value), i(index) {}

  bool operator<(const IndexedReal& other) const {
  return v < other.v;
  }
  
  double v;
  int i;

};

#endif	/* UTIL_H */

