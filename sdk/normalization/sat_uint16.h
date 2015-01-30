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

/// @file   sat_uint16.h
/// @brief  Saturating 16 bit integers.

#ifndef __SAT_UINT16_H_
#define __SAT_UINT16_H_

//! this will include stdint if possible.
#include "portability/affy-base-types.h"
//
#include <climits>
//

//! A more readable name.
#define SAT_UINT16_MAX USHRT_MAX

//! Convert an value to fit in the range of (0,SAT_UINT16_MAX)
#define SAT_CLAMP_VAL(x) ((in<0)?0:(in<SAT_UINT16_MAX)?(in):SAT_UINT16_MAX)


/// @brief          A class for saturating integers.
/// @param in          the value to make the sat_uint16 with.
/// @remark       This avoids overflow when preforming normalization with uint16
class sat_uint16 {
  uint16_t value;

 public:
  //! avoid a function call by using the macro
  //! Clamp does not assign a value, assign does.
  inline sat_uint16 clamp(int in) {
    return SAT_CLAMP_VAL(x);
  }
  inline void assign(int in) {
    value=SAT_CLAMP_VAL(in);
  }

  //! defaults to zero.
  sat_uint16()             { value=0;             }
  //! creation from a value
  sat_uint16(int   number) { assign(number);      }
  sat_uint16(uint64_t number) { assign((int)number); }
  sat_uint16(double number) { assign((int)number); }
  //! dont need to check range
  sat_uint16(const sat_uint16& number) { value=number.value; }

  //! Assignment from different types.
  sat_uint16& operator= (const double   number) { assign(int(number)); return *this; }
  sat_uint16& operator= (const int      number) { assign(int(number)); return *this; }
  sat_uint16& operator= (const uint64_t number) { assign(int(number)); return *this; }

  //! Basic math functions.
  sat_uint16 operator*(const double     number) { return clamp(int(value*number));  }
  sat_uint16 operator*(const int        number) { return clamp(value*number);       }
  sat_uint16 operator*(const sat_uint16 number) { return clamp(value*number.value); }
  sat_uint16 operator+(const double     number) { return clamp(int(value+number));  }
  sat_uint16 operator+(const int        number) { return clamp(value+number);       }
  sat_uint16 operator+(const sat_uint16 number) { return clamp(value+number.value); }
  sat_uint16 operator+(const uint64_t   number) { return clamp(int(value+number));  }
  sat_uint16 operator-(const double     number) { return clamp(int(value-number));  }
  sat_uint16 operator-(const int        number) { return clamp(value-number);       }
  sat_uint16 operator-(const sat_uint16 number) { return clamp(value-number.value); }
  sat_uint16 operator-(const uint64_t   number) { return clamp(int(value-number));  }
  sat_uint16 operator/(const double     number) { return clamp(int(value/number));  }
  sat_uint16 operator/(const int        number) { return clamp(value/number);       }

  //! Convert to an integer and compare.
  const int operator== (const int number)        const { return (number==value); }
  const int operator!= (const int number)        const { return (number!=value); }
  // these make it ambiguous -- arent clear they are needed.
  //const int operator== (const sat_uint16 number) const { return (int(number)==int(value); }
  //const int operator!= (const sat_uint16 number) const { return (number!=value); }

  const int operator<   (const int number) const { return (number< value);  }
  const int operator>   (const int number) const { return (number> value);  }
  const int operator<=  (const int number) const { return (number<=value);  }
  const int operator>=  (const int number) const { return (number>=value);  }

  //! almost standard inc/dec operators
  sat_uint16& operator++(int ignored) { if (value<SAT_UINT16_MAX) { value++; } ; return *this; }
  sat_uint16& operator--(int ignored) { if (0<value) { value--; }; return *this; }
  
  //
  sat_uint16& operator+=(const int number) { assign(value+number); return *this; }

  //! cast to different formats
  operator int()      const { return value;         }
  operator uint64_t() const { return value;         }
  operator float()    const { return float(value);  }
  operator double()   const { return double(value); }
};

#endif // __SAT_UINT16_H_
