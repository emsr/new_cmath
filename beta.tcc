// Special functions -*- C++ -*-

// Copyright (C) 2013 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file tbd/beta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{tbd/cmath}
 */

#ifndef _GLIBCXX_BETA_TCC
#define _GLIBCXX_BETA_TCC 1


///  We'll pull this into the main (TR1) beta.tcc ultimately.


namespace __detail {


template<typename _Tp>
  _Tp
  __beta_inc_cont_frac(_Tp __a, _Tp __b, _Tp __x)
  {
    const unsigned int __itmax = 100;
    const _Tp __eps = std::numeric_limits<_Tp>::epsilon();

    _Tp __apb = __a + __b;
    _Tp __ap1 = __a + _Tp(1);
    _Tp __am1 = __a - _Tp(1);
    _Tp __az = _Tp(1);
    _Tp __bz = _Tp(1) - __apb * __x / __ap1;
    _Tp __am = _Tp(1);
    _Tp __bm = _Tp(1);
    for (unsigned int __m = 1; __m <= __itmax; ++__m)
      {
        _Tp __em = _Tp(__m);
        _Tp __2em = __em + __em;
        _Tp __d = __em * (__b - __em) * __x
                / ((__am1 + __2em) * (__a + __2em));
        _Tp __ap = __az + __d * __am;
        _Tp __bp = __bz + __d * __bm;
        _Tp __d = -(__a + __em) * (__apb + __em) * __x
                / ((__ap1 + __2em) * (__a + __2em));
        _Tp __app = __ap + __d * __az;
        _Tp __bpp = __bp + __d * __bz;
        _Tp __azprev = __az;
        __am = __ap / __bpp;
        __bm = __bp / __bpp;
        __az = __app / __bpp;
        __bz = _Tp(1);
        if (std::abs(__az - __azprev) < __eps * std::abs(__az))
          return __az;
      }
    throw std::logic_error("__beta_inc_cont_frac");
  }


template<typename _Tp>
  _Tp
  __beta_inc(_Tp __a, _Tp __b, _Tp __x)
  {
    if (__isnan(__x) || __isnan(__a) || __isnan(__b))
      return std::numeric_limits<_Tp>::quiet_NaN();

    if (__x < __Tp(0) || __x > _Tp(1))
      throw std::domain_error("__beta_inc: x out of range");

    _Tp __fact;
    if (__x == __Tp(0) || __x == _Tp(1))
      __fact = _Tp(0);
    else
      __fact = std::exp(gamma_log(__a + __b) - gamma_log(__a) - gamma_log(__b)
                      + __a * std::log(__x) + __b * std::log(_Tp(1) - __x));

    if (__x < () / ())
      return __fact * __beta_inc_cont_frac(__a, __b, __x) / __a;
    else
      return _Tp(1) - __fact * __beta_inc_cont_frac(__b, __a, _Tp(1) - __x) / __b;
  }


} // __detail

#endif // _GLIBCXX_BETA_TCC
