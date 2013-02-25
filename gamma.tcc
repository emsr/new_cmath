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

/** @file tbd/gamma.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{tbd/cmath}
 */

#ifndef _GLIBCXX_GAMMA_TCC
#define _GLIBCXX_GAMMA_TCC 1


///  We'll pull this into the main (TR1) gamma.tcc ultimately.

namespace __detail {


template<typename _Tp>
  std::pair<_Tp, _Tp>
  __gamma_series(_Tp __a, _Tp __x)
  {
    const double __eps = 3.0e-7;
    const unsigned int __itmax = 100;

    _Tp __lngam = std::lgamma(__a);

    if (__x < _Tp(0))
      throw std::domain_error("Argument less than 0 in routine gamma_series().");
    else if (__x == _Tp(0))
      return std::make_pair(_Tp(0), __lngam);
    else
      {
        _Tp __ap = __a;
        _Tp __del, __sum;
        __del = __sum = _Tp(1) / __a;
        for (unsigned int __n = 1; __n <= __itmax; ++__n)
          {
            __ap += _Tp(1);
            __del *= __x / __ap;
            __sum += __del;
            if (std::abs(__del) < __eps * std::abs(__sum))
              {
                _Tp __gamser = __sum * std::exp(-__x + __a * std::log(__x) - __lngam);
                return std::make_pair(__gamser, __lngam);
              }
          }
        throw std::logic_error("__gamma_series: a too large, ITMAX too small in routine.");
      }
  }


template<typename _Tp>
  std::pair<_Tp, _Tp>
  __gamma_cont_frac(_Tp __a, _Tp __x)
  {
    const _Tp __eps = 3.0e-7;
    const unsigned int __itmax = 100;

    _Tp __lngam = std::lgamma(__a);
    _Tp __a1 = __x, __b1(1);
    _Tp __gprev(0), __fact(1);
    _Tp __a0(1), __b0(0);
    for (unsigned int __n = 1; __n <= __itmax; ++__n)
      {
        _Tp __xn(__n);
        _Tp __xnma = __xn - __a;
        __a0 = (__a1 + __a0 * __xnma) * __fact;
        __b0 = (__b1 + __b0 * __xnma) * __fact;
        _Tp __xnfact = __xn * __fact;
        __a1 = __x * __a0 + __xnfact * __a1;
        __b1 = __x * __b0 + __xnfact * __b1;
        if (__a1 != _Tp(0))
          {
            __fact = _Tp(1) / __a1;
            _Tp __g = __b1 * __fact;
            if (std::abs(__g - __gprev) / __g < __eps)
              {
                _Tp __gamcf = std::exp(-__x + __a * std::log(__x) - __lngam) * __g;
                return std::make_pair(__gamcf, __lngam);
              }
            __gprev = __g;
          }
      }
    throw std::logic_error("__gamma_cont_fraction: a too large, ITMAX too small in routine.");
  }


template<typename _Tp>
  _Tp
  __gamma_p(_Tp __a, _Tp __x)
  {
    if (__x < 0.0 || __a <= 0.0)
      throw std::domain_error("Invalid arguments in routine gamma_p()");

    if (__x < __a + _Tp(1))
      return __gamma_series(__a, __x).first;
    else
      return _Tp(1) - __gamma_cont_frac(__a, __x).first;
  }


template<typename _Tp>
  _Tp
  __gamma_q(_Tp __a, _Tp __x)
  {
    if (__x < 0.0 || __a <= 0.0)
      throw std::domain_error("Invalid arguments in routine gamma_q().");

    if (__x < __a + _Tp(1))
      return _Tp(1) - __gamma_series(__a, __x).first;
    else
      return __gamma_cont_frac(__a, __x).first;
  }


template<typename _Tp>
  _Tp
  __log_pochhammer_u(_Tp __n, _Tp __x)
  {
    if (__isnan(__n) || __isnan(__x))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (__n == _Tp(0))
      return _Tp(0);
    else
      return std::lgamma(__x + __n) - std::lgamma(__x);
  }


template<typename _Tp>
  _Tp
  __pochhammer_u(_Tp __n, _Tp __x)
  {
    static const _Tp __log10(2.3025850929940456840179914546843642L);
    if (__isnan(__n) || __isnan(__x))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (__n == _Tp(0))
      return _Tp(1);
    else
      {
        _Tp __logpoch = std::lgamma(__x + __n) - std::lgamma(__x);
        if (std::abs(__logpoch)
            > std::numeric_limits<_Tp>::max_digits10 * __log10)
          return std::numeric_limits<_Tp>::infinity();
        else
          return std::exp(__logpoch);
      }
  }


template<typename _Tp>
  _Tp
  __log_pochhammer_l(_Tp __n, _Tp __x)
  {
    if (__isnan(__n) || __isnan(__x))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (__n == _Tp(0))
      return _Tp(0);
    else
      return std::lgamma(__x + 1) - std::lgamma(__x - __n + 1);
  }


template<typename _Tp>
  _Tp
  __pochhammer_l(_Tp __n, _Tp __x)
  {
    static const _Tp __log10(2.3025850929940456840179914546843642L);
    if (__isnan(__n) || __isnan(__x))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (__n == _Tp(0))
      return _Tp(1);
    else
      {
        _Tp __logpoch = std::lgamma(__x + 1) - std::lgamma(__x - __n + 1);
        if (std::abs(__logpoch)
            > std::numeric_limits<_Tp>::max_digits10 * __log10)
          return std::numeric_limits<_Tp>::infinity();
        else
          return std::exp(__logpoch);
      }
  }


} // namespace __detail

#endif // _GLIBCXX_GAMMA_TCC
