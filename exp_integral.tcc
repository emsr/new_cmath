// Special functions -*- C++ -*-

// Copyright (C) 2006-2013 Free Software Foundation, Inc.
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

/** @file tr1/ell_integral.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{tr1/cmath}
 */

#ifndef _GLIBCXX_ELL_INTEGRAL_TCC
#define _GLIBCXX_ELL_INTEGRAL_TCC 1


namespace __detail {


template<typename _Tp>
  _Tp
  __logint(const _Tp __x)
  {
    if (__isnan(__x))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else if (std::abs(__x) == _Tp(1))
      return std::numeric_limits<_Tp>::infinity();
    else
      return std::tr1::expint(std::log(__x));
  }


template<typename _Tp>
  _Tp
  __coshint(const _Tp __x)
  {
    if (__isnan(__x))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else
      return (std::tr1::__detail::__expint_Ei(__x)
            - std::tr1::__detail::__expint_E1(__x)) / _Tp(2);
  }


template<typename _Tp>
  _Tp
  __sinhint(const _Tp __x)
  {
    if (__isnan(__x))
      return std::numeric_limits<_Tp>::quiet_NaN();
    else
      return (std::tr1::__detail::__expint_Ei(__x)
            + std::tr1::__detail::__expint_E1(__x)) / _Tp(2);
  }


} // namespace __detail

#endif // _GLIBCXX_ELL_INTEGRAL_TCC
