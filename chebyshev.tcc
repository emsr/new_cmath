namespace __detail {


template<typename _Tp>
  _Tp
  __chebyshev_recur(unsigned int __n, _Tp __x, _Tp __c0, _Tp __c1)
  {
    _Tp __c = _Tp(0);
    for (unsigned int __j = 1; __j < __n; ++__j)
    {
      __c = _Tp(2) * __x * __c1 - __c0;
      __c0 = __c1;
      __c1 = __c;
    }
    return __c;
  }

template<typename _Tp>
  _Tp
  __chebyshev_t(unsigned int __n, _Tp __x)
  {
    _Tp __t0 = _Tp(1);
    if (__n == 0)
      return __t0;

    _Tp __t1 = __x;
    if (__n == 1)
      return __t1;

    return __chebyshev_recur(__n, __x, __t0, __t1);
  }

template<typename _Tp>
  _Tp
  __chebyshev_u(unsigned int __n, _Tp __x)
  {
    _Tp __u0 = _Tp(1);
    if (__n == 0)
      return __u0;

    _Tp __u1 = _Tp(2) * __x;
    if (__n == 1)
      return __u1;

    return __chebyshev_recur(__n, __x, __u0, __u1);
  }

template<typename _Tp>
  _Tp
  __chebyshev_v(unsigned int __n, _Tp __x)
  {
    _Tp __v0 = _Tp(1);
    if (__n == 0)
      return __v0;

    _Tp __v1 = _Tp(2) * __x - _Tp(1);
    if (__n == 1)
      return __v1;

    return __chebyshev_recur(__n, __x, __v0, __v1);
  }

template<typename _Tp>
  _Tp
  __chebyshev_w(unsigned int __n, _Tp __x)
  {
    _Tp __w0 = _Tp(1);
    if (__n == 0)
      return __w0;

    _Tp __w1 = _Tp(2) * __x + _Tp(1);
    if (__n == 1)
      return __w1;

    return __chebyshev_recur(__n, __x, __w0, __w1);
  }


}
