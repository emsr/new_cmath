namespace __detail {

template<typename _Tp>
  _Tp
  __gegenbauer_poly(unsigned int __n, _Tp __alpha, _Tp __x)
  {
    _Tp __c0 = _Tp(1);
    if (__n == 0)
      return __c0;

    _Tp __c1 = _Tp(2) * __alpha * __x;
    if (__n == 1)
      return __c1;

    _Tp __c(0);
    for (unsigned int __nn = 2; __nn <= __n; ++__nn)
      {
        __c = (_Tp(2) * (_Tp(__nn) - _Tp(1) + __alpha) * __x * __c1
                    - (_Tp(__nn) - _Tp(2) + _Tp(2) * __alpha) * __c0) / _Tp(__nn);
        __c0 = __c1;
        __c1 = __c;
      }
    return __c;
  }

}
