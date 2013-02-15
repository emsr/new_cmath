namespace __detail {

template<typename _Tp>
  _Tp
  __gegenbauer_poly(unsigned int __n, _Tp __alpha, _Tp __x)
  {
    _Tp __c0 = _Tp(1);
    if (n == 0)
      return __c0;

    _Tp __c1 = _Tp(2) * alpha * x;
    if (n == 1)
      return __c1;

    const _Tp __ninv = _Tp(1) / __n;
    _Tp __c(0);
    for (unsigned int nn = 2; nn <= n; ++nn)
      {
        __c = __ninv * (_Tp(2) * (_Tp(nn) - _Tp(1) + alpha) * x * c1 - (_Tp(nn) - _Tp(2) + _Tp(2) * alpha) * c0);
        __c0 = __c1;
        __c1 = __c;
      }
    return __c;
  }

}
