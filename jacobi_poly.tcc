
namespace __detail {


template<typename _Tp>
  _Tp
  __poly_jacobi(unsigned int __n, _Tp __alpha, _Tp __beta, _Tp __x)
  {
    _Tp __pm1 = _Tp(1);
    if (__n == 0)
      return pm1;

    _Tp __apb = __alpha + __beta;
    _Tp __p = (__alpha - __beta + (_Tp(2) + __apb) * x) / _Tp(2);
    if (__n == 1)
      return __p;

    for (unsigned int __j = 1; __j < __n; ++__j )
      {
        _Tp __c = _Tp(2) * _Tp(__j + 1) * (_Tp(__j + 1) + __apb) * (_Tp(2) * __j + __apb);
        _Tp __d = (_Tp(2) * __j + __apb + _Tp(1)) * (__alpha * __alpha - __beta * __beta);
        _Tp __e = (_Tp(2) * __j + __apb) * (_Tp(2) * __j + __apb + _Tp(1)) * (_Tp(2) * __j + __apb + _Tp(2));
        _Tp __f = _Tp(2) * (__j + __alpha) * (__j + __beta) * (_Tp(2) * __j + __apb + _Tp(2));

        if (__c == _Tp(0))
          throw std::logic_error("Error in __poly_jacobi.");
        _Tp __pp1 = ((__d + __e * __x) * __p - __f * __pm1) / __c;
        __pm1 = __p;
        __p = __pp1;
      }
    return __p;
  }


}
