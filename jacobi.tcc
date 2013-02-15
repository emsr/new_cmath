#include <array>

namespace __detail {


template<typename _Tp>
  void
  __jacobi_cnsndn(_Tp __k, _Tp __u, _Tp& __cn, _Tp& __sn, _Tp& __dn)
  {
    _Tp __a, __b, __c, __d;
    int l, bo;

    const unsigned int __N = 20;
    const _Tp __CA = 0.00001;

    std::array<_Tp, __N> __m;
    std::array<_Tp, __N> __n;

    if (__k != _Tp(0))
      {
        bool __bo = (__k < _Tp(0));
        if (__bo)
          {
            __d = _Tp(1) - __k;
            __k /= -_Tp(1) / __d;
            __u *= (__d = std::sqrt(__d));
          }
        _Tp __a = _Tp(1);
        __dn = _Tp(1);
        unsigned int __l = __N;
        for (unsigned int __i = 0; __i < __N; ++__i)
          {
            __l = __i;
            __m[__i] = __a;
            __n[__i] = (__k = std::sqrt(__k));
            __c = 0.5 * (__a + __k);
            if (std::abs(__a - __k) <= __CA * __a)
              break;
            __k *= __a;
            __a = __c;
          }
        __u *= __c;
        __cn = std::cos(__u);
        __sn = std::sin(__u);
        if (__sn != _Tp(0))
          {
            __a = __cn / __sn;
            __c *= __a;
            for (unsigned int __ii = __l; __ii + 1 >= 1; --__ii)
              {
                _Tp __b = __m[__ii];
                __a *= __c;
                __c *= (__dn);
                __dn = (__n[__ii] + __a) / (__b + __a);
                __a = __c / __b;
              }
            __a = _Tp(1) / std::sqrt(__c * __c + _Tp(1));
            __sn = std::copysign(__a, __sn);
            __cn = __c * __sn;
          }
        if (__bo)
          {
            __a = __dn;
            __dn = __cn;
            __cn = __a;
            __sn /= __d;
          }
      }
    else
      {
        __cn = _Tp(1) / std::cosh(__u);
        __sn = std::tanh(__u);
        __dn = __cn;
      }
  }

}
