

#include <tr1/cmath>
#include <complex>
#include <stdexcept>


///
///  @brief This routine computes the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
///         integrals by continued fraction for positive argument.
///
template<typename _Tp>
  void
  __csint_cont_frac(_Tp __x, _Tp& __ci, _Tp& __si)
  {
    const unsigned int __max_iter = 100;
    const _Tp __eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();
    const _Tp __fp_min = std::numeric_limits<_Tp>::min();
    const _Tp __pi_2 = std::tr1::__detail::__numeric_constants<_Tp>::__pi_2();

    //  Evaluate Ci and Si by Lentz's modified method of continued fracions.
    std::complex<_Tp> __b(_Tp(1), __t);
    std::complex<_Tp> __c(_Tp(1) / __fp_min);
    std::complex<_Tp> __d(_Tp(1) / __b);
    std::complex<_Tp> __h(__d);
    unsigned int i = 2;
    while (true)
      {
        _Tp __a = -(i - 1) * (i - 1);
        __b += _Tp(2);
        __d = _Tp(1) / (__a * __d + __b);
        __c = __b + __a / __c;
        std::complex<_Tp> __del = __c * __d;
        __h *= __del;
        if (std::abs(__del.real() - _Tp(1))
          + std::abs(__del.imag()) < __eps)
          break;
        if (i > __max_iter)
          throw std::logic_error("Continued fraction evaluation failed in __cisi.");
        ++i;
      }
    __h *= std::polar(_Tp(1), -__t);
    __ci = -__h.real();
    __si = __pi_2 + __h.imag();

    return;
  }


///
///  @brief This routine computes the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
///         integrals by series summation for positive argument.
///
template<typename _Tp>
  void
  __csint_series(_Tp __x, _Tp& __ci, _Tp& __si)
  {
    const unsigned int __max_iter = 100;
    const _Tp __eps = _Tp(5) * std::numeric_limits<_Tp>::epsilon();
    const _Tp __fp_min = std::numeric_limits<_Tp>::min();
    const _Tp __gamma_e = std::tr1::__detail::__numeric_constants<_Tp>::__gamma_e();

    //  Evaluate Ci and Si by series simultaneously.
    _Tp __sumc(0), __sums(0);
    if (__t * __t < __fp_min)
      {
        //  Avoid underflow.
        __sumc = _Tp(0);
        __sums = __t;
      }
    else
      {
        //  Evaluate Si and Ci by series expansion.
        _Tp __sum(0);
        _Tp __sign(1), __fact(1);
        bool __odd = true;
        unsigned int __k = 1;
        while (true)
          {
            __fact *= __t / __k;
            _Tp __term = __fact / __k;
            __sum += __sign * __term;
            _Tp __err = __term / std::abs(__sum);
            if (__odd)
              {
                __sign = -__sign;
                __sums = __sum;
                __sum = __sumc;
              }
            else
              {
                __sumc = __sum;
                __sum = __sums;
              }
            if (__err < __eps)
              break;
            __odd = !__odd;
            ++__k;
            if (__k > __max_iter)
              throw std::logic_error("Series evaluation failed in __cisi.");
          }
      }
    __ci = __gamma_e + std::log(__t) + __sumc;
    __si = __sums;

    return;
  }


///
///  @brief This routine returns the cosine @f$ Ci(x) @f$ and sine @f$ Si(x) @f$
///         integrals as a pair.
///
///  The cosine integral is defined by:
///  @f[
///      Ci(x) = \gamma_E + \log(x) + \int_0^x dt \frac{\cos(t) - 1}{t}
///  @f]
///
///  The sine integral is defined by:
///  @f[
///      Si(x) = \int_0^x dt \frac{\sin(t)}{t}
///  @f]
///
template<typename _Tp>
  std::pair<_Tp, _Tp>
  __csint(_Tp __x)
  {
    _Tp __t = std::abs(__x);
    if (__t == _Tp(0))
      {
        __ci = -std::numeric_limits<_Tp>::infinity();
        __si = _Tp(0);
        return;
      }
    if (__t > _Tp(2))
      __csint_cont_frac(__x, __ci, __si);
    else
      __csint_series(__x, __ci, __si);

    if (__x < _Tp(0))
      __si = -__si;

    return std::make_pair(__ci, __si);
}


