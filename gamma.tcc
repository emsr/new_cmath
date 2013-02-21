

namespace __detail {


template<typename _Tp>
  _Tp
  __gamma_p(_Tp a, _Tp x)
  {
    if (x < 0.0 || a <= 0.0)
      throw std::domain_error("Invalid arguments in routine gamma_p()");

    if (x < (a + 1.0))
      return __gamma_series(a, x).first;
    else
      return 1.0 - __gamma_cont_frac(a, x).first;
  }


template<typename _Tp>
  _Tp
  __gamma_q(_Tp a, _Tp x)
  {
    if (x < 0.0 || a <= 0.0)
      throw std::domain_error("Invalid arguments in routine gamma_q().");

    if (x < (a + 1.0))
      return 1.0 - __gamma_series(a, x).first;
    else
      return __gamma_cont_frac(a, x).first;
  }


template<typename _Tp>
  std::pair<_Tp, _Tp>
  __gamma_series(_Tp a, _Tp x)
  {
    const double EPS = 3.0e-7;
    const unsigned int ITMAX = 100;

    _Tp lngam = ln_gamma(a);

    if (x < _Tp(0))
      throw std::domain_error("Argument less than 0 in routine gamma_series().");
    else if (x == _Tp(0))
      return std::make_pair(_Tp(0), lngam);
    else
      {
        _Tp ap = a;
        _Tp del, sum;
        del = sum = _Tp(1) / a;
        for (unsigned int n = 1; n <= ITMAX; n++)
          {
            ap += _Tp(1);
            del *= x / ap;
            sum += del;
            if (std::abs(del) < EPS * std::abs(sum))
              {
                _Tp gamser = sum * std::exp(-x + a * std::log(x) - lngam);
                return std::make_pair(gamser, lngam);
              }
          }
        throw std::logic_error("a too large, ITMAX too small in routine gamma_series().");
      }
  }


template<typename _Tp>
  std::pair<_Tp, _Tp>
  __gamma_cont_frac(_Tp a, _Tp x)
  {
    const _Tp EPS = 3.0e-7;
    const unsigned int ITMAX = 100;

    _Tp lngam = ln_gamma(a);
    _Tp a1 = x;
    _Tp gold(0), fact(1), b1(1);
    _Tp a0(1), b0(0);
    for (unsigned int n = 1; n <= ITMAX; ++n)
      {
        _Tp an(n);
        _Tp ana = an - a;
        _Tp a0 = (a1 + a0 * ana) * fact;
        _Tp b0 = (b1 + b0 * ana) * fact;
        _Tp anf = an * fact;
        a1 = x * a0 + anf * a1;
        _Tp b1 = x * b0 + anf * b1;
        if (a1 != _Tp(0))
          {

            fact = _Tp(1) / a1;
            _Tp g = b1 * fact;
            if (std::abs(g - gold) / g < EPS)
              {

                _Tp gamcf = std::exp(-x + a * std::log(x) - lngam) * g;
                return std::make_pair(gamcf, lngam);
              }
            gold = g;
          }
      }
    throw std::logic_error("a too large, ITMAX too small in routine gamma_cont_fraction.");
  }


}
