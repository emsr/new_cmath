

double
gamma_series(double a, double x)
{
  const double EPS = 3.0 * std::numeric_limits<double>::epsilon();
  const int ITMAX = 100;

  if (x < 0.0)
    throw ("Argument less than 0 in routine gamma_series().");
  else if (x == 0.0)
    return 0.0;
  else
    {
      double ap = a;
      double del = sum = 1.0 / a;
      double sum = 0.0;
      for (int n = 1; n <= ITMAX; n++)
        {
          ap += 1.0;
          del *= x / ap;
          sum += del;
          if (std::abs(del) < EPS * std::abs(sum))
            return sum * exp(-x + a*std::log(x) - std::lgamma(a));
        }
      throw ("Argument a too large, ITMAX too small in routine gamma_series().");
    }
}


double
gamma_cont_fract(double a, double x)
{
  const double EPS = 3.0 * std::numeric_limits<double>::epsilon();
  const int ITMAX = 100;

  double fact = 1.0, gold = 0.0;
  double a0 = 1.0, a1 = x;
  double b0 = 0.0, b1 = 1.0;
  for (int n = 1; n <= ITMAX; ++n)
    {
      double an = static_cast<double>(n);
      double ana = an - a;
      a0 = (a1 + a0 * ana) * fact;
      b0 = (b1 + b0 * ana) * fact;
      double anf = an * fact;
      a1 = x * a0 + anf * a1;
      b1 = x * b0 + anf * b1;
      if (a1 != 0.0)
        {
          fact = 1.0 / a1;
          double g = b1 * fact;
          if  (std::abs(g - gold) / g < EPS)
            return std::exp(-x + a * std::log(x) - std::lgamma(a)) * g;
          gold = g;
        }
    }
  throw ("Argument a too large, ITMAX too small in routine gamma_series().");
}


double
gamma_p(double a, double x)
{
  if  (x < 0.0 || a <= 0.0)
    throw ("Invalid arguments in routine gamma_p()");

  if  (x < a + 1.0)
    return gamma_series(a, x);
  else
    return 1.0 - gamma_cont_fract(a, x);
}


double
gamma_q(double a, double x)
{
  if (x < 0.0 || a <= 0.0)
    throw ("Invalid arguments in routine gamma_q().");

  if (x < a + 1.0)
    return 1.0 - gamma_series(a, x);
  else
    return gamma_cont_fract(a, x);
}
