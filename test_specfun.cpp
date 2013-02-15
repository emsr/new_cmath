// ../bin/bin/g++ -o test_chebyshev test_chebyshev.cpp


#include <iostream>
#include <iomanip>


#include "chebyshev.tcc"
#include "trig_integral.tcc"

double ci(double x) { double c, s; __cisi(x, c, s); return c; }
double si(double x) { double c, s; __cisi(x, c, s); return s; }


int
main()
{
    std::cout << std::endl;
    for (int i = 0; i <= 200; ++i)
    {
        double x = 0.01 * (i - 100);
        std::cout << std::setw(12) << x
                  << std::setw(12) << chebyshev_t(0, x)
                  << std::setw(12) << chebyshev_t(1, x)
                  << std::setw(12) << chebyshev_t(2, x)
                  << std::setw(12) << chebyshev_t(3, x)
                  << std::setw(12) << chebyshev_t(4, x)
                  << std::setw(12) << chebyshev_t(5, x)
                  << std::endl;
    }

    std::cout << std::endl;
    for (int i = 0; i <= 200; ++i)
    {
        double x = 0.01 * (i - 100);
        std::cout << std::setw(12) << x
                  << std::setw(12) << chebyshev_u(0, x)
                  << std::setw(12) << chebyshev_u(1, x)
                  << std::setw(12) << chebyshev_u(2, x)
                  << std::setw(12) << chebyshev_u(3, x)
                  << std::setw(12) << chebyshev_u(4, x)
                  << std::setw(12) << chebyshev_u(5, x)
                  << std::endl;
    }

    std::cout << std::endl;
    for (int i = 0; i <= 200; ++i)
    {
        double x = 0.01 * (i - 100);
        std::cout << std::setw(12) << x
                  << std::setw(12) << chebyshev_v(0, x)
                  << std::setw(12) << chebyshev_v(1, x)
                  << std::setw(12) << chebyshev_v(2, x)
                  << std::setw(12) << chebyshev_v(3, x)
                  << std::setw(12) << chebyshev_v(4, x)
                  << std::setw(12) << chebyshev_v(5, x)
                  << std::endl;
    }

    std::cout << std::endl;
    for (int i = 0; i <= 200; ++i)
    {
        double x = 0.01 * (i - 100);
        std::cout << std::setw(12) << x
                  << std::setw(12) << chebyshev_w(0, x)
                  << std::setw(12) << chebyshev_w(1, x)
                  << std::setw(12) << chebyshev_w(2, x)
                  << std::setw(12) << chebyshev_w(3, x)
                  << std::setw(12) << chebyshev_w(4, x)
                  << std::setw(12) << chebyshev_w(5, x)
                  << std::endl;
    }

    std::cout << std::endl;
    for (int i = 0; i <= 200; ++i)
    {
        double x = 0.01 * (i - 100);
        std::cout << std::setw(12) << x
                  << std::setw(12) << ci(x)
                  << std::setw(12) << si(x)
                  << std::endl;
    }
}

