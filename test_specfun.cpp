// ../bin/bin/g++ -std=c++11 -o test_specfun test_specfun.cpp


#include <iostream>
#include <iomanip>


#include "cmath"


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
                  << std::setw(12) << cosint(x)
                  << std::setw(12) << sinint(x)
                  << std::endl;
    }
}

