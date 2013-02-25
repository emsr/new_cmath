
all: test_specfun

test_specfun: test_specfun.cpp cmath *.tcc
	../bin/bin/g++ -std=c++11 -o test_specfun test_specfun.cpp
