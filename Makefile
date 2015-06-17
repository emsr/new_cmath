
all: test_specfun

test_specfun: test_specfun.cpp cmath *.tcc
	/home/ed/bin/bin/g++ -std=c++11 -o test_specfun test_specfun.cpp

test:
	LD_LIBRARY_PATH=/home/ed/bin/lib64:$$LD_LIBRARY_PATH ./test_specfun > test_specfun.txt
