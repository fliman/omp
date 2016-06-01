
CXXFLAGS = -g -O0 

omptest:omptest.o 
	g++ $< -o $@

omptest.o:omptest.cpp omp.hpp
	g++ ${CXXFLAGS} $< -c -o $@



clean:
	rm omptest *.o
