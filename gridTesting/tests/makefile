# CC=g++-10
# std=c++2a
# CFLAGS=-I.
CPP=g++
# CPP=CC
# CPPflags=-O3 -std=c++2a -fopenmp
# CPPflags=-O3 -std=c++2a 
CPPflags=-std=c++2a -fopenmp

grid_test.x: grid_test.cpp ../grid.hpp ../Utils.hpp ../vec3.hpp
	${CPP} ${CPPflags} -o grid_test.x grid_test.cpp

# test: ColliderSingleCore.cpp 
# 	${CPP} ${CPPflags} -o testSC.x ColliderSingleCore.cpp

clean:
	rm -f *.o
	rm -f *.x
