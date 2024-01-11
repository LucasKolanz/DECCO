# CC=g++-10
# std=c++2a
# CFLAGS=-I.
CPP=g++
# CPP=CC
# CPPflags=-O3 -std=c++2a -fopenmp
# CPPflags=-O3 -std=c++2a 
CPPflags=-std=c++2a 

randomGen.x: randomGen.cpp ../Utils.hpp
	${CPP} ${CPPflags} -o randomGen.x randomGen.cpp


getRand.x: getRand.cpp ../Utils.hpp
	${CPP} ${CPPflags} -o getRand.x getRand.cpp

# test: ColliderSingleCore.cpp 
# 	${CPP} ${CPPflags} -o testSC.x ColliderSingleCore.cpp

clean:
	rm -f *.o
	rm -f *.x
