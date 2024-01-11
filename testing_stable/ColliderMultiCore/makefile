# CC=g++-10
# std=c++2a
# CFLAGS=-I.
CPP=g++
# CPP=CC
# CPPflags=-O3 -std=c++2a -fopenmp
# CPPflags=-O3 -std=c++2a 
CPPflags=-std=c++2a 

ColliderMultiCore.x: ColliderMultiCore.cpp ../ball_group.hpp
	${CPP} ${CPPflags} -o ColliderMultiCore.x ColliderMultiCore.cpp

# test: ColliderSingleCore.cpp 
# 	${CPP} ${CPPflags} -o testSC.x ColliderSingleCore.cpp

clean:
	rm -f *.o
	rm -f *.x
