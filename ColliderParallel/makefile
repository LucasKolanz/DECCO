# CC=g++-10
# std=c++2a
# CFLAGS=-I.

ColliderSingletoParallel.o: ColliderSingletoParallel.cpp ../dust_const.hpp ../ball_group.hpp ../grid_group_parallel.hpp
	g++-10 -std=c++2a -o ColliderSingletoParallel.o ColliderSingletoParallel.cpp -fopenmp


test: ColliderSingletoParallel.cpp 
	g++-10 -std=c++2a -o testSC.o ColliderSingletoParallel.cpp

clean:
	rm -f *.o
