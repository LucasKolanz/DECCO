# CC=g++-10
# std=c++2a
# CFLAGS=-I.

ColliderSingleCore.o: ColliderSingleCore.cpp ../dust_const.hpp ../ball_group.hpp ../grid_group.hpp
	g++-10 -std=c++2a -o ColliderSingleCore.o ColliderSingleCore.cpp


test: ColliderSingleCore.cpp 
	g++-10 -std=c++2a -o testSC.o ColliderSingleCore.cpp

clean:
	rm -f *.o
