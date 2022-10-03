# CC=g++-10
# std=c++2a
# CFLAGS=-I.

grid_test.o: grid_test.cpp grid_group.hpp
	g++-10 -std=c++2a -o grid_test.o grid_test.cpp


# test: ColliderSingleCore.cpp 
# 	g++-10 -std=c++2a -o testSC.o ColliderSingleCore.cpp

clean:
	rm -f *.o
