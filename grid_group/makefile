# CC=g++-10
# std=c++2a
# CFLAGS=-I.

grid_test.o: grid_test.cpp grid_group.hpp
	g++ -O3 -std=c++2a -o grid_test.o grid_test.cpp 

par: grid_test_parallel.cpp grid_group_parallel.hpp
	g++ -O3 -std=c++2a -o grid_test_parallel.o grid_test_parallel.cpp -fopenmp

# test: ColliderSingleCore.cpp 
# 	g++-10 -std=c++2a -o testSC.o ColliderSingleCore.cpp

clean:
	rm -f *.o
