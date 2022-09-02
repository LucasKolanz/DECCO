# CC=g++-10
# std=c++2a
# CFLAGS=-I.

graph_test.o: graph_test.cpp graph_group.hpp
	g++-10 -std=c++2a -o graph_test.o graph_test.cpp


# test: ColliderSingleCore.cpp 
# 	g++-10 -std=c++2a -o testSC.o ColliderSingleCore.cpp

clean:
	rm -f *.o
