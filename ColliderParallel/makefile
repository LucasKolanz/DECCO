# CC=g++-10
# std=c++2a
# CFLAGS=-I.

ColliderSingletoParallel.o: ColliderSingletoParallel.cpp ../dust_const.hpp ../ball_group.hpp ../grid_group/grid_group_parallel.hpp
	g++-10 -O3 -ftree-vectorize -std=c++2a -o ColliderSingletoParallel.o ColliderSingletoParallel.cpp -fopenmp


par: ColliderParallel.cpp ../vec3.hpp initializations.hpp ../json/single_include/nlohmann/json.hpp
	g++ -O3 -ftree-vectorize -std=c++2a -o ColliderParallel.o ColliderParallel.cpp


test: parallel_test.cpp ../ball_group.hpp ../vec3.hpp ../timing/timing.hpp ../grid_group/grid_group.hpp
	g++ -O3 -ftree-vectorize -std=c++2a -o parallel_test.o parallel_test.cpp -fopenmp

clean:
	rm -f *.o
