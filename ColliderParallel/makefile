# CC=g++-10
# std=c++2a
# CFLAGS=-I.



ColliderParallel.o: ColliderParallel.cpp ../vec3.hpp ballSim.hpp ../json/single_include/nlohmann/json.hpp
	g++ -O3 -ftree-vectorize -std=c++2a -o ColliderParallel.o ColliderParallel.cpp -ltbb

stp: ColliderSingletoParallel.cpp ../dust_const.hpp ../ball_group.hpp ../grid_group/grid_group_parallel.hpp
	g++-10 -O3 -ftree-vectorize -std=c++2a -o ColliderSingletoParallel.o ColliderSingletoParallel.cpp -fopenmp

test: parallel_test.cpp ../ball_group.hpp ../vec3.hpp ../timing/timing.hpp ../grid_group/grid_group.hpp
	g++ -O3 -ftree-vectorize -std=c++2a -o parallel_test.o parallel_test.cpp -fopenmp

clean:
	rm -f *.o
