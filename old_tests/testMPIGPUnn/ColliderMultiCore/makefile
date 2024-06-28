# CC=g++-10
# std=c++2a
# CFLAGS=-I.
# CPP=CC
CPP=CC
# CPP=CC
# CPPflags=-O3 -std=c++2a -fopenmp
# CPPflags=-O3 -std=c++2a 
# CPPflags=-std=c++2a -fopenmp -O2 -DMPI_ENABLE -Mlarge_arrays -mp=gpu -gpu=cc80 -Minfo=mp,accel -Minline
CPPflags=-std=c++2a -fopenmp -O2 -Mlarge_arrays -mp=gpu -gpu=cc80 -Minfo=mp,accel -Minline

ColliderMultiCore.x: ColliderMultiCore.cpp ../dust_const_init.hpp ../Utils.hpp ball_group_multi_core.hpp ../../timing/timing.hpp
	${CPP} ${CPPflags} -o ColliderMultiCore.x ColliderMultiCore.cpp

# test: ColliderSingleCore.cpp 
# 	${CPP} ${CPPflags} -o testSC.x ColliderSingleCore.cpp

clean:
	rm -f *.o
	rm -f *.x
