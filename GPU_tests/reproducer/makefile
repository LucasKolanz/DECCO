
CPP=CC
CPPflags= -std=c++2a -O2 -mp=gpu -Mlarge_arrays -Minfo=mp,accel -Minline

works: works.cpp 
	${CPP} ${CPPflags} -Wl,-znoexecstack -o works.x works.cpp 


broken: broken.cpp 
	${CPP} ${CPPflags} -Wl,-znoexecstack -o broken.x broken.cpp 


clean:
	rm -f *.o
	rm -f *.x

