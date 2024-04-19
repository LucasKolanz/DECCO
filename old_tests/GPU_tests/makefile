# CC=g++-10
# std=c++2a
# CFLAGS=-I.
# CPP=CC
CPP=CC
# CPPflags= -std=c++2a -fopenmp -DPSEUDOPARTICLES=10 -DTHREADS=2 
# CPPflags= -std=c++2a -O2 -target-accel=nvidia80 -Mlarge_arrays# -Minfo=mp,accel -Minline#-fopenmp-targets=nvptx64 -Xopenmp-target=nvptx64 -march=host#-gpu=cc80
CPPflags= -std=c++2a -O2 -mp=gpu -acc=gpu -Mlarge_arrays -Minfo=accel# -Minline#-fopenmp-targets=nvptx64 -Xopenmp-target=nvptx64 -march=host#-gpu=cc80
ACCflags= -acc #-ta=tesla
 

GPUTestVec3.x: GPUTestVec3.cpp 
	${CPP} ${CPPflags} -Wl,-znoexecstack -o GPUTestVec3.x GPUTestVec3.cpp 
# 	clang++ -v -std=c++2a -O2 -fopenmp -fopenmp-targets=nvptx64 -o GPUTestVec3.x GPUTestVec3.cpp 
# 	-fopenmp-offload=sm_80 


acc: acc.cpp 
	${CPP} ${CPPflags} ${ACCflags} -Wl,-znoexecstack -o acc.x acc.cpp

compAcc: GPUTestACC.cpp 
	${CPP} ${CPPflags} ${ACCflags} -Wl,-znoexecstack -o GPUTestACC.x GPUTestACC.cpp

dataAcc: acc_noData.cpp 
	${CPP} ${CPPflags} ${ACCflags} -Wl,-znoexecstack -o acc_noData.x acc_noData.cpp

square: acc_sqma.cpp 
	${CPP} ${CPPflags} ${ACCflags} -Wl,-znoexecstack -o acc_sqma.x acc_sqma.cpp

update: acc_update.cpp 
	${CPP} ${CPPflags} ${ACCflags} -Wl,-znoexecstack -o acc_update.x acc_update.cpp

fast: acc_fastHybrid.cpp 
	${CPP} ${CPPflags} ${ACCflags} -Wl,-znoexecstack -o acc_fastHybrid.x acc_fastHybrid.cpp

# vec3: GPUTestVec3.cpp 
# 	${CPP} ${CPPflags} -o GPUTestVec3.x GPUTestVec3.cpp 

# breakup: GPUTestBu.cpp 
# 	${CPP} ${CPPflags} -o GPUTestBu.x GPUTestBu.cpp 

# 1: GPUTestv1.cpp 
# 	${CPP} ${CPPflags} -o GPUTestv1.x GPUTestv1.cpp 

clean:
	rm -f *.o
	rm -f *.x

# srun -N 1 -n 1 -c 4 -C gpu --cpu-bind=cores numactl --interleave=all ./GPUTest.x
#salloc --nodes 1 --qos interactive --time 01:00:00 --constraint  gpu --gpus 1 --account=m2651