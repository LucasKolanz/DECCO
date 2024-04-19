import os
# import subprocess

def main():
	# threads = [1,2,4,8,16,32,64]
	psize = [25000,26000,27000,28000,29000,30000]

	with open("timing2.csv",'a') as fp:
		fp.write("particles,time\n")

	for p in psize:
		os.system("CC -std=c++2a -O2 -fopenmp -mp=gpu -Mlarge_arrays -Minfo=mp,accel -Minline -o GPUTestVec3.x GPUTestVec3.cpp -DPSEUDOPARTICLES={}".format(p))
		os.system("srun -N 1 -n 1 -c 4 -C gpu --cpu-bind=cores numactl --interleave=all ./GPUTestVec3.x")

if __name__ == '__main__':
	main()