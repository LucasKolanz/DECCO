The error popps up when trying to use class methods from vec3.hpp in an
OpenMP target section that is in a different class. This is the case
even if the methods in vec3.hpp are wrapped in #pragma omp declare target.




To reproduce this error you need to:
	module load nvhpc

The version that is broken is made with:
	make broken

The version that works is made with:
	make works

I have been running this through an interactive queue with these options:
	salloc --nodes 1 --qos interactive --time 01:00:00 --constraint  gpu --gpus 1 --account=m2651

And using srun like this:
	srun -n 1 -N 1 -c 2 -G 1 {broken/works}.x