import os
import json
import multiprocessing as mp
import subprocess
import numpy as np

def matMul(x,y):
	print("Starting matmul with size {}".format(x.shape))
	z = np.zeros(x.shape)
	for i in range(x.shape[0]):
		for j in range(y.shape[1]):
			for k in range(x.shape[1]):
				z[i,j] = x[i][k]*y[k][j]
	print(z)
	print(z.shape)
	print("Finished matmul with size {}".format(x.shape))

if __name__ == '__main__':

	runs_at_once = 2
	jobs=8

	e = 200

	a = np.full(shape=(e,e),fill_value=54.0,dtype=np.float64)
	b = np.full(shape=(e,e),fill_value=45.0,dtype=np.float64)
	c = np.full(shape=(int(e/3),int(e/3)),fill_value=54.0,dtype=np.float64)
	d = np.full(shape=(int(e/3),int(e/3)),fill_value=45.0,dtype=np.float64)


	func_input = [(a,b),(c,d)]
	print("Starting Pool")
	with mp.Pool(processes=runs_at_once) as pool:
		for i in range(0, jobs):
			input_data = func_input[i % len(func_input)]
			pool.apply_async(matMul, input_data)

		pool.close()
		pool.join()
	

