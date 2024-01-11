import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kolanzl/Desktop/SpaceLab")
import utils as u
from multiprocessing import Pool

def pool_fract_dim(folder):
	o3dv = u.o3doctree(folder,overwrite_data=False)
	o3dv.make_tree()
	return o3dv.calc_fractal_dimension(show_graph=False)
	# print(folder.split('/')[-2])
	# return float(folder.split('/')[-2].strip('T_'))

if __name__ == '__main__':
	max_processes = 4
	data_prefolder = '/home/kolanzl/Desktop/SpaceLab/jobs/tempVariance_attempt'
	temps = [3,10,30,100,300,1000]
	# temps = [3]
	# attempts = [2,3]
	# attempts = [2]
	attempts = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]

	fract_dims = np.zeros((len(temps)*len(attempts)))

	data_folders = []
	for atti,attempt in enumerate(attempts):
		for tempi,temp in enumerate(temps):
			df = data_prefolder + str(attempt) + '/' + 'T_' + str(temp) + '/'
			data_folders.append(df)
			
	for i in range(0,len(data_folders),max_processes):
		with Pool(processes=max_processes) as pool:
			fract_dims[i:i+max_processes] = pool.map(pool_fract_dim,data_folders[i:i+max_processes])
		# pool.join()

	fract_dims = fract_dims.reshape((len(attempts),len(temps)))

	fig, ax = plt.subplots(1,1)
	yerr = np.std(fract_dims,axis=0)/np.sqrt(len(attempts))
	ax.errorbar(temps,np.average(fract_dims,axis=0),yerr=yerr)
	ax.set_xscale('log')
	ax.set_xlabel('Temp (K)')
	ax.set_ylabel('Avg Fractal Dimension')
	ax.set_title('Average Fractal Dim over {} sims'.format(len(attempts)))
	plt.savefig("TotalFractDim.png")
	plt.show()