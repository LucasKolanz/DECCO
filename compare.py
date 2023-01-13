import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import csv
sys.path.append("/home/kolanzl/Desktop/SpaceLab")
import utils as u
import porosity_FD as p

if __name__ == '__main__':
	data_prefolder = []
	data_prefolder.append('/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/comparison_test')
	data_prefolder.append('/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/comparison_test')

	# temps = [10]
	temps = [1000]
	Nums = [100]
	# Nums = [300]
	# attempts = [19]
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
	# attempts = [i for i in range(1,5)]
	attempts = [1]
	data = []


	for folder in data_prefolder:
		std_dev = []
		std_err = []
		# attempts = [attempts[3]]
		porositiesabc = np.zeros((len(Nums),len(temps),len(attempts)),dtype=np.float64)
		porositiesKBM = np.zeros((len(Nums),len(temps),len(attempts)),dtype=np.float64) 
		FD_data = np.zeros((len(Nums),len(temps),len(attempts)),dtype=np.float64)
		porositiesabc[:] = np.nan
		porositiesKBM[:] = np.nan
		FD_data[:] = np.nan
		for i,temp in enumerate(temps):
			for n,N in enumerate(Nums):
				if N == 300:
					a = attempts300
				else:
					a = attempts
				for j,attempt in enumerate(a):
					
					data_folder = folder + str(attempt) + '/' + 'N_' + str(N) + '/T_' + str(temp) + '/'
					print(data_folder)
					count = 0
					for root_dir, cur_dir, files in os.walk(data_folder):
					    count += len(files)
					if count/3 > N:
						porositiesabc[n,i,j] = p.porosity_measure1(data_folder)
						porositiesKBM[n,i,j] = p.porosity_measure2(data_folder)

						o3dv = u.o3doctree(data_folder,overwrite_data=False)
						o3dv.make_tree()
						FD_data[n,i,j] = o3dv.calc_fractal_dimension(show_graph=False)
					else:
						porositiesabc[n,i,j] = np.nan
						porositiesKBM[n,i,j] = np.nan
						FD_data[n,i,j] = np.nan
		print("For Data from {}".format(folder))
		print("Porosity (abc):\t{}".format(porositiesabc))
		print("Porosity (KBM):\t{}".format(porositiesKBM))
		print("Fractal Dimension:\t{}".format(FD_data))
