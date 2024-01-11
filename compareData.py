import os
import numpy as np
import sys
import h5py
relative_path = "../"
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'
sys.path.append(project_path)
# import utils as u
# import porosity_FD as p

def main():
	folder1 = "/home/lucas/Desktop/SpaceLab/jobs/test1/N_5/T_3/"
	folder2 = "/home/lucas/Desktop/SpaceLab_copy/SpaceLab/jobs/test1/N_5/T_3/"

	N = 5
	temp = 100
	show_FD_plots = False
	body = "2_R6e-05_v4e-01_cor0.63_mu0.1_rho2.25_k5e+00_Ha5e-12_dt4e-10_"
	for ind in range(N):
		file1 = f"{ind}_"

		if ind == 0:
			file2 = f"{body}"
		else:
			file2 = f"{ind}_{body}"

		print(file1)
		print(file2)

		f = h5py.File(folder1+file1+'data.h5','r')
		simData1 = np.array(f['/simData'][:])
		simData2 = np.loadtxt(folder2+file2+"simData.csv",delimiter=",",dtype=np.float64,skiprows=1).reshape(simData1.shape)

		constData1 = np.array(f['/constants'][:])
		constData2 = np.loadtxt(folder2+file2+"constants.csv",delimiter=",",dtype=np.float64,skiprows=0).reshape(constData1.shape)

		# print(constData1)
		# print(constData2)

		energyData1 = np.array(f['/energy'][:])
		energyData2 = np.loadtxt(folder2+file2+"energy.csv",delimiter=",",dtype=np.float64,skiprows=1).reshape(energyData1.shape)

		print("===================TESTING simData===================")
		err = False
		for i in range(simData1.size):
			if np.average([simData1[i],simData2[i]]) != 0.0:
				compare = abs((simData1[i]-simData2[i])/np.average([simData1[i],simData2[i]]))
				if (compare < 1e-5):
					continue
				else:
					err = True
					print(f"ERROR: {compare}")
		if not err:
			print("No errors")
		print("===================TEST Finished===================")
		print("===================TESTING constants===================")
		err = False
		for i in range(constData1.size):
			if np.average([constData1[i],constData2[i]]) != 0.0:
				compare = abs((constData1[i]-constData2[i])/np.average([constData1[i],constData2[i]]))
				if (compare < 1e-5):
					continue
				else:
					err = True
					print(f"ERROR: {compare}")
		if not err:
			print("No errors")
		print("===================TEST Finished===================")
		print("===================TESTING energy===================")
		err = False
		for i in range(energyData1.size):
			if np.average([energyData1[i],energyData2[i]]) != 0.0:
				compare = abs((energyData1[i]-energyData2[i])/np.average([energyData1[i],energyData2[i]]))
				if (compare < 1e-5):
					continue
				else:
					err = True
					print(f"ERROR: {compare}")
		if not err:
			print("No errors")
		print("===================TEST Finished===================")

if __name__ == '__main__':
	main()
