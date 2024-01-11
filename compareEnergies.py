import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab")
import utils as u
import porosity_FD as p

def main():
	base = os.getcwd() + "/jobs/"
	# folder1 = base + "multiCoreTest4/"
	# folder1 = base + "singleCoreComparison/"
	# folder1 = base + "multiCoreTest1/"
	# folder2 = base + "singleCoreComparison2/"
	# folder2 = base + "multiCoreTest9/"
	# folder2 = base + "OpenMpParallel1/"
	# folder2 = base + "vdwImprovement1/"
	# folder2 = base + "openMPallLoops1/"
	# folder2 = base + "piplineImprovements1/"
	# folder2 = base + "dynamicAffinityTests_th2_1/affinity_0-1/"
	# folder2 = base + "affinityTests_th2_1/affinity_0-1/"
	# folder2 = base + "affinityTests_th2_1/affinity_0-8/"
	# folder2 = base + "affinityTests_th8_1/affinity_0-1-2-3-4-5-6-7-8/"
	# folder2 = base + "singleCoreComparison_COPY7/"
	# folder1 = "/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab/jobs/accuracyTest11/N_10/T_100/"
	# folder2 = "/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab/jobs/accuracyTest15/N_10/T_100/"
	folder1 = "/global/homes/l/lpkolanz/SpaceLab/jobs/full2400Collide1/thread_32/"
	folder2 = "/pscratch/sd/l/lpkolanz/SpaceLab/testHybrid/jobs/fullCompHybrid1/node_1/"
	folder2 = "/pscratch/sd/l/lpkolanz/SpaceLab/testSqMat/jobs/fullCompSqMa1/node_1/"

	# max_ind = -1
	# for file in os.listdir(folder1):
	# 	if file[-4:] == ".csv":
	# 		try:
	# 			ind = int(file.split('_')[0])
	# 			if ind > max_ind:
	# 				max_ind = ind
	# 				body = file.split('_')[1:-1]
	# 		except ValueError:
	# 			continue

	# file1 = str(max_ind)+'_'+'_'.join(body)+"_simData.csv"
	# # file1 = str(max_ind)+'_'+'_'.join(body)+"_simData_COPY.csv"
	# # # file1 = "9_simData.csv"

	max_ind = -1
	for file in os.listdir(folder2):
		if file[-4:] == ".csv":
			try:
				ind = int(file.split('_')[0])
				if ind > max_ind:
					max_ind = ind
					body = file.split('_')[1:-1]
			except ValueError:
				continue

	# file2 = str(max_ind)+'_'+'_'.join(body)+"_simData.csv"
	# file2 = "9_simData.csv"

	

	N = 5
	temp = 100
	show_FD_plots = False
	allEnergy = []
	for ind in [1]:
		f1 = "{}_{}_simData.csv".format(ind,'_'.join(body))
		f2 = "{}_{}_simData.csv".format(ind,'_'.join(body))
		# f2 = "{}_simData.csv".format(ind)
		KE=[]
		PE=[]
		Etot=[]
		p=[]
		L=[]
		COM=[]
		for data_folder in [folder1,folder2]:
			print("========================================")
			print(data_folder)
			count = 0
			# for root_dir, cur_dir, files in os.walk(data_folder):
			#     count += len(files)
			# if count/3 > N:
			energy = u.get_last_line_energy(data_folder,ind)
			allEnergy.append(energy)
			COM.append(u.COM(data_folder,ind))
			PE.append(energy[1])
			KE.append(energy[2])
			Etot.append(energy[3])
			p.append(energy[4])
			L.append(energy[5])


		print("================Comparing {}================".format(ind))
		print("Data 1: {}".format(allEnergy[0]))
		print("Data 2: {}".format(allEnergy[1]))
		print("COM difference : {}".format(np.diff(COM,axis=0)))	
		print("PE difference  : {}".format(np.diff(PE)))	
		print("KE difference  : {}".format(np.diff(KE)))	
		print("Etot difference: {}".format(np.diff(Etot)))	
		print("p difference   : {}".format(np.diff(p)))	
		print("L difference   : {}".format(np.diff(L)))	
		print("====================END====================")

		# data1 = np.loadtxt(folder1+f1,delimiter=",",dtype=np.float64,skiprows=1)[-1]
		# data2 = np.loadtxt(folder2+f2,delimiter=",",dtype=np.float64,skiprows=1)[-1]
		# fig, ax = plt.subplots(1,1,figsize=(15,7))
		# ax.plot(np.arange(len(data1)),np.subtract(data2,data1))
		# plt.show()

if __name__ == '__main__':
	main()
