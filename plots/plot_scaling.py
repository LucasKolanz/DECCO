import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab")
import utils as u
import porosity_FD as p

def main():
	base = os.getcwd() + "/jobs/"
	# base = "/global/homes/l/lpkolanz/SpaceLab/jobs/"
	# folder1 = base + "forceTest1/"
	# folder1 = base + "singleCoreComparison/"
	# folder2 = base + "smallerDt3/"
	# folder2 = base + "smallerDt2/"
	# folder2 = base + "singleCoreComparison3/"
	# folder2 = base + "multiCoreTest4/"
	# folder2 = base + "multiCoreTest8/"
	# folder2 = base + "singleCoreComparison_COPY7/"
	# folder1 = "/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab/jobs/accuracyTest11/N_10/T_100/"
	# folder2 = "/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab/jobs/accuracyTest15/N_10/T_100/"


	


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

	# max_ind = -1
	# for file in os.listdir(folder2):
	# 	if file[-4:] == ".csv":
	# 		try:
	# 			ind = int(file.split('_')[0])
	# 			if ind > max_ind:
	# 				max_ind = ind
	# 				body = file.split('_')[1:-1]
	# 		except ValueError:
	# 			continue

	# file2 = str(max_ind)+'_'+'_'.join(body)+"_simData.csv"
	# file2 = "9_simData.csv"

	inds = np.arange(1,20)
	threads = [1,2,4,8,16,32,64]
	# folders = ["strongScaleGrowth1","weakScaleGrowth1"]
	folders = ["strongScaleGrowth1"]#,"weakScaleGrowth1"]
	# inds = np.arange(1,3)

	times = np.zeros((len(folders),len(threads)),dtype=np.float64)
	times[:,:] = np.nan

	
	temp = 100
	for f_i,folder in enumerate(folders):
		# lowest_index = 100
		th = threads
		# if folder == folders[1]:
		# 	th = th[:-1]
		for t_i,t in enumerate(th):
			timeFile = base + folder + "/thread_{}/time.csv".format(t,t)
			print("========================================")
			print(timeFile)
			try:
				with open(timeFile,'r') as tF:
					lines = tF.readlines()
				# if t == th[0]:
				# 	i = -1
				# 	while "ball,update time" in lines[i]:
				# 		i -= 1
				# 	lowest_index = lines[i].split(',')
				# 	print(lowest_index)
				# else:
				# 	continue
				try:
					time = float(lines[-1].split(',')[1][:-1])
				except:
					continue
				# lowest_index = int(lines[-1].split(',')[0])
				times[f_i,t_i] = time
			except FileNotFoundError:
				continue

	speedups = np.copy(times)
	# speedups[0,:] = speedups[0,0]/speedups[0,:]
	# speedups[1,:] = speedups[1,0]/speedups[1,:]

	print(speedups)			


	for f_i,folder in enumerate(folders):
		if folder == folders[0]:
			inds = threads
			title = "Strong Scaling of sim_one_step"
		else:
			inds = threads[:-1]
			title = "Weak Scaling of sim_one_step"
		# print(inds)
		# print(speedups[f_i,:len(inds)])
		fig, ax = plt.subplots(1,1,figsize=(15,7))
		plt.rcParams.update({'font.size': 20})
		ax.tick_params(axis='x',labelsize=20)
		ax.tick_params(axis='y',labelsize=20)
		ax.loglog(inds,speedups[f_i,:len(inds)])
		# ax.plot(inds,,label='multiCoreTest7')
		# ax.set_title(title)
		ax.set_xlabel("Number of Threads",fontsize=20)
		ax.set_ylabel("Time (s)",fontsize=20)
		# ax.legend()
		plt.tight_layout()
		plt.savefig("figures/{}Scaling.png".format(title.split(' ')[0]))


	plt.show()

if __name__ == '__main__':
	main()
