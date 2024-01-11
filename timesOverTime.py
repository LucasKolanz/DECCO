import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/lpkolanz/Desktop/SpaceLab_branch/SpaceLab")
import utils as u
import porosity_FD as p

def main():
	base = os.getcwd() + "/jobs/"
	# folder1 = base + "forceTest1/"
	folders = []
	folders.append(base + "singleCoreComparison/")
	folders.append(base + "piplineImprovements1/")
	folders.append(base + "OpenMpParallel1/")
	folders.append(base + "pipeAndOpenmp1/")
	# folders.append(base + "openMPallLoops1/")
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

	max_ind = -1
	for file in os.listdir(folders[0]):
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

	inds = np.arange(1,20)
	# inds = np.arange(1,3)

	times = np.zeros((len(folders),inds.size),dtype=np.float64)
	times[:,:] = np.nan
	temp = 100
	show_FD_plots = False
	# for ind_i,ind in enumerate(inds):
	f1 = "sim_errors.txt"
	f2 = "sim_err.log"

	# f2 = "{}_simData.csv".format(ind)
	for fold_i,data_folder in enumerate(folders):
		print("========================================")
		if os.path.isfile(data_folder+f1):
			f = f1 
		elif os.path.isfile(data_folder+f2):
			f = f2 
		else:
			print("ERROR: sim error file not found in folder: {}".format(data_folder))
			exit(0)
		
		with open(data_folder+f) as fp:
			count = 0
			while 1:
				line = fp.readline()
				if count >= inds.size or not line:
					break
				if ("Simulation complete!" in line):
					fp.readline() #read line with ball info
					fp.readline() #read line with simulation time info
					timeLine = fp.readline() #read line with run time info
					time = float(timeLine.split(' ')[-2])
					times[fold_i,count] = time
					count += 1


	# realname1=""
	# if folder1.split('/')[-2] == "multiCoreTest7":
	# 	realname1 = "backwardsLoop"
	# elif folder1.split('/')[-2] == "multiCoreTest8" or folder1.split('/')[-2] == "multiCoreTest9":
	# 	realname1 = "parallel"
	# else:
	# 	realname1 = folder1.split('/')[-2]

	# realname2=""
	# if folder2.split('/')[-2] == "multiCoreTest7":
	# 	realname2 = "backwardsLoop"
	# elif folder2.split('/')[-2] == "multiCoreTest8" or folder2.split('/')[-2] == "multiCoreTest9":
	# 	realname2 = "parallel"
	# else:
	# 	realname2 = folder2.split('/')[-2]

	fig, ax = plt.subplots(1,1,figsize=(15,7))
	plt.rcParams.update({'font.size': 20})
	# plt.rcParams.update({'xtick.labelsize': 20})
	# plt.rcParams.update({'ytick.labelsize': 20})
	# plt.rcParams.update({'axes.labelsize': 20})
	for f,folder in enumerate(folders):
		realname=""
		if folder.split('/')[-2] == "singleCoreComparison":
			realname = "Single core"
			mark = '-.'
		elif folder.split('/')[-2] == "OpenMpParallel1":# or folder.split('/')[-2] == "multiCoreTest9":
			realname = "OpenMP for"
			mark = '-'
		elif folder.split('/')[-2] == "piplineImprovements1":# or folder.split('/')[-2] == "multiCoreTest9":
			realname = "Pipeline improvements"
			mark = '--'
		elif folder.split('/')[-2] == "pipeAndOpenmp1":# or folder.split('/')[-2] == "multiCoreTest9":
			realname = "OpenMP for and pipeline improvements"
			mark = ':'
		else:
			realname = folder.split('/')[-2]
		ax.plot(inds,times[f],label=realname,linestyle=mark)

	ax.tick_params(axis='x',labelsize=20)
	ax.tick_params(axis='y',labelsize=20)
	# ax.set_yticklabels(fontsize=20)

	# plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
	# plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
	# plt.rc('ytick', labelsize=20)
	title = "Execution time vs ball number"
	# ax.set_title(title)
	ax.set_xlabel("Balls added",fontsize=20)
	ax.set_ylabel("Time (s)",fontsize=20)
	ax.legend()
	plt.savefig("figures/exTimeVsBall.png")

	# for m,method in enumerate([porositiesabc,porositiesKBM,contacts,FD_data]):
	# 	fig, ax = plt.subplots(1,1,figsize=(15,7))
	# 	ax.plot(inds,method[0]-method[1],label='single - {}'.format(realname))
	# 	# ax.plot(inds,,label='multiCoreTest7')
	# 	title = ""
	# 	if m == 0:
	# 		title = "Porositiesabc_difference"
	# 	elif m == 1:
	# 		title = "PorositiesKBM_difference"
	# 	elif m == 2:
	# 		title = "Contacts_difference"
	# 	elif m == 3:
	# 		title = "FD_difference"
	# 	ax.set_title(title)
	# 	ax.set_xlabel("ball")
	# 	ax.set_ylabel("Value")
	# 	ax.legend()
	# 	plt.savefig("figures/{}".format(realname)+title+"OverTime.png")


	plt.show()

if __name__ == '__main__':
	main()
