import os
import numpy as np
import matplotlib.pyplot as plt
import sys


relative_path = "."
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
import utils as u
import gen_data as gd


def main():
	base = os.getcwd() + "/jobs/"
	# folder1 = base + "multiCoreTest4/"
	# folder1 = base + "singleCoreComparison/"
	# folder1 = base + "multiCoreTest1/"
	# folder2 = base + "singleCoreComparison2/"
	# folder2 = base + "multiCoreTest7/"
	# folder2 = base + "singleCoreComparison_COPY7/"
	folder1 = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPAWELDTEST_10/M_3/N_300/T_1000/"
	folder2 = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPAWELDTEST_11/M_3/N_300/T_1000/"



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

	

	N = 6
	temp = 100
	show_FD_plots = False
	for ind in [6]:
		f1 = "{}_{}_simData.csv".format(ind,'_'.join(body))
		f2 = "{}_{}_simData.csv".format(ind,'_'.join(body))
		# f2 = "{}_simData.csv".format(ind)
		porositiesabc=[]
		porositiesKBM=[]
		contacts=[]
		FD_data=[]
		for data_folder in [folder1,folder2]:
			print("========================================")
			print(data_folder)
			count = 0
			# for root_dir, cur_dir, files in os.walk(data_folder):
			#     count += len(files)
			# if count/3 > N:
			porositiesabc.append(gd.calc_porosity_abc(data_folder,ind))
			porositiesKBM.append(gd.calc_porosity_KBM(data_folder,ind))
			# contacts.append(p.number_of_contacts(data_folder,ind))
			# if not np.isnan(porositiesabc[-1]):
			# 	o3dv = u.o3doctree(data_folder,overwrite_data=True,index=ind,Temp=temp)
			# 	o3dv.make_tree()
			# 	FD_data.append(o3dv.calc_fractal_dimension(show_graph=show_FD_plots))
		print("================Comparing {}================".format(ind))
		print("Porosity abc difference: {}".format(np.diff(porositiesabc)))	
		print("Porosity KBM difference: {}".format(np.diff(porositiesKBM)))	
		# print("contacts difference    : {}".format(np.diff(contacts)))	
		# print("Fract dim difference   : {}".format(np.diff(FD_data)))	
		print("====================END====================".format(ind))

		# data1 = np.loadtxt(folder1+f1,delimiter=",",dtype=np.float64,skiprows=1)[-1]
		# data2 = np.loadtxt(folder2+f2,delimiter=",",dtype=np.float64,skiprows=1)[-1]
		# fig, ax = plt.subplots(1,1,figsize=(15,7))
		# ax.plot(np.arange(len(data1)),np.subtract(data2,data1))
		# plt.show()

if __name__ == '__main__':
	main()
