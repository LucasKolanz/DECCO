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
	folder1 = base + "singleCoreComparison/"
	folder2 = base + "smallerDt3/"
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

	inds = np.arange(1,20)
	# inds = np.arange(1,3)

	porositiesabc=np.zeros((2,inds.size))
	porositiesKBM=np.zeros((2,inds.size))
	contacts=np.zeros((2,inds.size))
	FD_data=np.zeros((2,inds.size))

	temp = 100
	show_FD_plots = False
	for ind_i,ind in enumerate(inds):
		f1 = "{}_{}_simData.csv".format(ind,'_'.join(body))
		f2 = "{}_{}_simData.csv".format(ind,'_'.join(body))
		# f2 = "{}_simData.csv".format(ind)
		for fold_i,data_folder in enumerate([folder1,folder2]):
			print("========================================")
			print(data_folder)
			count = 0
			# for root_dir, cur_dir, files in os.walk(data_folder):
			#     count += len(files)
			# if count/3 > N:
			porositiesabc[fold_i,ind_i] = p.porosity_measure1(data_folder,ind)
			porositiesKBM[fold_i,ind_i] = p.porosity_measure2(data_folder,ind)
			contacts[fold_i,ind_i] = p.number_of_contacts(data_folder,ind)
			if not np.isnan(porositiesabc[fold_i,ind_i]):
				o3dv = u.o3doctree(data_folder,overwrite_data=True,index=ind,Temp=temp)
				o3dv.make_tree()
				FD_data[fold_i,ind_i] = o3dv.calc_fractal_dimension(show_graph=show_FD_plots)


	# fig, ax = plt.subplots(1,1,figsize=(15,7))
	# ax.plot(inds,porositiesabc[0],label='por abc')
	# ax.plot(inds,porositiesKBM[0],label='por KBM')
	# ax.plot(inds,contacts[0],label='contacts')
	# ax.plot(inds,FD_data[0],label='FD')
	# ax.set_title("folder: {}".format(folder1))
	# ax.set_xlabel("ball")
	# ax.set_ylabel("Value")
	# ax.legend()

	# fig1, ax1 = plt.subplots(1,1,figsize=(15,7))
	# ax1.plot(inds,porositiesabc[1],label='por abc')
	# ax1.plot(inds,porositiesKBM[1],label='por KBM')
	# ax1.plot(inds,contacts[1],label='contacts')
	# ax1.plot(inds,FD_data[1],label='FD')
	# ax1.set_title("folder: {}".format(folder2))
	# ax1.set_xlabel("ball")
	# ax1.set_ylabel("Value")
	# ax1.legend()

	# print(porositiesabc)
	# print(np.subtract(porositiesabc[0],porositiesabc[1]))
	# fig2, ax2 = plt.subplots(1,1,figsize=(15,7))
	# ax2.plot(inds,np.subtract(porositiesabc[0],porositiesabc[1])/porositiesabc[0],label='por abc')
	# ax2.plot(inds,np.subtract(porositiesKBM[0],porositiesKBM[1])/porositiesKBM[0],label='por KBM')
	# ax2.plot(inds,np.subtract(contacts[0],contacts[1])/contacts[0],label='contacts')
	# ax2.plot(inds,np.subtract(FD_data[0],FD_data[1])/FD_data[0],label='FD')
	# ax2.set_title("folder1: {}\nfolder2: {}".format(folder1,folder2))
	# ax2.set_xlabel("ball")
	# ax2.set_ylabel("Value")
	# ax2.legend()

	realname1=""
	if folder1.split('/')[-2] == "multiCoreTest7":
		realname1 = "backwardsLoop"
	elif folder1.split('/')[-2] == "multiCoreTest8" or folder1.split('/')[-2] == "multiCoreTest9":
		realname1 = "parallel"
	else:
		realname1 = folder1.split('/')[-2]

	realname2=""
	if folder2.split('/')[-2] == "multiCoreTest7":
		realname2 = "backwardsLoop"
	elif folder2.split('/')[-2] == "multiCoreTest8" or folder2.split('/')[-2] == "multiCoreTest9":
		realname2 = "parallel"
	else:
		realname2 = folder2.split('/')[-2]

	for m,method in enumerate([porositiesabc,porositiesKBM,contacts,FD_data]):
		fig, ax = plt.subplots(1,1,figsize=(15,7))
		ax.plot(inds,method[0],label=realname1)
		ax.plot(inds,method[1],label=realname2)
		title = ""
		if m == 0:
			title = "Porositiesabc"
		elif m == 1:
			title = "PorositiesKBM"
		elif m == 2:
			title = "Contacts"
		elif m == 3:
			title = "FD"
		ax.set_title(title)
		ax.set_xlabel("ball")
		ax.set_ylabel("Value")
		ax.legend()
		plt.savefig("figures/{}-{}".format(realname1,realname2)+title+"OverTime.png")

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
