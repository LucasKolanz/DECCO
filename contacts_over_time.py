import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import os
import csv

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
# import utils as u
import porosity_FD as pFD



def main():
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	data_directory = input_json["data_directory"]

	data_prefolder = data_directory + 'jobs/longer_jobs'
	data_prefolder = data_directory + 'jobsNovus/const'
	data_prefolder = data_directory + 'jobsCosine/lognorm'
	dataset_name = data_prefolder.split("/")[-1]

	attempts = [i for i in range(30)]
	attempts = [1]
	attempts300 = attempts


	Nums = [30,100,300]
	Nums = [30]

	temps = [3,10,30,100,300,1000]
	temps = [1000]
	

	# numConts = np.full(shape=(len(attempts),len(Nums),len(temps),),fill_value=np.nan)
	numConts = []

	newData = True
	savfile = data_directory+"/data/numContOverTime-a_1-N_30-T_1000.npy"

	rows = 51

	if newData:
		for T_i,T in enumerate(temps):
			for n_i,N in enumerate(Nums):
				if N == 300:
					a = attempts300
				else:
					a = attempts
				a = attempts
				for a_i,attempt in enumerate(a):
					data_folder = data_prefolder + str(attempt) + '/' + 'N_' + str(N) + '/T_' + str(T) + '/'
					# print(pFD.number_of_contacts(data_folder,29,5000))
					# print(pFD.number_of_contacts(data_folder,29,-1))
					for i in list(range(rows)):
						numConts.append(pFD.number_of_contacts(data_folder,29,i))
						# print(pFD.number_of_contacts(data_folder,29,-1))
						# print(pFD.number_of_contacts(data_folder,29,-2))

		np.save(savfile,np.array(numConts))
	else:
		numConts = np.load(savfile)
	#plot num of cont for all writes
	fig,ax = plt.subplots()

	# print(list(range(rows)))
	xdata = [i for i in list(range(rows))]
	# print(len(xdata))
	# print(numConts[:].shape)

	print(np.sum(np.where(numConts == numConts[-1],1,0)))

	# end = 500
	ax.plot(xdata[:end],numConts[:end])
	ax.set_title('Num Contacts vs time for N=30, T=3')
	ax.set_xlabel('Writes (time)')
	ax.set_ylabel('Average Num Contacts')
	# # ax.set_legend(['Rabc','RKBM'])
	# plt.errorbar(temps,)
	# ax.set_xscale('log')


	# fig.legend()
	plt.show()


if __name__ == '__main__':
	main()