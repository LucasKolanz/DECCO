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

	data_prefolder = data_directory + 'jobs/longer_job'
	data_prefolder = data_directory + 'jobsOld/tempVarianceRand_attempt'
	data_prefolder = data_directory + 'jobsCosine/lognorm'
	data_prefolder = data_directory + 'jobsNovus/const'
	dataset_name = data_prefolder.split("/")[-1]

	attempts = [i for i in range(5)]
	# attempts = [0]
	attempts300 = attempts


	Nums = [30,100,300]
	# Nums = [30]

	temps = [3,10,30,100,300,1000]
	# temps = [1000]
	

	# numConts = np.full(shape=(len(attempts),len(Nums),len(temps),),fill_value=np.nan)

	newData = True

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
					numConts = []
					angMom = []
					savfilenc = data_directory+f"/data/numContOverTime-{dataset_name}-a_{attempt}-N_{N}-T_{T}.npy"
					savfileam = data_directory+f"/data/angMomOverTime-{dataset_name}-a_{attempt}-N_{N}-T_{T}.npy"
					data_folder = data_prefolder + str(attempt) + '/' + 'N_' + str(N) + '/T_' + str(T) + '/'
					print(data_folder)
					# print(pFD.number_of_contacts(data_folder,29,5000))
					# print(pFD.number_of_contacts(data_folder,29,-1))
					if os.path.exists(data_folder):
						for i in list(range(rows)):
							numConts.append(pFD.number_of_contacts(data_folder,N-3,i))
							angMom.append(pFD.angular_momentum(data_folder,N-3,i))
							# print(pFD.number_of_contacts(data_folder,29,-1))
							# print(pFD.number_of_contacts(data_folder,29,-2))
					else:
						for i in list(range(rows)):
							numConts.append(np.nan)
							angMom.append(np.nan)
						

					#plot num of cont for all writes
					fig,ax = plt.subplots()

					xdata = [i for i in list(range(rows))]

					end = len(xdata)
					ax.plot(xdata[:end],numConts[:end],label='Contacts')
					ax.set_xlabel('Writes (time)')
					ax.set_ylabel('Average Num Contacts')

					ax2 = ax.twinx()
					ax2.plot(xdata[:end],angMom[:end],color='g',label='Ang Mom')
					# # ax.set_legend(['Rabc','RKBM'])
					ax2.set_ylabel('Angular Momentum')
					# plt.errorbar(temps,)
					# ax.set_xscale('log')

					plt.title(f'Num Contacts/ ang mom vs time for {dataset_name} a={attempt}, N={N}, T={T}')
					fig.legend()
					fig.tight_layout() 
					plt.savefig(data_directory+f"data/figures/contacts_and_angmom/{dataset_name}/numContAngMom-{dataset_name}-a_{attempt}-N_{N}-T_{T}.png")
					# plt.show()
	# 	np.save(savfilenc,np.array(numConts))
	# 	np.save(savfileam,np.array(angMom))
	# else:
	# 	numConts = np.load(savfilenc)
	# 	angMom = np.load(savfileam)


if __name__ == '__main__':
	main()