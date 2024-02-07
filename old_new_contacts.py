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
	data_prefolder = data_directory + 'jobsNovus/const'
	data_prefolder = data_directory + 'jobsCosine/lognorm'
	data_prefolder = data_directory + 'jobsCosine/lognorm_relax'
	data_prefolder = data_directory + 'jobsNovus/const_relax'
	dataset_name = data_prefolder.split("/")[-1]

	attempts = [i for i in range(30)]
	# attempts = [0]
	attempts300 = attempts


	Nums = [30,100,300]
	Nums = [30]

	temps = [3,10,30,100,300,1000]
	# temps = [1000]
	

	numConts = np.full(shape=(2,len(Nums),len(temps),len(attempts)),fill_value=np.nan)
	# angMom = np.full(shape=(len(attempts),len(Nums),len(temps),2),fill_value=np.nan)

	newData = True

	rows = 51

	for T_i,T in enumerate(temps):
		for n_i,N in enumerate(Nums):
			if N == 300:
				a = attempts300
			else:
				a = attempts
			a = attempts

			for a_i,attempt in enumerate(a):
				# savfilenc = data_directory+f"/data/numContOverTime-{dataset_name}-a_{attempt}-N_{N}-T_{T}.npy"
				# savfileam = data_directory+f"/data/angMomOverTime-{dataset_name}-a_{attempt}-N_{N}-T_{T}.npy"
				data_folder = data_prefolder + str(attempt) + '/' + 'N_' + str(N) + '/T_' + str(T) + '/'
				print(data_folder)
				# print(pFD.number_of_contacts(data_folder,29,5000))
				# print(pFD.number_of_contacts(data_folder,29,-1))
				if os.path.exists(data_folder):
					t1 = []
					# t2 = []
					for i in list(range(rows)):
						t1.append(pFD.number_of_contacts(data_folder,data_index=N-3,line=i,relax=True))
						# t2.append()
					numConts[0,n_i,T_i,a_i] = np.max(t1) # new relaxed data
					numConts[1,n_i,T_i,a_i] = pFD.number_of_contacts(data_folder,data_index=N-3,line=-1,relax=False) # old data

						# print(pFD.number_of_contacts(data_folder,29,-1))
						# print(pFD.number_of_contacts(data_folder,29,-2))
					#plot num of cont for all writes


	fig,ax = plt.subplots()
	print(np.sum(np.where(np.isnan(numConts),1,0)))
	ax.plot(temps,np.nanmean(numConts[0,0,:,:],axis=1),label='Contacts relaxed')
	ax.set_xlabel('Writes (time)')
	ax.set_ylabel('Average Num Contacts')

	# ax2 = ax.twinx()
	ax.plot(temps,np.mean(numConts[1,0,:,:],axis=1),label='Contacts old')
	# # ax.set_legend(['Rabc','RKBM'])
	# ax2.set_ylabel('Angular Momentum')
	# plt.errorbar(temps,)
	# ax.set_xscale('log')
	ax.set_xscale('log')


	plt.title(f'Num Contacts old vs new')
	# Adjusting legend position
	fig.legend(loc='lower left')#, bbox_to_anchor=(1, 1))

	# Adjust layout to make space for the legend
	plt.subplots_adjust(right=0.75)  # Adjust this value as needed to fit your legend

	fig.tight_layout() 
	plt.savefig(data_directory+f"data/figures/contacts_and_angmom/{dataset_name}/oldVsNewNumConts.png")
	plt.show()
	# plt.close(fig)


if __name__ == '__main__':
	main()