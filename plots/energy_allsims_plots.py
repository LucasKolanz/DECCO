import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
	
	# for i,c in enumerate([1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6]):
	# for i,c in enumerate([1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29]):
	fig, ax = plt.subplots(1,1,gridspec_kw = {'wspace':0, 'hspace':0})
	# ax = ax.flatten()
	# plots = [1.1,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.3,1.35,1.5,1.6]
	# plots = [0.5,0.6,0.7,0.8,0.9,1.1,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22]
	# plots = [1]

	sims = [2]
	title = "sims"
	tot_data = np.array([])
	tot_x_data = np.array([])

	path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/step_back3/N_3/T_1/"
	path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/step_back4/N_3/T_1/"
	
	for j,sim in enumerate(sims):
		if sim == 0:
			fileprefix = "1_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
		else:
			fileprefix = "{}_1_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_".format(sim)

		datafile = path + fileprefix + "energy.csv"

		# header = np.genfromtxt(datafile,delimiter=',',dtype=str)[0]
		data = np.loadtxt(datafile,delimiter=',',dtype=np.float64,skiprows=1)
		data = data.transpose()

		title += " {}".format(sim)
		tot_data.resize(tot_data.shape[0] + data.shape[1])
		tot_x_data.resize(tot_x_data.shape[0] + data.shape[1])
		tot_data[-data.shape[1]:] = data[5] #This is angmom data
		if j == 0:
			tot_x_data[-data.shape[1]:] = data[0] #This is angmom data
		else:
			tot_x_data[-data.shape[1]:] = data[0]+tot_x_data[data.shape[1]*(j)-1] #This is angmom data

	start = 0
	end = -1

	title += ' Ang Mom'
	# xdata = np.linspace(0,tot_data,tot_data.shape[0]) #dummy x data. x is actually time
	ax.plot(tot_x_data[start:end],tot_data[start:end])
	ax.set_title(title)
	ax.set_xlabel('time')
	ax.set_ylabel('Ang Mom')

	contactfile = path + "contact.csv"
	contact = np.loadtxt(contactfile,delimiter=',',dtype=str)
	for i,s in enumerate(sims):
		ax.axvline(float(contact[s][1])+(i)*tot_x_data[data.shape[1]-1],color='g')
	# ax.legend()

	plt.subplots_adjust(wspace=0, hspace=0)
	plt.tight_layout()
	plt.show()

if __name__ == '__main__':
	main()