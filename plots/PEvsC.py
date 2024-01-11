import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
	
	# for i,c in enumerate([1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6]):
	# for i,c in enumerate([1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29]):
	fig, ax = plt.subplots(1,2,gridspec_kw = {'wspace':0, 'hspace':0})
	# ax = ax.flatten()
	plots = [1.1,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.3,1.35,1.5,1.6]
	plots = [0.5,0.6,0.7,0.8,0.9,1.1,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22]
	c_vals = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.7,1.8,1.9,2.0]
	PEdata = []
	KEdata = []
	for i,c in enumerate(c_vals):
		path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/cuttoff_test/c_{}/".format(c)
		fileprefix = "1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"

		datafile = path + fileprefix + "energy.csv"

		# header = np.genfromtxt(datafile,delimiter=',',dtype=str)[0]
		data = np.loadtxt(datafile,delimiter=',',dtype=np.float64,skiprows=1)
		print(np.average(data[int(len(data)/2):][2]), c)
		print(np.average(data[int(len(data)/2):][1]), c)
		# exit(0)
		PEdata.append(data[-1][1])
		KEdata.append(data[-1][2])
		# print(data.shape)
	
	# print(c)
	# print(PEdata)
	ax[0].plot(c_vals,PEdata,label='PE vs c')
	ax[0].axvline(np.sqrt(2))
	ax[0].axvline(1/np.sqrt(2))
	ax[0].legend()

	ax[1].plot(c_vals,KEdata,label='KE vs c')
	ax[1].axvline(np.sqrt(2))
	ax[1].axvline(1/np.sqrt(2))
	ax[1].legend()

	plt.subplots_adjust(wspace=0, hspace=0)
	plt.tight_layout()
	plt.show()

if __name__ == '__main__':
	main()