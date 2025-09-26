import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
	
	# for i,c in enumerate([1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6]):
	# for i,c in enumerate([1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29]):
	fig, ax = plt.subplots(4,5,gridspec_kw = {'wspace':0, 'hspace':0})
	ax = ax.flatten()
	plots = [1.1,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.3,1.35,1.5,1.6]
	plots = [0.5,0.6,0.7,0.8,0.9,1.1,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22]
	plots = [0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,5.0]
	plots = [1.0]
	for i,c in enumerate(plots):
		path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/cuttoff_test/c_{}/".format(c)
		path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/step_back1/N_2/T_1/"
		fileprefix = "1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"

		datafile = path + fileprefix + "energy.csv"

		# header = np.genfromtxt(datafile,delimiter=',',dtype=str)[0]
		data = np.loadtxt(datafile,delimiter=',',dtype=np.float64,skiprows=1)
		data = data.transpose()
		# print(data.shape)
		xdata = np.linspace(0,1,data.shape[1]) #dummy x data. x is actually time

		# fig, ax = plt.subplots(1,1,figsize=(15,7))

		start = 0
		end = -1
		# print(xdata[start:end].shape)

		# ax[i].set_title('ang mom, cuttoff = {}'.format(c))
		ax[i].plot(xdata[start:end],data[5][start:end],label="c={}".format(c))
		ax[i].legend()

	plt.subplots_adjust(wspace=0, hspace=0)
	plt.tight_layout()
	plt.show()

if __name__ == '__main__':
	main()