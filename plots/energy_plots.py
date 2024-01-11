import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
	# path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/"
	# fileprefix = "1_2_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
	path = []
	fileprefix = []
	#These show progression of problem from a single simulation
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# fileprefix.append("2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("2_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# filepreot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))

	#These show 1) sliding and rolling friction 2) just sliding 3) just rolling
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls22/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls16/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls40/N_2/T_1/")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt3e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	
	h_val = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
	# for i in [1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6]:
	for i in h_val:
		# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/cuttoff_test/c_{}/".format(i))
		path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/hminTests1/N_2/hmin_{}/".format(i))
		fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")

	# print(data)
	# print(data[1])
	# # exit(0)
	fig, ax = plt.subplots(2,4,figsize=(15,7))

	ax = ax.flatten()
	for i,p in enumerate(path):

		datafile = p + fileprefix[i] + "energy.csv"
		header = np.loadtxt(datafile,delimiter=',',dtype=str)[0]
		data = np.loadtxt(datafile,delimiter=',',dtype=np.float64,skiprows=1)
		data = data.transpose()
		time = data[0]
		PE = data[1]
		KE = data[2]
		.
		.
		.

		xdata = np.linspace(0,1,data.shape[1]) #dummy x data. x is actually time


		# ax[0,0].plot(data[0],data[1],label="PE")
		# ax[0,1].plot(data[0],data[2],label="KE")
		# ax[0,2].plot(data[0],data[3],label="PE+KE")

		# ax[1,0].plot(data[0],data[1],label="PE")
		# ax[1,0].plot(data[0],data[2],label="KE")
		# ax[1,0].plot(data[0],data[3],label="PE+KE")
		
		# ax[1,1].plot(data[0],data[4],label="mom")
		start = 0
		end = -1
		print(xdata[start:end].shape)

		ax[i].plot(xdata[start:end],data[5][start:end],label="ang mom hmin={}".format(h_val[i]))
		# ax.set_title('')
		# ax[1,2].plot(data[0][start:end],data[5][start:end],label="ang mom")

		# ax[0,0].legend()
		# ax[0,1].legend()
		# ax[0,2].legend()
		# ax[1,0].legend()
		# ax[1,1].legend()
		# ax[1,2].legend()
		ax[i].legend()

		plt.tight_layout()
	plt.show()

if __name__ == '__main__':
	main()