import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from os.path import exists

def main():
	# path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/"
	# fileprefix = "1_2_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
	path = []
	fileprefix = []
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls41/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls42/N_2/T_1/")
	path.append("skip")
	path.append("skip")
	path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls40/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls39/N_2/T_1/")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	
	new_data = False
	xdata = []

	# print(data[1])
	# # exit(0)
	for pind,p in enumerate(path):
		if p == "skip":
			continue
		disp_balls = [0,1,2]
		
		sav_slid = p + 'sliding_b3.csv'
		sav_roll = p + 'rolling_b3.csv'

		has_slid = exists(sav_slid)
		has_roll = exists(sav_roll)

		if new_data or (not has_slid or not has_roll):
			datafile = p + fileprefix[pind] + "fricB3Data.csv"
			# distdatafile = p + fileprefix[pind] + "distB3Data.csv"
			print("Reading in data from {}".format(datafile))
			data = np.transpose(np.loadtxt(datafile, delimiter=',',skiprows=1))
			# distdata = np.transpose(np.loadtxt(distdatafile,delimiter=',',skiprows=1))

			rows = data.shape[1]
			cols = data.shape[0]

			print("Starting data analysis")

			rolling_mag = np.zeros((int(cols/6),rows))
			sliding_mag = np.zeros((int(cols/6),rows))
			xdata = np.linspace(0,1,rows) #dummy x data. x is actually time

			for i in range(0,int(cols/6)):
				rolling_mag[i,:] = [np.linalg.norm(data[i*6:i*6+3,j]) for j in range(rows)]
				sliding_mag[i,:] = [np.linalg.norm(data[i*6+3:i*6+6,j]) for j in range(rows)]

			print("Starting data save")

			np.savetxt(sav_slid,sliding_mag,delimiter=',')
			np.savetxt(sav_roll,rolling_mag,delimiter=',')

			
		else:
			print("Getting premade data")
			distdatafile = p + fileprefix[pind] + "distB3Data.csv"
			rolling_mag = np.loadtxt(sav_roll,delimiter=',')
			sliding_mag = np.loadtxt(sav_slid,delimiter=',')
			# distdata = np.transpose(np.loadtxt(distdatafile,delimiter=',',skiprows=1))
			xdata = np.linspace(0,1,rolling_mag.shape[1]) #dummy x data. x is actually time
		
		energydatafile = p + fileprefix[pind] + "energy.csv"
		energy_data = np.loadtxt(energydatafile,delimiter=',',dtype=np.float64,skiprows=1)
		energy_data = energy_data.transpose()
		
		print("Starting plots")

		fig, ax = plt.subplots(len(disp_balls)*3,1,figsize=(5,10))
		# fig, ax = plt.subplots(len(disp_balls)*3,1,figsize=(15,10))
		if pind == 0:
			startind = np.where(xdata>=0)[0][0]
			endind = np.where(xdata>=0.99)[0][0]
			fig.suptitle("Distance w.r.t. ball {} end zoom NO SLIDING FRIC".format(disp_balls[-1]+1))
		elif pind == 1:
			startind = np.where(xdata>=0)[0][0]
			endind = np.where(xdata>=0.99)[0][0]
			# startind = np.where(xdata>=0.047826)[0][0]
			# endind = np.where(xdata>=0.047827)[0][0]
			fig.suptitle("Distance w.r.t. ball {} end zoom NO ROLLING FRIC".format(disp_balls[-1]+1))
		else:
			startind = np.where(xdata>=0.0)[0][0]
			endind = np.where(xdata>=0.9880001)[0][0]
			# startind = np.where(xdata>=0.047826)[0][0]
			# endind = np.where(xdata>=0.047827)[0][0]
			fig.suptitle("Distance w.r.t. ball {} end zoom ROLL AND SLID".format(disp_balls[-1]+1))
		for i in disp_balls:
			# ax[i*3].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
			# ax[i*3].grid(visible=True, which='both', color='0.65', linestyle='-')
			# ax[i*3].minorticks_on()
			# ax[i*3].legend()
			# ax[i*3+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
			# ax[i*3+1].grid(visible=True, which='both', color='0.65', linestyle='-')
			# ax[i*3+1].minorticks_on()
			# ax[i*3+1].legend()
			ax[i*3+2].plot(energy_data[0][startind:endind],energy_data[5][startind:endind],label='tot ang mom'.format(i))
			ax[i*3+2].grid(visible=True, which='both', color='0.65', linestyle='-')
			ax[i*3+2].minorticks_on()
			ax[i*3+2].legend()
			# if i == 1:
				# print(distdata[i][startind:endind])
			# ax[i*3+2].set_ylim((1.45e-5,1.5e-5))
		plt.tight_layout()
		plt.savefig(p + "slidRollenergyB3_allb_zoom_rand.png".format(
i))

	plt.show()

if __name__ == '__main__':
	main()