import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from os.path import exists
from timeit import default_timer as timer

def count_rows(file_path):
    line_count = 0
    with open(file_path, 'rb') as file:
        buffer = file.raw.read(1024 * 1024)
        while buffer:
            line_count += buffer.count(b'\n')
            buffer = file.raw.read(1024 * 1024)
    return line_count

def count_cols(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().rstrip()
        comma_count = first_line.count(',')
    return comma_count+1

def main():
	# path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/"
	# fileprefix = "1_2_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
	path = []
	fileprefix = []
	path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls41/N_2/T_1/")
	path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls42/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls40/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls39/N_2/T_1/")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	
	new_data = False
	xdata = []

	# print(data[1])
	# # exit(0)
	for pind,p in enumerate(path):
		disp_balls = [0,1,2]
		
		sav_slid = p + 'sliding_b3.csv'
		sav_roll = p + 'rolling_b3.csv'

		has_slid = exists(sav_slid)
		has_roll = exists(sav_roll)

		if new_data or (not has_slid or not has_roll):
			datafile = p + fileprefix[pind] + "fricB3Data.csv"
			distdatafile = p + fileprefix[pind] + "distB3Data.csv"
			print("Reading in data from {}".format(datafile))
			data = np.transpose(np.loadtxt(datafile, delimiter=',',skiprows=1))
			distdata = np.transpose(np.loadtxt(distdatafile,delimiter=',',skiprows=1))

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
			# rolling_mag = np.loadtxt(sav_roll,delimiter=',')
			# sliding_mag = np.loadtxt(sav_slid,delimiter=',')
			distdata = np.transpose(np.loadtxt(distdatafile,delimiter=',',skiprows=1))
			xdata = np.linspace(0,1,distdata.shape[1]) #dummy x data. x is actually time

		
		print("Starting plots")

		# fig, ax = plt.subplots(len(disp_balls)*3,1,figsize=(15,10))
		# if pind == 0:
		# 	startind = np.where(xdata>=0.016)[0][0]
		# 	endind = np.where(xdata>=0.02)[0][0]
		# 	fig.suptitle("rolling, sliding, and dist w.r.t. ball {} impact 0 differnt sizes".format(disp_balls[-1]+1))
		# else:
		# 	startind = np.where(xdata>=0.0216)[0][0]
		# 	endind = np.where(xdata>=0.02195)[0][0]
		# 	fig.suptitle("rolling, sliding, and dist w.r.t. ball {} impact 0 same size".format(disp_balls[-1]+1))
		# for i in disp_balls:
		# 	ax[i*3].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
		# 	ax[i*3].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3].minorticks_on()
		# 	ax[i*3].legend()
		# 	ax[i*3+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
		# 	ax[i*3+1].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+1].minorticks_on()
		# 	ax[i*3+1].legend()
		# 	ax[i*3+2].plot(xdata[startind:endind],distdata[i][startind:endind],label='ball {} dist'.format(i))
		# 	ax[i*3+2].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+2].minorticks_on()
		# 	ax[i*3+2].legend()
		# plt.tight_layout()
		# plt.savefig(p + "slidRollDistB3_allb_zoom_impact0.png".format(i))

		# fig, ax = plt.subplots(len(disp_balls)*3,1,figsize=(15,10))
		# if pind == 0:
		# 	startind = np.where(xdata>=0.036)[0][0]
		# 	endind = np.where(xdata>=0.038)[0][0]
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} impact 1 different sizes".format(disp_balls[-1]+1))
		# else:
		# 	startind = np.where(xdata>=0.0225)[0][0]
		# 	endind = np.where(xdata>=0.0228)[0][0]
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} impact 1 same size".format(disp_balls[-1]+1))
		# for i in disp_balls:
		# 	ax[i*3].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
		# 	ax[i*3].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3].minorticks_on()
		# 	ax[i*3].legend()
		# 	ax[i*3+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
		# 	ax[i*3+1].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+1].minorticks_on()
		# 	ax[i*3+1].legend()
		# 	ax[i*3+2].plot(xdata[startind:endind],distdata[i][startind:endind],label='ball {} dist'.format(i))
		# 	ax[i*3+2].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+2].minorticks_on()
		# 	ax[i*3+2].legend()
		# plt.tight_layout()
		# plt.savefig(p + "slidRollDistB3_allb_zoom_impact1.png".format(i))

		# fig, ax = plt.subplots(len(disp_balls)*3,1,figsize=(15,10))
		# if pind == 0:
		# 	startind = np.where(xdata>=0.044)[0][0]
		# 	endind = np.where(xdata>=0.046)[0][0]
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} impact 2 different sizes".format(disp_balls[-1]+1))
		# else:
		# 	startind = np.where(xdata>=0.04774)[0][0]
		# 	endind = np.where(xdata>=0.04791)[0][0]
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} impact 2 same size".format(disp_balls[-1]+1))
		# for i in disp_balls:
		# 	ax[i*3].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
		# 	ax[i*3].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3].minorticks_on()
		# 	ax[i*3].legend()
		# 	ax[i*3+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
		# 	ax[i*3+1].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+1].minorticks_on()
		# 	ax[i*3+1].legend()
		# 	ax[i*3+2].plot(xdata[startind:endind],distdata[i][startind:endind],label='ball {} dist'.format(i))
		# 	ax[i*3+2].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+2].minorticks_on()
		# 	ax[i*3+2].legend()
		# plt.tight_layout()
		# plt.savefig(p + "slidRollDistB3_allb_zoom_impact2.png".format(i))


		# fig, ax = plt.subplots(len(disp_balls)*3,1,figsize=(15,10))
		# if pind == 0:
		# 	startind = np.where(xdata>=0.00)[0][0]
		# 	endind = np.where(xdata>=0.07)[0][0]
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} all impacts div 0 protec normal dir".format(disp_balls[-1]+1))
		# else:
		# 	startind = np.where(xdata>=0.00)[0][0]
		# 	endind = np.where(xdata>=0.07)[0][0]
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} all impacts normal run".format(disp_balls[-1]+1))
		# for i in disp_balls:
		# 	ax[i*3].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
		# 	ax[i*3].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3].minorticks_on()
		# 	ax[i*3].legend()
		# 	ax[i*3+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
		# 	ax[i*3+1].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+1].minorticks_on()
		# 	ax[i*3+1].legend()
		# 	ax[i*3+2].plot(xdata[startind:endind],distdata[i][startind:endind],label='ball {} dist'.format(i))
		# 	ax[i*3+2].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+2].minorticks_on()
		# 	ax[i*3+2].legend()
		# plt.tight_layout()
		# plt.savefig(p + "slidRollDistB3_allb_allimpacts.png".format(i))


		# fig, ax = plt.subplots(len(disp_balls)*3,1,figsize=(15,10))
		# if pind == 0:
		# 	startind = np.where(xdata>=0.98802)[0][0]
		# 	endind = np.where(xdata>=0.98808)[0][0]
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} rand zoom SAME SIZE".format(disp_balls[-1]+1))
		# elif pind == 1:
		# 	startind = np.where(xdata>=0.98801)[0][0]
		# 	endind = np.where(xdata>=0.98807)[0][0]
		# 	# startind = np.where(xdata>=0.047826)[0][0]
		# 	# endind = np.where(xdata>=0.047827)[0][0]
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} rand zoom DIFF SIZE".format(disp_balls[-1]+1))
		# else:
		# 	startind = np.where(xdata>=0.98801)[0][0]
		# 	endind = np.where(xdata>=0.98807)[0][0]
		# 	# startind = np.where(xdata>=0.047826)[0][0]
		# 	# endind = np.where(xdata>=0.047827)[0][0]
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} rand zoom DIFF SIZE dt/2".format(disp_balls[-1]+1))
		# for i in disp_balls:
		# 	ax[i*3].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
		# 	ax[i*3].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3].minorticks_on()
		# 	ax[i*3].legend()
		# 	ax[i*3+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
		# 	ax[i*3+1].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+1].minorticks_on()
		# 	ax[i*3+1].legend()
		# 	ax[i*3+2].plot(xdata[startind:endind],distdata[i][startind:endind],label='ball {} dist'.format(i))
		# 	ax[i*3+2].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*3+2].minorticks_on()
		# 	ax[i*3+2].legend()
		# 	# if i == 1:
		# 		# print(distdata[i][startind:endind])
		# 	# ax[i*3+2].set_ylim((1.45e-5,1.5e-5))
		# plt.tight_layout()
		# plt.savefig(p + "slidRollDistB3_allb_zoom_rand.png".format(i))

		fig, ax = plt.subplots(len(disp_balls),1,figsize=(15,10))
		if pind == 0:
			startind = np.where(xdata>=0.92)[0][0]
			endind = np.where(xdata>=0.925)[0][0]
			fig.suptitle("Distance w.r.t. ball {} end zoom NO SLIDING FRIC".format(disp_balls[-1]+1))
		elif pind == 1:
			startind = np.where(xdata>=0.92)[0][0]
			endind = np.where(xdata>=0.9202)[0][0]
			# startind = np.where(xdata>=0.047826)[0][0]
			# endind = np.where(xdata>=0.047827)[0][0]
			fig.suptitle("Distance w.r.t. ball {} end zoom NO ROLLING FRIC".format(disp_balls[-1]+1))
		else:
			startind = np.where(xdata>=0.92)[0][0]
			endind = np.where(xdata>=0.925)[0][0]
			# startind = np.where(xdata>=0.047826)[0][0]
			# endind = np.where(xdata>=0.047827)[0][0]
			fig.suptitle("Distance w.r.t. ball {} end zoom ROLL AND SLID".format(disp_balls[-1]+1))
		for i in disp_balls:
			ax[i].plot(xdata[startind:endind],distdata[i][startind:endind],label='ball {} dist'.format(i))
			ax[i].grid(visible=True, which='both', color='0.65', linestyle='-')
			ax[i].minorticks_on()
			ax[i].legend()
			# if i == 1:
				# print(distdata[i][startind:endind])
			# ax[i*3+2].set_ylim((1.45e-5,1.5e-5))
		plt.tight_layout()
		plt.savefig(p + "slidRollDistB3_allb_zoom_rand.png".format(i))


	plt.show()

if __name__ == '__main__':
	main()