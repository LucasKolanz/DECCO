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
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls21/N_2/T_1_copy/")
	path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls29/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls21/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls15/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls16/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# fileprefix.append("2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("third_1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("2_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("9_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")

	# disp_balls = [0]
	new_data = False
	xdata = []

	# print(data[1])
	# # exit(0)
	for pind,p in enumerate(path):
		disp_balls = [0,1]
		if pind == 1:
			disp_balls.append(2)
		
		sav_slid = p + 'sliding_b3.csv'
		sav_roll = p + 'rolling_b3.csv'

		has_slid = exists(sav_slid)
		has_roll = exists(sav_roll)

		if new_data or (not has_slid or not has_roll):
			datafile = p + fileprefix[pind] + "fricB3Data.csv"
			print("Reading in data from {}".format(datafile))
			data = np.transpose(np.loadtxt(datafile, delimiter=',',skiprows=1))
			# data = data[:,5000000:]
			# print(data.shape)
			rows = data.shape[1]
			cols = data.shape[0]

			print("Starting data analysis")

			rolling_mag = np.zeros((int(cols/6),rows))
			sliding_mag = np.zeros((int(cols/6),rows))
			xdata = np.linspace(0,1,rows) #dummy x data. x is actually time
			# print(data)
			# print(data.shape)
			for i in range(0,int(cols/6)):
				rolling_mag[i,:] = [np.linalg.norm(data[i*6:i*6+3,j]) for j in range(rows)]
				# print(rolling_mag[i,:])
				sliding_mag[i,:] = [np.linalg.norm(data[i*6+3:i*6+6,j]) for j in range(rows)]

			# print(rolling_mag)

			# exit(0)

			print("Starting data save")
			# dat = [sliding_mag,rolling_mag]
			# for i,sav in enumerate([sav_slid,sav_roll]):
			np.savetxt(sav_slid,sliding_mag,delimiter=',')
			np.savetxt(sav_roll,rolling_mag,delimiter=',')
		else:
			print("Getting premade data")
			rolling_mag = np.loadtxt(sav_roll,delimiter=',')
			sliding_mag = np.loadtxt(sav_slid,delimiter=',')
			xdata = np.linspace(0,1,rolling_mag.shape[1]) #dummy x data. x is actually time

		
		# ax = plt.axes(projection='3d')
		print("Starting plots")
		# fig, ax = plt.subplots(3,1,figsize=(15,10))
		# startind = np.where(xdata>=0.0475)[0][0]
		# endind = np.where(xdata>=0.0476)[0][0]
		# for i in disp_balls:
		# 	ax[i].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {}'.format(i))
		# 	# ax[i].plot(xdata,rolling_mag[i],label='ball {}'.format(i))
		# 	ax[i].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i].minorticks_on()
		# 	ax[i].legend()
		# fig.suptitle('rolling on B3')
		# plt.tight_layout()
		# plt.savefig(p + "RollFricB3_allb.png".format(i))

		# fig, ax = plt.subplots(3,1,figsize=(15,10))
		# # startind = np.where(xdata>=0.05)[0][0]
		# for i in disp_balls:
		# 	# ax[i].plot(xdata,sliding_mag[i],label='ball {}'.format(i))
		# 	# ax[i].plot(xdata[startind:],sliding_mag[i][startind:],label='ball {}'.format(i))
		# 	ax[i].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {}'.format(i))
		# 	ax[i].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i].minorticks_on()
		# 	ax[i].legend()
		# 	# ax[i].set_xlim((0.5,1))
		# fig.suptitle('sliding on B3')
		# plt.tight_layout()
		# plt.savefig(p + "slidFricB3_allb.png".format(i))

		# fig, ax = plt.subplots(len(disp_balls)*2,1,figsize=(15,10))
		# if pind == 1:
		# 	startind = np.where(xdata>=0.01708)[0][0]
		# 	endind = np.where(xdata>=0.01725)[0][0]
		# else:
		# 	startind = np.where(xdata>=0.0184)[0][0]
		# 	endind = np.where(xdata>=0.0186)[0][0]
		# for i in disp_balls:
		# 	# ax[i*2].plot(xdata,rolling_mag[i],label='ball {} roll'.format(i))
		# 	# ax[i*2].plot(xdata[startind:],rolling_mag[i][startind:],label='ball {} roll'.format(i))
		# 	ax[i*2].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
		# 	ax[i*2].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*2].minorticks_on()
		# 	ax[i*2].legend()
		# 	ax[i*2+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
		# 	# ax[i*2+1].plot(xdata,sliding_mag[i],label='ball {} slide'.format(i))
		# 	# ax[i*2+1].plot(xdata[startind:],sliding_mag[i][startind:],label='ball {} slide'.format(i))
		# 	ax[i*2+1].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*2+1].minorticks_on()
		# 	ax[i*2+1].legend()
		# 	# ax[i].set_xlim((0.5,1))
		# fig.suptitle("rolling and sliding w.r.t. ball {} impact 0".format(disp_balls[-1]+1))
		# plt.tight_layout()
		# plt.savefig(p + "slidRollFricB3_allb_zoom_impact0.png".format(i))

		# fig, ax = plt.subplots(len(disp_balls)*2,1,figsize=(15,10))
		# if pind == 1:
		# 	startind = np.where(xdata>=0.0439)[0][0]
		# 	endind = np.where(xdata>=0.04405)[0][0]
		# else:
		# 	startind = np.where(xdata>=0.0333)[0][0]
		# 	endind = np.where(xdata>=0.03355)[0][0]
		# for i in disp_balls:
		# 	# ax[i*2].plot(xdata,rolling_mag[i],label='ball {} roll'.format(i))
		# 	# ax[i*2].plot(xdata[startind:],rolling_mag[i][startind:],label='ball {} roll'.format(i))
		# 	ax[i*2].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
		# 	ax[i*2].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*2].minorticks_on()
		# 	ax[i*2].legend()
		# 	ax[i*2+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
		# 	# ax[i*2+1].plot(xdata,sliding_mag[i],label='ball {} slide'.format(i))
		# 	# ax[i*2+1].plot(xdata[startind:],sliding_mag[i][startind:],label='ball {} slide'.format(i))
		# 	ax[i*2+1].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 	ax[i*2+1].minorticks_on()
		# 	ax[i*2+1].legend()
		# 	# ax[i].set_xlim((0.5,1))
		# fig.suptitle("rolling and sliding w.r.t. ball {} impact 1".format(disp_balls[-1]+1))
		# plt.tight_layout()
		# plt.savefig(p + "slidRollFricB3_allb_zoom_impact1.png".format(i))

		# if len(disp_balls) > 2:
		# 	fig, ax = plt.subplots(len(disp_balls)*2,1,figsize=(15,10))
		# 	startind = np.where(xdata>=0.04715)[0][0]
		# 	endind = np.where(xdata>=0.0473)[0][0]
		# 	for i in disp_balls:
		# 		# ax[i*2].plot(xdata,rolling_mag[i],label='ball {} roll'.format(i))
		# 		# ax[i*2].plot(xdata[startind:],rolling_mag[i][startind:],label='ball {} roll'.format(i))
		# 		ax[i*2].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
		# 		ax[i*2].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 		ax[i*2].minorticks_on()
		# 		ax[i*2].legend()
		# 		ax[i*2+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
		# 		# ax[i*2+1].plot(xdata,sliding_mag[i],label='ball {} slide'.format(i))
		# 		# ax[i*2+1].plot(xdata[startind:],sliding_mag[i][startind:],label='ball {} slide'.format(i))
		# 		ax[i*2+1].grid(visible=True, which='both', color='0.65', linestyle='-')
		# 		ax[i*2+1].minorticks_on()
		# 		ax[i*2+1].legend()
		# 		# ax[i].set_xlim((0.5,1))
		# 	fig.suptitle("rolling and sliding w.r.t. ball {} impact 2".format(disp_balls[-1]+1))
		# 	plt.tight_layout()
		# 	plt.savefig(p + "slidRollFricB3_allb_zoom_impact2.png".format(i))

		fig, ax = plt.subplots(len(disp_balls)*2,1,figsize=(15,10))
		if pind == 1:
			startind = np.where(xdata>=0)[0][0]
			endind = np.where(xdata>=0.05)[0][0]
		else:
			startind = np.where(xdata>=0)[0][0]
			endind = np.where(xdata>=0.05)[0][0]
		for i in disp_balls:
			# ax[i*2].plot(xdata,rolling_mag[i],label='ball {} roll'.format(i))
			# ax[i*2].plot(xdata[startind:],rolling_mag[i][startind:],label='ball {} roll'.format(i))
			ax[i*2].plot(xdata[startind:endind],rolling_mag[i][startind:endind],label='ball {} roll'.format(i))
			# ax[i*2].plot(xdata,rolling_mag[i],label='ball {} roll'.format(i))
			ax[i*2].grid(visible=True, which='both', color='0.65', linestyle='-')
			ax[i*2].minorticks_on()
			ax[i*2].legend()
			ax[i*2+1].plot(xdata[startind:endind],sliding_mag[i][startind:endind],label='ball {} slide'.format(i))
			# ax[i*2+1].plot(xdata,sliding_mag[i],label='ball {} slide'.format(i))
			# ax[i*2+1].plot(xdata[startind:],sliding_mag[i][startind:],label='ball {} slide'.format(i))
			ax[i*2+1].grid(visible=True, which='both', color='0.65', linestyle='-')
			ax[i*2+1].minorticks_on()
			ax[i*2+1].legend()
			# ax[i].set_xlim((0.5,1))
		fig.suptitle("rolling and sliding w.r.t. ball {} all impacts".format(disp_balls[-1]+1))
		plt.tight_layout()
		plt.savefig(p + "slidRollFricB3_allb_allimpacts.png".format(i))

		# for i in disp_balls:
			# fig, ax = plt.subplots(2,1,figsize=(15,10))
			# ax = ax.flatten()
			# ax[0].plot(xdata, rolling_mag[i], label='ball {}'.format(i))
			# ax[0].set_title('rolling friction magnitude')
			# ax[0].set_xlabel('time (t/tmax)')
			# ax[0].set_ylabel('friction')
			# ax[0].grid(visible=True, which='both', color='0.65', linestyle='-')
			# ax[0].minorticks_on()
			# ax[0].legend()

			# ax[1].plot(xdata, sliding_mag[i], label='ball {}'.format(i))
			# ax[1].set_title('sliding friction magnitude')
			# ax[1].set_xlabel('time (t/tmax)')
			# ax[1].set_ylabel('friction')
			# ax[1].legend()
			# ax[1].grid(visible=True, which='both', color='0.65', linestyle='-')
			# ax[1].minorticks_on()
			# fig.suptitle(p[-9:])
			# plt.tight_layout()
			# plt.savefig(p + "slidRollFric_b-{}.png".format(i))

				
			# fig, ax = plt.subplots(2,1,figsize=(15,10))
			# ax = ax.flatten()
			# ax[0].plot(xdata, rolling_mag[i], label='ball {}'.format(i))
			# ax[0].set_title('rolling friction magnitude')
			# ax[0].set_xlabel('time (t/tmax)')
			# ax[0].set_ylabel('friction')
			# ax[0].grid(visible=True, which='both', color='0.65', linestyle='-')
			# ax[0].minorticks_on()
			# ax[0].legend()

			# ax[1].plot(xdata, sliding_mag[i], label='ball {}'.format(i))
			# ax[1].set_title('sliding friction magnitude')
			# ax[1].set_xlabel('time (t/tmax)')
			# ax[1].set_ylabel('friction')
			# if i == 3:
			# 	ax[1].set_ylim((0.00014,0.000165))
			# elif i == 1:
			# 	ax[1].set_ylim ((0.00016,0.00024))
			# ax[1].legend()		
			# ax[1].grid(visible=True, which='both', color='0.65', linestyle='-')
			# ax[1].minorticks_on()
			# fig.suptitle(p[-9:])
			# plt.tight_layout()
			# # plt.savefig(p + "slidRollFric_oscil_b-{}.png".format(i))

			# fig, ax = plt.subplots(2,1,figsize=(15,10))
			# ax = ax.flatten()
			# ax[0].plot(xdata, rolling_mag[i], label='ball {}'.format(i))
			# ax[0].set_title('rolling friction magnitude')
			# ax[0].set_xlabel('time (t/tmax)')
			# ax[0].set_ylabel('friction')
			# ax[0].set_xlim((0,0.05))
			# ax[0].grid(visible=True, which='both', color='0.65', linestyle='-')
			# ax[0].minorticks_on()
			# ax[0].legend()

			# ax[1].plot(xdata, sliding_mag[i], label='ball {}'.format(i))
			# ax[1].set_title('sliding friction magnitude')
			# ax[1].set_xlabel('time (t/tmax)')
			# ax[1].set_ylabel('friction')
			# ax[1].set_xlim((0,0.05))
			# ax[1].legend()		
			# ax[1].grid(visible=True, which='both', color='0.65', linestyle='-')
			# ax[1].minorticks_on()
			# fig.suptitle(p[-9:])
			# plt.tight_layout()
			# plt.savefig(p + "slidRollFric_impact_b-{}.png".format(i))

	plt.show()

if __name__ == '__main__':
	main()