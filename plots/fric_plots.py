import numpy as np
import matplotlib.pyplot as plt
import sys
from os.path import exists

def main():
	# path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/"
	# fileprefix = "1_2_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
	path = []
	fileprefix = []
	path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls20/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls14/N_2/T_1_copy/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls15/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls16/N_2/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# fileprefix.append("2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	fileprefix.append("third_1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("quarter_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("2_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("9_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")

	disp_balls = [0,1,2,3]
	# disp_balls = [0]
	new_data = False
	xdata = []

	# print(data[1])
	# # exit(0)
	for pind,p in enumerate(path):
		for b in disp_balls:
			sav_slid = p + 'sliding_b-{}.csv'.format(b)
			sav_roll = p + 'rolling_b-{}.csv'.format(b)

			has_slid = exists(sav_slid)
			has_roll = exists(sav_roll)

			if new_data or (not has_slid or not has_roll):
				datafile = p + fileprefix[pind] + "fricData.csv"
				print("Getting data from: {}{}".format(p,fileprefix[pind]))
				header = np.genfromtxt(datafile,delimiter=',',dtype=str)[0]
				data = np.genfromtxt(datafile,delimiter=',',dtype=np.float64,skip_header=1)
				datashape0,datashape1 = data.shape

				rolling = np.zeros((3,int(datashape1/6),datashape0))
				sliding = np.zeros((3,int(datashape1/6),datashape0))
				print("Starting data analysis")

				j = 0
				k = 0 
				for i in range(0,datashape1,3):
					if i % 2 == 0:
						rolling[:,j,:] = np.transpose(data[:,i:i+3])
						j += 1
					else:
						sliding[:,k,:] = np.transpose(data[:,i:i+3])
						k += 1
				data = []

				rolling_mag = np.zeros((int(datashape1/6),datashape0))
				sliding_mag = np.zeros((int(datashape1/6),datashape0))
				xdata = np.linspace(0,1,datashape0) #dummy x data. x is actually time
				for i in range(rolling.shape[1]):
					for j in range(rolling.shape[2]):
						rolling_mag[i,j] = np.linalg.norm(rolling[:,i,j])
						sliding_mag[i,j] = np.linalg.norm(sliding[:,i,j])

				rolling = []
				sliding = []

				
				print("Starting data save")
				# dat = [sliding_mag,rolling_mag]
				# for i,sav in enumerate([sav_slid,sav_roll]):
				np.savetxt(sav_slid,sliding_mag,delimiter=',')
				np.savetxt(sav_roll,rolling_mag,delimiter=',')
			else:
				print("Getting premade data")
				rolling_mag = np.genfromtxt(sav_roll,delimiter=',')
				sliding_mag = np.genfromtxt(sav_slid,delimiter=',')
				xdata = np.linspace(0,1,rolling_mag.shape[1]) #dummy x data. x is actually time

		
		# ax = plt.axes(projection='3d')
		print("Starting plots")
		fig, ax = plt.subplots(4,1,figsize=(15,10))
		for i in disp_balls:
			ax[i].plot(xdata,rolling_mag[i],label='ball {}'.format(i))
			ax[i].grid(visible=True, which='both', color='0.65', linestyle='-')
			ax[i].minorticks_on()
			ax[i].legend()
		fig.suptitle('rolling')
		plt.tight_layout()
		plt.savefig(p + "RollFric_allb.png".format(i))

		fig, ax = plt.subplots(4,1,figsize=(15,10))
		startind = np.where(xdata>=0.05)[0][0]
		for i in disp_balls:
			ax[i].plot(xdata[startind:],sliding_mag[i][startind:],label='ball {}'.format(i))
			ax[i].grid(visible=True, which='both', color='0.65', linestyle='-')
			ax[i].minorticks_on()
			ax[i].legend()
			# ax[i].set_xlim((0.5,1))
		fig.suptitle('sliding')
		plt.tight_layout()
		plt.savefig(p + "slidFric_allb.png".format(i))

		fig, ax = plt.subplots(8,1,figsize=(15,10))
		startind = np.where(xdata>=0.05)[0][0]
		for i in disp_balls:
			ax[i*2].plot(xdata[startind:],rolling_mag[i][startind:],label='ball {} roll'.format(i))
			ax[i*2].grid(visible=True, which='both', color='0.65', linestyle='-')
			ax[i*2].minorticks_on()
			ax[i*2].legend()
			ax[i*2+1].plot(xdata[startind:],sliding_mag[i][startind:],label='ball {} slide'.format(i))
			ax[i*2+1].grid(visible=True, which='both', color='0.65', linestyle='-')
			ax[i*2+1].minorticks_on()
			ax[i*2+1].legend()
			# ax[i].set_xlim((0.5,1))
		fig.suptitle(p[-9:])
		plt.tight_layout()
		plt.savefig(p + "slidRollFric_allb.png".format(i))

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