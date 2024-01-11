import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
	# path = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/"
	# fileprefix = "1_2_R1e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_"
	path = []
	fileprefix = []
	path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls11/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# path.append("/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/small_balls6/N_10/T_1/")
	# fileprefix.append("2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	fileprefix.append("1_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("2_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")
	# fileprefix.append("9_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_")

	
	# print(data)
	# print(data[1])
	# # exit(0)
	for pind,p in enumerate(path):
		datafile = p + fileprefix[pind] + "simData.csv"

		header = np.genfromtxt(datafile,delimiter=',',dtype=str)[0]
		data = np.genfromtxt(datafile,delimiter=',',dtype=np.float64,skip_header=1)
		# data = data.transpose()
		# print(data[0].size)
		# print(header)

		pos = np.zeros((3,int(data[0].size/11),data.shape[0]))

		for i in range(int(data[0].size/11)):
			for j in [0,1,2]:
				pos[j,i,:] = data[:,i*11+j]




		# fig, ax = plt.subplots(2,3)
		# ax = plt.axes(projection='3d')
		# for i in range(int(data[0].size/11)):
		# 	ax.plot3D(pos[0,i,:],pos[1,i,:],pos[2,i,:])
		# 	ax.set_xlim3d(-.85e-5,-0.75e-5)
		# 	# ax.set_xlim3d(-1e-5,1e-5)
		# 	ax.set_ylim3d(-9e-6,1e-5)
		# 	ax.set_zlim3d(-8e-6,-6e-7)
		# 	ax.set_xlabel('x')
		# 	ax.set_ylabel('y')
		# 	ax.set_zlabel('z')
		# ax.legend()


		# ax = plt.axes(projection='3d')
		maxind = int(data[0].size/11)-1
		startind = int(pos.shape[2]/10)
		# i = 0
		# print(np.sqrt((pos[0,i,:]-pos[0,maxind,:])+(pos[1,i,:]-pos[1,maxind,:])+(pos[2,i,:]-pos[2,maxind,:])))
		# exit(0)
		title = ""
		xdata = np.linspace(0,1,pos.shape[2]) #dummy x data. x is actually time
		for j in [0,1,2]:
			fig, ax = plt.subplots(3,1)
			for i in range(int(data[0].size/11)-1):
				dists = np.sqrt((pos[0,i,:]-pos[0,maxind,:])**2+(pos[1,i,:]-pos[1,maxind,:])**2+(pos[2,i,:]-pos[2,maxind,:])**2)
				if j == 0:
					ax[i].set_title('dist between balls 3 and {} (all)'.format(i))
					ax[i].plot(xdata, dists)
					title = "distB3all.png"
				elif j == 1:
					ax[i].set_title('dist between balls 3 and {} (front)'.format(i))
					ax[i].plot(xdata[0:startind], dists[0:startind])
					title = "distB3front.png"
				else:
					ax[i].set_title('dist between balls 3 and {} (tail)'.format(i))
					ax[i].plot(xdata[startind:], dists[startind:])
					title = "distB3tail.png"

				ax[i].set_xlabel('time (t/tmax)')
				ax[i].set_ylabel('dist (cm)')

			plt.savefig(p+fileprefix[pind]+title)
			plt.tight_layout()
	plt.show()

if __name__ == '__main__':
	main()