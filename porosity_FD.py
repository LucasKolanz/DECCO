import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import csv
sys.path.append("/home/kolanzl/Desktop/SpaceLab")
import utils as u

'''
porosity definitions 1-3 from:
MODELING POROUS DUST GRAINS WITH BALLISTIC AGGREGATES. I.
GEOMETRY AND OPTICAL PROPERTIES
'''

#this function taken from 
#https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
def translate_to_cofm(mass, data):
    # Position of centre of mass in original coordinates
    cofm = sum(mass * data) / (mass*data.shape[0])
    # Transform to CofM coordinates and return
    data -= cofm
    return data

#this function taken from 
#https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
def get_inertia_matrix(mass, data):
	# Moment of intertia tensor
	
	#should in general translate to center of mass
	#but data is already there
	# data = translate_to_cofm(mass, data)

	x, y, z = data.T

	Ixx = np.sum(mass * (y**2 + z**2))
	Iyy = np.sum(mass * (x**2 + z**2))
	Izz = np.sum(mass * (x**2 + y**2))
	Ixy = -np.sum(mass * x * y)
	Iyz = -np.sum(mass * y * z)
	Ixz = -np.sum(mass * x * z)
	I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
	# print(I)
	return I

#this function taken from 
#https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
def get_principal_moi(mass,data):
	I = get_inertia_matrix(mass,data)
	Ip = np.linalg.eigvals(I)
	# Sort and convert principal moments of inertia to SI (kg.m2)
	Ip.sort()
	return Ip

def porosity_measure1(data_folder):
	data,radius,mass,moi = u.get_data(data_folder)
	num_balls = data.shape[0]

	effective_radius = radius*np.power(num_balls,1/3) 
		
	principal_moi = get_principal_moi(mass,data)
	alphai = principal_moi/(0.4*num_balls*mass*effective_radius**2)
	
	a = effective_radius * np.sqrt(alphai[1] + alphai[2] - alphai[0])
	b = effective_radius * np.sqrt(alphai[2] + alphai[0] - alphai[1])
	c = effective_radius * np.sqrt(alphai[0] + alphai[1] - alphai[2])
	
	# Rabc = np.power(a*b*c,1/3)
	porosity = 1-(effective_radius**3/(a*b*c))

	return porosity

def porosity_measure2(data_folder):
	data,radius,mass,moi = u.get_data(data_folder)
	num_balls = data.shape[0]

	effective_radius = radius*np.power(num_balls,1/3) 
		
	principal_moi = get_principal_moi(mass,data)
	alphai = principal_moi/(0.4*num_balls*mass*effective_radius**2)

	RKBM = np.sqrt(np.sum(alphai)/3) * effective_radius

	porosity = 1-np.power((effective_radius/RKBM),3)
	return porosity

if __name__ == '__main__':
	data_prefolder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/tempVariance_attempt'
	data_prefolder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/tempVarianceRand_attempt'

	# temps = [10]
	temps = [3,10,30,100,300,1000]
	Nums = [30,100,300]
	# attempts = [19]
	attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
	attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
	attempts = [i for i in range(1,29)]
	data = []


	new_data = False


	if new_data:
		std_dev = []
		std_err = []
		# attempts = [attempts[3]]
		porositiesabc = np.zeros((len(Nums),len(temps),len(attempts)),dtype=np.float64)
		porositiesKBM = np.zeros((len(Nums),len(temps),len(attempts)),dtype=np.float64) 
		FD_data = np.zeros((len(Nums),len(temps),len(attempts)),dtype=np.float64)
		for i,temp in enumerate(temps):
			for n,N in enumerate(Nums):
				for j,attempt in enumerate(attempts):
					if N == 300 and attempt != 1:
						break
						# data_folder = data_prefolder + str(attempt) + '/' + 'T_' + str(temp) + '/'
					# else:
					data_folder = data_prefolder + str(attempt) + '/' + 'N_' + str(N) + '/T_' + str(temp) + '/'

					count = 0
					for root_dir, cur_dir, files in os.walk(data_folder):
					    count += len(files)
					if count/3 > N:
						porositiesabc[n,i,j] = porosity_measure1(data_folder)
						porositiesKBM[n,i,j] = porosity_measure2(data_folder)

						o3dv = u.o3doctree(data_folder,overwrite_data=False)
						o3dv.make_tree()
						FD_data[n,i,j] = o3dv.calc_fractal_dimension(show_graph=False)
					else:
						porositiesabc[n,i,j] = np.nan
						porositiesKBM[n,i,j] = np.nan
						FD_data[n,i,j] = np.nan
		# porositiesabc = np.array(porositiesabc,dtype=np.float64)
		# porositiesKBM = np.array(porositiesKBM,dtype=np.float64)
		# print(porositiesabc.shape)
		# print(porositiesabc.shape)

		# for i in range(len(attempts)):
		# 	pors = np.array([porositiesabc[:,i],porositiesKBM[:,i],np.array(porositiesabc[:,i])/np.array(porositiesKBM[:,i])])
		# 	plt.plot(temps,pors.T)
		# 	plt.title('Porosity run {}'.format(i))
		# 	plt.xlabel('Temperature in K')
		# 	plt.ylabel('Porosity')
		# 	plt.legend(['Rabc','RKBM','Rabc/RKBM'])
		# 	plt.xscale('log')
		# 	plt.show()

		porositiesabcavg = []
		porositiesKBMavg = []
		porositiesabcstd = []
		porositiesKBMstd = []
		FD_dataavg = []
		FD_datastd = []
		yerr_abc = []
		yerr_KBM = []
		yerr_FD = []
		for i,N in enumerate(Nums[:2]):
			porositiesabcavg.append(np.average(porositiesabc[i],axis=1))
			porositiesKBMavg.append(np.average(porositiesKBM[i],axis=1))
			porositiesabcstd.append(np.std(porositiesabc[i],axis=1))
			porositiesKBMstd.append(np.std(porositiesKBM[i],axis=1))
			FD_dataavg.append(np.average(FD_data[i],axis=1))
			FD_datastd.append(np.std(FD_data[i],axis=1))
			# print(porositiesabcavg[i])
			# print(porositiesKBMavg[i])

			# plotme = np.array([porositiesabcavg,porositiesKBMavg])
			yerr_abc.append(porositiesabcstd[i]/np.sqrt(len(attempts)))
			yerr_KBM.append(porositiesKBMstd[i]/np.sqrt(len(attempts)))
			yerr_FD.append(FD_datastd[i]/np.sqrt(len(attempts)))
			# print(yerr[0])



		#make table
		# fig, ax = plt.subplots()
		# fig.patch.set_visible(False)

		columns = ['Temperature']
		data = [temps]

		for i,N in enumerate(Nums[:2]):
			columns.append('N={} abc porositity'.format(N))
			data.append(porositiesabcavg[i])
			columns.append('N={} abc std err'.format(N)) 
			data.append(yerr_abc[i])
			columns.append('N={} KBM porositity'.format(N))
			data.append(porositiesKBMavg[i])
			columns.append('N={} KBM std err'.format(N))
			data.append(yerr_KBM[i])
			columns.append('N={} Fractal Dimension'.format(N))
			data.append(FD_dataavg[i])
			columns.append('N={} Fractal Dimension std err'.format(N))
			data.append(yerr_FD[i])

		columns.append('N={} abc porositity'.format(300))
		data.append(porositiesabc[2,:,0])
		columns.append('N={} KBM porositity'.format(300))
		data.append(porositiesKBM[2,:,0])
		columns.append('N={} Fractal Dimension'.format(300))
		data.append(FD_data[2,:,0])

		with open('data/averageData.csv','w') as file:
			write = csv.writer(file)
			write.writerow(columns)
			for row in np.transpose(np.array(data)):
				write.writerow(row)
	else:
		data = np.loadtxt('data/averageData.csv',delimiter=',',skiprows=1)




	# ax.axis('off')
	# ax.axis('tight')
	# ax.table(cellText=np.transpose(np.array(data)), colLabels=columns, loc='center')

	# fig.tight_layout()
	# plt.show()

	#plot each size separately
	# for i,N in enumerate(Nums):
	# 	fig,ax = plt.subplots()
	# 	# print(porositiesabc[i])
	# 	ax.errorbar(temps,porositiesabcavg[i],yerr=yerr_abc[i],label="Rabc",zorder=5)
	# 	ax.errorbar(temps,porositiesKBMavg[i],yerr=yerr_KBM[i],label="RKBM",zorder=5)
	# 	# ax.errorbar(temps,plotme[1],yerr=yerr1[1],label="RKBM",zorder=10)
		
	# 	ax.set_xlabel('Temperature in K')
	# 	ax.set_title('Porosity average over {} sims N={}'.format(len(attempts),N))
	# 	ax.set_ylabel('Porosity')
	# 	# ax.set_legend(['Rabc','RKBM'])
	# 	# plt.errorbar(temps,)
	# 	ax.set_xscale('log')

	# 	ax2 = ax.twinx()
	# 	# ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# 	# ax2.invert_yaxis()
	# 	ax2.set_ylabel('Avg Fractal Dimension')

	# 	fig.legend()
	# 	plt.savefig("figures/FractDimandPorosity_N{}.png".format(N))
	# 	plt.show()

	#plot all sizes together
	fig,ax = plt.subplots()
	ax2 = ax.twinx()
	# print(porositiesabc[i])
	styles = ['-','--','-.']
	colors = ['g','b','r']
	length = 6
	data = np.transpose(data)

	# length = data.shape[]
	for i,N in enumerate(Nums[:2]):
		ax.errorbar(temps,data[i*length+1],yerr=data[i*length+2],label="Rabc,N={}".format(N),\
			linestyle=styles[i],color="g",zorder=5)
		ax.errorbar(temps,data[i*length+3],yerr=data[i*length+4],label="RKBM,N={}".format(N),\
			linestyle=styles[i],color="b",zorder=5)
		ax2.errorbar(temps,data[i*length+5],yerr=data[i*length+6],label="FD, N={}".format(N),\
			linestyle=styles[i],color="r",zorder=5)

	# for i,N in enumerate(Nums[:2]):
	# 	ax.errorbar(temps,data[i*length],yerr=data[i*length+1],label="Rabc,N={}".format(N),\
	# 		linestyle=styles[i],color="g",zorder=5)
	# 	ax.errorbar(temps,data[i*length+2],yerr=yerr_KBM[i],label="RKBM,N={}".format(N),\
	# 		linestyle=styles[i],color="b",zorder=5)
	# 	ax2.errorbar(temps,FD_dataavg[i],yerr=yerr_FD[i],label="FD, N={}".format(N),\
	# 		linestyle=styles[i],color="r",zorder=5)


	ax.plot(temps,data[13],label='Rabc,N=300',linestyle=styles[2],\
			color='g',zorder=5)
	ax.plot(temps,data[14],label='RKBM,N=300',linestyle=styles[2],\
			color='b',zorder=5)
	ax2.plot(temps,data[15],label='FD, N=300',linestyle=styles[2],\
			color='r',zorder=5)

	# ax.errorbar(temps,plotme[1],yerr=yerr1[1],label="RKBM",zorder=10)
	
	ax.set_xlabel('Temperature in K')
	ax.set_title('Porosity average over {} sims'.format(len(attempts),N))
	ax.set_ylabel('Porosity')
	# ax.set_legend(['Rabc','RKBM'])
	# plt.errorbar(temps,)
	ax.set_xscale('log')

	# ax2 = ax.twinx()
	# ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# ax2.invert_yaxis()
	ax2.set_ylabel('Avg Fractal Dimension')

	fig.legend()
	plt.savefig("figures/FractDimandPorosity_allN.png")
	plt.show()


	# fig, ax = plt.