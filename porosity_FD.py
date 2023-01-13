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

# def dist(i,j,)

def number_of_contacts(data_folder):
	data,radius,mass,moi = u.get_data(data_folder)
	data = np.array(data)
	num_balls = data.shape[0]

	contacts = np.zeros((num_balls,num_balls),dtype=int)
	dist = lambda i,j: np.sqrt((data[i][0]-data[j][0])**2 + (data[i][1]-data[j][1])**2 + \
			(data[i][2]-data[j][2])**2)

	for i in range(num_balls):
		for j in range(num_balls):
			if i != j:
				contacts[i,j] = (dist(i,j) <= radius*2)
	
	return np.mean(np.sum(contacts,axis=1))

if __name__ == '__main__':
	data_prefolder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/tempVariance_attempt'
	data_prefolder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/tempVarianceRand_attempt'

	# temps = [10]
	temps = [3,10,30,100,300,1000]
	Nums = [30,100,300]
	# Nums = [300]
	# attempts = [19]
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
	# attempts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
	# attempts = [i for i in range(1,5)]
	attempts = [i for i in range(1,29)]
	attempts300 = [1,2,3,4,5,6,7,8,9,10]
	data = []

	new_data = False

	if new_data:
		std_dev = []
		std_err = []
		# attempts = [attempts[3]]
		porositiesabc = np.zeros((len(Nums),len(temps),len(attempts)),dtype=np.float64)
		porositiesKBM = np.zeros((len(Nums),len(temps),len(attempts)),dtype=np.float64) 
		FD_data = np.zeros((len(Nums),len(temps),len(attempts)),dtype=np.float64)
		contacts = np.zeros((len(Nums),len(temps),len(attempts)),dtype=float)
		porositiesabc[:] = np.nan
		porositiesKBM[:] = np.nan
		FD_data[:] = np.nan
		contacts[:] = np.nan
		for i,temp in enumerate(temps):
			for n,N in enumerate(Nums):
				if N == 300:
					a = attempts300
				else:
					a = attempts
				for j,attempt in enumerate(a):
					
					data_folder = data_prefolder + str(attempt) + '/' + 'N_' + str(N) + '/T_' + str(temp) + '/'
					count = 0
					for root_dir, cur_dir, files in os.walk(data_folder):
					    count += len(files)
					if count/3 > N:
						porositiesabc[n,i,j] = porosity_measure1(data_folder)
						porositiesKBM[n,i,j] = porosity_measure2(data_folder)
						contacts[n,i,j] = number_of_contacts(data_folder)

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
		contactsavg = []
		contactsstd = []
		FD_dataavg = []
		FD_datastd = []
		yerr_abc = []
		yerr_KBM = []
		yerr_FD = []
		yerr_ca = []
		for i,N in enumerate(Nums):
			if N == 300:
				a = attempts300
			else:
				a = attempts
			# porositiesabcavg.append(np.average(porositiesabc[i],axis=1))
			# porositiesKBMavg.append(np.average(porositiesKBM[i],axis=1))
			# FD_dataavg.append(np.average(FD_data[i],axis=1))
			# print('N={}\n{}'.format(N,np.nanmean(porositiesabc[i],axis=1)))
			porositiesabcavg.append(np.nanmean(porositiesabc[i],axis=1))
			porositiesKBMavg.append(np.nanmean(porositiesKBM[i],axis=1))
			FD_dataavg.append(np.nanmean(FD_data[i],axis=1))
			contactsavg.append(np.nanmean(contacts[i],axis=1))
			porositiesabcstd.append(np.nanstd(porositiesabc[i],axis=1))
			porositiesKBMstd.append(np.nanstd(porositiesKBM[i],axis=1))
			FD_datastd.append(np.nanstd(FD_data[i],axis=1))
			contactsstd.append(np.nanstd(contacts[i],axis=1))
			# print(porositiesabcavg[i])
			# print(porositiesKBMavg[i])

			# plotme = np.array([porositiesabcavg,porositiesKBMavg])
			yerr_abc.append(porositiesabcstd[i]/np.sqrt(len(a)))
			yerr_KBM.append(porositiesKBMstd[i]/np.sqrt(len(a)))
			yerr_FD.append(FD_datastd[i]/np.sqrt(len(a)))
			yerr_ca.append(contactsstd[i]/np.sqrt(len(a)))
			# print(yerr[0])



		#make table
		# fig, ax = plt.subplots()
		# fig.patch.set_visible(False)

		headers = ['Temperature']
		data = [temps]

		for i,N in enumerate(Nums):
			headers.append('N={} abc porosity'.format(N))
			data.append(porositiesabcavg[i])
			headers.append('N={} abc std err'.format(N)) 
			data.append(yerr_abc[i])
			headers.append('N={} KBM porosity'.format(N))
			data.append(porositiesKBMavg[i])
			headers.append('N={} KBM std err'.format(N))
			data.append(yerr_KBM[i])
			headers.append('N={} Fractal Dimension'.format(N))
			data.append(FD_dataavg[i])
			headers.append('N={} Fractal Dimension std err'.format(N))
			data.append(yerr_FD[i])
		for i,N in enumerate(Nums):
			headers.append('N={} average number of contacts'.format(N))
			data.append(contactsavg[i])
			headers.append('N={} average number of contacts std err'.format(N))
			data.append(yerr_ca[i])

		print(data)

		# headers.append('N={} abc porositity'.format(300))
		# data.append(porositiesabc[2,:,0])
		# headers.append('N={} KBM porositity'.format(300))
		# data.append(porositiesKBM[2,:,0])
		# headers.append('N={} Fractal Dimension'.format(300))
		# data.append(FD_data[2,:,0])
		sav = 'data/averageData.csv'
		with open(sav,'w') as file:
			write = csv.writer(file)
			write.writerow(headers)
			for row in np.transpose(np.array(data)):
				write.writerow(row)
			print("Data saved to {}".format(sav))
	else:
		headers = np.loadtxt('data/averageData.csv',delimiter=',',dtype=str)[0]
		data = np.loadtxt('data/averageData.csv',delimiter=',',skiprows=1,dtype=np.float64)
		data = np.transpose(data)


	styles = ['-','--','-.']
	colors = ['g','b','r']
	length = len(temps)


	# #plot each size separately for num contacts
	for i,N in enumerate(Nums):
		if N == 300:
			a = attempts300
		else:
			a = attempts
		fig,ax = plt.subplots()

		ax.errorbar(temps,data[len(data)-len(Nums)*2+2*i],yerr=data[len(data)-len(Nums)*2+2*i+1],\
					label="Avg # contacts N={}".format(N),color=colors[0],linestyle=styles[i],zorder=5)

		ax.set_xlabel('Temperature in K')
		ax.set_title('Averge number of contacts over {} sims N={}'.format(len(a),N))
		ax.set_ylabel('# Contacts')
		# ax.set_legend(['Rabc','RKBM'])
		# plt.errorbar(temps,)
		ax.set_xscale('log')

		# ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
		# ax2.invert_yaxis()

		fig.legend()
		plt.savefig("figures/avgContacts_N{}.png".format(N))
		plt.show()


	# #plot each size separately
	for i,N in enumerate(Nums):
		if N == 300:
			a = attempts300
		else:
			a = attempts
		fig,ax = plt.subplots()
		# print(porositiesabc[i])
		ax2 = ax.twinx()
		ax.errorbar(temps,data[i*length+1],yerr=data[i*length+2],\
					label="Rabc",color=colors[0],linestyle=styles[i],zorder=5)
		ax.errorbar(temps,data[i*length+3],yerr=data[i*length+4],\
					label="RKBM",color=colors[1],linestyle=styles[i],zorder=5)
		ax2.errorbar(temps,data[i*length+5],yerr=data[i*length+6],\
					label="FD",color=colors[2],linestyle=styles[i],zorder=5)
		
		ax.set_xlabel('Temperature in K')
		ax.set_title('Porosity average over {} sims N={}'.format(len(a),N))
		ax.set_ylabel('Porosity')
		# ax.set_legend(['Rabc','RKBM'])
		# plt.errorbar(temps,)
		ax.set_xscale('log')

		# ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
		# ax2.invert_yaxis()
		ax2.set_ylabel('Avg Fractal Dimension')

		fig.legend()
		plt.savefig("figures/FractDimandPorosity_N{}.png".format(N))
		plt.show()

	# plot each porosity measure separately
	for i,method in enumerate(["Rabc","RKBM","FD"]):
		# if Nums[i] == 300:
		# 	a = attempts300
		# else:
		# 	a = attempts
		
		fig,ax = plt.subplots()
		if i < 2:
			ax.set_ylabel('Porosity')
		else:
			# ax = ax.twinx()
			ax.set_ylabel('Avg Fractal Dimension')

		ax.errorbar(temps,data[2*i+1],yerr=data[2*i+2],\
				label="N={}".format(Nums[0]),color=colors[i],linestyle=styles[0],zorder=5)
		ax.errorbar(temps,data[2*i+length+1],yerr=data[2*i+length+2],\
				label="N={}".format(Nums[1]),color=colors[i],linestyle=styles[1],zorder=5)
		ax.errorbar(temps,data[2*i+length*2+1],yerr=data[2*i+length*2+2],\
				label="N={}".format(Nums[2]),color=colors[i],linestyle=styles[2],zorder=5)
		
		ax.set_xlabel('Temperature in K')
		ax.set_title('Porosity average with {}'.format(method))
		# ax.set_legend(['Rabc','RKBM'])
		# plt.errorbar(temps,)
		ax.set_xscale('log')

		# ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
		# ax2.invert_yaxis()

		fig.legend()
		plt.savefig("figures/FractDimandPorosity_{}.png".format(method))
		plt.show()

	#plot all sizes together
	fig,ax = plt.subplots()
	ax2 = ax.twinx()

	# length = data.shape[]
	for i,N in enumerate(Nums):
		# print('N={}\n{}'.format(N,data[i*length+1]))
		# print('N2={}\n{}'.format(N,data[i*length+2]))
		ax.errorbar(temps,data[i*length+1],yerr=data[i*length+2],label="Rabc,N={}".format(N),\
			linestyle=styles[i],color="g",zorder=5)
		ax.errorbar(temps,data[i*length+3],yerr=data[i*length+4],label="RKBM,N={}".format(N),\
			linestyle=styles[i],color="b",zorder=5)
		ax2.errorbar(temps,data[i*length+5],yerr=data[i*length+6],label="FD, N={}".format(N),\
			linestyle=styles[i],color="r",zorder=5)
	
	ax.set_xlabel('Temperature in K')
	ax.set_title('Average Porosity')
	ax.set_ylabel('Porosity')
	# ax.set_legend(['Rabc','RKBM'])
	# plt.errorbar(temps,)
	ax.set_xscale('log')

	# ax2 = ax.twinx()
	# ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# ax2.invert_yaxis()
	ax2.set_ylabel('Avg Fractal Dimension')

	fig.legend()
	plt.savefig("figures/FractDimandPorosity_all.png")
	plt.show()





	# fig, ax = plt.