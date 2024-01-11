import numpy as np
import matplotlib.pyplot as plt
import sys
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
	data_prefolder = '/home/kolanzl/Desktop/SpaceLab/jobs/tempVariance_attempt'

	temps = [3,10,30,100,300,1000]
	attempts = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
	std_dev = []
	std_err = []
	# attempts = [attempts[3]]
	porositiesabc = np.zeros((len(temps),len(attempts)),dtype=np.float64)
	porositiesKBM = np.zeros((len(temps),len(attempts)),dtype=np.float64) 
	FD_data = np.zeros((len(temps),len(attempts)),dtype=np.float64)
	for i,temp in enumerate(temps):
		for j,attempt in enumerate(attempts):
			data_folder = data_prefolder + str(attempt) + '/' + 'T_' + str(temp) + '/'
			porositiesabc[i,j] = porosity_measure1(data_folder)
			porositiesKBM[i,j] = porosity_measure2(data_folder)

			o3dv = u.o3doctree(data_folder,overwrite_data=False)
			o3dv.make_tree()
			FD_data[i,j] = o3dv.calc_fractal_dimension(show_graph=False)

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


	porositiesabcavg = np.average(porositiesabc,axis=1)
	porositiesKBMavg = np.average(porositiesKBM,axis=1)
	porositiesabcstd = np.std(porositiesabc,axis=1)
	porositiesKBMstd = np.std(porositiesKBM,axis=1)
	FD_dataavg = np.average(FD_data,axis=1)
	FD_datastd = np.std(FD_data,axis=1)

	plotme = np.array([porositiesabcavg,porositiesKBMavg])
	yerr1 = np.array([porositiesabcstd,porositiesKBMstd])/np.sqrt(len(attempts))
	yerr2 = FD_datastd/np.sqrt(len(attempts))
	# print(yerr[0])

	fig,ax = plt.subplots()
	ax.errorbar(temps,plotme[0],yerr=yerr1[0],label="Rabc",zorder=5)
	ax.errorbar(temps,plotme[1],yerr=yerr1[1],label="RKBM",zorder=10)
	# plt.errorbar(x, y + 3, yerr=yerr, label='both limits (default)')
	ax.set_xlabel('Temperature in K')
	ax.set_title('Porosity average over {} sims'.format(len(attempts)))
	ax.set_ylabel('Porosity')
	# ax.set_legend(['Rabc','RKBM'])
	# plt.errorbar(temps,)
	ax.set_xscale('log')

	ax2 = ax.twinx()
	ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	ax2.invert_yaxis()
	ax2.set_ylabel('Avg Fractal Dimension')

	fig.legend()
	plt.savefig("FractDimandPorosity.png")
	plt.show()
	#still need to do lots of other porosity measures


	'''
	Add standard deviation and standard error (std error / sqrt(N))
	'''
