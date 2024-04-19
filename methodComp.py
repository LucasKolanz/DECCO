import numpy as np
import matplotlib.pyplot as plt
import os 
import sys
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u
# import utils_old as u

def S(sigma):
	return np.sum(1/sigma**2)

def Si(i,sigma):
	return np.sum(i/(sigma**2))

def Sii(x,y,sigma):
	return np.sum(x*y/(sigma**2))

def main():

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	data_prefolder = path + 'jobsOld/tempVarianceRand_attempt'
	data_prefolder = path + 'jobsNovus/const'
	data_prefolder = path + 'jobsCosine/lognorm'
	data_prefolder = path + 'jobsNovus/const_relax'
	data_prefolder = path + 'jobsCosine/lognorm_relax'

	dataset_name = data_prefolder.split("/")[-1]

	sav = path+'data/{}_averageData.csv'.format(dataset_name)
	# figure_folder = 'figuresCompare/'
	figure_folder = path+'data/figures/'



	temps = [3,10,30,100,300,1000]
	x = np.log(temps)
	# # temps = [1000]
	# Nums = [30,100,300]
	# Nums = [300]
	
	
	# attempts = [i for i in range(30)]
	# attempts = [i for i in range(2)]
	# attempts = [18]

	data = np.loadtxt(sav,delimiter=' ',dtype=np.float64)
	print(f"data.shape: {data.shape}")
	print(data)
	properties = int(data.shape[0]/3) #18 #number of columns for every agg size
	print(properties)



	# data = np.transpose(data)
	# print(data.shape)
	# print(data[0])
	# exit(0)

	# x = np.log(data[0])
	y_data_unordered = data[:]

	# headers = []
	# order = [0,1,2,3,4,5,18,19,6,7,8,9,10,11,20,21,12,13,14,15,16,17,22,23]
	#This order should be Rabc value, Rabc uncert, RKBM value, RKBM uncert, FD value, FD uncert, NC value, NC uncert
	order30 = [0,1,3,4,6,7,9,10]
	order100 = [i+properties for i in order30]
	order300 = [i+2*properties for i in order30]

	order = order30
	order.extend(order100)
	order.extend(order300)


	y_data = np.full(shape=(len(order),len(temps)),fill_value=np.nan,dtype=np.float64)


	for i in range(len(order)):
		y_data[i] = y_data_unordered[order[i]]
		# headers.append(headers_unordered[order[i]+1])


	slope = []
	uncertainty = []
	slope_sigma = []

	for i in range(0,len(y_data),2):
		y = y_data[i]
		sigma = y_data[i+1]

		delta = S(sigma)*Sii(x,x,sigma)-(Si(x,sigma))**2
		slope.append(np.abs((S(sigma)*Sii(x,y,sigma)-Si(x,sigma)*Si(y,sigma))/delta))
		uncertainty.append(np.average(sigma))
		slope_sigma.append(np.sqrt(S(sigma)/delta))


	fig,ax = plt.subplots(figsize=(10,5))
	dummy_x_data = [0,.2,.4,.6]
	for N in [30,100,300]:
		
		
		rang = []
		if N == 30:
			rang = [0,4]
			shift = -0.035
		elif N == 100:
			rang = [4,8]
			shift = 0
		else:
			rang = [8,12]
			shift = 0.035
		
		x_data = [i+shift for i in dummy_x_data]
					# xerr=uncertainty[rang[0]:rang[1]],


		if np.sum(np.isnan(slope[rang[0]:rang[1]])) != rang[1]-rang[0]:
			ax.errorbar(x=x_data,y=slope[rang[0]:rang[1]],
						yerr=slope_sigma[rang[0]:rang[1]],
						fmt='o',linewidth=2, capsize=6,label="N={}".format(N))
		if N == 300:
			print(slope[rang[0]:rang[1]])

	plt.axvline(x = 0.1, color = 'black')
	plt.axvline(x = 0.3, color = 'black')
	plt.axvline(x = 0.5, color = 'black')

	ax.xaxis.set(ticks=dummy_x_data,
			ticklabels=['Rabc','RKBM','FD','NC'])

	ax.set_ylabel('Slope')

	fig.legend()

	plt.savefig(figure_folder+f'/{dataset_name}_methodComp.png')

	plt.show()

		

if __name__ == '__main__':
	main()