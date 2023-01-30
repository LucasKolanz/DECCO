import numpy as np
import matplotlib.pyplot as plt

def S(sigma):
	return np.sum(1/sigma**2)

def Si(i,sigma):
	return np.sum(i/(sigma**2))

def Sii(x,y,sigma):
	return np.sum(x*y/(sigma**2))

def main():
	project_path = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/'
	headers_unordered = np.loadtxt(project_path+"data/averageData.csv",delimiter=',',dtype=str)[0]
	data = np.loadtxt(project_path+"data/averageData.csv",delimiter=',',dtype=np.float64,skiprows=1)
	data = np.transpose(data)
	# print(headers)
	# print(headers.shape)
	# print(np.transpose(data))
	# print(np.transpose(data).shape)

	x = np.log(data[0])
	y_data_unordered = data[1:]

	y_data = np.empty_like(y_data_unordered,dtype=np.float64)
	headers = []
	order = [0,1,2,3,4,5,18,19,6,7,8,9,10,11,20,21,12,13,14,15,16,17,22,23]

	for i in range(y_data.shape[0]):
		y_data[i] = y_data_unordered[order[i]]
		headers.append(headers_unordered[order[i]+1])


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
		ax.errorbar(x=x_data,y=slope[rang[0]:rang[1]],
					yerr=slope_sigma[rang[0]:rang[1]],
					xerr=uncertainty[rang[0]:rang[1]],
					fmt='o',linewidth=2, capsize=6,label="N={}".format(N))

	plt.axvline(x = 0.1, color = 'black')
	plt.axvline(x = 0.3, color = 'black')
	plt.axvline(x = 0.5, color = 'black')

	ax.xaxis.set(ticks=dummy_x_data,
			ticklabels=['Rabc','RKBM','FD','NC'])

	ax.set_ylabel('Slope')

	fig.legend()

	plt.savefig(project_path+'figures/methodComp.png')

	plt.show()

		

if __name__ == '__main__':
	main()