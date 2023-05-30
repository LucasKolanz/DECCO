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
	Nums = [30,100,300]
	method = [r'$R_{abc}$',r'$R_{KBM}$',r'$Fractal\ Dimension$',r'$Number\ of\ Contacts$']

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

	# print(headers)
	# print(len(headers))
	# exit(0)

	# fig,ax = plt.subplots(figsize=(10,5))
	plt.rcParams.update({'font.size': 10}) 
	# plt.rcParams["font.family"] = "Times New Roman"
	fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(10,5))
	# plt.rcParams["mathtext.fontset"] = True
	# dummy_x_data = [0,.2,.4,.6]
	limits = ()
	for n,N in enumerate(Nums):
		for m,M in enumerate(method):
			
			index = m + 4*n
			ax[m].errorbar(x=N,y=slope[index],
						yerr=slope_sigma[index],
						xerr=uncertainty[index],
						fmt='o',linewidth=2, capsize=6)
			# if M == method[-1] and N == Nums[-1]:
			# 	limits = ax[m].get_ylim()

	for i,a in enumerate(ax):
		if a in ax[1:]:
			a.get_yaxis().set_visible(False)

		a.text(i+0.5,0.95,method[i],horizontalalignment='center',
	 			verticalalignment='center',transform=ax[0].transAxes)
		a.set_xscale('log')
		# a.set_ylim(limits)
		a.set_xlim((20,450))
		a.set_ylim((-0.001,0.035))
	fig.subplots_adjust(wspace=0.0)
	# plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
	fig.supxlabel("Aggregate size")
	fig.supylabel("Slope")
	# plt.axvline(x = 0.1, color = 'black')
	# plt.axvline(x = 0.3, color = 'black')
	# plt.axvline(x = 0.5, color = 'black')

	# ax.xaxis.set(ticks=dummy_x_data,
			# ticklabels=['Rabc','RKBM','FD','NC'])

	# ax.set_ylabel('Slope')

	# fig.legend()

	plt.savefig(project_path+'figures/methodComp.png')

	plt.show()


	# print(headers)
	# print(y_data)

	# Pabc30 = y_data[0]
	# PKBM30 = y_data[2]
	# Pabc100 = y_data[8]
	# PKBM100 = y_data[10]
	# Pabc300 = y_data[16]
	# PKBM300 = y_data[18]

	# print(Pabc30/PKBM30)
	# print(Pabc100/PKBM100)
	# print(Pabc300/PKBM300)
	# p30 = Pabc30/PKBM30
	# p100 = Pabc100/PKBM100
	# p300 = Pabc300/PKBM300

	# avgPdifftot = (np.sum(p30) + np.sum(p100) + np.sum(p300))/(len(p30)+len(p100)+len(p300)) 
	# avgPdiff30 = (np.sum(p30))/(len(p30)) 
	# avgPdiff100 = (np.sum(p100))/(len(p100)) 
	# avgPdiff300 = (np.sum(p300))/(len(p300)) 
	# print("Overall, on average, PKBM is {:.3f}% larger than Pabc".format(100*(1-avgPdifftot)))
	# print("For N=30, on average, PKBM is {:.3f}% larger than Pabc".format(100*(1-avgPdiff30)))
	# print("For N=100, on average, PKBM is {:.3f}% larger than Pabc".format(100*(1-avgPdiff100)))
	# print("For N=300, on average, PKBM is {:.3f}% larger than Pabc".format(100*(1-avgPdiff300)))

	# print("")

	# print(PKBM30/Pabc30)
	# print(PKBM100/Pabc100)
	# print(PKBM300/Pabc300)

		

if __name__ == '__main__':
	main()