"""
This file was originally written for SpaceLab/DECCO to do data processing

Author: Lucas Kolanz

This file plots the data produced by porosity_FD.py. These are, Porosity abc, Porosity KBM, 
average number of contacts, fractal dimension, bulk density (1-Porosity_KBM), and the final angular momentum. This data is then averaged 
over attempts and saved. The saved data contains the average, the uncertainty, and the number of attempts included in the average.

"""
import matplotlib.pyplot as plt

properties = 18 #number of columns to save for every Num






if __name__ == '__main__':
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	data_prefolder = path + 'jobsOld/tempVarianceRand_attempt'
	data_prefolder = path + 'jobsCosine/lognorm'
	data_prefolder = path + 'jobsNovus/const'
	data_prefolder = path + 'jobsCosine/constMinHmin'
	data_prefolder = path + 'jobsNovus/const_relax'
	data_prefolder = path + 'jobsCosine/lognorm_relax'

	# data_prefolder = path + 'jobs/BAPA'

	dataset_name = data_prefolder.split("/")[-1]

	sav = path+'data/{}_averageData.csv'.format(dataset_name)
	# figure_folder = 'figuresCompare/'
	figure_folder = path+'data/figures/'


	temps = [3,10,30,100,300,1000]
	# temps = [3,10]
	Nums = [30,100,300]
	Nums = [297]
	
	
	attempts = [i for i in range(30)]
	# attempts = [i for i in range(10)]
	# attempts = [i for i in range(2)]
	attempts = [7]



	
	porositiesabcavg = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	porositiesabcstd = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_abc = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	ABC_numruns = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	
	porositiesKBMavg = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	porositiesKBMstd = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_KBM = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	KBM_numruns = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	
	contactsavg = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	contactsstd = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_ca = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	contacts_numruns = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	
	FD_dataavg = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	FD_datastd = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_FD = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	FD_numruns = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	
	angmomavg = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	angmomstd = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_angmom = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	angmom_numruns = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)

	BD_dataavg = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	BD_datastd = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	yerr_BD = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	BD_numruns = np.full(shape=(len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)





	relax = False
	if dataset_name.split('_')[-1] == "relax":
		relax = True

	print(f"relax: {relax}")





	#Do you want to see plots of the data as they are made?
	show_plots = True
	#Do you want to save the plots once they are made?
	save_plots = False
	#Do you want the number of runs next to each point on the plots
	#so you know how many more runs need to finish
	include_totals = True




	# headers = np.loadtxt(sav,delimiter=',',dtype=str)[0]
	print(f"Opening data {sav}")
	data = np.loadtxt(sav,delimiter=' ',skiprows=0,dtype=np.float64)
	# print(data)
	# print(data.shape)
	# print(porositiesabcavg.shape)

	#TODO make this compatible with a single temp
	for i,N in enumerate(Nums):
		# if N == 300:
		# 	print(data[i*properties+9,:])
		# 	print(data[i*properties+10,:])
		# 	exit(0)
		porositiesabcavg[i] = data[i*properties,:]
		yerr_abc[i] = data[i*properties+1,:]
		ABC_numruns[i] = data[i*properties+2,:]

		porositiesKBMavg[i] = data[i*properties+3,:]
		yerr_KBM[i] = data[i*properties+4,:]
		KBM_numruns[i] = data[i*properties+5,:]

		FD_dataavg[i] = data[i*properties+6,:]
		yerr_FD[i] = data[i*properties+7,:]
		FD_numruns[i] = data[i*properties+8,:]

		contactsavg[i] = data[i*properties+9,:]
		yerr_ca[i] = data[i*properties+10,:]
		contacts_numruns[i] = data[i*properties+11,:]

		angmomavg[i] = data[i*properties+12,:]
		yerr_angmom[i] = data[i*properties+13,:]
		angmom_numruns[i] = data[i*properties+14,:]

		BD_dataavg[i] = data[i*properties+15,:]
		yerr_BD[i] = data[i*properties+16,:]
		BD_numruns[i] = data[i*properties+17,:]


	
	print("======================Starting figures======================")
	print(data.shape)
	print("Data has {} nan values".format(np.count_nonzero(np.isnan(data))))
	

	styles = ['-','--','-.','--.']
	colors = ['g','b','r','orange','black','red']
	length = len(temps)

	

	# print(porositiesabcavg)
	# print(yerr_abc)

	# #plot porosity vs size for all temps 
	# fig,ax = plt.subplots()
	# for i,t in enumerate(temps):
	# 	# print(porositiesabcavg)
	# 	# print(shape(porositiesabcavg))
	# 	# exit(0)
	# 	ax.errorbar(Nums,porositiesabcavg[:,i],yerr=yerr_abc[:,i],\
	# 				label="t={}K".format(t),zorder=5)
	# 				# label="t={}K".format(t),color=colors[0],linestyle=styles[i],zorder=5)

	# ax.set_xlabel('Number of particles')
	# # ax.set_title('Measure of porosity vs agg size')
	# ax.set_ylabel('Rabc Porosity')
	# # ax.set_legend(['Rabc','RKBM'])
	# # plt.errorbar(temps,)
	# ax.set_xscale('log')

	# # ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# # ax2.invert_yaxis()

	# fig.legend()
	# plt.savefig("figures/abcVsNum.png")
	# if show_plots:
	# 	plt.show()
	# plt.close("all")
	# # exit(0)

	# # #plot porosity vs size for all temps 
	# fig,ax = plt.subplots()
	# for i,t in enumerate(temps):
	# 	ax.errorbar(Nums,porositiesKBMavg[:,i],yerr=yerr_KBM[:,i],\
	# 				label="t={}K".format(t),zorder=5)
	# 				# label="t={}K".format(t),color=colors[0],linestyle=styles[i],zorder=5)

	# ax.set_xlabel('Number of particles')
	# # ax.set_title('Measure of porosity vs agg size')
	# ax.set_ylabel('RKBM Porosity')
	# # ax.set_legend(['Rabc','RKBM'])
	# # plt.errorbar(temps,)
	# ax.set_xscale('log')

	# # ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# # ax2.invert_yaxis()

	# fig.legend()
	# plt.savefig("figures/KBMVsNum.png")
	# if show_plots:
	# 	plt.show()
	# plt.close("all")

	# # #plot FD vs size for all temps 
	# fig,ax = plt.subplots()
	# for i,t in enumerate(temps):
	# 	ax.errorbar(Nums,FD_dataavg[:,i],yerr=yerr_FD[:,i],\
	# 				label="t={}K".format(t),zorder=5)
	# 				# label="t={}K".format(t),color=colors[0],linestyle=styles[i],zorder=5)

	# ax.set_xlabel('Number of particles')
	# # ax.set_title('Measure of Fractal Dimension vs agg size')
	# ax.set_ylabel('FD')
	# # ax.set_legend(['Rabc','RKBM'])
	# # plt.errorbar(temps,)
	# ax.set_xscale('log')

	# # ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# # ax2.invert_yaxis()

	# fig.legend()
	# plt.savefig("figures/FDVsNum.png")
	# if show_plots:
	# 	plt.show()
	# plt.close("all")

	# # print(contactsavg[:,0])
	# # print(yerr_ca[:,0])

	# # #plot num contacts vs size for all temps 
	# fig,ax = plt.subplots()
	# for i,t in enumerate(temps):
	# 	ax.errorbar(Nums,contactsavg[:,i],yerr=yerr_ca[:,i],\
	# 				label="t={}K".format(t),zorder=5)
	# 				# label="t={}K".format(t),color=colors[0],linestyle=styles[i],zorder=5)

	# ax.set_xlabel('Number of particles')
	# # ax.set_title('Averge contacts vs agg size')
	# ax.set_ylabel('Number of contacts')
	# # ax.set_legend(['Rabc','RKBM'])
	# # plt.errorbar(temps,)
	# ax.set_xscale('log')

	# # ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# # ax2.invert_yaxis()

	# fig.legend()
	# plt.savefig("figures/ContactsVsNum.png")
	# if show_plots:
	# 	plt.show()
	# plt.close("all")

	# # #plot each size separately
	# for i,N in enumerate(Nums):
	# 	if N == 300:
	# 		a = attempts300
	# 	else:
	# 		a = attempts
	# 	fig,ax = plt.subplots()
	# 	# print(porositiesabc[i])
	# 	ax2 = ax.twinx()
	# 	ax.errorbar(temps,data[i*length+1],yerr=data[i*length+2],\
	# 				label="Rabc",color=colors[0],linestyle=styles[i],zorder=5)
	# 	ax.errorbar(temps,data[i*length+3],yerr=data[i*length+4],\
	# 				label="RKBM",color=colors[1],linestyle=styles[i],zorder=5)
	# 	ax2.errorbar(temps,data[i*length+5],yerr=data[i*length+6],\
	# 				label="FD",color=colors[2],linestyle=styles[i],zorder=5)
		
	# 	ax.set_xlabel('Temperature in K')
	# 	# ax.set_title('Porosity average over {} sims N={}'.format(len(a),N))
	# 	ax.set_ylabel('Porosity')
	# 	# ax.set_legend(['Rabc','RKBM'])
	# 	# plt.errorbar(temps,)
	# 	ax.set_xscale('log')

	# 	# ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# 	# ax2.invert_yaxis()
	# 	ax2.set_ylabel('Avg Fractal Dimension')

	# 	fig.legend()
	# 	plt.savefig("figures/FractDimandPorosity_N{}.png".format(N))
	# 	if show_plots:
	# 		plt.show()
	#	plt.close("all")
	plt.rcParams.update({
	    'font.size': 18,
	    'text.usetex': True,
	    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	})
	# plot each porosity measure separately
	# plt.rcParams.update({'font.size': 15})
	for i,method in enumerate(["Rabc","RKBM","FD","# Contacts","Ang mom","Bulk Density"]):
		# if Nums[i] == 300:
		# 	a = attempts300
		# else:
		# 	a = attempts
		
		fig,ax = plt.subplots()
		if i == 0:
			ax.set_ylabel(r'$\bm{\mathcal{P}_{abc}}$')
		elif i == 1:
			ax.set_ylabel(r'$\bm{\mathcal{P}_{KBM}}$')
		elif i == 2:
			# ax = ax.twinx()
			ax.set_ylabel(r'$\bm{\mathcal{D}_{f}}$')
		elif i == 3:
			ax.set_ylabel(r'$\bm{\mathcal{ANC}}$')
		elif i == 4:	
			ax.set_ylabel('Total Angular Momentum')
		elif i == 5:
			ax.set_ylabel('Avg Bulk Density')


	# for i,N in enumerate(Nums):
	# 			headers.append('N={} abc porosity'.format(N))
	# 			data[i*6,:] = porositiesabcavg[i]
	# 			headers.append('N={} abc std err'.format(N)) 
	# 			data[i*6+1,:] = yerr_abc[i]
	# 			headers.append('N={} KBM porosity'.format(N))
	# 			data[i*6+2,:] = porositiesKBMavg[i]
	# 			headers.append('N={} KBM std err'.format(N))
	# 			data[i*6+3,:] = yerr_KBM[i]
	# 			headers.append('N={} Fractal Dimension'.format(N))
	# 			data[i*6+4,:] = FD_dataavg[i]
	# 			headers.append('N={} Fractal Dimension std err'.format(N))
	# 			data[i*6+5,:] = yerr_FD[i]
	# 		for i,N in enumerate(Nums):
	# 			headers.append('N={} average number of contacts'.format(N))
	# 			data[6*2 + i*2,:] = contactsavg[i]
	# 			headers.append('N={} average number of contacts std err'.format(N))
	# 			data[6*2 + 1 + i*2,:] = yerr_ca[i]

		# Plots with points
		# for j,N in enumerate(Nums):
		# 	ax.errorbar(temps,data[3*i+properties*j],yerr=data[3*i+1+properties*j],\
		# 			label="N={}".format(N),color=colors[i],\
		# 			linestyle=styles[j],marker='.',markersize=10,zorder=5)
		# 	print(f"N={N} number of completed sims:")
		# 	print(data[3*i+2+properties*j])
		# 	# ax.errorbar(temps,data[2*i+length+1],yerr=data[2*i+length+2],\
		# 	# 		label="N={}".format(Nums[1]),color=colors[i],\
		# 	# 		linestyle=styles[1],marker='.',markersize=10,zorder=5)
		# 	# ax.errorbar(temps,data[2*i+length*2+1],yerr=data[2*i+length*2+2],\
		# 	# 		label="N={}".format(Nums[2]),color=colors[i],\
		# 	# 		linestyle=styles[2],marker='.',markersize=10,zorder=5)
		# # else:
		# # 	for j,N in enumerate(Nums):
		# # 		ax.errorbar(temps,data[len(data)-len(Nums)*2+2*j],yerr=data[len(data)-len(Nums)*2+2*j+1],\
		# # 				label="N={}".format(N),color='orange',\
		# # 				linestyle=styles[j],marker='.',markersize=10,zorder=5)
		# bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		# # print(bbox.width, bbox.height)
		# ax.set_xlabel('Temperature in K')
		# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
		# # ax.set_legend(['Rabc','RKBM'])
		# # plt.errorbar(temps,)
		# ax.set_xscale('log')

		# # ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
		# # ax2.invert_yaxis()
		# # if i == 1:
		# if True:
		# 	fig.legend(loc='upper right',bbox_to_anchor=(0.98, 0.97))
		# plt.tight_layout()
		# plt.savefig("{}{}_FractDimandPorosity_{}.png".format(figure_folder,dataset_name,method.replace(" ","")))
		# if show_plots:
		# 	plt.show()
		# 	print()

		#Plots with numbers of sims
		for j,N in enumerate(Nums):
			ax.errorbar(temps,data[3*i+properties*j],yerr=data[3*i+1+properties*j],\
					label="N={}".format(N),color=colors[i],\
					linestyle=styles[j],marker='.',markersize=10,zorder=5)

			if include_totals:
				for k, txt in enumerate(data[3*i+2+properties*j]):
					ax.annotate("{:0.0f}".format(txt), (temps[k], data[3*i+properties*j][k]))

		bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		ax.set_xlabel('Temperature in K')
		# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
		ax.set_xscale('log')
		if i == 1:
			fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
		plt.tight_layout()
		if save_plots:
			plt.savefig("{}{}_{}_avgPlot.png".format(figure_folder,dataset_name,method.replace(" ","")))
		if show_plots:
			plt.show()

	plt.close("all")


	# # #plot each size for num contacts
	# fig,ax = plt.subplots()
	# for i,N in enumerate(Nums):
	# 	# if N == 300:
	# 	# 	a = attempts300
	# 	# else:
	# 	# 	a = attempts

	# 	ax.errorbar(temps,data[len(data)-len(Nums)*2+2*i],yerr=data[len(data)-len(Nums)*2+2*i+1],\
	# 				label="N={}".format(N),color='orange',linestyle=styles[i],zorder=5)

	# ax.set_xlabel('Temperature in K')
	# # ax.set_title('Averge number of contacts over {} sims N={}'.format(len(a),N))
	# ax.set_ylabel('# Contacts')
	# # ax.set_legend(['Rabc','RKBM'])
	# # plt.errorbar(temps,)
	# ax.set_xscale('log')

	# # ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# # ax2.invert_yaxis()

	# # fig.legend()
	# plt.savefig("figures/avgContacts.png")
	# if show_plots:
	# 	plt.show()
	# plt.close("all")
	# plt.rcParams.update({'font.size': 10})


	#plot all sizes together
	# fig,ax = plt.subplots()
	# ax2 = ax.twinx()

	# # length = data.shape[]
	# for i,N in enumerate(Nums):
	# 	# print('N={}\n{}'.format(N,data[i*length+1]))
	# 	# print('N2={}\n{}'.format(N,data[i*length+2]))
	# 	ax.errorbar(temps,data[i*length+1],yerr=data[i*length+2],label="Rabc,N={}".format(N),\
	# 		linestyle=styles[i],color="g",zorder=5)
	# 	ax.errorbar(temps,data[i*length+3],yerr=data[i*length+4],label="RKBM,N={}".format(N),\
	# 		linestyle=styles[i],color="b",zorder=5)
	# 	ax2.errorbar(temps,data[i*length+5],yerr=data[i*length+6],label="FD, N={}".format(N),\
	# 		linestyle=styles[i],color="r",zorder=5)
	
	# ax.set_xlabel('Temperature in K')
	# # ax.set_title('Average Porosity')
	# ax.set_ylabel('Porosity')
	# # ax.set_legend(['Rabc','RKBM'])
	# # plt.errorbar(temps,)
	# ax.set_xscale('log')

	# # ax2 = ax.twinx()
	# # ax2.errorbar(temps,FD_dataavg,yerr=yerr2,label="Frac Dim",color='r',zorder=0)
	# # ax2.invert_yaxis()
	# ax2.set_ylabel('Avg Fractal Dimension')

	# fig.legend()
	# plt.savefig("figures/FractDimandPorosity_all.png")
	# if show_plots:
	# 	plt.show()





	# fig, ax = plt.

