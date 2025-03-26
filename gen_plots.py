"""
This file was originally written for SpaceLab/DECCO to do data processing

Author: Lucas Kolanz

This file plots the data produced by porosity_FD.py. These are, Porosity abc, Porosity KBM, 
average number of contacts, fractal dimension, bulk density (1-Porosity_KBM), and the final angular momentum. This data is then averaged 
over attempts and saved. The saved data contains the average, the uncertainty, and the number of attempts included in the average.

"""
import os
import sys
import json
import matplotlib.pyplot as plt
import numpy as np

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
# import utils as u

import gen_data as gd

styles = ['-','--','-.',':']
colors = ['g','b','orange','r']

def Tanaka(sizes,initRg,temp):
	Kb = 1.380649e-16; #in erg/K
	rho0 = 2.25 #g/cm^3
	B = 0.4
	Pkbm_c = 0.55

	m = rho0*(4/3)*np.pi*(1e-5)**3
	rhomax = np.sqrt(2)*np.pi*rho0/6
	avgvels = np.sqrt(8*Kb*temp/(np.pi*m))
	Eimp = 0.5*m*avgvels**2
	rhoc = (1-Pkbm_c)*rho0
	Eroll = (1/(4*np.pi))*(3*B*m*avgvels**2*rho0**2)*(rhoc**(-1)-rhomax**(-1))**2##Eimp*1000 # I have no idea what this should be for us
	initV = (4*np.pi/3)*(5/3)**(3/2)*initRg**3
	
	print(f"Eimp: {Eimp}")
	print(f"Eroll: {Eroll}")
	print(f"avgvels: {avgvels}")

	Vols = [initV]
	for i,size in enumerate(sizes[:-1]):
		# print(1/np.sqrt((Vols[i]/(size-1) - m/rhomax)**(-2) +(3*B*Eimp*rho0**2)/(2*np.pi*size*Eroll*m**2)))
		Vfin = size*((1/np.sqrt((Vols[i]/(size-1) - m/rhomax)**(-2) +\
						(3*B*Eimp*rho0**2)/(2*np.pi*size*Eroll*m**2))) + m/rhomax)

		Vols.append(Vfin)


	return 1-np.array(sizes)*(m/rho0)/np.array(Vols)


def label_from_header(header):

	if header == gd.data_headers[0]:
		return r'$\bm{\mathcal{P}_{abc}}$'
	elif header == gd.data_headers[1]:
		return r'$\bm{\mathcal{P}_{KBM}}$'
	elif header == gd.data_headers[2]:
		return r'$\bm{\mathcal{ANC}}$'
	elif header == gd.data_headers[3]:
		return r'$\bm{\mathcal{D}_{f}}$'
	else:
		return ""


def gen_relax_vs_tense_seqstick_plots(distribution,show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	tense_data_prefolder = path + 'jobs/SeqStickLognorm_'
	relax_data_prefolder = path + 'jobs/SeqStickLognormrelax_'

	dataset_name = tense_data_prefolder.split("/")[-1]

	figure_folder = path+'data/figures/'


	temps = [1000]
	# temps = [3,10]
	Nums = [300]

	
	
	attempts = [i for i in range(30)]


	# requested_data_headers = gd.data_headers[:2]
	requested_data_headers = gd.data_headers[:2] + [gd.data_headers[-1]]



	tense_raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	relax_raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	for a_i,a in enumerate(attempts):
		for n_i,n in enumerate(Nums):
			for t_i,t in enumerate(temps):
				for data_prefolder in [tense_data_prefolder,relax_data_prefolder]:
					rel = ""
					if data_prefolder == tense_data_prefolder:
						raw_data = tense_raw_data
					elif data_prefolder == relax_data_prefolder:
						raw_data = relax_raw_data
						rel="relax_"

					folder = f"{data_prefolder}{a}/N_{n}/"
					print(folder)

					if os.path.exists(folder+f"{rel}job_data.csv"):
						with open(folder+f"{rel}job_data.csv",'r') as fp:
							existing_data = fp.readlines()

						existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
						#even though the data can have other sizes in it, 
						#we only want the data of size n
						if n not in existing_sizes:
							print(f"ERROR: Data of size {n} does not exist for {folder}.")
							continue
						index = existing_sizes.index(n)*4
						existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
						existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
						print(existing_headers_for_size)
						print(existing_values_for_size)

						for h_i,header in enumerate(requested_data_headers):
							if header in existing_headers_for_size:
								print(raw_data.size)
								print(existing_headers_for_size.index(header))
								raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]

	tense_avg_data = np.nanmean(tense_raw_data,axis=1)
	tense_std_data = np.nanstd(tense_raw_data,axis=1)
	tense_num_data = np.count_nonzero(~np.isnan(tense_raw_data),axis=1)
	tense_err_data = tense_std_data/np.sqrt(tense_num_data)

	relax_avg_data = np.nanmean(relax_raw_data,axis=1)
	relax_std_data = np.nanstd(relax_raw_data,axis=1)
	relax_num_data = np.count_nonzero(~np.isnan(relax_raw_data),axis=1)
	relax_err_data = relax_std_data/np.sqrt(relax_num_data)


	
	print("======================Starting figures======================")
	# print(data.shape)
	print("Tense data has {} nan values".format(np.count_nonzero(np.isnan(tense_avg_data))))
	print("Relax data has {} nan values".format(np.count_nonzero(np.isnan(relax_avg_data))))
	
	print("For tense data:")
	for h_i,header in enumerate(requested_data_headers):
		print(f"\t{distribution} {header}: {tense_avg_data[h_i]} +- {tense_err_data[h_i]}")

	print("For relax data:")
	for h_i,header in enumerate(requested_data_headers):
		print(f"\t{distribution} {header}: {relax_avg_data[h_i]} +- {relax_err_data[h_i]}")


def gen_relax_vs_tense_BPCA_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	tense_data_prefolder = path + 'jobsCosine/lognorm_'
	relax_data_prefolder = path + 'jobsCosine/lognormrelax_'

	dataset_name = tense_data_prefolder.split("/")[-1]

	figure_folder = path+'data/figures/'


	temps = [3,10,30,100,300,1000]
	# temps = [3,10]
	Nums = [30,100,300]

	
	
	attempts = [i for i in range(30)]


	# requested_data_headers = gd.data_headers[:2]
	requested_data_headers = gd.data_headers[:2] + [gd.data_headers[-1]]



	tense_raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	relax_raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	for a_i,a in enumerate(attempts):
		for n_i,n in enumerate(Nums):
			for t_i,t in enumerate(temps):
				for data_prefolder in [tense_data_prefolder,relax_data_prefolder]:
					rel = ""
					if data_prefolder == tense_data_prefolder:
						raw_data = tense_raw_data
					elif data_prefolder == relax_data_prefolder:
						raw_data = relax_raw_data
						rel="relax_"

					folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
					print(folder)

					if os.path.exists(folder+f"{rel}job_data.csv"):
						with open(folder+f"{rel}job_data.csv",'r') as fp:
							existing_data = fp.readlines()

						existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
						#even though the data can have other sizes in it, 
						#we only want the data of size n
						if n not in existing_sizes:
							print(f"ERROR: Data of size {n} does not exist for {folder}.")
							continue
						index = existing_sizes.index(n)*4
						existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
						existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
						
						for h_i,header in enumerate(gd.data_headers):
							if header in existing_headers_for_size:

								raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]

	tense_avg_data = np.nanmean(tense_raw_data,axis=1)
	tense_std_data = np.nanstd(tense_raw_data,axis=1)
	tense_num_data = np.count_nonzero(~np.isnan(tense_raw_data),axis=1)
	tense_err_data = tense_std_data/np.sqrt(tense_num_data)

	relax_avg_data = np.nanmean(relax_raw_data,axis=1)
	relax_std_data = np.nanstd(relax_raw_data,axis=1)
	relax_num_data = np.count_nonzero(~np.isnan(relax_raw_data),axis=1)
	relax_err_data = relax_std_data/np.sqrt(relax_num_data)


	
	print("======================Starting figures======================")
	# print(data.shape)
	print("Tense data has {} nan values".format(np.count_nonzero(np.isnan(tense_avg_data))))
	print("RElax data has {} nan values".format(np.count_nonzero(np.isnan(relax_avg_data))))
	

	
	length = len(temps)


	#	plt.close("all")
	plt.rcParams.update({
	    'font.size': 18,
	    'text.usetex': True,
	    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	})

	#Plot metric vs M for all metrics and all N and temps
	for h_i,header in enumerate(requested_data_headers):
		for n_i,n in enumerate(Nums):


			fig,ax = plt.subplots()

			print(tense_avg_data[h_i,n_i,:])
			print(relax_avg_data[h_i,n_i,:])

			ax.errorbar(temps,tense_avg_data[h_i,n_i,:],yerr=tense_err_data[h_i,n_i,:],\
					label=f"tense N={n}",\
					linestyle=styles[h_i],marker='.',markersize=10,zorder=5)
			ax.errorbar(temps,relax_avg_data[h_i,n_i,:],yerr=relax_err_data[h_i,n_i,:],\
					label=f"relax N={n}",\
					linestyle=styles[h_i],marker='.',markersize=10,zorder=5)

			# if include_totals:
			# 	for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
			# 		ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

			bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
			ax.set_xlabel('Temp K')
			ax.set_ylabel(header)
			# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
			ax.set_xscale('log')
			# if i == 1:
			fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
			plt.tight_layout()
			if save_plots:
				plt.savefig("{}{}_{}_tenseVsRelax.png".format(figure_folder,dataset_name,header))
			if show_plots:
				plt.show() 

def gen_BAPA_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	data_prefolder = path + 'jobs/BAPA_'

	dataset_name = data_prefolder.split("/")[-1]

	figure_folder = path+'data/figures/'


	temps = [1000]
	# temps = [3,10]
	Nums = [300]
	M = [1,20,30,50,60,100]
	
	
	attempts = [i for i in range(20)]


	# data_shape = (len(M),len(Nums),len(temps))
	


	raw_data = np.full(shape=(len(gd.data_headers),len(attempts),len(M),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	for a_i,a in enumerate(attempts):
		for m_i,m in enumerate(M):
			for n_i,n in enumerate(Nums):
				for t_i,t in enumerate(temps):
					folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
					if os.path.exists(folder+"job_data.csv"):
						with open(folder+"job_data.csv",'r') as fp:
							existing_data = fp.readlines()

						existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
						#even though the data can have other sizes in it, 
						#we only want the data of size n
						if n not in existing_sizes:
							print(f"ERROR: Data of size {n} does not exist for {folder}.")
							continue
						index = existing_sizes.index(n)*4
						existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
						existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
						
						for h_i,header in enumerate(gd.data_headers):
							if header in existing_headers_for_size:
								raw_data[h_i,a_i,m_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]

	avg_data = np.nanmean(raw_data,axis=1)
	std_data = np.nanstd(raw_data,axis=1)
	num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
	err_data = std_data/np.sqrt(num_data)

	
	print("======================Starting figures======================")
	# print(data.shape)
	print("Data has {} nan values".format(np.count_nonzero(np.isnan(avg_data))))
	


	length = len(temps)


	#	plt.close("all")
	plt.rcParams.update({
	    'font.size': 18,
	    'text.usetex': True,
	    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	})

	#Plot metric vs M for all metrics and all N and temps
	for h_i,header in enumerate(gd.data_headers):
		for n_i,n in enumerate(Nums):
			for t_i,t in enumerate(temps):

				fig,ax = plt.subplots()


				ax.errorbar(M,avg_data[h_i,:,n_i,t_i],yerr=err_data[h_i,:,n_i,t_i],\
						label=f"N={n},T={t}",color=colors[h_i],\
						linestyle=styles[h_i],marker='.',markersize=10,zorder=5)

				if include_totals:
					for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
						ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

				bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
				ax.set_xlabel('Fragment size')
				ax.set_ylabel(header)
				# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
				# ax.set_xscale('log')
				# if i == 1:
				fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
				plt.tight_layout()
				if save_plots:
					plt.savefig("{}{}_{}_avgPlot.png".format(figure_folder,dataset_name,header))
				if show_plots:
					plt.show() 



def gen_BPCA_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	data_prefolder = path + 'jobsNovus/constrelax_'
	data_prefolder = path + 'jobsCosine/lognormrelax_'

	dataset_name = data_prefolder.split("/")[-1].strip("_")
	figure_folder = path+f'data/figures/BPCA_per_step/{dataset_name}/'

	if save_plots and not os.path.exists(figure_folder):
		os.makedirs(figure_folder)


	temps = [1000]
	# temps = [3,10]
	Nums = [300]
	M = [1,20,30,50,60,100]
	
	
	attempts = [i for i in range(20)]


	# data_shape = (len(M),len(Nums),len(temps))
	


	raw_data = np.full(shape=(len(gd.data_headers),len(attempts),len(M),len(Nums),len(temps)),fill_value=np.nan,dtype=np.float64)
	for a_i,a in enumerate(attempts):
		for m_i,m in enumerate(M):
			for n_i,n in enumerate(Nums):
				for t_i,t in enumerate(temps):
					folder = f"{data_prefolder}{a}/M_{m}/N_{n}/T_{t}/"
					if os.path.exists(folder+"job_data.csv"):
						with open(folder+"job_data.csv",'r') as fp:
							existing_data = fp.readlines()

						existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
						#even though the data can have other sizes in it, 
						#we only want the data of size n
						if n not in existing_sizes:
							print(f"ERROR: Data of size {n} does not exist for {folder}.")
							continue
						index = existing_sizes.index(n)*4
						existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
						existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
						
						for h_i,header in enumerate(gd.data_headers):
							if header in existing_headers_for_size:
								raw_data[h_i,a_i,m_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]

	avg_data = np.nanmean(raw_data,axis=1)
	std_data = np.nanstd(raw_data,axis=1)
	num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
	err_data = std_data/np.sqrt(num_data)



	
	print("======================Starting figures======================")
	# print(data.shape)
	print("Data has {} nan values".format(np.count_nonzero(np.isnan(avg_data))))
	

	# styles = ['-','--','-.',':']
	# # styles = ['-','--','-.','--.']
	# colors = ['g','b','r','orange','black','red']
	length = len(temps)


	#	plt.close("all")
	plt.rcParams.update({
	    'font.size': 18,
	    'text.usetex': True,
	    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	})

	#Plot metric vs M for all metrics and all N and temps
	for h_i,header in enumerate(gd.data_headers):
		for n_i,n in enumerate(Nums):
			for t_i,t in enumerate(temps):

				fig,ax = plt.subplots()


				ax.errorbar(M,avg_data[h_i,:,n_i,t_i],yerr=err_data[h_i,:,n_i,t_i],\
						label=f"N={n},T={t}",color=colors[h_i],\
						linestyle=styles[h_i],marker='.',markersize=10,zorder=5)

				if include_totals:
					for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
						ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

				bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
				ax.set_xlabel('Fragment size')
				ax.set_ylabel(header)
				# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
				# ax.set_xscale('log')
				# if i == 1:
				fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
				plt.tight_layout()
				if save_plots:
					plt.savefig("{}{}_{}_avgPlot.png".format(figure_folder,dataset_name,header))
				if show_plots:
					plt.show() 



def gen_BPCA_vs_time_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	data_prefolder = path + 'jobsNovus/const_'
	data_prefolder = path + 'jobsCosine/lognorm_'

	dataset_name = data_prefolder.split("/")[-1].strip("_")
	figure_folder = path+f'data/figures/BPCA_per_step/{dataset_name}/'

	if save_plots and not os.path.exists(figure_folder):
		os.makedirs(figure_folder)


	temp = 1000
	# temps = [3,10]
	n = 300

	attempts = [i for i in range(30)]
	
	sizes=list(range(30,301))


	# requested_data_headers = gd.data_headers[:2]
	requested_data_headers = gd.data_headers[:2] + [gd.data_headers[-1]]

	


	for a_i,a in enumerate(attempts):
		raw_data = np.full(shape=(len(requested_data_headers),len(sizes)),fill_value=np.nan,dtype=np.float64)
		folder = f"{data_prefolder}{a}/N_{n}/T_{temp}/"
		print(folder)
		for s_i,size in enumerate(sizes):
			if os.path.exists(folder+"job_data.csv"):
				with open(folder+"job_data.csv",'r') as fp:
					existing_data = fp.readlines()

				existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
				#even though the data can have other sizes in it, 
				#we only want the data of size n
				if size not in existing_sizes:
					print(f"ERROR: Data of size {n} does not exist for {folder}.")
					continue
				index = existing_sizes.index(size)*4
				existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
				existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
				
				for h_i,header in enumerate(gd.data_headers):
					if header in existing_headers_for_size:
						raw_data[h_i,s_i] = existing_values_for_size[existing_headers_for_size.index(header)]




	
		print("======================Starting figures======================")
		# print(data.shape)
		print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
		



		#	plt.close("all")
		plt.rcParams.update({
		    'font.size': 18,
		    'text.usetex': True,
		    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
		})

		#Plot metric vs M for all metrics and all N and temps
		for h_i,header in enumerate(requested_data_headers):

			fig,ax = plt.subplots()

			ax.plot(sizes,raw_data[h_i,:],\
					# label=f"{header} N={n},T={temp}",\
					color=colors[h_i],\
					linestyle=styles[h_i],\
					marker='.',markersize=10,zorder=5)

			# if include_totals:
			# 	for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
			# 		ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

			bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
			ax.set_xlabel('aggregate size (number of particles)')



			ax.set_ylabel(label_from_header(header))
			# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
			# ax.set_xscale('log')
			# if i == 1:
			# fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
			plt.tight_layout()
			if save_plots:
				plt.savefig("{}{}_{}_a-{}_t-{}_overtime.png".format(figure_folder,dataset_name,header,a,temp))
			if show_plots:
				plt.show() 
			plt.close()


def gen_BPCA_vs_time_avg_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	# data_prefolder = path + 'jobsCosine/lognorm_'

	# dataset_name = data_prefolder.split("/")[-1]
	# figure_folder = path+'data/figures/BPCA_per_step/'

	data_prefolder = path + 'jobsNovus/const_'
	data_prefolder = path + 'jobsCosine/lognorm_'

	dataset_name = data_prefolder.split("/")[-1].strip("_")
	figure_folder = path+f'data/figures/BPCA_per_step/{dataset_name}/'

	if save_plots and not os.path.exists(figure_folder):
		os.makedirs(figure_folder)


	temps = [1000]
	# temps = [3,10]
	n = 300

	attempts = [i for i in range(30)]
	
	sizes=list(range(30,301))


	# requested_data_headers = gd.data_headers[:2] + [gd.data_headers[-1]]
	requested_data_headers = [gd.data_headers[1]]

	


	raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(temps),len(sizes)),fill_value=np.nan,dtype=np.float64)
	print(raw_data.shape)
	for a_i,a in enumerate(attempts):
		for t_i,temp in enumerate(temps):
			folder = f"{data_prefolder}{a}/N_{n}/T_{temp}/"
			# print(folder)
			for s_i,size in enumerate(sizes):
				if os.path.exists(folder+"job_data.csv"):
					with open(folder+"job_data.csv",'r') as fp:
						existing_data = fp.readlines()

					existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
					#even though the data can have other sizes in it, 
					#we only want the data of size n
					if size not in existing_sizes:
						print(f"ERROR: Data of size {n} does not exist for {folder}.")
						continue
					index = existing_sizes.index(size)*4
					existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
					existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
					
					for h_i,header in enumerate(requested_data_headers):
						if header in existing_headers_for_size:
							raw_data[h_i,a_i,t_i,s_i] = existing_values_for_size[existing_headers_for_size.index(header)]
				else:
					pass
					# print(f"file doesn't exist: {folder}job_data.csv")

	avg_data = np.nanmean(raw_data,axis=1)
	std_data = np.nanstd(raw_data,axis=1)
	num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
	err_data = std_data/np.sqrt(num_data)

	print(avg_data.shape)
	
	print("======================Starting figures======================")
	# print(data.shape)
	print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
	



	#	plt.close("all")
	plt.rcParams.update({
	    'font.size': 18,
	    'text.usetex': True,
	    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	})

	for h_i,header in enumerate(requested_data_headers):
		#Plot metric vs M for all metrics and all N and temps
		# for a_i,attempt in enumerate(attempts): 
		fig,ax = plt.subplots()
		for t_i,temp in enumerate(temps):


			ax.errorbar(sizes,avg_data[h_i,t_i,:],yerr=err_data[h_i],\
					label=f"{header} N={n},T={temp}",color=colors[t_i],\
					linestyle=styles[h_i],marker='.',markersize=10,zorder=5)


			r0 = 1e-5
			sizes=list(range(30,301))
			Vtot = sizes[0]*(4*np.pi/3)*(r0)**3
			reff = (3*Vtot/(4*np.pi))**(1/3)

			data = Tanaka(sizes,(reff/(1-avg_data[h_i,t_i,0])**(1/3))/(5/3)**(1/2),temp)

			ax.plot(sizes,data,\
					label=f"Tanaka prediction T={temp}",color=colors[t_i+1],\
					linestyle=styles[h_i],marker='*',markersize=10,zorder=5)

		# if include_totals:
		# 	for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
		# 		ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

		bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		ax.set_xlabel('aggregate size')
		ax.set_ylabel(header)
		# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
		# ax.set_xscale('log')
		# if i == 1:
		fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
		plt.tight_layout()
		if save_plots:
			plt.savefig("{}{}_{}_avgovertime.png".format(figure_folder,dataset_name,header))
		if show_plots:
			plt.show() 

	# for h_i,header in enumerate(requested_data_headers):
	# 	#Plot metric vs M for all metrics and all N and temps
	# 	fig,ax = plt.subplots()
	# 	for t_i,temp in enumerate(temps):


	# 		ax.errorbar(sizes,avg_data[h_i,t_i,:],yerr=err_data[h_i],\
	# 				label=f"{header} N={n},T={temp}",color=colors[t_i],\
	# 				linestyle=styles[h_i],marker='.',markersize=10,zorder=5)


	# 	# if include_totals:
	# 	# 	for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
	# 	# 		ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

	# 	bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
	# 	ax.set_xlabel('aggregate size')
	# 	ax.set_ylabel(header)
	# 	# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
	# 	# ax.set_xscale('log')
	# 	# if i == 1:
	# 	fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
	# 	plt.tight_layout()
	# 	if save_plots:
	# 		plt.savefig("{}{}_{}_avgovertime.png".format(figure_folder,dataset_name,header))
	# 	if show_plots:
	# 		plt.show()


def gen_seqstick_plots(distribution):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	if distribution == "lognormal":
		data_prefolder = path + 'jobs/SeqStickLognorm_'
		data_prefolder = path + 'jobs/SeqStickLognormrelax_'
	elif distribution == "constant":
		data_prefolder = path + 'jobs/SeqStickConst_'
		data_prefolder = path + 'jobs/SeqStickConstrelax_'
	else:
		print("Distribution not recognized")
		exit(-1)

	dataset_name = data_prefolder.split("/")[-1]
	figure_folder = path+'data/figures/'



	n = 300

	attempts = [i for i in range(30)]
	
	size=300

	relax = ("relax" in data_prefolder)
	print(f"relax: {relax}")
	rel = ""
	if relax:
		rel = "relax_"


	requested_data_headers = gd.data_headers[:2] + [gd.data_headers[-1]]


	raw_data = np.full(shape=(len(requested_data_headers),len(attempts)),fill_value=np.nan,dtype=np.float64)
	print(f"raw_data shape: {raw_data.shape}")

	for a_i,attempt in enumerate(attempts):
		folder = f"{data_prefolder}{attempt}/N_{n}/"
		# print(folder)
		if os.path.exists(folder+f"{rel}job_data.csv"):
			with open(folder+f"{rel}job_data.csv",'r') as fp:
				existing_data = fp.readlines()

			existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
			#even though the data can have other sizes in it, 
			#we only want the data of size n
			if size not in existing_sizes:
				print(f"ERROR: Data of size {n} does not exist for {folder}.")
				continue
			index = existing_sizes.index(size)*4
			existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
			existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")

			# print(existing_values_for_size)
			
			for h_i,header in enumerate(requested_data_headers):
				if header in existing_headers_for_size:
					raw_data[h_i,a_i] = existing_values_for_size[existing_headers_for_size.index(header)]
				else:
					print(f"Header {header} doesnt exist for dir {folder}")
		else:
			print(f"NO DATA FILE FOR FOLDER: {folder}")


	print(raw_data)
	avg_data = np.nanmean(raw_data,axis=1)
	std_data = np.nanstd(raw_data,axis=1)
	num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
	err_data = std_data/np.sqrt(num_data)

	print(f"avg_data shape: {avg_data.shape}")

	for h_i,header in enumerate(requested_data_headers):
		print(f"{distribution} {header}: {avg_data[h_i]} +- {err_data[h_i]} for {num_data[h_i]} data points.")


def gen_BPCA_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	data_prefolders = []
	data_prefolders.append(path + 'jobsCosine/lognormrelax_')
	data_prefolders.append(path + 'jobsNovus/constrelax_')

	for data_prefolder in data_prefolders:
		dataset_name = data_prefolder.split("/")[-1].strip("_")
		figure_folder = path+f'data/figures/{dataset_name}/'

		if save_plots and not os.path.exists(figure_folder):
			os.makedirs(figure_folder)


		temps = [3,10,30,100,300,1000]
		N = [30,100,300]
		attempts = [i for i in range(30)]

		
		


		# data_file = "test_job_data.csv" #without centering #mean mass
		# data_file = "test_maxnc_job_data.csv" #max nc
		# data_file = "DELETE_job_data.csv" #with centering 
		data_file = "job_data.csv" #with centering 
		# data_file = "nonrelax_job_data.csv" 


		data_file = "job_data.csv"
		# data_file = "nonrelax_job_data.csv" #This nonrelax data follows the Df figure in paper



		bool_headers = [1,1,1,1]
		# requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
		requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

		relax = not ("nonrelax" in data_file)
		print(f"relax: {relax}")
		rel = ""
		if relax:
			rel = "relax_"

		raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
		for a_i,a in enumerate(attempts):
			for n_i,n in enumerate(N):
				size = n
				for t_i,t in enumerate(temps):
					folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
					full_path_data_file = folder+f"{rel}{data_file}"
					if os.path.exists(full_path_data_file):
						with open(full_path_data_file,'r') as fp:
							existing_data = fp.readlines()

						existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
						#even though the data can have other sizes in it, 
						#we only want the data of size n
						if size not in existing_sizes:
							print(f"ERROR: Data of size {n} does not exist for {folder}.")
							continue
						index = existing_sizes.index(size)*4
						existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
						existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
						
						for h_i,header in enumerate(requested_data_headers):
							if header in existing_headers_for_size:
								raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]
					else:
						print(f"DNE: {full_path_data_file}")

		avg_data = np.nanmean(raw_data,axis=1)
		std_data = np.nanstd(raw_data,axis=1)
		num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
		err_data = std_data/np.sqrt(num_data)

		print(avg_data)


		# print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
		# print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

		
		print("======================Starting figures======================")
		# print(data.shape)
		for h_i,header in enumerate(requested_data_headers):
			print(f"Header {header} has {np.count_nonzero(np.isnan(raw_data[h_i]))} nan values")
		

		# styles = ['-','--','-.',':']
		# # styles = ['-','--','-.','--.']
		# colors = ['g','b','r','orange','black','red']


		#	plt.close("all")
		plt.rcParams.update({
		    'font.size': 18,
		    'text.usetex': True,
		    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
		})

		#Plot metric vs M for all metrics and all N and temps
		for h_i,header in enumerate(requested_data_headers):

			fig,ax = plt.subplots()

			for n_i,n in enumerate(N):
				# print(avg_data[h_i,n_i,:])
				ax.errorbar(temps,avg_data[h_i,n_i,:],yerr=err_data[h_i,n_i,:],\
						label=f"N={n}",\
						color=colors[h_i],\
						linestyle=styles[n_i],\
						marker='.',markersize=10,zorder=5)

				if include_totals:
					for txt_i, txt in enumerate(num_data[h_i,n_i,:]):
						ax.annotate("{:0.0f}".format(txt), (temps[txt_i], avg_data[h_i,n_i,txt_i]))

			bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
			ax.set_xlabel('Temperature in K')



			ax.set_ylabel(label_from_header(header))
			# ax.set_title(f'{dataset_name} relax: {relax}')
			ax.set_xscale('log')
			if header == requested_data_headers[1]:
				fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
			plt.tight_layout()
			if save_plots:
				plt.savefig("{}{}_{}_overtemp.png".format(figure_folder,dataset_name,header))
			if show_plots:
				plt.show() 
			plt.close()


def gen_BPCA_ratio_bugbetter_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]


	temps = [3,10,30,100,300,1000]
	N = [30,100,300]
	attempts = [i for i in range(30)]

	
	data_files = []
	data_files.append("job_data.csv")
	data_files.append("test_job_data.csv") #without centering #mean mass



	bool_headers = [1,1,0,0]
	# requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
	requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

	# data_prefolders = []
	# data_prefolders.append(path + 'jobsNovus/constrelax_')
	# data_prefolders.append(path + 'jobsCosine/lognormrelax_')

	data_prefolder = path + 'jobsNovus/constrelax_'
	# data_prefolder = path + 'jobsCosine/lognormrelax_'
	
	avg_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	std_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	num_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	err_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	for d_i,data_file in enumerate(data_files):
		dataset_name = data_prefolder.split("/")[-1].strip("_")
		figure_folder = path+f'data/figures/{dataset_name}/'

		if save_plots and not os.path.exists(figure_folder):
			os.makedirs(figure_folder)

		relax = ("relax" in data_prefolder)

		print(f"relax: {relax}")
		rel = ""
		if relax:
			rel = "relax_"

		raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
		for a_i,a in enumerate(attempts):
			for n_i,n in enumerate(N):
				size = n
				for t_i,t in enumerate(temps):
					folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
					full_data_path = folder+f"{rel}{data_file}"
					if os.path.exists(full_data_path):
						print(f"opening {full_data_path}")
						with open(full_data_path,'r') as fp:
							existing_data = fp.readlines()

						existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
						#even though the data can have other sizes in it, 
						#we only want the data of size n
						if size not in existing_sizes:
							print(f"ERROR: Data of size {n} does not exist for {folder}.")
							continue
						index = existing_sizes.index(size)*4
						existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
						existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
						
						for h_i,header in enumerate(requested_data_headers):
							if header in existing_headers_for_size:
								raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]
					else:
						print(f"DNE: {full_data_path}")


		avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
		std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
		num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
		err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


	ratio_data = avg_data[0]/avg_data[1]
	ratio_errs = ratio_data*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)


	# print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
	# print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

	
	print("======================Starting figures======================")
	# print(data.shape)
	print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
	

	


	#	plt.close("all")
	plt.rcParams.update({
	    'font.size': 18,
	    'text.usetex': True,
	    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	})

	#Plot metric vs M for all metrics and all N and temps
	for h_i,header in enumerate(requested_data_headers):

		fig,ax = plt.subplots()

		for n_i,n in enumerate(N):
			# print(avg_data[h_i,n_i,:])
			# ax.plot(temps,ratio_data[h_i,n_i,:],\
			ax.errorbar(temps,ratio_data[h_i,n_i,:],yerr=ratio_errs[h_i,n_i,:],\
					label=f"N={n}",\
					color=colors[h_i],\
					linestyle=styles[n_i],\
					marker='.',markersize=10,zorder=5)

		# if include_totals:
		# 	for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
		# 		ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

		bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		ax.set_xlabel('Temperature in K')

		ax.axhline(1)

		ax.set_ylabel(label_from_header(header))
		# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
		ax.set_xscale('log')
		if header == requested_data_headers[-1]:
			fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
		plt.tight_layout()
		if save_plots:
			plt.savefig("{}{}_{}_ratiobettertobugovertemp.png".format(figure_folder,dataset_name,header))
		if show_plots:
			plt.show() 
		plt.close()



def gen_BPCA_ratio_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]


	temps = [3,10,30,100,300,1000]
	N = [30,100,300]
	attempts = [i for i in range(30)]

	

	data_file = "nonrelax_job_data.csv" #This nonrelax data follows the Df figure in paper
	data_file = "job_data.csv"



	bool_headers = [1,1,1,1]
	# requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
	requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

	data_prefolders = []
	data_prefolders.append(path + 'jobsNovus/constrelax_')
	data_prefolders.append(path + 'jobsCosine/lognormrelax_')
	
	avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	std_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	num_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	err_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	for d_i,data_prefolder in enumerate(data_prefolders):
		dataset_name = data_prefolder.split("/")[-1].strip("_")
		figure_folder = path+f'data/figures/'

		if save_plots and not os.path.exists(figure_folder):
			os.makedirs(figure_folder)

		relax = not ("nonrelax" in data_file)

		print(f"relax: {relax}")
		rel = ""
		if relax:
			rel = "relax_"

		raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
		for a_i,a in enumerate(attempts):
			for n_i,n in enumerate(N):
				size = n
				for t_i,t in enumerate(temps):
					folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
					full_data_path = folder+f"{rel}{data_file}"
					if os.path.exists(full_data_path):
						print(f"opening {full_data_path}")
						with open(full_data_path,'r') as fp:
							existing_data = fp.readlines()

						existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
						#even though the data can have other sizes in it, 
						#we only want the data of size n
						if size not in existing_sizes:
							print(f"ERROR: Data of size {n} does not exist for {folder}.")
							continue
						index = existing_sizes.index(size)*4
						existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
						existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
						
						for h_i,header in enumerate(requested_data_headers):
							if header in existing_headers_for_size:
								raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]
					else:
						print(f"DNE: {full_data_path}")


		avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
		std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
		num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
		err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


	ratio_data = avg_data[0]/avg_data[1]
	ratio_errs = ratio_data*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)


	# print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
	# print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

	
	print("======================Starting figures======================")
	for h_i,header in enumerate(requested_data_headers):
			print(f"Header {header} has {np.count_nonzero(np.isnan(raw_data[h_i]))} nan values")
	

	


	#	plt.close("all")
	plt.rcParams.update({
	    'font.size': 18,
	    'text.usetex': True,
	    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	})

	#Plot metric vs M for all metrics and all N and temps
	for h_i,header in enumerate(requested_data_headers):

		fig,ax = plt.subplots()

		for n_i,n in enumerate(N):
			# print(avg_data[h_i,n_i,:])
			# ax.plot(temps,ratio_data[h_i,n_i,:],\
			ax.errorbar(temps,ratio_data[h_i,n_i,:],yerr=ratio_errs[h_i,n_i,:],\
					label=f"N={n}",\
					color=colors[h_i],\
					linestyle=styles[n_i],\
					marker='.',markersize=10,zorder=5)

		# if include_totals:
		# 	for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
		# 		ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

		bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		ax.set_xlabel('Temperature in K')

		ax.axhline(1)

		ax.set_ylabel(label_from_header(header))
		# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
		ax.set_xscale('log')
		if header == requested_data_headers[-1]:
			fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
		plt.tight_layout()
		if save_plots:
			plt.savefig("{}{}_ratioovertemp.png".format(figure_folder,header))
		if show_plots:
			plt.show() 
		plt.close()

def gen_BPCA_double_ratio_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]


	temps = [3,10,30,100,300,1000]
	N = [30,100,300]
	attempts = [i for i in range(30)]

	
	data_files = []
	data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
	data_files.append("job_data.csv")



	bool_headers = [1,1,0,1]
	# requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
	requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

	data_prefolders = []
	data_prefolders.append(path + 'jobsNovus/constrelax_')
	data_prefolders.append(path + 'jobsCosine/lognormrelax_')
	
	ratio_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	ratio_errs = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	double_ratio_data = np.full(shape=(len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	for df_i,data_file in enumerate(data_files):
		avg_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
		std_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
		num_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
		err_data = np.full(shape=(len(data_prefolders),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
		for d_i,data_prefolder in enumerate(data_prefolders):
			dataset_name = data_prefolder.split("/")[-1].strip("_")
			figure_folder = path+f'data/figures/{dataset_name}/'

			if save_plots and not os.path.exists(figure_folder):
				os.makedirs(figure_folder)

			relax = False
			if data_file == "job_data.csv":
				relax = True
			print(f"relax: {relax}")
			rel = ""
			if relax:
				rel = "relax_"

			raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
			for a_i,a in enumerate(attempts):
				for n_i,n in enumerate(N):
					size = n
					for t_i,t in enumerate(temps):
						folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
						full_data_path = folder+f"{rel}{data_file}"
						if os.path.exists(full_data_path):
							print(f"opening {full_data_path}")
							with open(full_data_path,'r') as fp:
								existing_data = fp.readlines()

							existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
							#even though the data can have other sizes in it, 
							#we only want the data of size n
							if size not in existing_sizes:
								print(f"ERROR: Data of size {n} does not exist for {folder}.")
								continue
							index = existing_sizes.index(size)*4
							existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
							existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
							
							for h_i,header in enumerate(requested_data_headers):
								if header in existing_headers_for_size:
									raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]
						else:
							print(f"DNE: {full_data_path}")


			avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
			std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
			num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
			err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


		ratio_data[df_i] = avg_data[0]/avg_data[1]
		ratio_errs[df_i] = ratio_data[df_i]*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)


	double_ratio_data = ratio_data[0]/ratio_data[1]
	double_ratio_errs = double_ratio_data*np.sqrt((ratio_errs[0]/ratio_data[0])**2+(ratio_errs[1]/ratio_data[1])**2)
	# print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
	# print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

	
	print("======================Starting figures======================")
	# print(data.shape)
	print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
	


	#	plt.close("all")
	plt.rcParams.update({
	    'font.size': 18,
	    'text.usetex': True,
	    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	})

	#Plot metric vs M for all metrics and all N and temps
	for h_i,header in enumerate(requested_data_headers):

		fig,ax = plt.subplots()

		for n_i,n in enumerate(N):
			# print(avg_data[h_i,n_i,:])

			# ax.plot(temps,double_ratio_data[h_i,n_i,:],\
			ax.errorbar(temps,double_ratio_data[h_i,n_i,:],yerr=double_ratio_errs[h_i,n_i,:],\
					label=f"N={n}",\
					color=colors[h_i],\
					linestyle=styles[n_i],\
					marker='.',markersize=10,zorder=5)

		# if include_totals:
		# 	for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
		# 		ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

		bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		ax.set_xlabel('Temperature in K')

		ax.axhline(1)

		ax.set_ylabel(label_from_header(header))
		# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
		ax.set_xscale('log')
		if header == requested_data_headers[-1]:
			fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
		plt.tight_layout()
		if save_plots:
			plt.savefig("{}{}_{}_doubleratioovertemp.png".format(figure_folder,dataset_name,header))
		if show_plots:
			plt.show() 
		plt.close()



def gen_BPCA_ratio_nonreltorel_vs_temp_plots(show_plots=True,save_plots=False,include_totals=False):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]


	temps = [3,10,30,100,300,1000]
	N = [30,100,300]
	attempts = [i for i in range(30)]

	
	data_files = []
	data_files.append("nonrelax_job_data.csv") #This nonrelax data follows the Df figure in paper
	data_files.append("job_data.csv")


	bool_headers = [1,1,1,1]
	# requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
	requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

	data_prefolder = path + 'jobsCosine/lognormrelax_'
	data_prefolder = path + 'jobsNovus/constrelax_'
	
	avg_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	std_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	num_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	err_data = np.full(shape=(len(data_files),len(requested_data_headers),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
	for d_i,data_file in enumerate(data_files):
		dataset_name = data_prefolder.split("/")[-1].strip("_")
		figure_folder = path+f'data/figures/{dataset_name}/'

		if save_plots and not os.path.exists(figure_folder):
			os.makedirs(figure_folder)

		relax = not ("nonrelax" in data_file)
		# relax = False
		print(f"relax: {relax}")
		rel = ""
		if relax:
			rel = "relax_"

		raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
		for a_i,a in enumerate(attempts):
			for n_i,n in enumerate(N):
				size = n
				for t_i,t in enumerate(temps):
					folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
					full_data_path = folder+f"{rel}{data_file}"
					if os.path.exists(full_data_path):
						print(f"opening {full_data_path}")
						with open(full_data_path,'r') as fp:
							existing_data = fp.readlines()

						existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
						#even though the data can have other sizes in it, 
						#we only want the data of size n
						if size not in existing_sizes:
							print(f"ERROR: Data of size {n} does not exist for {folder}.")
							continue
						index = existing_sizes.index(size)*4
						existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
						existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
						
						for h_i,header in enumerate(requested_data_headers):
							if header in existing_headers_for_size:
								raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]
					else:
						print(f"DNE: {full_data_path}")


		avg_data[d_i,:,:,:] = np.nanmean(raw_data,axis=1)
		std_data[d_i,:,:,:] = np.nanstd(raw_data,axis=1)
		num_data[d_i,:,:,:] = np.count_nonzero(~np.isnan(raw_data),axis=1)
		err_data[d_i,:,:,:] = std_data[d_i,:,:,:]/np.sqrt(num_data[d_i,:,:,:])


	ratio_data = avg_data[0]/avg_data[1]
	ratio_errs = ratio_data*np.sqrt((err_data[0]/avg_data[0])**2+(err_data[1]/avg_data[1])**2)


	# print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
	# print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

	
	print("======================Starting figures======================")
	# print(data.shape)
	print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
	

	


	#	plt.close("all")
	plt.rcParams.update({
	    'font.size': 18,
	    'text.usetex': True,
	    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	})

	#Plot metric vs M for all metrics and all N and temps
	for h_i,header in enumerate(requested_data_headers):

		fig,ax = plt.subplots()

		for n_i,n in enumerate(N):
			# print(avg_data[h_i,n_i,:])
			# ax.plot(temps,ratio_data[h_i,n_i,:],\
			ax.errorbar(temps,ratio_data[h_i,n_i,:],yerr=ratio_errs[h_i,n_i,:],\
					label=f"N={n}",\
					color=colors[h_i],\
					linestyle=styles[n_i],\
					marker='.',markersize=10,zorder=5)

		# if include_totals:
		# 	for txt_i, txt in enumerate(num_data[h_i,:,n_i,t_i]):
		# 		ax.annotate("{:0.0f}".format(txt), (M[txt_i], avg_data[h_i,txt_i,n_i,t_i]))

		bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		ax.set_xlabel('Temperature in K')

		ax.axhline(1)

		ax.set_ylabel(label_from_header(header))
		# ax.set_title('{} {} vs Temp'.format(dataset_name,method))
		ax.set_xscale('log')
		if header == requested_data_headers[-1]:
			fig.legend(loc='upper right',bbox_to_anchor=(0.97, 0.96))
		plt.tight_layout()
		if save_plots:
			plt.savefig("{}{}_{}_rationonreltorelovertemp.png".format(figure_folder,dataset_name,header))
		if show_plots:
			plt.show() 
		plt.close()

def gen_BPCA_temp_sensitivity_plots(show_plots=True,save_plots=False,include_totals=False):
	def S(sigma):
		return np.sum(1/sigma**2)

	def Si(i,sigma):
		return np.sum(i/(sigma**2))

	def Sii(x,y,sigma):
		return np.sum(x*y/(sigma**2))

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]


	temps = [3,10,30,100,300,1000]
	x = np.log(temps)
	N = [30,100,300]
	attempts = [i for i in range(30)]

	
	data_file = "test_RKBMs_job_data.csv" #print both ways of RKBM
	data_file = "test_job_data.csv" #without centering #mean mass
	data_file = "test_maxnc_job_data.csv" #max nc
	data_file = "job_data.csv" #with centering


	bool_headers = [1,1,1,1]
	# requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
	requested_data_headers = [gd.data_headers[i] for i in range(len(gd.data_headers)) if bool_headers[i]]

	requested_data_headers = requested_data_headers[:2] + [requested_data_headers[3]] + [requested_data_headers[2]]

	data_prefolders = []
	data_prefolders.append(path + 'jobsNovus/constrelax_')
	data_prefolders.append(path + 'jobsCosine/lognormrelax_')
	
	for data_prefolder in data_prefolders:
		dataset_name = data_prefolder.split("/")[-1].strip("_")
		figure_folder = path+f'data/figures/{dataset_name}/'

		if save_plots and not os.path.exists(figure_folder):
			os.makedirs(figure_folder)

		relax = not ("nonrelax" in data_file)
		print(f"relax: {relax}")
		rel = ""
		if relax:
			rel = "relax_"

		slope_data = np.full(shape=(4,len(N)),fill_value=0,dtype=np.float64)
		slope_sigma_data = np.full(shape=(4,len(N)),fill_value=0,dtype=np.float64)
		raw_data = np.full(shape=(len(requested_data_headers),len(attempts),len(N),len(temps)),fill_value=np.nan,dtype=np.float64)
		for a_i,a in enumerate(attempts):
			for n_i,n in enumerate(N):
				size = n
				for t_i,t in enumerate(temps):
					folder = f"{data_prefolder}{a}/N_{n}/T_{t}/"
					if os.path.exists(folder+f"{rel}{data_file}"):
						full_data_path = folder+f"{rel}{data_file}"
						print(f"opening {full_data_path}")
						with open(full_data_path,'r') as fp:
							existing_data = fp.readlines()

						existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]
						#even though the data can have other sizes in it, 
						#we only want the data of size n
						if size not in existing_sizes:
							print(f"ERROR: Data of size {n} does not exist for {folder}.")
							continue
						index = existing_sizes.index(size)*4
						existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
						existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")
						
						for h_i,header in enumerate(requested_data_headers):
							if header in existing_headers_for_size:
								raw_data[h_i,a_i,n_i,t_i] = existing_values_for_size[existing_headers_for_size.index(header)]


		avg_data = np.nanmean(raw_data,axis=1)
		std_data = np.nanstd(raw_data,axis=1)
		num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
		err_data = std_data/np.sqrt(num_data)

		for h_i,header in enumerate(requested_data_headers):
			for n_i,n in enumerate(N):
				y = avg_data[h_i,n_i,:]
				sigma = err_data[h_i,n_i,:]

				delta = S(sigma)*Sii(x,x,sigma)-(Si(x,sigma))**2
				slope_data[h_i,n_i] = np.abs((S(sigma)*Sii(x,y,sigma)-Si(x,sigma)*Si(y,sigma))/delta)
				slope_sigma_data[h_i,n_i] = np.sqrt(S(sigma)/delta)
		# print(f"{requested_data_headers[0]}: {avg_data[0,2,0]} +- {err_data[0,2,0]}")
		# print(f"{requested_data_headers[1]}: {avg_data[1,2,0]} +- {err_data[1,2,0]}")

		
		print("======================Starting figures======================")
		# print(data.shape)
		print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
		

		#	plt.close("all")
		plt.rcParams.update({
		    'font.size': 18,
		    'text.usetex': True,
		    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
		})

		fig,ax = plt.subplots(figsize=(10,5))
		dummy_x_data = [0,.2,.4,.6]


		for n_i,n in enumerate(N):
			rang = []
			if n == 30:
				rang = [0,4]
				shift = -0.035
			elif n == 100:
				rang = [4,8]
				shift = 0
			else:
				rang = [8,12]
				shift = 0.035
			x_data = [i+shift for i in dummy_x_data]

			ax.errorbar(x=x_data,y=slope_data[:,n_i],yerr=slope_sigma_data[:,n_i],\
					fmt='o',linewidth=2, capsize=6,label=f"N={n}")


		plt.axvline(x = 0.1, color = 'black')
		plt.axvline(x = 0.3, color = 'black')
		plt.axvline(x = 0.5, color = 'black')


		ticklabels = [label_from_header(i) for i in requested_data_headers]
		ticklabels.extend([""]*(len(dummy_x_data)-len(requested_data_headers)))

		ax.xaxis.set(ticks=dummy_x_data,
				ticklabels=ticklabels)

		plt.hlines(0,xmin=-0.0685,xmax=0.6685,linestyle="--",color='black')

		ax.set_ylabel('Sensitivity to Temperature')
		ax.set_ylim(-0.0010770496039509136, 0.02782241381226227)
		ax.set_xlim(-0.0685, 0.6685)
		
		fig.legend(loc='lower left',bbox_to_anchor=(0.785,0.65))
		plt.tight_layout()
		if save_plots:
			plt.savefig("{}{}_methodComp.png".format(figure_folder,dataset_name))
		if show_plots:
			plt.show() 
		plt.close()



if __name__ == '__main__':
	#Do you want to see plots of the data as they are made?
	show_plots = True
	#Do you want to save the plots once they are made?
	save_plots = True
	#Do you want the number of runs next to each point on the plots
	#so you know how many more runs need to finish
	include_totals = False




	# gen_BAPA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_BPCA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_BPCA_vs_time_avg_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_BPCA_vs_time_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_seqstick_plots(distribution="lognormal")
	# gen_seqstick_plots(distribution="constant")


	# gen_relax_vs_tense_BPCA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_relax_vs_tense_seqstick_plots(distribution="lognormal",show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)

	# gen_BPCA_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_BPCA_ratio_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	gen_BPCA_temp_sensitivity_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_BPCA_ratio_bugbetter_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_BPCA_double_ratio_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_BPCA_ratio_nonreltorel_vs_temp_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)



	# T3 = 0.43736817467052647
	# T1000 = 0.3927614742047937
	# temp = 3
	# r0 = 1e-5

	# sizes=list(range(30,3001))
	# Vtot = sizes[0]*(4*np.pi/3)*(r0)**3
	# reff = (3*Vtot/(4*np.pi))**(1/3)

	# data = Tanaka(sizes,(reff/(1-T3)**(1/3))/(5/3)**(1/2),temp)
	
	# fig,ax = plt.subplots()
	# ax.plot(sizes,data,temp,\
	# 		label=f"Tanaka prediction T={temp}")

	# ax.set_ylim(0.42,0.46)

	
	# plt.show() 
	# plt.close()