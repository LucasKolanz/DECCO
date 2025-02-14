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
	

	styles = ['-','--','-.',':']
	# styles = ['-','--','-.','--.']
	colors = ['g','b','r','orange','black','red']
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

	data_prefolder = path + 'jobsCosine/lognorm_'

	dataset_name = data_prefolder.split("/")[-1]
	figure_folder = path+'data/figures/BPCA_per_step/'


	temp = 3
	# temps = [3,10]
	n = 300

	attempts = [i for i in range(30)]
	
	sizes=list(range(30,301))


	requested_data_headers = gd.data_headers[:2]
	


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
		

		styles = ['-','--','-.',':']
		# styles = ['-','--','-.','--.']
		colors = ['g','b','r','orange','black','red']


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
					label=f"{header} N={n},T={temp}",color=colors[h_i],\
					linestyle=styles[h_i],marker='.',markersize=10,zorder=5)

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
				plt.savefig("{}{}_{}_a-{}_t-{}.png".format(figure_folder,dataset_name,header,a,temp))
			if show_plots:
				plt.show() 


def gen_seqstick_plots(distribution):
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	if distribution == "lognormal":
		data_prefolder = path + 'jobs/SeqStickLognorm_'
	elif distribution == "constant":
		data_prefolder = path + 'jobs/SeqStickConst_'
	else:
		print("Distribution not recognized")
		exit(-1)

	dataset_name = data_prefolder.split("/")[-1]
	figure_folder = path+'data/figures/'



	n = 300

	attempts = [i for i in range(30)]
	
	size=300


	requested_data_headers = gd.data_headers[:2]
	


	raw_data = np.full(shape=(len(requested_data_headers),len(attempts)),fill_value=np.nan,dtype=np.float64)
	for a_i,attempt in enumerate(attempts):
		folder = f"{data_prefolder}{attempt}/N_{n}/"
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
					raw_data[h_i,a_i] = existing_values_for_size[existing_headers_for_size.index(header)]



	avg_data = np.nanmean(raw_data,axis=1)
	std_data = np.nanstd(raw_data,axis=1)
	num_data = np.count_nonzero(~np.isnan(raw_data),axis=1)
	err_data = std_data/np.sqrt(num_data)


	for h_i,header in enumerate(requested_data_headers):
		print(f"{distribution} {header}: {avg_data[h_i]} +- {err_data[h_i]}")

	# print("======================Starting figures======================")
	# # print(data.shape)
	# print("Data has {} nan values".format(np.count_nonzero(np.isnan(raw_data))))
	

	# styles = ['-','--','-.',':']
	# # styles = ['-','--','-.','--.']
	# colors = ['g','b','r','orange','black','red']


	# #	plt.close("all")
	# plt.rcParams.update({
	#     'font.size': 18,
	#     'text.usetex': True,
	#     'text.latex.preamble': r'\usepackage{amsmath} \usepackage{bm}'
	# })

	# #Plot metric vs M for all metrics and all N and temps
	# for h_i,header in enumerate(requested_data_headers):

	# 	fig,ax = plt.subplots()

	# 	ax.errorbar([size],[avg_data[h_i]],yerr=[err_data[h_i]],\
	# 			label=f"{header} N={n}",color=colors[h_i],\
	# 			linestyle=styles[h_i],marker='.',markersize=10,zorder=5)

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
	# 		plt.savefig("{}{}_{}_seqStickConst.png".format(figure_folder,dataset_name,header))
	# 	if show_plots:
	# 		plt.show() 


if __name__ == '__main__':
	#Do you want to see plots of the data as they are made?
	show_plots = True
	#Do you want to save the plots once they are made?
	save_plots = True
	#Do you want the number of runs next to each point on the plots
	#so you know how many more runs need to finish
	include_totals = True




	# gen_BAPA_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	# gen_BPCA_vs_time_plots(show_plots=show_plots,save_plots=save_plots,include_totals=include_totals)
	gen_seqstick_plots(distribution="lognormal")
	gen_seqstick_plots(distribution="constant")