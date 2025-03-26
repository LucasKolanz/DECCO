"""
This file was originally written for SpaceLab/DECCO to do data processing

Author: Lucas Kolanz


"""





import sys
import glob
import os
import json
import numpy as np

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u




def get_last_velocity(data_folder,size,relax=False):
	data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	






if __name__ == '__main__':
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	# data_prefolder = path + 'jobsNovus/const_relax'
	# data_prefolder = path + 'jobsCosine/lognormrelax_*'


	size = 300
	data_folders = []

	# data_folders = path + 'jobs/BAPA_*'

	data_folders.append(path + 'jobsCosine/lognorm_*/')
	data_folders.append(path + 'jobsNovus/const_*/')
	# data_folders = [path + f'jobsCosine/lognormrelax_*/N_{size}/*']
	# data_folders = [path + 'jobsNovus/const_*/N_300/T_1000/']
	# data_folders = data_folders + [path + 'jobsNovus/const_*/N_300/T_3/']
	# data_folders = data_folders + [path + f'jobsNovus/constrelax_*/N_{size}/*']
	# data_folders = path + 'jobsNovus/constrelax_*/N_300/T_*/'

	possible_dirs = []
	for data_folder in data_folders:
		possible_dirs.extend(u.get_directores_containing(data_folder,["timing.txt"]))




	#list of the functions that calculate the data you want
	#if adding to this list, name your function calc_*header_name*
	#where *header_name* is what you want this data to be called
	#in the header of the job_data.csv file. The function should 
	#only take in the directory the data is in and the size of which 
	#to calculate as an input.
	#It should return a single data value.
	requested_data_functions = data_functions[:2]# + [data_functions[-1]]
	requested_data_headers = data_headers[:2]# + [data_headers[-1]]


	# requested_data_functions = data_functions[:]
	# requested_data_headers = data_headers[:]


	#list of intermediate sizes to calculate data for.
	# requested_sizes = list(range(30,301))
	requested_sizes = [size]


	# possible_dirs = ['/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_20/N_300/T_1000/']
	for directory in possible_dirs:
		relax = ("relax" in directory)
		print(f"relax: {relax}")
		rel = ""
		if relax:
			rel = "relax_"
		print(f"Generating data for: {directory}")
		if os.path.exists(directory+"timing.txt"):
			existing_data = []
			if os.path.exists(directory+f"{rel}job_data.csv"):
				with open(directory+f"{rel}job_data.csv",'r') as fp:
					existing_data = fp.readlines()
			existing_sizes = [int(i.split('=')[1].strip("\n\t ")) for i in existing_data if i[:2] == "N="]

			#This will be all the data we want to write in the end
			lines = []
			
			#contains both the sizes we want and the sizes we already have
			all_sizes = sorted(set(existing_sizes+requested_sizes))

			for size in all_sizes:

				existing_headers_for_size = []
				existing_values_for_size = []
				if size in existing_sizes:
					#calculate the index of the existing size
					index = existing_sizes.index(size)*4
					existing_headers_for_size = existing_data[index+1].strip("\n\t ").split(",")
					existing_values_for_size = existing_data[index+2].strip("\n\t ").split(",")

				#If this size wasn't requested then we don't want to overwrite it
				if size in requested_sizes:
					headers,values = calc_from_size(size,directory,existing_headers_for_size,existing_values_for_size,requested_data_headers,requested_data_functions,relax)
				else:
					headers = existing_headers_for_size
					values = existing_values_for_size


				lines.append(f"N={size}\n")
				lines.append(','.join(headers)+'\n')
				lines.append(','.join(values)+'\n')
				lines.append('\n')


			with open(directory+f"{rel}job_data.csv",'w') as fp:
				fp.writelines(lines)

		else:
			print(f"Job is not finished in {directory}")


