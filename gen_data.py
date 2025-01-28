"""
This file was originally written for SpaceLab/DECCO to do data processing

Author: Lucas Kolanz

This file goes through all folders in a specified job set and calculates the specified data for this run
for the folders which have completed jobs in them.

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




def calc_porosity_abc(data_folder,size):
	data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=False)
	if data is None:
		return np.nan
	# num_balls = data.shape[0]

	effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 


	# effective_radius = radius*np.power(num_balls,1/3) 
		
	principal_moi = u.get_principal_moi(np.mean(mass),data)
	# principal_moi = get_principal_moi(mass,data)
	
	
	alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)
	# alphai = principal_moi/(0.4*num_balls*mass*effective_radius**2)
	
	a = effective_radius * np.sqrt(alphai[1] + alphai[2] - alphai[0])
	b = effective_radius * np.sqrt(alphai[2] + alphai[0] - alphai[1])
	c = effective_radius * np.sqrt(alphai[0] + alphai[1] - alphai[2])
	
	# Rabc = np.power(a*b*c,1/3)
	porosity = 1-(effective_radius**3/(a*b*c))


	return porosity

def calc_porosity_KBM(data_folder,size):
	data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=False)
	if data is None:
		return np.nan
	# num_balls = data.shape[0]

	effective_radius = np.power(np.sum(np.power(radius,3)),1/3)  
	# effective_radius = radius*np.power(num_balls,1/3) 
		
	principal_moi = u.get_principal_moi(np.mean(mass),data)
	# principal_moi = get_principal_moi(mass,data)

	alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)
	# alphai = principal_moi/(0.4*num_balls*mass*effective_radius**2)

	RKBM = np.sqrt(np.sum(alphai)/3) * effective_radius


	porosity = 1-np.power((effective_radius/RKBM),3)
	return porosity

def calc_number_of_contacts(data_folder,size):
	data,radius,mass,moi = u.get_data(data_folder,data_index=size,linenum=-1,relax=False)

	if data is None:
		return np.nan
	data = np.array(data)
	num_balls = data.shape[0]

	contacts = np.zeros((num_balls,num_balls),dtype=int)
	dist = lambda i,j: np.sqrt((data[i][0]-data[j][0])**2 + (data[i][1]-data[j][1])**2 + \
			(data[i][2]-data[j][2])**2)

	for i in range(num_balls):
		for j in range(num_balls):
			if i != j:
				contacts[i,j] = (dist(i,j) <= (radius[i]+radius[j]))

	return np.mean(np.sum(contacts,axis=1))

def calc_fractal_dimension(data_folder,size):
	overwrite_octree_data = False
	show_FD_plots = False
	o3dv = u.o3doctree(data_folder,overwrite_data=overwrite_octree_data,index=size,Temp=-1,relax=False)
	o3dv.make_tree()
	return o3dv.calc_fractal_dimension(show_graph=show_FD_plots)





def calc_from_size(size,directory,existing_headers,existing_values,requested_headers,requested_functions):
	headers = []
	values = []

	for h,req_header in enumerate(requested_headers):
		if req_header in existing_headers:
			index = existing_headers.index(req_header)
			headers.append(req_header)
			values.append(existing_values[index])
		else:
			headers.append(req_header)
			values.append(str(requested_functions[h](directory,size)))

	return headers,values

data_functions = [calc_porosity_abc,calc_porosity_KBM,calc_number_of_contacts,calc_fractal_dimension]
data_headers = [i.__name__[5:] for i in data_functions]

if __name__ == '__main__':
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	data_prefolder = path + 'jobsNovus/const_relax'
	data_prefolder = path + 'jobsCosine/lognorm_relax'
	data_prefolder = path + 'jobs/BAPA_*'

	possible_dirs = u.get_directores_containing(data_prefolder,["timing.txt"])


	#list of the functions that calculate the data you want
	#if adding to this list, name your function calc_*header_name*
	#where *header_name* is what you want this data to be called
	#in the header of the job_data.csv file. The function should 
	#only take in the directory the data is in and the size of which 
	#to calculate as an input.
	#It should return a single data value.
	requested_data_functions = data_functions
	requested_data_headers = data_headers


	#list of intermediate sizes to calculate data for.
	requested_sizes = [300]

	# possible_dirs = ['/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobs/BAPA_0/M_20/N_300/T_1000/']
	for directory in possible_dirs:

		if os.path.exists(directory+"timing.txt"):
			existing_data = []
			if os.path.exists(directory+"job_data.csv"):
				with open(directory+"job_data.csv",'r') as fp:
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
					headers,values = calc_from_size(size,directory,existing_headers_for_size,existing_values_for_size,requested_data_headers,requested_data_functions)
				else:
					headers = existing_headers_for_size
					values = existing_values_for_size


				lines.append(f"N={size}\n")
				lines.append(','.join(headers)+'\n')
				lines.append(','.join(values)+'\n')
				lines.append('\n')


			with open(directory+"job_data.csv",'w') as fp:
				fp.writelines(lines)

		else:
			print(f"Job is not finished in {directory}")

