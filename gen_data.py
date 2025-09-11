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





def calc_porosity_abc(data_folder,size,relax=False):
	data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	if data is None:
		return np.nan
	# num_balls = data.shape[0]

	vol = [(4*np.pi/3)*r**3 for r in radius]

	effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 

		
	# principal_moi = u.get_principal_moi(vol,data)
	principal_moi = u.get_principal_moi(mass,data)
	
	
	alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)
	# alphai = principal_moi/(0.4*np.sum(vol)*effective_radius**2)
	
	a = effective_radius * np.sqrt(alphai[1] + alphai[2] - alphai[0])
	b = effective_radius * np.sqrt(alphai[2] + alphai[0] - alphai[1])
	c = effective_radius * np.sqrt(alphai[0] + alphai[1] - alphai[2])
	
	# Rabc = np.power(a*b*c,1/3)
	porosity = 1-(effective_radius**3/(a*b*c))

	return porosity

def calc_porosity_KBM(data_folder,size,relax=False):
	data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	if data is None:
		return np.nan
	# num_balls = data.shape[0]

	effective_radius = np.power(np.sum(np.power(radius,3)),1/3)   
		
	principal_moi = u.get_principal_moi(mass,data)
	# principal_moi = u.get_principal_moi(np.mean(mass),data)

	alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)

	RKBM = np.sqrt(np.sum(alphai)/3) * effective_radius

	# total_mass = np.sum(mass)
	# cofm = np.einsum('i,ij->j', mass, data) / total_mass

	# Rgyr = np.sqrt(np.sum(mass*(np.linalg.norm(data,axis=1)-np.linalg.norm(cofm))**2)/total_mass)

	porosity = 1-np.power((effective_radius/RKBM),3)
	# print(f"Rgyr way: {1-np.power((effective_radius/((5/3)**(1/2)*Rgyr)),3)}")
	# print(f"other way: {porosity}")

	return porosity

def calc_number_of_contacts(data_folder,data_index=-1,relax=False):
	print(data_folder)
	line = 0
	max_nc = -1
	nc = number_of_contacts(data_folder,data_index,line,relax)
	while not np.isnan(nc):
		max_nc = max(max_nc,nc)
		line += 1 
		nc = number_of_contacts(data_folder,data_index,line,relax)
	return max_nc

def number_of_contacts(data_folder,size,line,relax=False):
	data,radius,mass,moi = u.get_data(data_folder,data_index=size,linenum=line,relax=relax)

	if data is None:
		return np.nan
	data = np.array(data)
	num_balls = data.shape[0]

	contacts = np.zeros((num_balls,num_balls),dtype=int)
	dist = lambda i,j: np.sqrt((data[i][0]-data[j][0])**2 + (data[i][1]-data[j][1])**2 + \
			(data[i][2]-data[j][2])**2)

	##IT FEELS LIKE THIS IS OVERCOUNTING BUT IT ISN'T.
	##EACH BALL FELLS ITS OWN CONTACTS.
	for i in range(num_balls):
		for j in range(num_balls):
			if i != j:
				contacts[i,j] = (dist(i,j) <= (radius[i]+radius[j]))

	return np.mean(np.sum(contacts,axis=1))

def calc_fractal_dimension(data_folder,size,relax=False):
	overwrite_octree_data = True
	show_FD_plots = False
	o3dv = u.o3doctree(data_folder,overwrite_data=overwrite_octree_data,index=size,Temp=-1,relax=relax)
	o3dv.make_tree()
	FD = o3dv.calc_fractal_dimension(show_graph=show_FD_plots)
	return FD



data_functions = [calc_porosity_abc,calc_porosity_KBM,calc_number_of_contacts,calc_fractal_dimension]
data_headers = [i.__name__[5:] for i in data_functions]

def calc_from_size(size,directory,existing_headers,existing_values,requested_headers,relax=False,overwrite=False):
	headers = []
	values = []

	for h,header in enumerate(data_headers):
		if not overwrite:
			#we already have this header include it
			if header in existing_headers:
				index = existing_headers.index(header)
				headers.append(header)
				values.append(existing_values[index])
			#we dont have this header and we want it, so calculate
			elif (header not in existing_headers and header in requested_headers):
				headers.append(header)
				values.append(str(data_functions[h](directory,size,relax)))
		else:
			#we already have this header and didn't ask for it, so just include it
			if header in existing_headers and header not in requested_headers:
				index = existing_headers.index(header)
				headers.append(header)
				values.append(existing_values[index])
			#if overwrite, we want to calculate reguardless, if header is requested
			elif header in requested_headers:
				headers.append(header)
				values.append(str(data_functions[h](directory,size,relax)))

	return headers,values


if __name__ == '__main__':
	# print(calc_porosity_abc("/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_0/N_300/T_1000/",300,True))
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]

	# data_prefolder = path + 'jobsNovus/const_relax'
	# data_prefolder = path + 'jobsCosine/lognormrelax_*'

	# data_file = "test_job_data.csv" #without centering #mean mass
	# data_file = "test_RKBMs_job_data.csv" #print both ways of RKBM
	# data_file = "test_maxnc_job_data.csv" #max nc
	# data_file = "DELETE_job_data.csv" #with centering 
	# data_file = "nonrelax_job_data.csv" #This nonrelax data follows the Df figure in paper
	data_file = "job_data.csv" #with centering




	N = [300]

	#list of the functions that calculate the data you want
	#if adding to this list, name your function calc_*header_name*
	#where *header_name* is what you want this data to be called
	#in the header of the job_data.csv file. The function should 
	#only take in the directory the data is in and the size of which 
	#to calculate as an input.
	#It should return a single data value.
	bool_headers = [1,1,0,0]
	# requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
	requested_data_headers = [data_headers[i] for i in range(len(data_headers)) if bool_headers[i]]

	overwritedata = False

	for n_i,n in enumerate(N):
	
		data_folders = []
		data_folders = [path + 'jobs/BAPA_*']
		# data_folders = [path + 'jobs/BAPA_0/M_60/*']
		# data_folders = [path + 'jobs/SeqStickConst_*/']
		# data_folders = [path + 'jobsCosine/lognorm_*/N_300/T_*/']
		# data_folders = [path + 'jobsCosine/lognorm_*/N_300/T_*/']
		# data_folders = [path + 'jobsNovus/const_*/N_300/T_1000/']
		# data_folders = data_folders + [path + 'jobsNovus/const_*/N_300/T_3/']
		# data_folders.append(path + f'jobsNovus/constrelax_*/N_{n}/*')
		# data_folders = [path + 'jobsNovus/constrelax_*/N_300/T_*/']

		# data_folders.append(path + 'jobs/SeqStickLognormrelax*/')
		# data_folders.append(path + 'jobs/SeqStickConstrelax*/')
		
		# data_folders.append(path + f'jobsNovus/constrelax_*/N_{n}/*')
		# data_folders.append(path + f'jobsCosine/lognormrelax_*/N_{n}/*')

		possible_dirs = []
		for data_folder in data_folders:
			possible_dirs.extend(u.get_directores_containing(data_folder,["timing.txt"]))


		#list of intermediate sizes to calculate data for.

		# requested_sizes = list(range(30,301))
		requested_sizes = [n]


		for directory in possible_dirs:
			relax = ("relax" in directory)
			# relax = False
			# print(f"relax: {relax}")
			rel = ""
			if relax:
				rel = "relax_"
			# directory = "/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsNovus/constrelax_10/N_30/T_1000/"
			print(f"Generating data for: {directory}{rel}{data_file}")

			if os.path.exists(directory+"timing.txt"):
				existing_data = []
				if os.path.exists(directory+f"{rel}{data_file}"):
					with open(directory+f"{rel}{data_file}",'r') as fp:
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
					#Or if overwrite is true then overwrite anyway
					if overwritedata or size in requested_sizes:
						headers,values = calc_from_size(size,directory,existing_headers_for_size,existing_values_for_size,requested_data_headers,relax,overwritedata)
						print(headers)
						print(values)
					else:
						headers = existing_headers_for_size
						values = existing_values_for_size


					lines.append(f"N={size}\n")
					lines.append(','.join(headers)+'\n')
					lines.append(','.join(values)+'\n')
					lines.append('\n')


				with open(directory+f"{rel}{data_file}",'w') as fp:
					fp.writelines(lines)

			else:
				print(f"Job is not finished in {directory}")
			



