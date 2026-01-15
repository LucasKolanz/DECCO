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
import csv
import numpy as np

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

sys.path.append(project_path+"utilities/")
# sys.path.append("/home/kolanzl/Desktop/SpaceLab/")
import utils as u

# Path for the single global log file
GLOBAL_LOG = f"{project_path}../SpaceLab_data/logs/"

#from Suyama 2012
#Just geometric cross section, no normalization or nothing. 
def calc_geometric_cross_section(data_folder,size,relax=False,**kwargs):
	makeVisual = kwargs.get("makeVisual", False)

	pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	if pos is None:
		return np.nan

	if makeVisual:
		orientations = 1
	else:
		orientations = 30
	
	gcs = 0
	direction = (0,0,1)
	for _ in range(orientations):
		if not makeVisual:
			direction = u.make_rand_vector(3)
			gcs += u.calc_geometric_cross_section(pos,radius,direction=direction,write=makeVisual,directory=data_folder)
		else:
			gcs += u.calc_geometric_cross_section(pos,radius,direction=direction,mesh_factor=1.0,write=makeVisual,directory=data_folder)
	gcs /= orientations
	
	# normalization = np.pi*np.sum(np.power(radius,2))

	if makeVisual:
		R = np.sqrt(gcs/np.pi)
		COM = u.calcCOM(mass,pos)
		u.writeGCSRadius(R,COM,data_folder)

	return gcs


def calc_porosity_ch(data_folder,size,relax=False,**kwargs):
	makeVisual = kwargs.get("makeVisual", False)

	pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	if pos is None:
		return np.nan
	
	effective_radius_cubed = np.sum(np.power(radius,3))

	hull_vol = u.calc_hull_volume(pos,radius,write=makeVisual,directory=data_folder)

	#if the MVEE failed, return NAN
	if np.isnan(hull_vol):
		return np.nan

	porosity = 1-((4.0/3.0)*np.pi*effective_radius_cubed)/(hull_vol)

	return porosity

def calc_porosity_fee(data_folder, size, relax=False,**kwargs):
	makeVisual = kwargs.get("makeVisual", False)

	# Load particle data
	pos, radius, mass, moi = u.get_data(data_folder, data_index=size, relax=relax)
	if pos is None:
		return np.nan

	# Compute reference volume
	effective_radius_cubed = np.sum(radius**3)

	# Compute ellipsoid + QC
	status, max_q, c, R, axes = u.calc_fully_enclosed_ellipsoid(pos, radius, write=makeVisual,directory=data_folder)

	# --- Logging ---
	logfile = GLOBAL_LOG+"ellipsoid_log.csv"
	log_exists = os.path.isfile(logfile)

	row = {
		"data_folder": data_folder,
		"size": size,
		"status": status,
		"max_q": max_q
	}

	with open(logfile, "a", newline="") as f:
		writer = csv.DictWriter(f, fieldnames=row.keys())
		if not log_exists:
			writer.writeheader()
		writer.writerow(row)

	if status != "ok":
		return np.nan

	# Compute porosity
	porosity = 1 - effective_radius_cubed / (axes[0] * axes[1] * axes[2])

	return porosity


def calc_porosity_fes(data_folder,size,relax=False,**kwargs):
	makeVisual = kwargs.get("makeVisual", False)

	pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	if pos is None:
		return np.nan
		
	effective_radius_cubed = np.sum(np.power(radius,3))

	_,full_sphere_radius = u.calc_fully_enclosed_sphere(pos,radius,write=makeVisual,directory=data_folder)

	porosity = 1-effective_radius_cubed/np.power(full_sphere_radius,3)

	return porosity

def calc_asymmetry_parameter(data_folder,size,relax=False,**kwargs):
	data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	if data is None:
		return np.nan
	# num_balls = data.shape[0]

	vol = [(4*np.pi/3)*r**3 for r in radius]

	effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 
		
	principal_moi = u.get_principal_moi(mass,data)[::-1]
	alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)

	# return alphai[0]/np.sqrt(alphai[1]*alphai[2]) #From Draine Sensitivity of Polarization to Grain Shape. I. Convex Shapes
	return np.sqrt(alphai[0]/(alphai[1]+alphai[2]-alphai[0])) #From Draine Sensitivity of Polarization to Grain Shape. II. Aggregates


def calc_stretch_parameter(data_folder,size,relax=False,**kwargs):
	data,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	if data is None:
		return np.nan
	# num_balls = data.shape[0]

	vol = [(4*np.pi/3)*r**3 for r in radius]

	effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 

		
	principal_moi = u.get_principal_moi(mass,data)[::-1]
	alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)

	return alphai[1]/np.sqrt(alphai[0]*alphai[2])

def calc_porosity_abc(data_folder,size,relax=False, **kwargs):
	makeVisual = kwargs.get("makeVisual", False)

	pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	if pos is None:
		return np.nan

	effective_radius = np.power(np.sum(np.power(radius,3)),1/3) 

	a,b,c = u.calc_equivalent_ellipsoid_principal_axes(effective_radius,pos,mass,data_folder,write=makeVisual)
		
	porosity = 1-(effective_radius**3/(a*b*c))

	return porosity

def calc_porosity_KBM(data_folder,size,relax=False,**kwargs):
	makeVisual = kwargs.get("makeVisual", False)

	pos,radius,mass,moi = u.get_data(data_folder,data_index=size,relax=relax)
	if pos is None:
		return np.nan

	effective_radius = np.power(np.sum(np.power(radius,3)),1/3)   

	RKBM = u.calc_gyration_radius(effective_radius,pos,mass,data_folder,write=makeVisual)

	porosity = 1-np.power((effective_radius/RKBM),3)

	return porosity

def calc_number_of_contacts(data_folder,data_index=-1,relax=False,**kwargs):
	makeVisual = kwargs.get("makeVisual", False)

	return u.calc_max_number_of_contacts(data_folder,data_index,relax,makeVisual)

	



def calc_fractal_dimension(data_folder,size,relax=False,**kwargs):
	makeVisual = kwargs.get("makeVisual", False)

	overwrite_octree_data = True
	show_FD_plots = False
	o3dv = u.o3doctree(data_folder,overwrite_data=overwrite_octree_data,index=size,Temp=-1,relax=relax)
	o3dv.make_tree()
	FD = o3dv.calc_fractal_dimension(show_graph=show_FD_plots)
	return FD



data_functions = [calc_porosity_abc,calc_porosity_KBM, \
					calc_number_of_contacts,calc_fractal_dimension, \
					calc_asymmetry_parameter,calc_stretch_parameter, \
					calc_porosity_fee,calc_porosity_fes, \
					calc_porosity_ch,calc_geometric_cross_section]
data_headers = [i.__name__[5:] for i in data_functions]

def calc_from_size(size,directory,existing_headers,existing_values,requested_headers,relax=False,overwrite=False,makeVisual=False):
	headers = []
	values = []

	for h,header in enumerate(data_headers):
		if not overwrite:
			#we already have this header include it
			if header in existing_headers:
				index = existing_headers.index(header)
				headers.append(header)
				#if the value is nan then we need to recalculate
				if existing_values[index] == "nan":
					existing_values[index] = str(data_functions[h](directory,size,relax,makeVisual=makeVisual))
				values.append(existing_values[index])
			#we dont have this header and we want it, so calculate
			elif (header not in existing_headers and header in requested_headers):
				headers.append(header)
				val = np.nan
				for _ in range(10):
					val = data_functions[h](directory,size,relax,makeVisual=makeVisual)
					if not np.isnan(val):
						break
				else:
					print(f"Non nan value not calculated for directory: {directory}")
				values.append(str(val))
		else:
			#we already have this header and didn't ask for it, so just include it
			if header in existing_headers and header not in requested_headers:
				index = existing_headers.index(header)
				headers.append(header)
				values.append(existing_values[index])
			#if overwrite, we want to calculate reguardless, if header is requested
			elif header in requested_headers:
				headers.append(header)
				val = np.nan
				for _ in range(10):
					val = data_functions[h](directory,size,relax,makeVisual=makeVisual)
					if not np.isnan(val):
						break
				else:
					print(f"Non nan value not calculated for directory: {directory}")
				values.append(str(val))

	return headers,values


if __name__ == '__main__':
	# print(calc_porosity_abc("/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_0/N_300/T_1000/",300,True))
	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)
	
	path = input_json["data_directory"]




	data_file = "job_data.csv" 




	N = [300]
	# N = [30,100,300]

	#list of the functions that calculate the data you want
	#if adding to this list, name your function calc_*header_name*
	#where *header_name* is what you want this data to be called
	#in the header of the job_data.csv file. The function should 
	#only take in the directory the data is in and the size of which 
	#to calculate as an input.
	#It should return a single data value.
	# bool_headers = [1,1,0,0,1,0,0,0]
	bool_headers = [0,0,1,0,0,0,0,0,0,0]
	bool_headers = [1,1,0,0,0,0,1,1,1,1]
	# requested_data_functions = [data_functions[i] for i in range(len(data_functions)) if bool_headers[i]]
	requested_data_headers = [data_headers[i] for i in range(len(data_headers)) if bool_headers[i]]

	overwritedata = False
	makeVisual = False

	for n_i,n in enumerate(N):
	
		data_folders = []
		# data_folders = [path + 'jobs/AsymBAPA_*']
		# data_folders = [path + 'jobs/BAPA_*']
		# data_folders = [path + 'jobs/BAPA_*/M_1/*']
		# data_folders = [path + 'jobs/constrollingfric*']
		# data_folders = [path + 'jobs/BAPA_0/M_60/*']
		data_folders = [path + 'jobs/SeqStickConst_*/']
		data_folders = [path + 'jobs/SeqStickLognorm_*/']
		# data_folders = [path + 'jobsCosine/lognorm_*/N_300/T_*/']
		# data_folders = [path + 'jobsCosine/lognorm_*/N_300/T_*/']
		# data_folders = [path + 'jobsNovus/const_*/N_300/T_1000/']
		# data_folders = data_folders + [path + 'jobsNovus/const_*/N_300/T_3/']
		# data_folders.append(path + f'jobsNovus/constrelax_*/N_{n}/*')
		# data_folders = [path + 'jobsNovus/constrelax_*/N_300/T_*/']

		# data_folders.append(path + 'jobs/SeqStickLognormrelax*/')
		# data_folders.append(path + 'jobs/SeqStickConstrelax*/')
		# data_folders.append(path + f'jobs/constrollingfricrelax*/N_{n}/*')
		
		# data_folders.append(path + f'jobsNovus/const_*/N_{n}/*')
		# data_folders.append(path + f'jobsCosine/lognorm_*/N_{n}/*')
		# data_folders.append(path + f'jobsNovus/constrelax_*/N_{n}/*')
		# data_folders.append(path + f'jobsCosine/lognormrelax_*/N_{n}/*')
		# data_folders.append('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsNovus/constrelax_4/N_300/T_30/')
		# data_folders.append('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_12/N_300/T_3/')
		# data_folders.append('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_7/N_300/T_1000/')
		# data_folders.append('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_19/N_300/T_1000/')
		# data_folders.append('/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/jobsCosine/lognormrelax_20/N_300/T_1000/')

		possible_dirs = []
		for data_folder in data_folders:
			possible_dirs.extend(u.get_directores_containing(data_folder,["timing.txt"]))

		#list of intermediate sizes to calculate data for.

		# requested_sizes = list(range(30,301))
		# requested_sizes = [[u.find_max_index(directory)] for directory in possible_dirs]
		requested_sizes = [[n] for directory in possible_dirs]

		for d_i,directory in enumerate(possible_dirs):
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
				all_sizes = sorted(set(existing_sizes+requested_sizes[d_i]))
				# all_sizes = [300]
				# print(all_sizes)
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
					if overwritedata or size in requested_sizes[d_i]:
						headers,values = calc_from_size(size,directory,existing_headers_for_size,existing_values_for_size,requested_data_headers,relax,overwritedata,makeVisual)
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
				# print(lines)

			else:
				print(f"Job is not finished in {directory}")
			



