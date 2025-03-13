#TODO: move the getting data stuff into the other python helper file and keep this one for calculations
#	   and such


import numpy as np
import random
# from scipy.spatial.transform import Rotation as R
import os,glob
import sys
from pathlib import Path
import json
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import h5py	
#from treelib import Node, Tree
import time
from itertools import combinations
# cwd = os.getcwd()
# os.system("cd /home/kolanzl/Open3D/build")
# sys.path.append("/home/kolanzl/Open3D/build/lib/python_package/open3d")
import open3d as o3d
##include <pybind11/stl.h>`? Or <pybind11/complex.h>,
# <pybind11/functional.h>, <pybind11/chrono.h>

data_columns = 11

#next three functions translated from matlab code from
#https://blogs.mathworks.com/cleve/files/menger.m
# def menger(level):
# 	V = [[-3,-3,-3],[-3,-3,3],[-3,3,-3],[-3,3,3],[3,-3,-3],[3,-3,3],[3,3,-3],[3,3,3]]
# 	V = np.array(V) 
# 	V = sponge(V,level)
# 	return V

# def sponge(V,level):

# 	if level > 0:
# 		V = V/3
# 		for x in [-2,0,2]:
# 			for y in [-2,0,2]:
# 				for z in [-2,0,2]:
# 					if np.sum(np.array([x,y,z])) > 0:
# 						sponge(V)
# 	# else:
# 		# cube(V)			
# 	return V	

# def cube(V):
	# return

#this function taken from 
#https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
def translate_to_cofm(mass, data):
	# Position of centre of mass in original coordinates
	# cofm = sum(mass * data) / (mass*data.shape[0])
	if mass is not list:
		mass = np.full(len(data),fill_value=mass)
	total_mass = np.sum(mass)
	cofm = np.einsum('i,ij->j', mass, data) / total_mass  # Shape: (particles, 3)
	# Transform to CofM coordinates and return
	data -= cofm
	return data

#this function taken from 
#https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
def get_inertia_matrix(mass, data):
	# Moment of intertia tensor
	
	#should in general translate to center of mass
	#but data is already there
	data = translate_to_cofm(mass, data)

	x, y, z = data.T

	Ixx = np.sum(mass * (y**2 + z**2))
	Iyy = np.sum(mass * (x**2 + z**2))
	Izz = np.sum(mass * (x**2 + y**2))
	Ixy = -np.sum(mass * x * y)
	Iyz = -np.sum(mass * y * z)
	Ixz = -np.sum(mass * x * z)
	I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
	# print(I)
	return I

#this function taken from 
#https://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
def get_principal_moi(mass,data):
	I = get_inertia_matrix(mass,data)
	Ip = np.linalg.eigvals(I)
	# Sort and convert principal moments of inertia to SI (kg.m2)
	Ip.sort()
	return Ip

def find_max_index(folder):
	files = os.listdir(folder)
	max_index = -1
	for file in files:
		if file.split("_")[0].isnumeric():
			if file.endswith("simData.csv") or file.endswith("constants.csv") or file.endswith("energy.csv") or file.endswith("data.h5"):
				# max_index = max(index_from_file(file),max_index)
				max_index = max(int(file.split("_")[0]),max_index)
	return max_index

def get_directores(root_dir):
	directories_with_file = []
	
	for dirs in glob.glob(root_dir):
		for current_dir, subdirs, files in os.walk(dirs):
			# Check if it's a deepest directory (no subdirectories)
			if not subdirs:
				directories_with_file.append(current_dir+'/')
	
	return directories_with_file

def get_directores_containing(root_dir,necessary_files):
	directories_with_file = []
	
	for dirs in glob.glob(root_dir):
		for current_dir, subdirs, files in os.walk(dirs):
			# Check if it's a deepest directory (no subdirectories)
			if not subdirs:
				has_necessary_files = True
				for necessary_file in necessary_files:
					if necessary_file not in files:
						has_necessary_files = False
				if has_necessary_files:
					if current_dir[-1] != '/':
						current_dir += '/'
					directories_with_file.append(current_dir)
	return directories_with_file
  

#returns all indices in a directory if {index}_checkpoint.txt exists, or 
#if the data file exits and timing.txt exists, or we just want to incldue unfinished ones too.
def get_all_indices(directory,checkpoint=False):
	indices = []
	is_finished = os.path.exists(directory+"timing.txt") 

	for file in os.listdir(directory):
		if file.split("_")[0].isnumeric():
			#if this index has checkpointed, or this is the data file and the simulation has finished
			if file.endswith("checkpoint.txt") or (not checkpoint and ((file.endswith("simData.csv") or (file.endswith("data.h5"))) and is_finished)):
				index = int(file.split("_")[0])
				#make sure this index isn't already in indices
				if index not in indices:
					indices.append(index)

	return list(sorted(set(indices)))

def index_from_file(file):
	file_split = file.split("_")
	if not file_split[1].isnumeric():
		return 0
	else:
		return int(file_split[0])
	exit(0)

def calc_rotational_kinetic_energy(positions, velocities, masses):
	masses = masses[:positions.shape[1]]
	 # Total mass of the aggregate
	total_mass = np.sum(masses)

	# Compute center of mass positions and velocities at each timestep

	center_of_mass_positions = np.einsum('i,tij->tj', masses, positions) / total_mass  # Shape: (timesteps, 3)
	center_of_mass_velocities = np.einsum('i,tij->tj', masses, velocities) / total_mass  # Shape: (timesteps, 3)

	# Compute relative positions and velocities (r_i - r_cm and v_i - v_cm)
	relative_positions = positions - center_of_mass_positions[:, np.newaxis, :]  # Shape: (timesteps, balls, 3)
	relative_velocities = velocities - center_of_mass_velocities[:, np.newaxis, :]  # Shape: (timesteps, balls, 3)

	# Compute rotational kinetic energy at each timestep
	# T_rot = 0.5 * sum_i m_i * |v_i'|^2
	kinetic_energies = 0.5 * masses[np.newaxis, :, np.newaxis] * relative_velocities**2  # Shape: (timesteps, balls, 3)
	rotational_kinetic_energy = np.sum(kinetic_energies, axis=(1,2))  # Shape: (timesteps,)

	return rotational_kinetic_energy

def plot_rotational_kinetic_energy(rotational_kinetic_energy, horizontal_line_value):
	import matplotlib.pyplot as plt
	"""
	Plots the rotational kinetic energy of an aggregate of spheres versus time,
	calculated about the center of mass.

	Parameters:
	positions (numpy.ndarray): Array of shape (timesteps, balls, 3) containing position data.
	velocities (numpy.ndarray): Array of shape (timesteps, balls, 3) containing velocity data.
	horizontal_line_value (float): Value at which to plot a horizontal reference line.
	masses (numpy.ndarray, optional): Array of shape (balls,) containing masses of the spheres.
									  If None, all masses are assumed to be equal and set to 1.

	The units of time are arbitrary; each timestep is considered as one unit of time.
	"""
	timesteps, = rotational_kinetic_energy.shape

	# Generate time array
	time = np.arange(timesteps)

	# Plot rotational kinetic energy versus time
	plt.figure(figsize=(10, 6))
	plt.plot(time, rotational_kinetic_energy, label='Rotational Kinetic Energy')
	plt.axhline(y=horizontal_line_value, color='red', linestyle='--', label='Reference Value')
	plt.xlabel('Time (units)')
	plt.ylabel('Rotational Kinetic Energy (units)')
	plt.title('Rotational Kinetic Energy vs. Time (About Center of Mass)')
	plt.legend()
	plt.grid(True)
	plt.show()

def plot_3d_dots(dots):
	import matplotlib.pyplot as plt
	"""
	Plots 3D dots in space using matplotlib.
	
	Parameters:
	dots (list of tuples or numpy array of shape (n, 3)): List of 3D coordinates.
	"""
	# Unpack x, y, z coordinates
	x, y, z = zip(*dots)
	
	# Create a 3D plot
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	
	# Plot points
	ax.scatter(x, y, z, c='b', marker='o')
	
	# Set labels
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	
	# Show plot
	plt.show()

def plot_3d_spheres(dots, radii):
	import matplotlib.pyplot as plt
	"""
	Plots 3D spheres in space using matplotlib.
	
	Parameters:
	dots (list of tuples or numpy array of shape (n, 3)): List of 3D coordinates for sphere centers.
	radii (list of floats or numpy array of shape (n,)): Radii for each sphere.
	"""
	# Create a 3D plot
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	
	# Plot each sphere
	for (x, y, z), radius in zip(dots, radii):
		# Create a sphere
		u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
		sphere_x = x + radius * np.cos(u) * np.sin(v)
		sphere_y = y + radius * np.sin(u) * np.sin(v)
		sphere_z = z + radius * np.cos(v)
		
		# Plot sphere surface
		ax.plot_surface(sphere_x, sphere_y, sphere_z, color='b', alpha=0.5)
	
	# Set labels
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	
	# Show plot
	plt.show()

def plot(verts,center,radius):
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')

	for i in range(len(verts)):
		print(verts[i],i)
		ax.scatter(verts[i][0],verts[i][1],verts[i][2],marker='*',color='b')
	ax.scatter(center[0],center[1],center[2],marker='.',color='r')
	ax.set_xlabel('X (row)')
	ax.set_ylabel('Y (col)')
	ax.set_zlabel('Z (dep)')

	# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	# x = np.cos(u)*np.sin(v) - center[0]
	# y = np.sin(u)*np.sin(v) - center[1]
	# z = np.cos(v) - center[2]
	# ax.plot_wireframe(x, y, z, color="r")
	plt.show()


#
def get_data_file(data_folder,data_index=-1,relax=False): #Works with csv or h5
	old = False
	file_suffix = ""
	rel = ''
	if relax:
		rel = "RELAX"
	files = os.listdir(data_folder)
	for file in files:
		if file.endswith(f"{rel}simData.csv"):
			file_suffix = f"_{rel}simData.csv"
			if file.count("_") > 1 and not relax: #relax are never old
				old = True
		if file.endswith(f"{rel}data.h5"):
			file_suffix = f"_{rel}data.h5"
	try:
		file_indicies = np.array([file.split('_')[0] for file in files\
					if file.endswith(file_suffix)],dtype=np.int64)
 
	except:
		# print("ERROR: ") 
		files = [file for file in files if file.endswith(file_suffix)]
		files = [file for file in files if '_' in file]
		file_indicies = np.array([int(file.split('_')[0]) for file in files],dtype=np.int64)

	if data_index == -1:
		index = np.max(file_indicies)
	else:
		index = data_index

	if old and data_index == 0:
		data_file = [file for file in files \
					if file.endswith(file_suffix)]
		data_file = [file for file in data_file if file.split('_')[1][0] == "R"]
	else:
		data_file = [file for file in files \
					if file.endswith(file_suffix) and file.startswith(str(index)+"_")]


	if len(data_file) == 1:
		return data_file[0]
	elif len(data_file) == 2 and len(set(file_indicies)) == 1:
		if (len(data_file[0]) > len(data_file[1])):
			return data_file[0] 
		else:
			return data_file[1]
	else:

		data_file = [file for file in files \
				if file.endswith(f"{rel}simData.csv") and file.startswith(str(index)+'_2')]
		if len(data_file) == 1:
			return data_file[0]
		elif len(data_file) == 2:
			if len(data_file[0]) > len(data_file[1]):
				return data_file[0]
			else:
				return data_file[1]
		print(f"data file of index {index} in folder '{data_folder}' not found.")
		print("Now exiting.")
		exit(-1)

def get_energy_file(data_folder,data_index=-1,relax=False):
	files = os.listdir(data_folder)
	rel = ""
	if relax:
		rel = "RELAX"
	try:
		file_indicies = np.array([file.split('_')[0] for file in files\
					if file.endswith(f"{rel}energy.csv")],dtype=np.int64)
	except: 
		files = [file for file in files if file.endswith(f'{rel}energy.csv')]
		files = [file for file in files if '_' in file]
		file_indicies = np.array([int(file.split('_')[0]) for file in files],dtype=np.int64)
	# 	file_indicies = 

	if data_index == -1:
		index = np.max(file_indicies)
	else:
		index = data_index

	# print(file_indicies)
	# print(np.max(file_indicies))

	# print("index: {}".format(index))

	data_file = [file for file in files \
				if (file.endswith(f"{rel}energy.csv") or file.endswith(f"{rel}data.h5")) and file.startswith(str(index))]
	# print(data_file)

	if len(data_file) == 1 or len(set(data_file)) == 1:
		return data_file[0]
	elif len(data_file) == 2 and len(set(file_indicies)) == 1:
		if (len(data_file[0]) > len(data_file[1])):
			return data_file[0] 
		else:
			return data_file[1]
	else:
		data_file = [file for file in files \
				if file.endswith(f"{rel}energy.csv") and file.startswith(str(index)+'_2')]
		if len(data_file) == 1:
			return data_file[0]
		elif len(data_file) == 2:
			if len(data_file[0]) > len(data_file[1]):
				return data_file[0]
			else:
				return data_file[1]
		print("energy file in folder '{}' not found.".format(data_folder))
		print("Now exiting.")
		exit(-1)

def get_line_h5data_from_file(file,linenum=-1):
	width = -1
	with h5py.File(file, 'r') as file:
		# data = file['/simData'][:]
		metadata = file['/simData'].attrs['metadata']
		for meta in metadata.split("\n"):
			md = meta.split(": ")
			# print(meta.split(": "))
			if md[0] == "row width":
				width = int(md[1])
		dat = file['/simData'][:]

		totlines = int(dat.shape[0]/width)
		
		if linenum < 0:
			start = width*(totlines+linenum)
		else:
			start = width*(linenum)
		stop = start + width

		data = np.array(dat)[start:stop]

	return data

def get_line_h5_energy_data_from_file(file,linenum=-1):
	width = -1
	with h5py.File(file, 'r') as file:
		# data = file['/simData'][:]
		metadata = file['/energy'].attrs['metadata']
		for meta in metadata.split("\n"):
			md = meta.split(": ")
			# print(meta.split(": "))
			if md[0] == "row width":
				width = int(md[1])
		dat = file['/energy'][:]

		totlines = int(dat.shape[0]/width)
		
		if linenum < 0:
			start = width*(totlines+linenum)
		else:
			start = width*(linenum)
		stop = start + width

		data = np.array(dat)[start:stop]

	return data

def find_files(folder, pattern):
	"""
	Finds all files in the given folder that match the glob pattern.

	Parameters:
	folder (str): The path to the folder.
	pattern (str): The glob pattern to match files.

	Returns:
	list: A list of file paths matching the pattern.
	"""
	p = Path(folder)
	files = [str(f) for f in p.glob(pattern) if f.is_file()]
	numbers = [int(file.split('/')[-1].split('_')[0]) for file in files]

	return [file for _, file in sorted(zip(numbers, files))]


def get_all_line_data(data_folder,data_index=-1,linenum=-1,relax=False): #Works with csv and h5
	csv_data = False
	h5_data = False
	data_file = get_data_file(data_folder,data_index,relax=relax)
	if data_file.endswith(".csv"):
		csv_data = True
	elif data_file.endswith(".h5"):
		h5_data = True
	if csv_data:
		try:
			data = np.loadtxt(data_folder + data_file,skiprows=1,dtype=float,delimiter=',')
			if data.ndim > 1:
				data = data[linenum]
		except Exception as e:
			with open(data_folder + data_file) as f:
				for line in f:
					pass
				last_line = line
			data = np.array([last_line.split(',')],dtype=np.float64)
			# print(data)
			print(f"WARNING while getting data index {data_index} in folder: {data_folder}")
			print(e)
			return None
	elif h5_data:
		data = get_line_h5data_from_file(data_folder+data_file,linenum)
	else:
		print("ERROR: datatype not recognized by utils.py: {data_file}")

	return data

def get_last_line_data(data_folder,data_index=-1,relax=False): #Works with csv and h5
	return get_line_data(data_folder,data_index,-1,relax=relax)

def get_line_data(data_folder,data_index=-1,linenum=-1,relax=False): #Works with csv and h5
	data = get_all_line_data(data_folder,data_index,linenum,relax=relax)
	return format_pos(data)

def get_last_line_energy(data_folder,data_index=-1,relax=False):
	energy_file = get_energy_file(data_folder,data_index,relax=relax)
	if energy_file.endswith(".csv"):
		try:
			energy = np.loadtxt(data_folder + energy_file,skiprows=1,dtype=float,delimiter=',')
			if energy.ndim > 1:
				energy = energy[-1]
			print(energy)
		except Exception as e:
			with open(data_folder + energy_file) as f:
				for line in f:
					pass
				last_line = line
			energy = np.array([last_line.split(',')],dtype=np.float64)
			print("ERROR CAUGHT getting energy in folder: {}".format(data_folder))
			print(e)
			# return None
	elif energy_file.endswith(".h5"):
		energy = get_line_h5_energy_data_from_file(data_folder + energy_file,-1)
	else:
		print(f"ERROR: file extension not recognized for file '{energy_file}'")
	# print("DATA LEN: {} for file {}{}".format(data.size,data_folder,data_file))
	# print("FOR {} Balls".format(data.size/11))
	return energy

def get_radii(data_folder,data_index=-1,relax=False):
	return get_constants(data_folder,data_index,relax)[0]

def get_masses(data_folder,data_index=-1,relax=False):
	return get_constants(data_folder,data_index,relax)[1]

def get_moi(data_folder,data_index=-1,relax=False):
	return get_constants(data_folder,data_index,relax)[2]


#Writes a single line to the given simData file
#pos, vel, and w should be in the form [[0x,0y,0z],[1x,1y,1z], . . .] where the number is the particle index
def write_simData(data_file,pos,w,vel):
	num_particles = len(pos)
	if data_file.endswith(".csv"):
		if not os.path.exists(data_file):
			with open(data_file, "w") as fp:
				for i in range(num_particles):
					header_sentence = f"x{i},y{i},z{i},wx{i},wy{i},wz{i},wmag{i},vx{i},vy{i},vz{i},bound{i}"
					if i < num_particles-1:
						fp.write(f"{header_sentence},")
					else:
						fp.write(header_sentence)

		with open(data_file, "a") as fp:
			for i in range(len(pos)):
				sentence = f"{pos[i][0]},{pos[i][1]},{pos[i][2]},{w[i][0]},{w[i][1]},{w[i][2]},{np.linalg.norm(w[i])},{vel[i][0]},{vel[i][1]},{vel[i][2]},0"
				if i == 0:
					fp.write(f"\n{sentence},")
				elif i < num_particles-1:
					fp.write(f"{sentence},")
				else:
					fp.write(f"{sentence}\n")

	elif data_file.endswith(".h5"):
		print("Write simData with h5 not implimented yet.")

#Writes a single line to the given simData file
def write_energy(data_file,time, PE, KE, E, p, L):
	if data_file.endswith(".csv"):
		if not os.path.exists(data_file):
			with open(data_file, "w") as fp:
				fp.write("Time,PE,KE,E,p,L\n")

		with open(data_file, "a") as fp:
			fp.write(f"{time},{PE},{KE},{E},{p},{L}\n")

	elif data_file.endswith(".h5"):
		print("Write simData with h5 not implimented yet.")

def write_constants(data_file, radii, mass, moi):
	if data_file.endswith(".csv"):
		with open(data_file, "w") as fp:
			for i in range(len(radii)):
				fp.write(f"{radii[i]},{mass[i]},{moi[i]}\n")

	elif data_file.endswith(".h5"):
		print("Write constants with h5 not implimented yet.")

def get_constants(data_folder,data_index=-1,relax=False):#Works with csv and h5
	data_file = get_data_file(data_folder,data_index,relax=relax)
	data_file = data_file.replace('simData','constants')	
	if data_file.endswith(".csv"):
		data_constants = np.loadtxt(data_folder+data_file,skiprows=0,dtype=float,delimiter=',')
	elif data_file.endswith(".h5"):
		with h5py.File(data_folder+data_file, 'r') as file:
			consts = file['/constants'][:]
			data_constants = np.array(consts).reshape(int(len(consts)/3),3)
	else:
		print(f"ERROR: data file type not recognized by utils.py: {data_file}")
	return data_constants[:,0],data_constants[:,1],data_constants[:,2]

def get_all_constants(data_folder,data_index=-1,relax=False): #Works with csv and h5
	data_file = get_data_file(data_folder,data_index,relax=relax)
	data_file = data_file.replace('simData','constants')	
	if data_file.endswith(".csv"):
		data_constants = np.loadtxt(data_folder+data_file,skiprows=0,dtype=float,delimiter=',')
	elif data_file.endswith(".h5"):
		with h5py.File(data_folder+data_file, 'r') as file:
			consts = file['/constants'][:]
			data_constants = np.array(consts).reshape(int(len(consts)/3),3)
	else:
		print(f"ERROR: data file type not recognized by utils.py: {data_file}")
	return data_constants

#a line of data is in the following format
#x0,y0,z0,wx0,wy0,wz0,wmag0,vx0,vy0,vz0,bound0
def format_pos(data):
	if data is not None and not np.isscalar(data):
		data = np.reshape(data,(int(data.size/data_columns),data_columns))
		data = data[:,:3] #pos is first three of every row
	return data

def format_w(data):
	data = np.reshape(data,(int(data.size/data_columns),data_columns))
	data = data[:,3:6] #w is second three of every row
	return data

def format_vel(data):
	data = np.reshape(data,(int(data.size/data_columns),data_columns))
	data = data[:,7:10] #vel is after 3x pos, 3x w, 1x w mag, and is 3 long
	return data

def calcCOM(pos,mass):
	com = np.array([0,0,0],dtype=np.float64)
	mtot = 0

	for ball in range(pos.shape[0]):
		com += mass[ball]*pos[ball]
		mtot += mass[ball]

	return com/mtot

def COM(data_folder,data_index=-1,relax=False):
	data = get_last_line_data(data_folder,data_index,relax=relax)

	consts = get_all_constants(data_folder,data_index,relax=relax)
	com = np.array([0,0,0],dtype=np.float64)
	mtot = 0

	for ball in range(data.shape[0]):
		com += consts[ball][0]*data[ball]
		mtot += consts[ball][0]

	return com/mtot

def get_data(data_folder,data_index=-1,linenum=-1,relax=False): #Works with both csv and h5
	if data_folder == '/home/kolanzl/Desktop/bin/merger.csv':
		data = np.loadtxt(data_folder,delimiter=',')
		radius = 1
		mass = 1
		moi = 1
	else:
		# data_file = get_data_file(data_folder,data_index,relax=relax)
		radius,mass,moi = get_constants(data_folder,data_index,relax=relax)

		data = get_line_data(data_folder,data_index,linenum,relax=relax)

	return data,radius,mass,moi

def get_all_pos_data(data_file):
	csv_data = False
	h5_data = False
	if data_file.endswith(".csv"):
		csv_data = True
	elif data_file.endswith(".h5"):
		h5_data = True
	if csv_data:
		try:
			data = np.loadtxt(data_file,skiprows=1,dtype=float,delimiter=',')
			if data.ndim > 1:
				return np.array([format_pos(d) for d in data])
		except Exception as e:
			# with open(data_folder + data_file) as f:
			# 	for line in f:
			# 		pass
			# 	last_line = line
			# data = np.array([last_line.split(',')],dtype=np.float64)
			# print(data)
			print("WARNING while getting data in folder: {}".format(data_folder))
			print(e)
			return None
	elif h5_data:
		# data = get_line_h5data_from_file(data_folder+data_file,linenum)
		print("ERROR: h5 not implimented for get_all_pos_data function in utils.py")
	
	else:
		print("ERROR: datatype not recognized by utils.py: {data_file}")

	return None

def get_all_vel_data(data_file):
	csv_data = False
	h5_data = False
	if data_file.endswith(".csv"):
		csv_data = True
	elif data_file.endswith(".h5"):
		h5_data = True
	if csv_data:
		try:
			data = np.loadtxt(data_file,skiprows=1,dtype=float,delimiter=',')
			if data.ndim > 1:
				return np.array([format_vel(d) for d in data])
		except Exception as e:
			# with open(data_folder + data_file) as f:
			# 	for line in f:
			# 		pass
			# 	last_line = line
			# data = np.array([last_line.split(',')],dtype=np.float64)
			# print(data)
			print("WARNING while getting data in folder: {}".format(data_folder))
			print(e)
			return None
	elif h5_data:
		# data = get_line_h5data_from_file(data_folder+data_file,linenum)
		print("ERROR: h5 not implimented for get_all_pos_data function in utils.py")
	
	else:
		print("ERROR: datatype not recognized by utils.py: {data_file}")

	return None

def get_all_data(data_folder,data_index=-1,linenum=-1,relax=False): #Works with both csv and h5

	# data_file = get_data_file(data_folder,data_index)

	radius,mass,moi = get_constants(data_folder,data_index,relax=relax)

	data = get_all_line_data(data_folder,data_index,linenum,relax=relax)
	pos = format_pos(data)
	w = format_w(data)
	vel = format_vel(data)


	return pos,vel,w,radius,mass,moi


def get_data_range(data_folder,data_index=-1,relax=False):
	if data_folder == '/home/kolanzl/Desktop/bin/merger.csv':
		data = np.loadtxt(data_folder,delimiter=',')
		radius = 1
		mass = 1
		moi = 1
	else:
		data = get_last_line_data(data_folder,data_index,relax=relax)
		radius,m,moi = get_constants(data_folder,data_index,relax=relax)

	max_x = np.max(data[:,0]) + np.max(radius)
	min_x = np.min(data[:,0]) - np.max(radius)
	max_y = np.max(data[:,1]) + np.max(radius)
	min_y = np.min(data[:,1]) - np.max(radius)
	max_z = np.max(data[:,2]) + np.max(radius)
	min_z = np.min(data[:,2]) - np.max(radius)
	
	return max_x,min_x,max_y,min_y,max_z,min_z

#following functions taken from 
#http://www.open3d.org/docs/release/tutorial/geometry/voxelization.html#Voxel-carving
def xyz_spherical(xyz):
	x = xyz[0]
	y = xyz[1]
	z = xyz[2]
	r = np.sqrt(x * x + y * y + z * z)
	r_x = np.arccos(y / r)
	r_y = np.arctan2(z, x)
	return [r, r_x, r_y]


def get_rotation_matrix(r_x, r_y):
	rot_x = np.asarray([[1, 0, 0], [0, np.cos(r_x), -np.sin(r_x)],
						[0, np.sin(r_x), np.cos(r_x)]])
	rot_y = np.asarray([[np.cos(r_y), 0, np.sin(r_y)], [0, 1, 0],
						[-np.sin(r_y), 0, np.cos(r_y)]])
	return rot_y.dot(rot_x)


def get_extrinsic(xyz):
	rvec = xyz_spherical(xyz)
	r = get_rotation_matrix(rvec[1], rvec[2])
	t = np.asarray([0, 0, 2]).transpose()
	trans = np.eye(4)
	trans[:3, :3] = r
	trans[:3, 3] = t
	return trans


def preprocess(model):
	import open3d as o3d
	min_bound = model.get_min_bound()
	max_bound = model.get_max_bound()
	center = min_bound + (max_bound - min_bound) / 2.0
	scale = np.linalg.norm(max_bound - min_bound) / 2.0
	# scale = 1
	vertices = np.asarray(model.vertices)
	vertices -= center
	model.vertices = o3d.utility.Vector3dVector(vertices / scale)
	return model

def preprocess_pt(model):
	import open3d as o3d
	min_bound = model.get_min_bound()
	max_bound = model.get_max_bound()
	print("max: {}\tmin: {}".format(max_bound,min_bound))
	center = min_bound + (max_bound - min_bound) / 2.0
	scale = np.linalg.norm(max_bound - min_bound) / 2.0
	# scale = 1
	vertices = np.asarray(model.points)
	vertices -= center
	model.points = o3d.utility.Vector3dVector(vertices / scale)
	return model


def vox_carve(mesh,
				  cubic_size,
				  voxel_resolution,
				  w=300,
				  h=300,
				  use_depth=True,
				  surface_method='pointcloud'):
	# mesh.compute_vertex_normals()
	import open3d as o3d
	mesh.estimate_normals()
	camera_sphere = o3d.geometry.TriangleMesh.create_sphere()

	# setup dense voxel grid
	voxel_carving = o3d.geometry.VoxelGrid.create_dense(
		width=cubic_size,
		height=cubic_size,
		depth=cubic_size,
		voxel_size=cubic_size / voxel_resolution,
		origin=[-cubic_size / 2.0, -cubic_size / 2.0, -cubic_size / 2.0],
		color=[1.0, 0.7, 0.0])
	print("Vox size: {}".format(cubic_size / voxel_resolution))
	# rescale geometry
	camera_sphere = preprocess(camera_sphere)
	mesh = preprocess_pt(mesh)

	# setup visualizer to render depthmaps
	vis = o3d.visualization.Visualizer()
	vis.create_window(width=w, height=h, visible=False)
	vis.add_geometry(mesh)
	vis.get_render_option().mesh_show_back_face = True
	ctr = vis.get_view_control()
	param = ctr.convert_to_pinhole_camera_parameters()

	# carve voxel grid
	pcd_agg = o3d.geometry.PointCloud()
	centers_pts = np.zeros((len(camera_sphere.vertices), 3))
	for cid, xyz in enumerate(camera_sphere.vertices):
		# get new camera pose
		trans = get_extrinsic(xyz)
		param.extrinsic = trans
		c = np.linalg.inv(trans).dot(np.asarray([0, 0, 0, 1]).transpose())
		centers_pts[cid, :] = c[:3]
		ctr.convert_from_pinhole_camera_parameters(param)

		# capture depth image and make a point cloud
		vis.poll_events()
		vis.update_renderer()
		depth = vis.capture_depth_float_buffer(False)
		pcd_agg += o3d.geometry.PointCloud.create_from_depth_image(
			o3d.geometry.Image(depth),
			param.intrinsic,
			param.extrinsic,
			depth_scale=1)

		# depth map carving method
		if use_depth:
			voxel_carving.carve_depth_map(o3d.geometry.Image(depth), param)
		else:
			voxel_carving.carve_silhouette(o3d.geometry.Image(depth), param)
		# print("Carve view %03d/%03d" % (cid + 1, len(camera_sphere.vertices)))
	vis.destroy_window()

	# add voxel grid survace
	print('Surface voxel grid from %s' % surface_method)
	if surface_method == 'pointcloud':
		voxel_surface = o3d.geometry.VoxelGrid.create_from_point_cloud_within_bounds(
			pcd_agg,
			voxel_size=cubic_size / voxel_resolution,
			min_bound=(-cubic_size / 2, -cubic_size / 2, -cubic_size / 2),
			max_bound=(cubic_size / 2, cubic_size / 2, cubic_size / 2))
	elif surface_method == 'mesh':
		voxel_surface = o3d.geometry.VoxelGrid.create_from_triangle_mesh_within_bounds(
			mesh,
			voxel_size=cubic_size / voxel_resolution,
			min_bound=(-cubic_size / 2, -cubic_size / 2, -cubic_size / 2),
			max_bound=(cubic_size / 2, cubic_size / 2, cubic_size / 2))
	else:
		raise Exception('invalid surface method')
	voxel_carving_surface = voxel_surface + voxel_carving

	return voxel_carving_surface, voxel_carving, voxel_surface

def dist(pt1,pt2):
	return np.sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2 + (pt1[2]-pt2[2])**2)


class datamgr(object):
	"""docstring for datamgr"""
	# def __init__(self, data_folder,voxel_buffer=5,ppb=3000):
	def __init__(self, data_folder,index=-1,ppb=30000,Temp=-1,relax=False):
		super(datamgr, self).__init__()
		self.relax = relax
		self.data_folder = data_folder
		self.index = index
		if data_folder != '/home/kolanzl/Desktop/bin/merger.csv' and Temp < 0:
			for i in data_folder.split('/'):
				if i[:2] == "T_":
					self.Temp = float(i.split("_")[1])
					print(f"automatically setting Temp in datamgr to {self.Temp}K")

		else:
			self.Temp = Temp
		#how many points in single ball pointcloud shell
		self.ppb = ppb
		self.data,self.radius,self.mass,self.moi = get_data(self.data_folder,self.index,relax=self.relax)
		self.nBalls = self.data.shape[0]
		# self.buffer = voxel_buffer # how many extra voxels in each direction 
		self.data_range = get_data_range(self.data_folder,self.index,relax=self.relax)

	def shift_to_first_quad(self,data_range=None):
		if data_range is None:
			data_range = get_data_range(self.data_folder,self.index,relax=self.relax)
		# print("SHIFTED")

		self.data[:,0] -= data_range[1] 
		self.data[:,1] -= data_range[3] 
		self.data[:,2] -= data_range[5] 


	def vox_init(self,num_vox):
		data_abs_max = max(self.data_range,key=abs) 
		self.vox_size = (data_abs_max*2)/num_vox
		# self.vox_size = 
		# print(self.vox_per_radius)
		# self.vox_rep = np.zeros((num_vox+self.buffer*2,num_vox+self.buffer*2,num_vox+self.buffer*2))

	#Function written by chatGPT
	def rotation_matrix(v1, v2):
		"""
		Returns the rotation matrix between two vectors v1 and v2.
		Both v1 and v2 must be numpy arrays with the same shape.

		:param v1: First vector
		:param v2: Second vector
		:return: Rotation matrix
		"""
		v1 = np.array(v1)
		v2 = np.array(v2)
		if v1.shape != v2.shape:
			raise ValueError("Both vectors must have the same shape.")
		v1 = v1 / np.linalg.norm(v1)
		v2 = v2 / np.linalg.norm(v2)
		v = np.cross(v1, v2)
		s = np.linalg.norm(v)
		c = np.dot(v1, v2)
		vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
		rotation_matrix = np.eye(v1.shape[0]) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))
		return rotation_matrix

	def orient_data(self):
		max_lengsq = -1
		pt1 = []
		pt2 = []
		for i,p1 in enumerate(self.data):
			for j,p2 in enumerate(self.data):
				if i != j:
					lengsq = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2
					if max_lengsq < lengsq:
						max_lengsq = lengsq
						pt1 = p1
						pt2 = p2
		print(max_lengsq)
		print(pt1)
		print(pt2)


	def gen_whole_pt_cloud(self):
		# self.orient_data()
		# exit(0)
		self.shift_to_first_quad()

		accums = []
		for ind,pt in enumerate(self.data):
			radii = np.linspace(self.radius[ind]/100,self.radius[ind],100)

			accum = [self.ppb*(radius**2/self.radius[ind]**2) for radius in radii]
			accum = np.array(accum,dtype=int)
			accum = np.where(accum < 100, 100, accum)
			accums.append(accum)
		summm = 0
		return_array = np.zeros((np.sum(np.array(accums)),3))
		for ind,pt in enumerate(self.data):
			radii = np.linspace(self.radius[ind]/100,self.radius[ind],100)

			for i,radius in enumerate(radii):
				start_index = int(np.sum(accums[:ind])) + np.sum(accums[ind][:i])
				end_index = int(np.sum(accums[:ind])) + np.sum(accums[ind][:i+1])
				return_array[start_index:end_index] = self.gen_pt_cloud(pt,radius,accums[ind][i]) 
		return return_array


	#evenly spaced points on sphere code form:
	#http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
	def gen_pt_cloud(self,pt,radius,num_pts):
		goldenRatio = (1 + 5**0.5)/2

		return_array = np.zeros((num_pts,3))	
		i = np.arange(0, num_pts)

		theta = 2 * np.pi * i / goldenRatio
		phi = np.arccos(1 - 2*(i+0.5)/num_pts)

		# return_array = np.zeros((self.ppb,3))
		return_array[:,0] = radius*(np.cos(theta) * np.sin(phi)) + pt[0]
		return_array[:,1] = radius*(np.sin(theta) * np.sin(phi)) + pt[1]
		return_array[:,2] = radius*np.cos(phi) + pt[2]
		
		return return_array

	# def get_center_pt(self,ind):
	# 	return [ind[0]*self.radius+self.radius/2,ind[1]*self.radius+self.radius/2,ind[2]*self.radius+self.radius/2]


class o3doctree(object):
	"""docstring for o3doctree"""
	def __init__(self, data_folder=None,ppb=30000,verbose=False,overwrite_data=False, \
				visualize_pcd=False, visualize_octree=False, index=-1,Temp=-1,relax=False):
	# def __init__(self, data_folder, max_depth=8,ppb=600000,verbose=False):
		super(o3doctree, self).__init__()
		import open3d as o3d
		self.data_folder = data_folder
		self.ppb = ppb
		self.verbose = verbose
		self.overwrite_data = overwrite_data
		self.visualize_pcd = visualize_pcd
		self.visualize_octree = visualize_octree
		self.bestfitlen = 4
		self.index = index
		self.relax = relax
		self.Temp = -1
		if Temp > 0:
			self.Temp = Temp
		else:
			#extract temp from data_folder
			for i in data_folder.split('/'):
				if i[:2] == "T_":
					self.Temp = float(i.split("_")[1])
					print(f"automatically setting Temp in o3doctree to {self.Temp}K")
		if self.Temp < 0:
			print("WARNING: Could not set temp based on data_folder in o3doctree constructor.")
			self.Temp = -1


	def make_tree(self):

		self.dm = datamgr(self.data_folder,self.index,self.ppb,Temp=self.Temp,relax=self.relax)

		bounds = [self.dm.data_range[0]-self.dm.data_range[1],self.dm.data_range[2]-self.dm.data_range[3],self.dm.data_range[4]-self.dm.data_range[5]]
		
		max_bound = max(bounds) + np.max(self.dm.radius)*2
		# max_bound = max(bounds) + 2*self.dm.radius
		n=0
		rad = 999
		# while rad > self.dm.radius:
		while rad > np.min(self.dm.radius):
			n += 1
			rad = max_bound/(2**n)

		self.max_depth = n


		pcd_file = self.data_folder + "pointcloud_ppb-{}.pcd".format(self.ppb)
		oct_file = self.data_folder + "octree_ppb-{}".format(self.ppb)
		fractdim_data_file = self.data_folder + "fractdim_ppb-{}.csv".format(self.ppb)
		
		make_data = False
		try:
			print("Loading FD data for :{}".format(fractdim_data_file))
			d = np.loadtxt(fractdim_data_file, delimiter=',')
			# print(d)
			# print(np.sum(d[:,1]))
			if np.sum(d[:,1]) == 0:
				print("Data needs recomputing")
				os.remove(fractdim_data_file)
				make_data = True
			else:
				with open(fractdim_data_file,'r') as f:
					header = f.readline()
				self.octree_size = float(header.strip('\n').strip('# '))
				self.s_data = d[:,0]
				self.Ns_data = d[:,1]

			# with open(fractdim_data_file,'r') as f:
			# 	data = f.readlines()
			# 	print(data)
			# if len(data) == 0:
			# 	make_data = True
			# if np.sum(d[:,0]) == 0:
		# except IOError:
		except:
			# print("Computing data for :{}".format(fractdim_data_file))
			make_data = True

		octree = []
		if make_data or self.overwrite_data:
			octverb = ''
			print("Generating FD data for :{}".format(fractdim_data_file))
			if os.path.isfile(oct_file) and not self.overwrite_data:
				octstart = time.process_time()
				octverb = 'Getting'
				octree = o3d.io.read_octree(oct_file)
				if self.visualize_octree:
					self.show_octree(self.verbose)
			else:
				if self.verbose:
					pcdstart = time.process_time()
				
				pcdverb = ''
				# if os.path.isfile(pcd_file):
				pcd = []
				make_pcd_data = False
				try:
					with open(pcd_file,'r') as f:
						data = f.readlines()
						if len(data) == 0:
							make_pcd_data = True
							os.remove(pcd_file)
				# except IOError:
				except:
					make_pcd_data = True
				if not make_pcd_data and not self.overwrite_data:
					pcdverb = 'Getting'
					pcd = o3d.io.read_point_cloud(pcd_file)
					if self.visualize_pcd:
						self.show_pcd(pcd,self.verbose)
				else:
					pcdverb = 'Making'
					point_cloud = self.dm.gen_whole_pt_cloud()
					pcd = o3d.geometry.PointCloud()
					pcd.points = o3d.utility.Vector3dVector(point_cloud)
					pcd.colors = o3d.utility.Vector3dVector(np.random.uniform(0, 1, size=point_cloud.shape))
					#Saving the pointcloud doesn't really help
					# o3d.io.write_point_cloud(pcd_file, pcd)
					if self.visualize_pcd:
						self.show_pcd(pcd,self.verbose)
					# exit(0)
			
				if self.verbose:
					pcdend = time.process_time()
					print("{} pcd took {:.2f}s".format(pcdverb,pcdend-pcdstart))

				octverb = 'Making'
				octstart = time.process_time()

				octree = o3d.geometry.Octree(max_depth=self.max_depth)
				if self.visualize_octree:
					self.show_octree(octree,self.verbose)
				octree.convert_from_point_cloud(pcd, size_expand=0.01)

				#Until the documentation for open3d says what file extension works for octree data, this can't be saved
				# o3d.io.write_octree(oct_file, octree)


			if self.verbose:
				octend = time.process_time()
				print("{} octree took {:.2f}s".format(octverb, octend-octstart))
				start = time.process_time()
				print("Starting octree traversal")
		
			self.tree_info = []
			for i in range(self.max_depth):
				self.tree_info.append(0)
			# print(self.tree_info)
			if self.verbose:
				end = time.process_time()
				print("Traversing octree took {:.2f}s".format(end-start))
			octree.traverse(self.f_traverse)

			self.s_data = np.zeros((self.max_depth))
			self.Ns_data = np.zeros((self.max_depth))
			for i in range(self.max_depth):
				self.s_data[i] = (2**-(i+1))
				self.Ns_data[i] = self.tree_info[i]
			save_data = np.zeros((self.max_depth,2))
			save_data[:,0] = self.s_data
			save_data[:,1] = self.Ns_data
			np.savetxt(fractdim_data_file,save_data, delimiter=',',header=str(octree.size))
			self.octree_size = octree.size
		# else:
			
	#TODO  This function should find the orientation that minimizes 
	#	   the original fractal dimension (depth of 1)
	# def point_orientation(self,point_cloud):
	# 	# print(point_cloud)
	# 	# exit(0)
	# 	i = 0
	# 	best_i = 0
	# 	rotations = []
	# 	pcd = o3d.geometry.PointCloud()
	# 	while (i < 10):
	# 		xrot = random.uniform(0,360)
	# 		yrot = random.uniform(0,360)
	# 		zrot = random.uniform(0,360)
	# 		rotation_matrix = R.from_euler('xyz',[xrot,yrot,zrot],degrees=True).as_matrix()
	# 		rotations.append(rotation_matrix)
	# 		temp_point_cloud[:] = rotation_matrix @ point_cloud[:]
	# 		# temp_point_cloud = [rotation_matrix @ i for i in point_cloud]
	# 		pcd.points = o3d.utility.Vector3dVector(temp_point_cloud)
	# 		octree = o3d.geometry.Octree(max_depth=1)#check max_depth def
	# 		octree.convert_from_point_cloud(pcd, size_expand=0.01)#check size_expand def
	# 		self.tree_info = [0]
			
	# 		octree.traverse(self.f_traverse)
	# 		print(self.tree_info)
	# 		exit(0)
	# 		i+=1
	# 		# print(rotation_matrix.as_matrix())
	# 	exit(0)
	# 	return o3d.utility.Vector3dVector(point_cloud)

	# def add_menger_points(self,data):
	# 	dlen = data.shape
	# 	print(dlen)
	# 	exit(0)

	def test_menger_sponge(self):
		self.data_folder = '.'
		merger_file = '/home/kolanzl/Desktop/bin/merger.csv'
		self.dm = datamgr(merger_file,self.ppb)

		self.tree_info = []
		max_depth = 8
		self.s_data = np.zeros((max_depth))
		self.Ns_data = np.zeros((max_depth))
		for i in range(max_depth):
			self.tree_info.append(0)
		menger_points = np.loadtxt(merger_file,delimiter=',')
		# menger_points = self.add_menger_points(menger_points)

		pcd = o3d.geometry.PointCloud()
		pcd.points = o3d.utility.Vector3dVector(menger_points)
		pcd.colors = o3d.utility.Vector3dVector(np.random.uniform(0, 1, size=menger_points.shape))

		octree = o3d.geometry.Octree(max_depth=max_depth)
		octree.convert_from_point_cloud(pcd, size_expand=0.01)
		# o3d.visualization.draw_geometries([octree])
		# exit(0)

		octree.traverse(self.f_traverse)
		for i in range(max_depth):
			self.s_data[i] = (2**-(i+1))
			self.Ns_data[i] = self.tree_info[i]
		
		self.calc_fractal_dimension(True)
		# self.Ns_data[-1] -=
		print(np.log(self.Ns_data)/np.log(1/self.s_data))


	def show_octree(self,octree,verbose):
		if verbose:
			start = time.process_time()
		o3d.visualization.draw_geometries([octree])
		if verbose:
			end = time.process_time()
			print("Visualizing octree took {}".format(end-start))

	def show_pcd(self,pcd,verbose):
		if verbose:
			start = time.process_time()
		o3d.visualization.draw_geometries([pcd])
		if verbose:
			end = time.process_time()
			print("Visualizing pcd took {}".format(end-start))	

	def bestfit(self,x_data,y_data,length,min_rang=None,max_range=None):
		acceptable_indicies = []
		if min_rang is not None and max_range is not None:
			for i in range(len(x_data)):
				#Note that min_range/max_range is before taking the inverse (1/min_range etc)
				#so min_range is actually larger than max_range. 
				if x_data[i] < min_rang and x_data[i] > max_range:
					acceptable_indicies.append(i)

			x_data = x_data[acceptable_indicies]
			y_data = y_data[acceptable_indicies]
		if length > len(x_data):
			length = len(x_data)
		if length == 1:
			print("ERROR: Cannot fit to one data point")

		index_list = np.arange(0,len(x_data))
		rsq = 0
		rsq_index = 0 
		fit = []
		combos = combinations(index_list,length)
		combos = [list(c) for c in list(combos)]
		winning_combo = []

		for i,comb in enumerate(combos):
			# print(comb)
			ind_combo = np.array(sorted(comb))
			fit_x_pts = np.log(1/x_data[ind_combo])
			fit_y_pts = np.log(y_data[ind_combo])
			fit.append(np.polyfit(fit_x_pts,fit_y_pts,1))
			y_predict = np.array(fit_x_pts*fit[-1][0] + fit[-1][1])
			corr_matrix = np.corrcoef(fit_y_pts, y_predict)
			corr = corr_matrix[0,1]
			new_rsq = corr**2
			if 1-new_rsq < 1 - rsq:
				rsq = new_rsq
				rsq_index = i
				winning_combo = comb

		# print(x_data[winning_combo])
		return fit[rsq_index]


	def calc_fractal_dimension(self,show_graph=False):
		import matplotlib.pyplot as plt
		OIsize = self.octree_size
		S0 = 1
		
		# x_IO = np.zeros((self.max_depth))
		
		fract_dim_fit = self.bestfit(self.s_data,self.Ns_data,self.bestfitlen)
		# fract_dim_fit = self.bestfit(self.s_data,self.Ns_data,self.bestfitlen,1,self.dm.radius/OIsize)

		fig, ax = plt.subplots(2,1)
		ax.flatten()


		ax[0].set_title('D = {:.2f}, unit side length'.format(fract_dim_fit[0]))
		# ax[1].set_title('D = {:.2f} for T = {}'.format(fract_dim_fit[0],self.dm.Temp))
		ax[0].plot(1/self.s_data,self.Ns_data,marker='*',label='Frac dim data')
		ax[0].loglog(1/self.s_data,np.exp(np.log(1/self.s_data)*fract_dim_fit[0]+fract_dim_fit[1]),label='log(y) = {:.2f}*log(x) + {:.2f}'.format(fract_dim_fit[0],fract_dim_fit[1]))
		ax[0].set_xlabel('log(1/(Unit side lengths))')
		ax[0].set_ylabel('log(Number of boxes to enclose)')


		ax[1].plot(1/self.s_data,np.log(self.Ns_data)/np.log(1/self.s_data))
		ax[1].set_xscale('log')
		ax[1].set_xlabel('log(1/(Unit side lengths))')
		ax[1].set_ylabel('Fractal Dimension')
		
		avg = np.mean(np.log(self.Ns_data)/np.log(1/self.s_data))
		# print('y avg: {}'.format(avg))
		ax[1].axhline(y=np.mean(np.log(self.Ns_data)/np.log(1/self.s_data)), color='g')

		ax[0].axvline(1/1, color='g')
		ax[0].axvline(1/(np.mean(self.dm.radius)/OIsize), color='g')
		# ax[0].axvline(1/(self.dm.radius/OIsize), color='g')

		fig.legend()
		# plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
		# mode="expand", borderaxespad=0, ncol=3)
		fig.set_figheight(10)
		fig.set_figwidth(10)
		plt.tight_layout()
		plt.savefig(self.data_folder+"FractDim.png")
		if show_graph:
			plt.show()
		plt.close()
		# fig.close()

		return fract_dim_fit[0]
	

	#function adapted from:
	#http://www.open3d.org/docs/release/python_example/geometry/octree/index.html
	def f_traverse(self,node, node_info):

		if isinstance(node, o3d.geometry.OctreePointColorLeafNode):
			self.tree_info[node_info.depth-1] += 1
		if isinstance(node, o3d.geometry.OctreeInternalPointNode):
			for child in node.children:
				if child is not None:
					self.tree_info[node_info.depth-1] += 1
					break
		#if return True, f_traverse will stop prematurely
		early_stop = False
		return early_stop


#############################Novus Functions#############################
def get_squeue_output():
    try:
        # Run the squeue command and capture its output

        result = subprocess.run(['squeue', '-o', '"%.20u %.25j"'], capture_output=True, text=True)
        output = result.stdout
        return output
    except subprocess.CalledProcessError as e:
        # Handle any errors that occur during the command execution
        print(f"Error executing squeue: {e}")
        return None


def same_job(fullpath, job_name):

	fpsplit = fullpath.split('/')
	start_ind = fpsplit.index("SpaceLab_data") + 1

	fpattrs = re.split(r'\D+',"".join(fpsplit[start_ind:-1]))
	fpattrs = [int(i) for i in fpattrs if len(i) > 0]
	
	qattrs = re.split(r'\D+',job_name)
	qattrs = [int(i) for i in qattrs if len(i) > 0]

	if len(fpattrs) != len(qattrs):
		print("ERROR IN same_job")
		exit(0)

	for i in range(len(qattrs)):
		if fpattrs[i] != qattrs[i]:
			return False
	return True

def on_queue(fullpath):
	queue_out = get_squeue_output()
	for line in queue_out.split('\n')[1:]:
		line = line.strip('"').split()
		if len(line) > 0:
			if line[0] == "kolanzl" and line[1] != "interactive":
				if same_job(fullpath,line[1]):
					return True
	return False
