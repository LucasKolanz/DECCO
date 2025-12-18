#TODO: move the getting data stuff into the other python helper file and keep this one for calculations
#	   and such

import numpy as np
# import scipy as sp
# import cvxpy as cp
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
from itertools import product
import subprocess
import re
# cwd = os.getcwd()
# os.system("cd /home/kolanzl/Open3D/build")
# sys.path.append("/home/kolanzl/Open3D/build/lib/python_package/open3d")
# import open3d as o3d
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

def get_h5data_from_file(file):
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
		
		data = np.array(dat).reshape(-1, width)

	return data

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

def get_h5_energy_data_from_file(file):
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

		if width > 0:
			total = len(dat)
			dat = dat.reshape(int(total/width),width)
		
		data = np.array(dat)

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


#returns all pos,vel,w (from all timesteps) for a given folder and data index
#as arrays of shape: pos[timestep,ball,xyz]
def get_simData(data_folder,data_index=-1,relax=False): #Works with csv and h5

	data_file = get_data_file(data_folder,data_index,relax=relax)

	if data_file.endswith(".csv"):
		data = np.loadtxt(data_folder + data_file,skiprows=1,dtype=float,delimiter=',')

	elif data_file.endswith(".h5"):
		data = get_h5data_from_file(data_folder+data_file)
	else:
		print("ERROR: datatype not recognized by utils.py: {data_file}")

	steps = data.shape[0]
	balls = int(data.shape[1]/data_columns)

	#shape is pos[timestep,ball,xyz]
	pos = np.zeros((steps,balls,3),dtype=np.float64)
	vel = np.zeros((steps,balls,3),dtype=np.float64)
	w = np.zeros((steps,balls,3),dtype=np.float64)

	for step in range(steps):
		data_line = data[step,:]
		pos[step] = format_pos(data_line)
		vel[step] = format_vel(data_line)
		w[step] = format_w(data_line)

	return pos,vel,w

#returns all the energy data as seperat arrays from a given folder for a given index
def get_energy(data_folder,data_index=-1,relax=False):
	energy = []
	energy_file = get_energy_file(data_folder,data_index,relax=relax)
	if energy_file.endswith(".csv"):
		energy = np.loadtxt(data_folder + energy_file,skiprows=1,dtype=float,delimiter=',')
	elif energy_file.endswith(".h5"):
		energy = get_h5_energy_data_from_file(data_folder + energy_file)
	else:
		print(f"ERROR: file extension not recognized for file '{energy_file}'")

	print(energy)
	time = energy[:,0]
	PE = energy[:,1]
	KE = energy[:,2]
	E = energy[:,3]
	p = energy[:,4]
	L = energy[:,5]

	return time,PE,KE,E,p,L

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
#x0,y0,z0,wx0,wy0,wz0,wmag0,vx0,vy0,vz0,bound0, ... ,xN,yN,zN,wxN,wyN,wzN,wmagN,vxN,vyN,vzN,boundN
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

def calcCOM(mass, pos):
	mass = np.asarray(mass, dtype=np.float64)
	pos  = np.asarray(pos,  dtype=np.float64)

	if pos.shape[0] != mass.shape[0]:
		raise ValueError("mass and pos must have same length")
	if pos.shape[1] != 3:
		raise ValueError("pos must be an (N,3) array")

	mtot = mass.sum()
	com = np.sum(pos * mass[:, None], axis=0)
	return com / mtot

def COM(data_folder,data_index=-1,relax=False):
	data = get_last_line_data(data_folder,data_index,relax=relax)

	consts = get_all_constants(data_folder,data_index,relax=relax)
	com = np.array([0,0,0],dtype=np.float64)
	mtot = 0

	for ball in range(data.shape[0]):
		com += consts[ball][0]*data[ball]
		mtot += consts[ball][0]

	return com/mtot


#########################################################################
#Start functions for writing metric visualizations to file
#########################################################################
#writes bounding sphere to temp file for blender to read later
def write_gyration_radius_sphere(center,radius,directory=""):
	"""
	center : np.array shape (3,)
	radius : float
	"""

	if directory == "":
		import tempfile
		# Use system temp directory, e.g. /tmp on Linux
		directory = tempfile.gettempdir()
	path = os.path.join(directory, "gyration_radius_sphere.json")

	data = {
		"center": center.tolist(),
		"radius": float(radius)
	}

	# Write with high precision so Blender gets the exact geometry
	with open(path, "w") as f:
		json.dump(data, f, indent=2)

	print(f"[INFO] Gyration Radius sphere written to {path}")
	return path

#writes bounding sphere to temp file for blender to read later
def write_enclosing_sphere(center,radius,directory=""):
	"""
	center : np.array shape (3,)
	radius : float
	"""

	if directory == "":
		import tempfile
		# Use system temp directory, e.g. /tmp on Linux
		directory = tempfile.gettempdir()
	path = os.path.join(directory, "enclosing_sphere.json")

	data = {
		"center": center.tolist(),
		"radius": float(radius)
	}

	# Write with high precision so Blender gets the exact geometry
	with open(path, "w") as f:
		json.dump(data, f, indent=2)

	print(f"[INFO] Enclosing sphere written to {path}")
	return path

def write_convex_hull(points, hull, radii=None, directory=""):
	import tempfile

	points = np.asarray(points, float)

	# Map global indices -> compact local indices
	vmap = {idx: i for i, idx in enumerate(hull.vertices)}
	simplices_mapped = [
		[vmap[a], vmap[b], vmap[c]]
		for (a, b, c) in hull.simplices
	]

	if directory == "":
		import tempfile
		# Use system temp directory, e.g. /tmp on Linux
		directory = tempfile.gettempdir()
	path = os.path.join(directory, "convex_hull.json")

	data = {
		"points": points[hull.vertices].tolist(),   # only hull vertices
		"simplices": simplices_mapped,              # remapped faces
		"hull_vertices": hull.vertices.tolist()
	}

	if radii is not None:
		data["radii"] = np.asarray(radii, float).tolist()


	with open(path, "w") as f:
		json.dump(data, f, indent=2)

	print(f"[INFO] Convex hull written to {path}")
	return path

def write_equivalent_ellipsoid(center, axes, R,directory=""):
	"""
	Write the bounding ellipsoid parameters to a JSON file in directory.

	Parameters
	----------
	center : (3,) array-like
		Ellipsoid center.
	axes   : (3,) array-like
		Ellipsoid semi-axis lengths (a, b, c).
	R      : (3,3) array-like
		Rotation matrix whose columns are the principal axis directions.

	Creates /tmp/bounding_ellipsoid.json with high precision.
	"""

	center = np.asarray(center, dtype=float)
	axes   = np.asarray(axes,   dtype=float)
	R      = np.asarray(R,      dtype=float)

	if directory == "":
		import tempfile
		# Use system temp directory, e.g. /tmp on Linux
		directory = tempfile.gettempdir()
	path = os.path.join(directory, "equivalent_ellipsoid.json")

	data = {
		"center": center.tolist(),
		"axes": axes.tolist(),
		"R": R.tolist(),        # full 3×3 rotation matrix
	}

	with open(path, "w") as f:
		json.dump(data, f, indent=2)

	print(f"[INFO] Equivalent ellipsoid written to {path}")
	return path

def write_enclosing_ellipsoid(center, axes, R,directory=""):
	"""
	Write the bounding ellipsoid parameters to a JSON file in /tmp.

	Parameters
	----------
	center : (3,) array-like
		Ellipsoid center.
	axes   : (3,) array-like
		Ellipsoid semi-axis lengths (a, b, c).
	R      : (3,3) array-like
		Rotation matrix whose columns are the principal axis directions.

	Creates /tmp/bounding_ellipsoid.json with high precision.
	"""

	center = np.asarray(center, dtype=float)
	axes   = np.asarray(axes,   dtype=float)
	R      = np.asarray(R,      dtype=float)

	if directory == "":
		import tempfile
		# Use system temp directory, e.g. /tmp on Linux
		directory = tempfile.gettempdir()
	path = os.path.join(directory, "enclosing_ellipsoid.json")

	data = {
		"center": center.tolist(),
		"axes": axes.tolist(),
		"R": R.tolist(),        # full 3×3 rotation matrix
	}

	with open(path, "w") as f:
		json.dump(data, f, indent=2)

	print(f"[INFO] Bounding ellipsoid written to {path}")
	return path

def write_shadow_grid(xs, ys, shadow, k=None, delta=None,
							 proj_points=None, radii=None, directory=""):
	"""
	Write the 2D projected grid + shadow mask to a temporary JSON file
	so Blender can visualize it.

	Parameters
	----------
	xs : (Nx,)  array of x-coordinates of grid cell centers
	ys : (Ny,)  array of y-coordinates of grid cell centers
	shadow : (Ny, Nx) boolean array (True = cell is shadowed)
	k : (3,) projection/view direction (optional)
	delta : float, grid spacing (optional)
	proj_points : (N,2) projected monomer centers (optional)
	radii : (N,) monomer radii (optional)

	Returns
	-------
	path : str
		Path to the JSON file.
	"""

	xs = np.asarray(xs, float)
	ys = np.asarray(ys, float)
	shadow = np.asarray(shadow, bool)

	Ny, Nx = shadow.shape
	assert xs.shape[0] == Nx
	assert ys.shape[0] == Ny

	if directory == "":
		import tempfile
		# Use system temp directory, e.g. /tmp on Linux
		directory = tempfile.gettempdir()
	path = os.path.join(directory, "shadow_grid.json")

	data = {
		"xs": xs.tolist(),         # x grid centers
		"ys": ys.tolist(),         # y grid centers
		"shadow": shadow.tolist(), # 2D boolean mask
	}

	if k is not None:
		data["k"] = np.asarray(k, float).tolist()

	if delta is not None:
		data["delta"] = float(delta)

	if proj_points is not None:
		data["proj_points"] = np.asarray(proj_points, float).tolist()

	if radii is not None:
		data["radii"] = np.asarray(radii, float).tolist()

	with open(path, "w") as f:
		json.dump(data, f, indent=2)

	print(f"[INFO] Shadow grid written to {path}")
	return path

#########################################################################
#END functions for writing metric visualizations to file
#########################################################################

def number_of_contacts(pos,radii):
	# pos,radii,mass,moi = u.get_data(data_folder,data_index=size,linenum=line,relax=relax)

	if pos is None:
		return np.nan
	pos = np.array(pos)

	contacts = u.calc_contacts(pos,radii)

	# contacts = np.zeros((num_balls,num_balls),dtype=int)
	# dist = lambda i,j: np.sqrt((data[i][0]-data[j][0])**2 + (data[i][1]-data[j][1])**2 + \
	# 		(data[i][2]-data[j][2])**2)

	# ##IT FEELS LIKE THIS IS OVERCOUNTING BUT IT ISN'T.
	# ##EACH BALL FELLS ITS OWN CONTACTS.
	# for i in range(num_balls):
	# 	for j in range(num_balls):
	# 		if i != j:
	# 			contacts[i,j] = (dist(i,j) <= (radius[i]+radius[j]))

	return np.mean(np.sum(contacts,axis=1))

def calc_contacts(pos,radii):
	num_balls = pos.shape[0]
	contacts = np.zeros((num_balls,num_balls),dtype=int)
	dist = lambda i,j: np.sqrt((pos[i][0]-pos[j][0])**2 + (pos[i][1]-pos[j][1])**2 + \
			(pos[i][2]-pos[j][2])**2)

	##IT FEELS LIKE THIS IS OVERCOUNTING BUT IT ISN'T.
	##EACH BALL FELLS ITS OWN CONTACTS.
	for i in range(num_balls):
		for j in range(num_balls):
			if i != j:
				contacts[i,j] = (dist(i,j) <= (radii[i]+radii[j]))
	return contacts

#########################################################################
#Start functions for calculating metrics
#########################################################################

def calc_max_number_of_contacts(data_folder,data_index=-1,relax=False,makeVisual=False):

	line = 0
	max_nc = -1
	max_line = -1
	pos,radii,mass,moi = u.get_data(data_folder,data_index=data_index,linenum=line,relax=relax)
	nc = number_of_contacts(pos,radii)
	while not np.isnan(nc):
		# max_nc = max(max_nc,nc)
		line += 1 
		pos,radii,mass,moi = u.get_data(data_folder,data_index=data_index,linenum=line,relax=relax)
		nc = number_of_contacts(pos,radii)
		if max_nc < nc:
			max_nc = nc
			max_line = line

	if makeVisual:
		pos,radii,mass,moi = u.get_data(data_folder,data_index=data_index,linenum=max_line,relax=relax)

	return max_nc

#finds smallest needed sphere to enclose the aggregate 
#returns that center and the radius of the enclosing sphere
def calc_fully_enclosed_sphere(pos,radii,write=False,directory=""):
	
	"""
	pos   : (N,3) array of centers
	radii : (N,)  array of radii (can all be different)
	"""
	import cvxpy as cp

	pos = np.asarray(pos, dtype=float)
	radii = np.asarray(radii, dtype=float)
	assert pos.shape[0] == radii.shape[0]
	assert pos.shape[1] == 3

	N = pos.shape[0]

	# Optimization variables: center c and radius R
	c = cp.Variable(3)
	R = cp.Variable()

	constraints = [R >= 0]
	for i in range(N):
		# ||c - p_i||_2 + r_i <= R
		constraints.append(cp.norm(c - pos[i, :], 2) + radii[i] <= R)

	prob = cp.Problem(cp.Minimize(R), constraints)
	prob.solve(solver=cp.ECOS, verbose=False)

	if prob.status not in ("optimal", "optimal_inaccurate"):
		raise RuntimeError(f"ECOS failed: {prob.status}")

	c_opt = np.array(c.value, dtype=float)

	# Recompute R from c_opt to guarantee containment
	d = np.linalg.norm(pos - c_opt, axis=1) + radii
	R_opt = float(d.max())

	# Optional debug print:
	# print(f"R_opt: {R_opt}\nmax(dist+ri): {d.max()}\nslack:{R_opt - d.max()}", R_opt - d.max())

	if write:
		write_enclosing_sphere(c_opt, R_opt, directory)

	return c_opt, R_opt


def calc_fully_enclosed_ellipsoid(pos, radii, samples_per_sphere=8192,write=False, directory="", verbose=False):
	"""
	Minimum-volume enclosing ellipsoid of a union of spheres in R^3,
	approximated by sampling sphere surfaces and running MVEE.
	"""

	relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + "../"
	project_path = os.path.abspath(relative_path) + '/'
	sys.path.append(project_path + "utilities/mvee/")
	from mvee import mvee2

	pos = np.asarray(pos, dtype=float)
	radii = np.asarray(radii, dtype=float)
	m, d = pos.shape
	assert d == 3, "This helper is for 3D only."

	# ---- sample surfaces ----
	dirs = fibonacci_sphere(samples_per_sphere)  # (S, 3)

	samples_list = []
	for c, r in zip(pos, radii):
		pts = c + r * dirs
		samples_list.append(pts)

	samples = np.vstack(samples_list).astype(float)  # (m*S, 3)

	# ---- normalize for MVEE ----
	mean = samples.mean(axis=0)
	samples_centered = samples - mean
	scale = np.max(np.linalg.norm(samples_centered, axis=1))
	if scale < 1e-12:
		scale = 1.0

	samples_norm = samples_centered / scale
	X = samples_norm.T

	# ---- run MVEE ----
	# ret = mvee2(X, full_output=True, epsilon=1e-6)
	ret = mvee2(X, hybrid={}, full_output=True, epsilon=1e-6)

	center_n, R_n, axes_n, A_n = ellipse_from_mvee2_output(ret)

	status, max_q = check_ellipsoid_encloses_samples(samples_norm, center_n, A_n, verbose=verbose)

	# ---- undo normalization ----
	center = center_n * scale + mean
	axes   = axes_n * scale
	R      = R_n
	A      = A_n / (scale**2)   # optional: A in original coords if you want it

	# ---- validate numerically ----
	if verbose:
		print(f"Ellipsoid check: status={status}, max_q={max_q:.6e}")

	if write:
		write_enclosing_ellipsoid(center, axes, R, directory)
	return status, max_q, center, R, axes

def calc_gyration_radius(effective_radius,pos,mass,directory,write=False):
	principal_moi = u.get_principal_moi(mass,pos)[::-1]

	alphai = principal_moi/(0.4*np.sum(mass)*effective_radius**2)

	RKBM = np.sqrt(np.sum(alphai)/3) * effective_radius

	if makeVisual:
		COM = u.calcCOM(mass,pos)
		write_gyration_radius_sphere(COM,RKBM,data_folder)

	return RKBM

def calc_equivalent_ellipsoid_principal_axes(effective_radius,pos,mass,directory,write=False):
	principal_moi = get_principal_moi(mass,pos)[::-1]

	total_mass = np.sum(mass)
	alphai = principal_moi/(0.4*total_mass*effective_radius**2)
	
	a = effective_radius * np.sqrt(alphai[1] + alphai[2] - alphai[0])
	b = effective_radius * np.sqrt(alphai[2] + alphai[0] - alphai[1])
	c = effective_radius * np.sqrt(alphai[0] + alphai[1] - alphai[2])

	if write:
		center = calcCOM(mass,pos)
		moi_world = get_inertia_matrix(mass,pos)
		
		#parallel axis theorem incase COM isnt origin
		moi_com = moi_world-total_mass*(np.dot(center,center)*np.eye(3)-np.outer(center,center))

		#enforce symmetry in case of numeric noise
		moi_com = 0.5 * (moi_com + moi_com.T)
		
		#Eigendecomposition: principal moments and axes
		vals, vecs = np.linalg.eigh(moi_com)
		# vals = principal moments (I1, I2, I3)
		# vecs[:, i] = eigenvector for vals[i]
		
		#Sort for consistent axis order (e.g. ascending inertia)
		idx = np.argsort(vals)
		vals = vals[idx]
		vecs = vecs[:, idx]

		R = vecs
		# Ensure right-handed coordinate system (det = +1)
		if np.linalg.det(R) < 0:
			R[:, 0] *= -1  # flip one axis

		axes = [c,b,a]
		write_equivalent_ellipsoid(center, axes, R, directory)

	return a,b,c

def calc_hull_volume(pos,radii,samples_per_sphere=64,write=False,directory=""):
	"""
	makes a convec hull around the aggregate and calculates it's volume

	pos    : (m,3) centers of spheres
	radii  : (m,) radii of spheres
	samples_per_sphere : number of directions sampled on each sphere surface

	Returns:
		volume : volume enclosed by the convex hull
	"""
	import scipy as sp

	pos = np.asarray(pos, dtype=float)
	radii = np.asarray(radii, dtype=float)
	m, d = pos.shape
	assert d == 3, "This helper is for 3D only."

	# #Get points on spheres closest to MVEE
	# status,max_q,center,R,axes = calc_fully_enclosed_ellipsoid(pos,radii)
	# # exit(0)

	# # If the MVEE failed return nan so we can try again later
	# if status != "ok":
	# 	return np.nan

	# support_points = np.vstack([
	# 	sphere_point_relative_to_ellipsoid(pos[i], radii[i], center, axes, R)
	# 	for i in range(len(pos))
	# ])
	

	dirs = fibonacci_sphere(samples_per_sphere)  # (S, 3)
	samples_list = []
	for c, r in zip(pos, radii):
		pts = c + r * dirs
		samples_list.append(pts)

	support_points = np.vstack(samples_list).astype(float)  # (m*S, 3)

	# --- Build convex hull ---
	# print("calcing hull")
	hull = sp.spatial.ConvexHull(support_points)

	# sanity check for Blender
	if hull.simplices.size == 0:
		raise RuntimeError("Convex hull is degenerate (no faces). Cannot export.")

	if write:
		write_convex_hull(support_points, hull, directory=directory)

	return hull.volume

def calc_geometric_cross_section(
	pos,
	radii,
	direction=(0.0, 0.0, 1.0),
	r1=None,
	mesh_factor=0.0055,
	write = False,
	directory = ""
):
	"""
	Compute geometric cross section of an aggregate by mesh-counting the projected shadow.

	Parameters
	----------
	pos : (N, 3) array
		Monomer centers.
	radii : (N,) array or float
		Monomer radii. If float, same radius for all.
	r1 : float, optional
		Reference monomer radius used for mesh size. If None and radii is array,
		uses min(radii).
	k : (3,) array-like, optional
		Viewing direction (will be normalized).
	mesh_factor : float, optional
		Grid spacing as fraction of r1 (Suyama uses 0.0055).

	Returns
	-------
	sigma : float
		Geometric cross section for this viewing direction.
	"""

	pos = np.asarray(pos, dtype=float)
	N = pos.shape[0]
	if np.isscalar(radii):
		radii = np.full(N, float(radii))
	else:
		radii = np.asarray(radii, dtype=float)

	if r1 is None:
		r1 = float(np.min(radii))
	delta = mesh_factor * r1

	# Orthonormal basis (u, v) spanning plane ⟂ k
	k = np.asarray(direction, dtype=float)
	k /= np.linalg.norm(k)
	u, v = _orthonormal_basis_from_k(k)

	# Project centers into (X, Y) plane
	X = pos @ u   # shape (N,)
	Y = pos @ v   # shape (N,)

	r_max = np.max(radii)

	X_min = X.min() - r_max
	X_max = X.max() + r_max
	Y_min = Y.min() - r_max
	Y_max = Y.max() + r_max

	Nx = int(np.ceil((X_max - X_min) / delta))
	Ny = int(np.ceil((Y_max - Y_min) / delta))

	# Centers of grid cells
	xs = X_min + (np.arange(Nx) + 0.5) * delta
	ys = Y_min + (np.arange(Ny) + 0.5) * delta

	shadow = np.zeros((Ny, Nx), dtype=bool)

	# For each monomer, mark cells whose centers fall inside its projected circle
	for Xi, Yi, ri in zip(X, Y, radii):
		if ri <= 0:
			continue

		# Find index ranges potentially affected by this sphere
		ix_min = max(0, int(np.floor((Xi - ri - X_min) / delta)))
		ix_max = min(Nx - 1, int(np.floor((Xi + ri - X_min) / delta)))
		iy_min = max(0, int(np.floor((Yi - ri - Y_min) / delta)))
		iy_max = min(Ny - 1, int(np.floor((Yi + ri - Y_min) / delta)))

		if ix_min > ix_max or iy_min > iy_max:
			continue

		local_x = xs[ix_min:ix_max+1]
		local_y = ys[iy_min:iy_max+1]

		dx2 = (local_x[None, :] - Xi)**2   # shape (1, n_x_local)
		dy2 = (local_y[:, None] - Yi)**2   # shape (n_y_local, 1)

		mask_local = dx2 + dy2 <= ri**2
		shadow[iy_min:iy_max+1, ix_min:ix_max+1] |= mask_local

	N_shadow = np.count_nonzero(shadow)
	sigma = N_shadow * (delta**2)

	if write:
		write_shadow_grid(xs, ys, shadow, k=k, delta=delta,
								 proj_points=None, radii=None, directory=directory)
	return sigma

# def mvee_points(points, tol=0.0001):
# 	"""
# 	Finds the ellipse equation in "center form"
# 	(x-c).T * A * (x-c) = 1
# 	"""
# 	import numpy.linalg as la

# 	N, d = points.shape
# 	Q = np.column_stack((points, np.ones(N))).T
# 	err = tol+1.0
# 	u = np.ones(N)/N
# 	while err > tol:
# 		# assert u.sum() == 1 # invariant
# 		X = np.dot(np.dot(Q, np.diag(u)), Q.T)
# 		M = np.diag(np.dot(np.dot(Q.T, la.inv(X)), Q))
# 		jdx = np.argmax(M)
# 		step_size = (M[jdx]-d-1.0)/((d+1)*(M[jdx]-1.0))
# 		new_u = (1-step_size)*u
# 		new_u[jdx] += step_size
# 		err = la.norm(new_u-u)
# 		u = new_u
# 	c = np.dot(u, points)
# 	A = la.inv(np.dot(np.dot(points.T, np.diag(u)), points)
# 			   - np.multiply.outer(c, c))/d

# 	# Ensure A is symmetric (numerical noise)
# 	A = 0.5 * (A + A.T)

# 	# Eigen-decomposition
# 	lam, R = np.linalg.eigh(A)

# 	# Safety: remove tiny negative eigenvals due to numerical noise
# 	lam = np.maximum(lam, 1e-15)

# 	# Semi-axes
# 	axes = 1.0 / np.sqrt(lam)    # (a1, a2, a3)

# 	center = c                   # center of ellipsoid
# 	return center,R,axes,0.0



def fibonacci_sphere(n_points):
	"""
	Fibonacci sphere directions, roughly uniform on S^2.
	Returns (n_points, 3) unit vectors.
	"""
	n_points = int(n_points)
	k = np.arange(n_points, dtype=float)
	phi = np.pi * (1.0 + np.sqrt(5.0))  # golden angle

	z = 1.0 - 2.0 * (k + 0.5) / n_points
	r = np.sqrt(np.maximum(0.0, 1.0 - z*z))
	theta = phi * k

	x = r * np.cos(theta)
	y = r * np.sin(theta)

	return np.stack([x, y, z], axis=1)  # (n_points, 3)

def ellipse_from_mvee2_output(retvals):
	"""
	Extract center, rotation R, axes, and canonical A from mvee2() output.

	mvee2 returns:
		retvals["mat"] = B = H^{-1}

	The ellipsoid in normalized coordinates is:
		x^T B x <= n_dim

	We convert to the standard form:
		(x - c)^T A (x - c) <= 1
	with:
		A = B / n_dim
		axes_i = sqrt(n_dim / mu_i)
	where mu_i are eigenvalues of B.
	"""
	# ---- center ----
	if "c" in retvals:
		c = np.asarray(retvals["c"], float).flatten()
	else:
		c = np.zeros(retvals["mat"].shape[0], float)

	# ---- B = H^{-1} ----
	B = np.asarray(retvals["mat"], float)
	B = 0.5 * (B + B.T)      # symmetrize for safety

	n_dim = B.shape[0]

	# ---- eigendecomposition of B ----
	mu, R = np.linalg.eigh(B)   # B = R diag(mu) R^T
	mu = np.maximum(mu, 1e-15)  # avoid divide by zero

	# ---- canonical A and axes ----
	# A = B / n_dim
	# axes_i^2 = n_dim / mu_i => axes_i = sqrt(n_dim / mu_i)
	axes = np.sqrt(n_dim / mu)

	A = R @ np.diag(mu / n_dim) @ R.T

	return c, R, axes, A

def check_ellipsoid_encloses_samples(samples, center, A, verbose=False,
									 tol_ok=1e-3, tol_warn=1e-2):
	"""
	samples : (N, d) in the SAME normalized space as center, A
	center  : (d,)
	A       : (d, d) canonical ellipsoid matrix (x-c)^T A (x-c) <= 1
	"""
	samples = np.asarray(samples, float)
	center  = np.asarray(center, float)

	Xc = samples - center  # (N, d)

	# Correct quadratic form:
	q = np.sum((Xc @ A) * Xc, axis=1)
	max_q = float(np.max(q))

	if verbose:
		print(f"check: max_q = {max_q:.6e}")

	if max_q <= 1.0 + tol_ok:
		status = "ok"
	elif max_q <= 1.0 + tol_warn:
		status = "warn"
	else:
		status = "fail"

	return status, max_q




def analyze_enclosing_ellipsoid_for_balls(pos, radii, center, A,
										  n_dirs_check=64,
										  eps_list=(1e-6, 1e-5, 1e-4)):
	"""
	Perform a numerical error analysis for an enclosing ellipsoid of balls.

	Ellipsoid is given as:
		E = { x | (x - center)^T A (x - center) <= 1 }

	We:
	  1) Sample points on each ball surface and check max violation:
		   v = max_j ( (x_j - center)^T A (x_j - center) - 1 )
		 If v <= 0 (or very small +), the ellipsoid encloses all sampled points.

	  2) Define F(c) = max_j (x_j - c)^T A (x_j - c) over the same samples.
		 For small perturbations c' = c + eps * e_k, we check
			 F(c') - F(c).
		 If these are >= 0 for all small eps and coordinate directions,
		 that’s analogous to your sphere test and suggests c is locally optimal
		 for the sampled points.

	Parameters
	----------
	pos : (N,3) array of centers of balls
	radii : (N,) array of radii
	center : (3,) ellipsoid center
	A : (3,3) positive-definite matrix
	n_dirs_check : int, number of directions per sphere to sample
	eps_list : iterable of epsilons for the perturbation test

	Prints diagnostics and returns a dict with key metrics.
	"""
	pos = np.asarray(pos, dtype=float)
	radii = np.asarray(radii, dtype=float)
	center = np.asarray(center, dtype=float)
	A = np.asarray(A, dtype=float)

	assert pos.shape[0] == radii.shape[0]
	assert pos.shape[1] == 3
	assert center.shape == (3,)
	assert A.shape == (3, 3)

	N = pos.shape[0]

	# 1. Build a denser sample of surface points on all balls
	dirs = _fibonacci_sphere(n_dirs_check)  # (n_dirs, 3)

	pts = []
	ball_idx = []  # which ball each point comes from
	for i in range(N):
		p = pos[i]
		r = radii[i]
		pts_i = p[None, :] + r * dirs  # (n_dirs_check, 3)
		pts.append(pts_i)
		ball_idx.extend([i] * n_dirs_check)
	pts = np.vstack(pts)                 # (N * n_dirs_check, 3)
	ball_idx = np.array(ball_idx, dtype=int)

	# 2. Evaluate the quadratic form (x - c)^T A (x - c) on these points
	diff = pts - center[None, :]          # (M,3)
	# q_j = (x_j - c)^T A (x_j - c)
	q = np.einsum("ij,ij->i", diff @ A, diff)  # (M,)

	max_q = float(q.max())
	worst_idx = int(q.argmax())
	worst_ball = int(ball_idx[worst_idx])

	violation = max_q - 1.0  # >0 means outside if we treat ==1 as boundary

	print("=== Ellipsoid containment check (sampled points) ===")
	print(f"Max q = (x-c)^T A (x-c) over sampled points: {max_q:.6e}")
	print(f"Max violation max(q - 1): {violation:.6e}")
	print(f"Worst sample index: {worst_idx}, from ball index: {worst_ball}")
	print(f"Worst sample point: {pts[worst_idx]}")

	# 3. Define F(c) = max_j (x_j - c)^T A (x_j - c) on these sampled points
	def F(c_vec):
		c_vec = np.asarray(c_vec, dtype=float)
		d = pts - c_vec[None, :]
		q_local = np.einsum("ij,ij->i", d @ A, d)
		return float(q_local.max())

	F_c = F(center)
	print("\n=== Center optimality (perturbation) check ===")
	print(f"F(c_opt) = max_j (x_j - c_opt)^T A (x_j - c_opt) = {F_c:.6e}")

	# 4. Check F(c + eps * e_k) - F(c) for coordinate directions
	directions = np.eye(3)  # e_x, e_y, e_z
	results = {}
	for eps in eps_list:
		for k in range(3):
			d = directions[k]
			c_pert = center + eps * d
			F_pert = F(c_pert)
			delta = F_pert - F_c
			print(f"eps={eps:g}, dir=e_{k}:  F(c+eps*e_{k}) - F(c) = {delta:.6e}")
			results[(eps, k)] = delta

	return {
		"max_q": max_q,
		"violation": violation,
		"worst_sample_index": worst_idx,
		"worst_ball_index": worst_ball,
		"worst_point": pts[worst_idx],
		"F_c": F_c,
		"perturbation_deltas": results,
	}








def ellipsoid_eval(x, C, axes, R):
	y = R.T @ (x - C)
	return (y[0]/axes[0])**2 + (y[1]/axes[1])**2 + (y[2]/axes[2])**2

#c is center of sphere
#r is radius of sphere
#C is center of ellipsoid
#axes is pricipal axes of ellipsoid
#R is rotation matirx of ellipsoid
def sphere_point_relative_to_ellipsoid(c, r, C, axes, R):
	"""
	If sphere is fully inside ellipsoid → return closest sphere point to ellipsoid.
	If any part of sphere sticks outside → return the sphere point farthest outside.
	"""

	x_ell = project_point_onto_ellipsoid(c, C, axes, R)
	d = x_ell - c
	n = np.linalg.norm(d)
	if n == 0:
		# sphere centered at ellipsoid center → choose arbitrary direction
		d_unit = np.array([1,0,0])
	else:
		d_unit = d / n

	# Candidate points
	p_near = c + r * d_unit     # closest sphere point to ellipsoid
	p_far  = c - r * d_unit     # farthest sphere point from ellipsoid

	# Determine if sphere is fully inside
	if ellipsoid_eval(p_far, C, axes, R) <= 1.0:
		# whole sphere inside ellipsoid
		return p_near
	else:
		# sphere pokes out
		return p_far

#find the point on the ellipsoid closest to p (the center of a grain)
#C is the center of the ellipsoid, axes are the three semi-axes lengths of the ellipsoid, R is the rotation matrix of the ellipsoid
# def project_point_onto_ellipsoid(p, C, axes, R):
# 	from scipy.optimize import root_scalar
# 	import numpy as np

# 	p = np.asarray(p, float)
# 	C = np.asarray(C, float)
# 	axes = np.asarray(axes, float)
# 	R = np.asarray(R, float)

# 	#Transform into ellipsoid coordinates.
# 	#u is point p in these coordinates
# 	u = R.T @ (p - C)

# 	#now a2[0] = a^2, a2[1] = b^2, a2[2] = c^2
# 	a2 = axes**2

# 	# print(f"point: {p}")
# 	# print(f"center: {C}")

# 	#use newtons method to solve for lam, the lagrange multiplier
# 	def F(lam):
# 		return np.sum((u*u)*a2 / (a2 + lam)**2) - 1.0

# 	def Fprime(lam):
# 		return np.sum(-2*(u*u)*a2 / (a2 + lam)**3)

# 	# Newton always converges if p is not exactly the ellipsoid center.
# 	sol = root_scalar(F, fprime=Fprime, x0=0.0, method='newton')
# 	lam = sol.root

# 	#use lam to find the point we want
# 	x_local = u * a2 / (a2 + lam)

# 	# print(f"returning: {C + R @ x_local}")

# 	return C + R @ x_local

#find the point on the ellipsoid closest to p (the center of a grain)
#C is the center of the ellipsoid, axes are the three semi-axes lengths of the ellipsoid, R is the rotation matrix of the ellipsoid
def project_point_onto_ellipsoid(p, C, axes, R):
	from scipy.optimize import root_scalar
	import numpy as np

	p = np.asarray(p, float)
	C = np.asarray(C, float)
	axes = np.asarray(axes, float)
	R = np.asarray(R, float)

	# Transform to ellipsoid coordinates
	u = R.T @ (p - C)
	a2 = axes**2

	# Function
	def F(lam):
		return np.sum((u*u)*a2 / (a2 + lam)**2) - 1.0

	# Safe bracket
	#This choice guarantees F(lam_lo) < 0 and F(lam_high) > 0
	#Thus the solution has to be inside this range.
	lam_lo = -np.min(a2) * 0.9999999  # stays just inside valid region
	lam_hi =  np.max(a2) * 100       # safely positive

	# Use Brent's method (fast + guaranteed convergence)
	sol = root_scalar(F, bracket=(lam_lo, lam_hi), method='brentq')

	lam = sol.root

	# Closest point in ellipsoid local coords
	x_local = u * a2 / (a2 + lam)

	# Back to world
	return C + R @ x_local


def _orthonormal_basis_from_k(k):
	"""Given a unit vector k, return two unit vectors u, v that span the plane ⟂ k."""
	k = np.asarray(k, dtype=float)
	k /= np.linalg.norm(k)

	# Pick a reference vector not parallel to k
	if abs(k[2]) < 0.9:
		a = np.array([0.0, 0.0, 1.0])
	else:
		a = np.array([0.0, 1.0, 0.0])

	u = np.cross(k, a)
	u /= np.linalg.norm(u)
	v = np.cross(k, u)
	# v should already be unit if k, u are unit and orthogonal
	return u, v



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
		import open3d as o3d
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
		import open3d as o3d
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
		import open3d as o3d
		if verbose:
			start = time.process_time()
		o3d.visualization.draw_geometries([octree])
		if verbose:
			end = time.process_time()
			print("Visualizing octree took {}".format(end-start))

	def show_pcd(self,pcd,verbose):
		import open3d as o3d
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
		import open3d as o3d
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


############################Helpful Functions############################
def unroll(*lists):
	normalized = [lst if len(lst) > 0 else [None] for lst in lists]
	return list(product(*normalized))

#############################Novus Functions#############################
def get_squeue_output():
	try:
		# Run the squeue command and capture its output

		result = subprocess.run(['squeue', '-o', '%u %j'], capture_output=True, text=True)
		output = result.stdout
		return output
	except subprocess.CalledProcessError as e:
		# Handle any errors that occur during the command execution
		print(f"Error executing squeue: {e}")
		return None


# def same_job(fullpath, job_name):

# 	print(job_name)

# 	fpsplit = fullpath.split('/')
# 	start_ind = fpsplit.index("SpaceLab_data") + 1

# 	fpattrs = re.split(r'\D+',"".join(fpsplit[start_ind:-1]))
# 	fpattrs = [int(i) for i in fpattrs if len(i) > 0]
	
# 	qattrs = re.split(r'\D+',job_name)
# 	qattrs = [int(i) for i in qattrs if len(i) > 0]


# 	print(fpattrs)
# 	print(qattrs)

# 	if len(fpattrs) != len(qattrs):
# 		return False
# 		# print("ERROR IN same_job")
# 		# exit(0)

# 	for i in range(len(qattrs)):
# 		if fpattrs[i] != qattrs[i]:
# 			return False
	# return True

def on_queue(job_name):
	queue_out = get_squeue_output()
	for line in queue_out.split('\n')[1:]:
		line = line.strip('"').split()
		if len(line) > 0:
			if line[0] == "kolanzl" and line[1] != "interactive":
				# if same_job(fullpath,line[1]):
				if job_name == line[1]:
					return True
	return False


def rand_seed():
	import datetime

	random.seed(datetime.datetime.now().timestamp())
	# Generating a random integer from 0 to the maximum unsigned integer in C++
	# In C++, the maximum value for an unsigned int is typically 2^32 - 1
	max_unsigned_int_cpp = 2**32 - 1
	random_unsigned_int = random.randint(0, max_unsigned_int_cpp)
	return random_unsigned_int



def make_rand_vector(dims):
	vec = [random.gauss(0, 1) for i in range(dims)]
	mag = sum(x**2 for x in vec) ** .5
	return [x/mag for x in vec]
