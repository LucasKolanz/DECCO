"""
This file was originally written for SpaceLab/DECCO to check for errors in the output of csv files.

Author: Lucas Kolanz

This file checks the output in all folders with a specified pattern to check for several errors, defined below. Adding
errors to this file should be done by defining a function that, given a folder, will check if output files in that
folder have the error. This function should then be added to the for loop at the bottom of the file.

"""




##Check output of simulations for errors

#Error -1: sim_err.log has ERROR in its tail output.
#
#Error 0: Is the sim finished as indicated by the presence of timing.txt, but also lacking in "Simulation Complete!" in the output
#
#Error 1: number of balls is incorrect in at least 1 sim_data output. 
#
#Error 2: incorrect number of rows in simData
#
#Error 3: do the first two balls of the sim overlap eachother?
#
#Error 4: integer overflow in number of steps
#
#Error 5: Are any balls within a few radii of the center of mass? If not, balls might not be contacting the growing aggregate
#
#Error 6: Are any balls not touching the aggregate at the end? (Only applicable to BPCA growth as an error)
#			It seems to be possible for a ball to not be touching at the moment of the end of the simulation, but still
#			be a part of the aggregate if the aggregate isn't totally relaxed. To make sure this is the case and 
#			it isn't an actual error6, load the aggregate in Blender, select all balls and move them all into frame.
#			If a ball is outside the aggregate, you will then be able to tell. 


import os
import glob
import numpy as np
# import utils as u
import h5py
import sys
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

simData_properties = 11
# sys.path.append(project_path+"utilities/")
# from utils import index_from_file
# sys.path.append(project_path)
# import porosity_FD as pFD


def check_final_error(job_base,error,\
				N=[30,100,300],\
				Temps=[3,10,30,100,300,1000],\
				num_attempts=30):

	errors = []
	if isinstance(num_attempts,int):
		attempts = [i for i in range(num_attempts)]
	elif isinstance(num_attempts,list):
		attempts = num_attempts
	# N=[30]
	# Temps = [3]
	valid_count = 0

	for n in N:
		for Temp in Temps:
			for attempt in attempts:
				job = job_base.replace("$a$",str(attempt)).replace("$n$",str(n)).replace("$t$",str(Temp))
				if os.path.exists(job+"timing.txt"):
					file_base = get_file_base(job) 
							
					
					valid_count += 1
					output = error(job+str(n-1)+file_base + "simData.csv")
					if output > 0:
						errors.append(job)

	print(f"{len(errors)} errors, out of {valid_count} valid runs, out of {len(N)*len(attempts)*len(Temps)} runs.")
	return errors


def where_is_smallest_error2(job_base,error,\
				N=[30,100,300],\
				Temps=[3,10,30,100,300,1000],\
				num_attempts=30):

	errors = []
	output = []
	if isinstance(num_attempts,int):
		attempts = [i for i in range(num_attempts)]
	elif isinstance(num_attempts,list):
		attempts = num_attempts
	# N=[30]
	# Temps = [3]
	valid_count = 0

	for n in N:
		for Temp in Temps:
			for attempt in attempts:
				job = job_base.replace("$a$",str(attempt)).replace("$n$",str(n)).replace("$t$",str(Temp))
				if os.path.exists(job+"timing.txt"):
					valid_count += 1
					out = error(job)
					if out > 0:
						output.append(out)
						# print(f"{job} errored at {output}")
						errors.append(job)

	print(f"{len(errors)} errors, out of {valid_count} valid runs, out of {len(N)*len(attempts)*len(Temps)} runs.")
	return errors,output

#returns the last line of the file file_path
def tail(file_path,n):
	count = 0
	with open(file_path, 'rb') as f:
		try:  # catch OSError in case of a one line file 
			f.seek(-2, os.SEEK_END)
			while count <= n:
				temp = f.read(1)
				if temp == b'\n':
					count += 1
				if count < n:
					f.seek(-2, os.SEEK_CUR)
				# else:
				# 	f.seek(-1, os.SEEK_CUR)
		except OSError:
			print("OSERROR")
			f.seek(0)
		last_lines = f.read().decode()
	return last_lines


def ck_error2_by_file(file,verbose=False,relax=False):
	try:
		temp = np.loadtxt(file,skiprows=1,delimiter=',')
		if verbose:
			print(f"shape of file: {temp.shape}")
		x,y = temp.shape
		if relax:
			if not (((x-1)*1.0)/50.0).is_integer():
				if verbose:
					print(f"ERROR: The x dimension of the data is {x}, n in (n * 50)+1 is not an integer for file {file}")
				return True
		else:
			if x != 51:
				if verbose:
					print(f"ERROR: The x dimension of the data is {x}, not 51. for file {file}")
				return True
	except ValueError as E:
		if verbose:
			print(f"ERROR ValueError: {E}")
			print(f"for the file: {file}")
		return True
	return False

def index_from_file(file):
	file_split = file.split("_")

	old_name_scheme = (len(file_split)>2)

	if old_name_scheme and not file_split[1].isnumeric():
		return 0
	else:
		return int(file_split[0])


def error2_index(fullpath,verbose=False):
	directory = os.fsencode(fullpath)
	lowest_index = 9999999999
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("simData.csv"): 
			try:
				temp = np.loadtxt(fullpath+filename,skiprows=1,delimiter=',')
				x,y = temp.shape
				ind = index_from_file(filename)
				print(f"ind={ind}, x={x}, y={y}")
				if x != 51:
					if verbose:
						print(f"ERROR: x={x} for file={filename}")
					if ind < lowest_index:
						lowest_index = ind
			except ValueError as E:
				ind = index_from_file(filename)
				print(f"otherind={ind}, x={x}, y={y}")
				if verbose:
					print(f"ERROR: {E} for file={filename}")
				if ind < lowest_index:
					lowest_index = ind
	if lowest_index != 9999999999:
		return lowest_index
	return False

def where_did_error1_start(fullpath):
	directory = os.fsencode(fullpath)
	min_ind=9999999999
	min_balls = 9999999999
	test_file = ''
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("constants.csv") and filename[2] != "R": 
			index = int(filename.split('_')[0])
			with open(fullpath+filename, 'r') as fp: #number of lines in this file is the number of balls in sim
				for count, line in enumerate(fp):
					pass
			balls = count+1 #IDK why you need to add one but this method doesn't count every line, it misses one at the beginning or end
			if (balls != index+3):
				if (index < min_ind):
					min_balls = balls
					min_ind = index
					test_file = filename

	if len(test_file) == 0:
		print(f"Error1 NOT detected in {fullpath}")
	else:
		print(f"Error1 detected in {fullpath}")
		print(f"Min index {min_ind} has {min_balls} balls")

	
def check_error(job_base,error,\
				N=[30,100,300],\
				Temps=[3,10,30,100,300,1000],\
				num_attempts=30,relax=False):

	errors = []
	if isinstance(num_attempts,int):
		attempts = [i for i in range(num_attempts)]
	elif isinstance(num_attempts,list):
		attempts = num_attempts
	# N=[30]
	# Temps = [3]
	valid_count = 0

	for n in N:
		for Temp in Temps:
			for attempt in attempts:
				job = job_base.replace("$a$",str(attempt)).replace("$n$",str(n)).replace("$t$",str(Temp))
				# print(job)
				if os.path.exists(job):
					if os.path.exists(job+"timing.txt"):
						valid_count += 1
					output = error(job,relax=relax)
					if output > 0:
						errors.append(job)
				else:
					continue
					# print(f"{job} doesn't exist")

	print(f"{len(errors)} errors, out of {valid_count} valid runs, out of {len(N)*len(attempts)*len(Temps)} runs.")
	if valid_count == 0:
		print("WARNING: valid_count is zero.")
	return errors

def get_file_base(folder):
	directory = os.fsencode(folder)
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("simData.csv"):
			file_base = '_' + '_'.join(filename.split('/')[-1].split('_')[1:-1]) + '_'
			return file_base
	return "ERROR: NO FILE FOUND"; 

def dist(x1,y1,z1,x2,y2,z2):
	return np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

def errorn1(fullpath,relax=False):
	has_err = os.path.exists(fullpath+"sim_err.log")
	has_errors = os.path.exists(fullpath+"sim_errors.txt")

	error_file = ''
	if has_err:
		error_file = "sim_err.log"
	elif has_errors:
		error_file = "sim_errors.txt"
	else:
		print(f"NO ERROR FILE FOR {fullpath}")
		return False

	tail_out = tail(fullpath+error_file,10).split('\n')
	error = False
	for i in tail_out:
		if "ERROR" in i:
			error = True
			break 
	return error

def error0(fullpath,relax=False):
	num_files = len(glob.glob(fullpath+'*'))

	if os.path.exists(fullpath+"timing.txt") and num_files > 1:
		has_err = os.path.exists(fullpath+"sim_err.log")
		has_errors = os.path.exists(fullpath+"sim_errors.txt")

		error_file = ''
		if has_err:
			error_file = "sim_err.log"
		elif has_errors:
			error_file = "sim_errors.txt"
		else:
			print(f"NO ERROR FILE FOR {fullpath}")
			# return True

		tail_out = tail(fullpath+error_file,10).split('\n')
		for i in tail_out:
			if "Simulation complete!" in i or "Simulation already complete." in i:
				return False
		return True
	else:
		return False



#If fullpath has error1 in it, return 1, if not return 0
def error1(fullpath,relax=False):
	directory = os.fsencode(fullpath)
	max_ind=-1
	test_file = ''
	rel = ""
	if relax:
		rel = "RELAX"
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith(f"{rel}constants.csv"): 

			index = int(filename.split('_')[0])
			if (index > max_ind):
				max_ind = index
				test_file = filename
	# print(fullpath+test_file)

	if (not test_file): #if test_file is empty string the sim hasnt started yet
		return False

	with open(fullpath+test_file, 'r') as fp: #number of lines in this file is the number of balls in sim
		for count, line in enumerate(fp):
			pass
	balls = count+1 #IDK why you need to add one but this method doesn't count every line, it misses one at the beginning or end

	if balls == max_ind+3: ### THIS IS SPECIFIC TO BPCA GROWTH RUNS
		return False
	else:
		# print("balls in sim                : {}".format(balls))
		# print("balls that should be in sim : {}".format(max_ind+3))
		# print("initially specified N value : {}".format(correct_N))
		return True

#If the (number rows of constants)*11 != (simData columns)
#Then there was some kind of write error
def error2(fullpath,verbose=False,relax=False):
	if os.path.exists(fullpath+"timing.txt"):
		rel = ""
		if relax:
			rel = "RELAX"
		directory = os.fsencode(fullpath)
		for file in os.listdir(directory):
			filename = os.fsdecode(file)
			if filename.endswith(f"{rel}simData.csv"): 
				err = ck_error2_by_file(fullpath+filename,verbose,relax=relax)
				if err:
					return True
	return False

def error3(fullpath,relax=False):
	directory = os.fsencode(fullpath)
	rel = ""
	if relax:
		rel = "RELAX"
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith(f"{relax}simData.csv"):
			if filename.startswith("0_") or filename.split("_")[1][0] == "R":
				simData = np.loadtxt(fullpath+filename,skiprows=1,delimiter=',',dtype=np.float64)[0]
				constants = np.loadtxt(fullpath+filename.replace("simData.csv","constants.csv"),skiprows=0,delimiter=',',dtype=np.float64)
				radii = constants[:,0]
				if (dist(simData[0],simData[1],simData[2],simData[11],simData[12],simData[13]) <= radii[0]+radii[1]):
					return True
				else:
					return False
				
	return False

def error4(fullpath,relax=None):

	has_err = os.path.exists(fullpath+"sim_err.log")
	has_errors = os.path.exists(fullpath+"sim_errors.txt")

	error_file = ''
	if has_err:
		error_file = "sim_err.log"
	elif has_errors:
		error_file = "sim_errors.txt"
	else:
		print(f"NO ERROR FILE FOR {fullpath}")
		return False

	tail_out = tail(fullpath+error_file,10).split('\n')
	error = False
	for i in tail_out:
		if "ERROR: STEPS IS NEGATIVE" in i:
			error = True
			break 
	return error



def error5(fullpath,relax=False):
	num_files = len(glob.glob(fullpath+'*'))

	if os.path.exists(fullpath+"timing.txt") and num_files > 1:
		directory = os.fsencode(fullpath)
		max_index = -1
		max_filename = ""

		rel = ""
		if relax:
			rel = "RELAX"
		#find the highest index file
		for file in os.listdir(directory):
			filename = os.fsdecode(file)
			if filename.endswith(f"{rel}simData.csv"):
				if filename.startswith("0_") or (not relax and filename.split("_")[1][0] == "R"): #this is for index zero
					index = 0 
				else:
					index = int(filename.split("_")[0])
				
				if index > max_index:
					max_index = index
					max_filename = filename

		simData = np.loadtxt(fullpath+max_filename,skiprows=1,delimiter=',',dtype=np.float64)[-1]
		simData = simData.reshape(int(simData.shape[0]/simData_properties),simData_properties) #now it is simData[ball,property]
		pos = simData[:,:3]
		constants = np.loadtxt(fullpath+max_filename.replace("simData.csv","constants.csv"),skiprows=0,delimiter=',',dtype=np.float64)
		radii = constants[:,0]
		mass = constants[:,1]
		com = get_COM(pos,mass)
		# data,radii,mass,moi = u.get_data(fullpath,max_index,linenum=-1,relax=relax)
		# r_g = pFD.get_gyration_radius(fullpath,max_index)


		for i in range(radii.shape[0]):
			if (dist(pos[i,0],pos[i,1],pos[i,2],com[0],com[1],com[2]) <= 2*radii[i]): #if a ball's center is closer or equal to two of its radii to the aggregates center of mass. I take this to mean the aggregate is aggregating correctly and not shooting all over the place  
				return False
		return True
	return False


def are_spheres_connected(pos, radii):
	from collections import deque

	def distance(a, b):
		return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5

	def is_connected(index1, index2):
		return distance(pos[index1], pos[index2]) <= 2*radii[index1] + 2*radii[index2]

	n = len(pos)
	if n == 0:
		return False

	# Build the adjacency list
	adj_list = {i: [] for i in range(n)}
	for i in range(n):
		for j in range(i + 1, n):
			if is_connected(i, j):
			    adj_list[i].append(j)
			    adj_list[j].append(i)

	# BFS to check connectivity
	visited = [False] * n
	queue = deque([0])
	visited[0] = True
	while queue:
		node = queue.popleft()
		for neighbor in adj_list[node]:
			if not visited[neighbor]:
				visited[neighbor] = True
				queue.append(neighbor)

	return all(visited)


def error6(fullpath,relax=None):
	num_files = len(glob.glob(fullpath+'*'))

	if os.path.exists(fullpath+"timing.txt") and num_files > 1:

		directory = os.fsencode(fullpath)
		max_index = -1
		max_filename = ""

		rel = ""
		if relax:
			rel = "RELAX"
		#find the highest index file
		for file in os.listdir(directory):
			filename = os.fsdecode(file)
			if filename.endswith(f"{rel}simData.csv"):
				if filename.startswith("0_") or (not relax and filename.split("_")[1][0] == "R"): #this is for index zero
					index = 0 
				else:
					index = int(filename.split("_")[0])
				
				if index > max_index:
					max_index = index
					max_filename = filename

		simData = np.loadtxt(fullpath+max_filename,skiprows=1,delimiter=',',dtype=np.float64)[-1]
		simData = simData.reshape(int(simData.shape[0]/simData_properties),simData_properties) #now it is simData[ball,property]
		pos = simData[:,:3]
		constants = np.loadtxt(fullpath+max_filename.replace("simData.csv","constants.csv"),skiprows=0,delimiter=',',dtype=np.float64)
		radii = constants[:,0]
		# mass = constants[:,1]
		
		connected = are_spheres_connected(pos,radii)

		if connected:
			return False
		else:
			return True

	return False




def get_COM(pos,mass):
	
	com = np.array([0,0,0],dtype=np.float64)
	mtot = 0

	for ball in range(pos.shape[0]):
		com += mass[ball]*pos[ball]
		mtot += mass[ball]

	return com/mtot



def main():


	# folder = "/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/lognorm6/N_100/T_100/"
	# folder = "/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/lognorm24/N_300/T_3/"
	# ind = where_did_error1_start(folder)
	# ind = error2_index(folder,True)
	# print(ind)


	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)


	job = input_json["data_directory"] + 'jobsCosine/lognorm_relax$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobs/lognorm$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobsCosine/lognorm$a$/N_$n$/T_$t$/'
	

	attempts = [i for i in range(30)]
	# attempts = [0]

	N = [30,100,300]
	# N=[100]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	errorDic = {}

	relax = False

	if job.split("/")[-4].split("_")[-1].strip("$a$") == "relax":
		relax = True



	# for i,error in enumerate([error6]):
	for i,error in enumerate([errorn1,error0,error1,error2,error3,error4,error5,error6]):
		print(f"======================================{error.__name__}======================================")
		error_folders = check_error(job,error,N,Temps,attempts,relax=relax)
		for folder in error_folders:
			if folder in errorDic.keys():
				errorDic[folder].append(i)
			else:
				errorDic[folder] = [f"{error.__name__}"]

	if not errorDic:
		print("No Errors detected")
	else:
		for key in errorDic.keys():
			print(f"Errors in folder {key}")
			for error in errorDic[key]:
				print(f"\t{error}")





























	# errorgen_folders = check_error(job,error_general,N,Temps,attempts)
	# print(errorgen_folders)
	# error1_folders = check_error(job,error1,N,Temps,attempts)
	# print(error1_folders)

	# error2_folders = check_error(job,error2,N,Temps,attempts)
	# print(error2_folders)

	###
	### This section finds the index at which a folder had error2
	###
	# error2_folders,o = where_is_smallest_error2(job,error2_index,N,Temps,attempts)
	# print(error2_folders)
	# print(o)
	# zipped = zip(error2_folders,o)
	# for i in zipped:
	# 	print(f"Error started at index {i[1]} for folder {i[0]}")
	# mind = np.argmin(np.array(o))
	# print(f'min index is {mind} for job {error2_folders[0][mind]}')


	# folder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/erroredJobs/lognorm12/N_30/T_10/'
	# folder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/error2Test10/N_10/T_10/'
	# folder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/error2Test2/N_10/T_10/'

	# print(error2_index(folder,True))
	# print(error2(folder,True))

if __name__ == '__main__':
	main()