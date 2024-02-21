##Check output of simulations for errors

#Error -1: sim_err.log has ERROR in its tail output.
#
#Error 0: Is the sim finished as indicated by the presence of timing.txt, but also lacking in "Simulation Complete!" in the output
#
#Error 1: Are there any nan values in any of the data when timing.txt exists?
#
#Error 2: Are there any nan values in sims that aren't the highest index, when timing.txt doesnt exist
#
#Error 3: Metadata check. Is there any metadata associated with the largest index file if timing.txt DOESN'T exists?
#		  This is really more of a warning as it means the code won't restart from an ideal spot if stopped.
#
#Error 4: integer overflow in number of steps
#
#Error 5: did we get "Simulation complete!" within the last 10 lines of sim_error.log and timing.txt exists?
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
import json
import check_for_errors as cfe

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'



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
		except OSError as e:
			print(f"OSERROR {e} on tail line {n} for file {file_path}")
			f.seek(0)
		last_lines = f.read().decode()
	return last_lines

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
	if os.path.exists(fullpath+"timing.txt"):
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
			if "Simulation complete!" in i:
				error = False
				break 
		
	return False

#If fullpath has error1 in it, return 1, if not return 0
def error1(fullpath,relax=False):
	if os.path.exists(fullpath+"timing.txt"):
		# Loop through all files in the directory fullpath
		for filename in os.listdir(fullpath):
			if filename.endswith(".h5"):
				filepath = os.path.join(fullpath, filename)

				# Open the HDF5 file
				with h5py.File(filepath, 'r') as file:
					simData = np.array(file['/simData'][:])
					simDatanans = np.sum(np.where(np.isnan(simData),1,0))
					if simDatanans > 0:
						
						print(f"simData has {simDatanans} nans: " + filepath)
						return True

					energy = np.array(file['/energy'][:])
					energynans = np.sum(np.where(np.isnan(energy),1,0))
					if energynans > 0:
						print("energy: " + filepath)
						return True

					constants = np.array(file['/constants'][:])
					constantsnans = np.sum(np.where(np.isnan(constants),1,0))
					if constantsnans > 0:
						print("constants: " + filepath)
						return True
	return False
	
#If fullpath has error2 in it, return 1, if not return 0
def error2(fullpath,relax=False):
	if not os.path.exists(fullpath+"timing.txt"):
		max_ind = get_max_ind(fullpath)
		# Loop through all files in the directory fullpath
		for filename in os.listdir(fullpath):
			if filename.endswith(".h5") and not filename.startswith(f"{max_ind}_"):
				filepath = os.path.join(fullpath, filename)

				# Open the HDF5 file
				with h5py.File(filepath, 'r') as file:
					simData = np.array(file['/simData'][:])
					simDatanans = np.sum(np.where(np.isnan(simData),1,0))
					if simDatanans > 0:
						
						print(f"simData has {simDatanans} nans: " + filepath)
						return True

					energy = np.array(file['/energy'][:])
					energynans = np.sum(np.where(np.isnan(energy),1,0))
					if energynans > 0:
						print("energy: " + filepath)
						return True

					constants = np.array(file['/constants'][:])
					constantsnans = np.sum(np.where(np.isnan(constants),1,0))
					if constantsnans > 0:
						print("constants: " + filepath)
						return True
	return False

#If fullpath has error3 in it, return 1, if not return 0
def error3(fullpath,relax=False):
	if os.path.exists(fullpath+"timing.txt"):
		max_ind = get_max_ind(fullpath)
		file = fullpath+str(max_ind)+"_data.h5"
		with h5py.File(file, 'r') as f:
			dataset = f['constants']
			metadata = {attr: dataset.attrs[attr] for attr in dataset.attrs}
			if len(metadata) > 0:
				return False
			else:
				return True
	return False


def error4(fullpath,relax=False):
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
	if os.path.exists(fullpath+"timing.txt"):
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
		if "Simulation complete!" in tail_out:
			return False
		else:
			return True
	else:
		return False



def error6(fullpath,relax=None):
	if (os.path.exists(fullpath+"timing.txt")):
		directory = os.fsencode(fullpath)

		N = int(fullpath.split("/")[-3].split("_")[-1])
		rel = ""
		if relax:
			rel = "RELAX"
		#find the highest index file
		for file in os.listdir(directory):
			filename = os.fsdecode(file)
			if filename.endswith(f"{rel}data.h5"):

				index = int(filename.split("_")[0])
				
				if index == N-3:
								
					with h5py.File(fullpath+filename, 'r') as file:
						constants = np.array(file['/constants'][:])
						radii = constants[np.where(np.arange(constants.shape[0])%3==0)]
						num_spheres = radii.shape[0]

						simData_single_ball_width = 11

						simData = np.array(file['/simData'][-simData_single_ball_width*num_spheres:])
						simData = simData.reshape(num_spheres,simData_single_ball_width)
						pos = simData[:,:3]

						connected = cfe.are_spheres_connected(pos,radii)

						if connected:
							return False
						else:
							return True

	return False



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
	return errors

def get_max_ind(fullpath):
	max_ind = -1
	for filename in os.listdir(fullpath):
		if filename.endswith(".h5"):
			file_ind = int(filename.split("_")[0])
			if file_ind > max_ind:
				max_ind = file_ind
	return max_ind


def main():

	curr_folder = os.getcwd() + '/'

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	job = curr_folder + 'jobs/tempVarianceRand_attempt$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobs/lognorm$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobs/weakseed$a$/N_$n$/T_$t$/'
	job = curr_folder + 'erroredJobs/lognorm$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobsNovus/testError$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobs/const$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobsNovus/const$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobsNovus/const_relax$a$/N_$n$/T_$t$/'
	print(job)


	attempts = [i for i in range(30)]
	# attempts = [0]

	N = [30,100,300]
	# N=[100]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]

	errorDic = {}


	# for i,error in enumerate([errorn1,error0,error1,error2,error3,error4]):
	for i,error in enumerate([error6]):
		print(f"======================================{error.__name__}======================================")
		error_folders = check_error(job,error,N,Temps,attempts)
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


	# error_folders = check_error(job,error1,N,Temps,attempts)
	# error_folders.extend(check_error(job,error2,N,Temps,attempts))
	# error_folders.extend(check_error(job,error4,N,Temps,attempts))

	# print(error_folders)

	# errorgen_folders = check_error(job,error_general,N,Temps,attempts)
	# print(errorgen_folders)

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