##Check output of simulations for errors

#Error 1: Are there any nan values in any of the data?
#Error general: did we get "Simulation complete!" within the last 10 lines of sim_error.log

import os
import glob
import numpy as np
import utils as u
import h5py

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



#If the (number rows of constants)*11 != (simData columns)
#Then there was some kind of write error
def error2(fullpath,verbose=False):
	directory = os.fsencode(fullpath)
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("simData.csv"): 
			err = ck_error2_by_file(fullpath+filename,verbose)
			if err:
				return True
	return False



#If fullpath has error1 in it, return 1, if not return 0
def error1(fullpath):
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
					exit(0)
					return True

				energy = np.array(file['/energy'][:])
				energynans = np.sum(np.where(np.isnan(energy),1,0))
				if energynans > 0:
					print("energy: " + filepath)
					exit(0)
					return True

				constants = np.array(file['/constants'][:])
				constantsnans = np.sum(np.where(np.isnan(constants),1,0))
				if constantsnans > 0:
					print("constants: " + filepath)
					exit(0)
					return True
	return False


def error_general(fullpath):
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



def check_error(job_base,error,\
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
					valid_count += 1
					output = error(job)
					if output > 0:
						errors.append(job)

	print(f"{len(errors)} errors, out of {valid_count} valid runs, out of {len(N)*len(attempts)*len(Temps)} runs.")
	return errors

def get_file_base(folder):
	directory = os.fsencode(folder)
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("simData.csv"):
			file_base = '_' + '_'.join(filename.split('/')[-1].split('_')[1:-1]) + '_'
			return file_base
	return "ERROR: NO FILE FOUND"; 



def main():

	curr_folder = os.getcwd() + '/'

	job = curr_folder + 'jobs/tempVarianceRand_attempt$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobs/lognorm$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobs/weakseed$a$/N_$n$/T_$t$/'
	job = curr_folder + 'erroredJobs/lognorm$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobsNovus/testError$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobsNovus/const$a$/N_$n$/T_$t$/'



	attempts = [i for i in range(30)]
	# attempts = [0]

	N = [30,100,300]
	# N=[5]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]


	error1_folders = check_error(job,error1,N,Temps,attempts)
	print(error1_folders)

	errorgen_folders = check_error(job,error_general,N,Temps,attempts)
	print(errorgen_folders)

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