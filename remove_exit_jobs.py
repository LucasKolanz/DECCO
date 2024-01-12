##Check output of simulations for errors

#exit_error: Is there a timing.txt file and only 1 h5 file?


import os
import glob
import numpy as np
# import utils as u
import h5py






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
				
				output = error(job)
				if output > 0:
					if not os.path.exists(job+"timing.txt"):
						print(f"Sim started but not finished: {job}")
					else:
						valid_count += 1
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
	job = '/home/kolanzl/novus/kolanzl/SpaceLab_data/jobs/const$a$/N_$n$/T_$t$/'



	attempts = [i for i in range(30)]
	# attempts = [0]

	N = [30,100,300]
	# N=[5]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]


	remove_folders = check_error(job,error1,N,Temps,attempts)
	# print(remove_folders)

	for folder in remove_folders:
		command = f"rm {folder}*"
		print(command)
		# os.system(command)


if __name__ == '__main__':
	main()