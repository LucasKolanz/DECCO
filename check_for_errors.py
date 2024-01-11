##Check output of simulations for errors

#Error 1: restart failed, number of balls didn't carry over so there are an incorrect 
#			number of colums in at least 1 sim_data output
#Error 2: sim fails to write all the data (so there aren't the correct number of colums in simData)
#Error *: possible signed integer over flow in number of steps for a sim
#			Indicated by a specific sequence at the end of sim_errors.txt
#Error general: did we get "Simulation complete!" within the last 10 lines of sim_error.log

import os
import glob
import numpy as np
import utils as u

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

# #works like command line function tail
# def tail( file_path, lines=10 ):
# 	with open(file_path, 'rb') as f:
# 	    total_lines_wanted = lines

# 	    BLOCK_SIZE = 1024
# 	    f.seek(0, 2)
# 	    block_end_byte = f.tell()
# 	    lines_to_go = total_lines_wanted
# 	    block_number = -1
# 	    blocks = [] # blocks of size BLOCK_SIZE, in reverse order starting
# 	                # from the end of the file
# 	    while lines_to_go > 0 and block_end_byte > 0:
# 	        if (block_end_byte - BLOCK_SIZE > 0):
# 	            # read the last block we haven't yet read
# 	            f.seek(block_number*BLOCK_SIZE, 2)
# 	            blocks.append(f.read(BLOCK_SIZE))
# 	        else:
# 	            # file too small, start from begining
# 	            f.seek(0,0)
# 	            # only read what was not read
# 	            blocks.append(f.read(block_end_byte))
# 	        lines_found = blocks[-1].count('\n')
# 	        lines_to_go -= lines_found
# 	        block_end_byte -= BLOCK_SIZE
# 	        block_number -= 1
# 	    all_read_text = ''.join(reversed(blocks))
# 	    return '\n'.join(all_read_text.splitlines()[-total_lines_wanted:])

def ck_error2_by_file(file,verbose=True ):
	try:
		temp = np.loadtxt(file,skiprows=1,delimiter=',')
		if verbose:
			print(f"shape of file: {temp.shape}")
		x,y = temp.shape
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

#If the (number rows of constants)*11 != (simData columns)
#Then there was some kind of write error
def error2(fullpath,verbose=False):
	directory = os.fsencode(fullpath)
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("simData.csv"): 
			err = ck_error2_by_file(filename,verbose)
			if err:
				return True
	return False

def error2_index(fullpath,verbose=False):
	directory = os.fsencode(fullpath)
	lowest_index = 9999999999
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("simData.csv"): 
			try:
				temp = np.loadtxt(fullpath+filename,skiprows=1,delimiter=',')
				x,y = temp.shape
				if x != 51:
					ind = u.index_from_file(filename)
					if verbose:
						print(f"ERROR: x={x} for file={filename}")
					if ind < lowest_index:
						lowest_index = ind
			except ValueError as E:
				ind = u.index_from_file(filename)
				if verbose:
					print(f"ERROR: {E} for file={filename}")
				if ind < lowest_index:
					lowest_index = ind
	if lowest_index != 9999999999:
		return lowest_index
	return False


#If fullpath has error1 in it, return 1, if not return 0
def error1(fullpath):
	directory = os.fsencode(fullpath)
	max_ind=-1
	test_file = ''
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("constants.csv"): 

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

def main():

	curr_folder = os.getcwd() + '/'

	job = curr_folder + 'jobs/tempVarianceRand_attempt$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobs/lognorm$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobs/weakseed$a$/N_$n$/T_$t$/'
	job = curr_folder + 'jobsCosine/lognorm$a$/N_$n$/T_$t$/'
	job = curr_folder + 'erroredJobs/lognorm$a$/N_$n$/T_$t$/'



	attempts = [i for i in range(30)]
	# attempts = [0]

	N = [30,100,300]
	# N=[30]

	Temps = [3,10,30,100,300,1000]
	# Temps = [1000]

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
	folder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/error2Test10/N_10/T_10/'
	folder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/error2Test2/N_10/T_10/'

	print(error2_index(folder,True))
	print(error2(folder,True))

if __name__ == '__main__':
	main()