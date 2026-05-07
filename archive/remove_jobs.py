


import os
import glob
import numpy as np
# import utils as u
import h5py
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'


def get_file_base(folder):
	directory = os.fsencode(folder)
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("simData.csv"):
			file_base = '_' + '_'.join(filename.split('/')[-1].split('_')[1:-1]) + '_'
			return file_base
	return "ERROR: NO FILE FOUND"; 

def get_folders(job_base,\
				N=[30,100,300],\
				Temps=[3,10,30,100,300,1000],\
				num_attempts=30,relax=False):

	folders = []
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
				folders.append(job)
				# if os.path.exists(job):
				# 	if os.path.exists(job+"timing.txt"):
				# 		valid_count += 1
				# 		folders.append(job)
				# else:
				# 	continue
					# print(f"{job} doesn't exist")

	return folders

def main():

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)

	
	job = input_json["data_directory"] + 'jobs/const$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobsNovus/const_relax$a$/N_$n$/T_$t$/'
	job = input_json["data_directory"] + 'jobsCosine/lognorm_relax$a$/N_$n$/T_$t$/'
	print(job)

	move_to = input_json['data_directory'] + "erroredRelaxJobs/"

	relax = False
	if job.split('/')[-4].split("_")[-1].strip("$a$") == "relax":
		relax = True


	# attempts = [i for i in range(30)]
	# attempts = [1]

	# N = [30,100,300]
	# N=[30,100]

	# Temps = [3,10,30,100,300,1000]
	# Temps = [3]


	folders = get_folders(job,N,Temps,attempts,relax=relax)
	print(folders)

	for folder in folders:
		full_move_to = move_to + "/".join(folder.split("/")[-4:])
		if not os.path.exists(full_move_to):
			os.makedirs(full_move_to)
		command = f"mv {folder}* {full_move_to}." 
		# command = f"mv {full_move_to}* {folder}." 
		print(command)
		# os.system(command)


if __name__ == '__main__':
	main()