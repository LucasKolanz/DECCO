#Removes data files that are in the form {index}_{datatype}.csv
#from all leaf directories of a given root directory

import os
import glob
import numpy as np
import subprocess
# import check_for_errors as cfe

# import utils as u
import h5py
import json

relative_path = ""
relative_path = '/'.join(__file__.split('/')[:-1]) + '/' + relative_path
project_path = os.path.abspath(relative_path) + '/'

def get_new_data_files(fullpath):
	directory = os.fsencode(fullpath)
	num_olds = 0
	new_files = []
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("simData.csv"):
			if filename.count('_') > 1:
				num_olds += 1
			else:
				new_files.append(fullpath+filename)
				new_files.append(fullpath+filename.replace("simData.csv","constants.csv"))
				new_files.append(fullpath+filename.replace("simData.csv","energy.csv"))

	if len(new_files)/3 < num_olds:
		return new_files
	else:
		return []

def all_old(fullpath):
	directory = os.fsencode(fullpath)
	num_olds = 0
	num_news = 0
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		if filename.endswith("simData.csv"):
			if filename.count('_') > 1:
				num_olds += 1
			else:
				num_news += 1

	if num_olds == 0 and num_news > 0:
		print(num_news)
		return True
	else:
		return False

def new_files(job_base,\
				N=[30,100,300],\
				Temps=[3,10,30,100,300,1000],\
				num_attempts=30):

	errors = []
	if isinstance(num_attempts,int):
		attempts = [i for i in range(num_attempts)]
	elif isinstance(num_attempts,list):
		attempts = num_attempts

	files = []
	for n in N:
		for Temp in Temps:
			for attempt in attempts:
				job = job_base.replace("$a$",str(attempt)).replace("$n$",str(n)).replace("$t$",str(Temp))
				if os.path.exists(job):
					files.extend(get_new_data_files(job))
					# print(output)

	# print(f"{len(errors)} errors, out of {valid_count} valid runs, out of {len(N)*len(attempts)*len(Temps)} runs.")
	return files



def all_olds(job_base,\
				N=[30,100,300],\
				Temps=[3,10,30,100,300,1000],\
				num_attempts=30):

	if isinstance(num_attempts,int):
		attempts = [i for i in range(num_attempts)]
	elif isinstance(num_attempts,list):
		attempts = num_attempts

	files = []
	for n in N:
		for Temp in Temps:
			for attempt in attempts:
				job = job_base.replace("$a$",str(attempt)).replace("$n$",str(n)).replace("$t$",str(Temp))
				if os.path.exists(job):
					if all_old(job):
						files.append(job)

	# print(f"{len(errors)} errors, out of {valid_count} valid runs, out of {len(N)*len(attempts)*len(Temps)} runs.")
	return files

def main():

	with open(project_path+"default_files/default_input.json",'r') as fp:
		input_json = json.load(fp)



	job_folder = 'jobs/'###FOR COSINE

	job = input_json["data_directory"] + job_folder + 'lognorm$a$/N_$n$/T_$t$/'

	attempts = [i for i in range(30)]
	# attempts = [0,1,6]

	N = [30,100,300]
	# N=[300]

	Temps = [3,10,30,100,300,1000]
	# Temps = [3]


	all_old_file_folders = all_olds(job,N,Temps,attempts)
	print(all_old_file_folders)

	# remove_files = new_files(job,N,Temps,attempts) #out of the errored jobs, which ones do not have error 4

	# if len(remove_files) == 0:
	# 	print("No moves found")

	for folder in all_old_file_folders:

		command = f"rm {folder}*"
		# os.system(command)
		print(command)
		# exit(0)


if __name__ == '__main__':
	main()